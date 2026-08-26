#!/usr/bin/env Rscript
# =============================================================================
# data-raw/liver_1.8/01_build_network_1.8.R
# -----------------------------------------------------------------------------
# Construit data/liver_1.8_network.rda = liver_1.7 dont la COUCHE PATHWAY est
# restreinte aux pathways Reactome hepatotoxiques. Pur elagage de noeuds :
# liver_1.8 == sous-graphe induit de liver_1.7 sur les noeuds gardes (tous les
# autres noeuds/aretes/attributs sont inchanges).
#
# Source de l'elagage = mapping pathway-level PUBLIE de Tsai et al. 2025
# (Toxicol Sci, doi:10.1093/toxsci/kfaf036), Supplemental File 1 =
# ressources/toxsci-25-0020-File008.xlsx. Assignation pathway -> Key
# Characteristic par 2 scientifiques + consensus. On garde les pathways
# assignes a l'un des 11 termes umbrella des 12 KC des HEPATOTOXICANTS
# (Rusyn et al. 2021 ; domaine "HEP" du crosswalk Tsai).
#
# Noms de pathway -> R-HSA via ressources/ReactomePathways.txt.
# Dependance dev : readxl (parse xlsx) ; non ajoutee aux Imports du package.
#
# Executer depuis la racine du package CENetwork.
# =============================================================================

suppressMessages({library(readxl); library(igraph); library(tidyverse)})

PKG   <- "/Users/antoine/Documents/DEV/CENetwork"
RES   <- file.path(PKG, "data-raw/liver_1.8/ressources")
F008  <- file.path(RES, "toxsci-25-0020-File008.xlsx")
NAMEF <- file.path(RES, "ReactomePathways.txt")
NET17 <- file.path(PKG, "data/liver_1.7_network.rda")
OUT   <- file.path(PKG, "data/liver_1.8_network.rda")
TSV   <- file.path(PKG, "data-raw/liver_1.8/hepatotox_pathways.tsv")

# 11 termes umbrella des KC hepatotoxicants (Rusyn 2021 ; KC2+KC3 fusionnent)
HEP_UMBRELLA <- c(
  "IS METABOLICALLY ACTIVATED OR REACTIVE",                        # KC1
  "ALTERS CELL PROLIFERATION AND CELL DEATH",                      # KC2, KC3
  "DISRUPTS TRANSPORT FUNCTION",                                   # KC4
  "INDUCES OXIDATIVE STRESS",                                      # KC5
  "CAUSES INFLAMMATION",                                           # KC6
  "CAUSES MITOCHONDRIAL DYSFUNCTION",                              # KC7
  "ACTIVATES STRESS SIGNALING PATHWAYS",                           # KC8
  "CAUSES CHOLESTASIS",                                            # KC9
  "DISRUPTS CELLULAR CYTOSKELETON",                                # KC10
  "CAUSES LIVER FIBROSIS",                                         # KC11
  "DISRUPTS LIPID AND PROTEIN METABOLISM OR CAUSES DYSLIPIDEMIA")  # KC12

# ---- 1. Set de pathways hepatotox depuis File008 (matrice pathway x KC) ------
# Format : 1 feuille par racine top-level Reactome ; header des 34 termes KC en
# LIGNE 2 ; noms de pathway en col 1 des la ligne 3 ; valeurs 0/1.
sheets <- setdiff(excel_sheets(F008), "TopPage")
parse_sheet <- function(s) {
  m    <- suppressMessages(read_excel(F008, sheet = s, col_names = FALSE))
  kc   <- as.character(unlist(m[2, ]))
  keep <- which(!is.na(kc) & kc != "" & kc != "group")
  body <- m[-(1:2), , drop = FALSE]
  pw   <- as.character(body[[1]])
  map_dfr(keep, function(j) tibble(
    pathway_name = pw, kc_term = kc[j],
    v = suppressWarnings(as.numeric(as.character(body[[j]]))))) %>%
    filter(!is.na(v), v == 1)
}
hep <- map_dfr(sheets, parse_sheet) %>% filter(kc_term %in% HEP_UMBRELLA)

nm <- read_tsv(NAMEF, col_names = c("id", "name", "sp"), show_col_types = FALSE) %>%
  filter(sp == "Homo sapiens") %>% distinct(name, .keep_all = TRUE)
hep <- hep %>% left_join(nm, by = c("pathway_name" = "name"))
n_miss <- hep %>% filter(is.na(id)) %>% distinct(pathway_name) %>% nrow()
hep_ids <- hep %>% filter(!is.na(id)) %>% distinct(id) %>% pull(id)
message(sprintf("Pathways hepatotox (Tsai File008, R-HSA resolus) : %d (%d noms non resolus, ignores)",
                length(hep_ids), n_miss))

# table de tracabilite (pathway_id, kc concatenes)
hep %>% filter(!is.na(id)) %>%
  group_by(pathway_id = id, pathway_name) %>%
  summarise(kc = paste(sort(unique(kc_term)), collapse = ";"), .groups = "drop") %>%
  arrange(kc, pathway_name) %>% write_tsv(TSV)

# ---- 2. Elagage : liver_1.8 = liver_1.7 prive des pathways non hepatotox -----
load(NET17)  # -> liver_1.7_network
g   <- liver_1.7_network
pw  <- which(V(g)$type == "pathway")
drop <- V(g)$name[pw[!(V(g)$pathway_id[pw] %in% hep_ids)]]
liver_1.8_network <- delete_vertices(g, drop)

message(sprintf("Noeuds %d -> %d | Aretes %d -> %d | Pathways %d -> %d",
                vcount(g), vcount(liver_1.8_network), ecount(g), ecount(liver_1.8_network),
                length(pw), sum(V(g)$type == "pathway") - length(drop)))

save(liver_1.8_network, file = OUT)
message("Ecrit : ", OUT, "\n        ", TSV)
