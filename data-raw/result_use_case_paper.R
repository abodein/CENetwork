# Generate result_use_case_paper_1.8.Rda  (network liver_1.8, pruned pathway layer)
#
# Runs the full diffusion + permutation pipeline for two hepatotoxic drugs:
#   - Acetaminophen (APAP) anchored on DB00316
#   - Valproic acid  (VPA)  anchored on DB00313
#
# Output (saved to <project>/papier/result_use_case_paper_1.8.Rda):
#   res_apap, res_report_aceta, enrich_apap
#   res_vpa,  res_report_valpro, enrich_vpa
#
# This file replaces the original produce_report_use_cases.R and adds:
#   - mandatory_nodes as a dedicated argument (no longer mixed into signature_vids)
#   - do_permutation = TRUE with n_permut = 1000
#   - compute_diffusion_enrichment() for GO and pathway enrichment stats
#
# Run once from the package root:
#   source("data-raw/result_use_case_paper.R")
#
# Runtime: ~30–60 min depending on hardware (1000 permutations × 2 drugs).
# Use n_permut = 100 for a quick test.

devtools::load_all()

data(liver_1.8_network)
data(liver_1.8_rwr_closest_dfr)
data(signature_maison)

target_layers <- c("drug/compound", "pathway", "side_effect", "GO")
n_permut      <- 100
seed_permut   <- 123

path_out <- file.path("~", "Documents", "Projects", "LO_cosmetic_europe",
                      "papier", "result_use_case_paper_1.8.Rda")

# --- Valproic acid -----------------------------------------------------------

message("=== Valproic acid ===")
sig_vpa <- signature_maison$`valproic acid_all_all`

res_vpa <- get_route(
  network         = liver_1.8_network,
  closest_dfr     = liver_1.8_rwr_closest_dfr,
  signature_vids  = sig_vpa,
  target_type     = target_layers,
  mandatory_nodes = "DB00313",
  do_permutation  = TRUE,
  n_permut        = n_permut,
  seed_permut     = seed_permut
)

enrich_vpa <- compute_diffusion_enrichment(
  route_result    = res_vpa,
  target_types    = c("GO", "pathway"),
  min_visit_count = 1L,
  padj_method     = "BH"
)

res_report_valpro <- report(res_vpa, complete_network = liver_1.8_network)

message(sprintf("VPA: %d nodes, %d edges",
                igraph::vcount(res_vpa$network),
                igraph::ecount(res_vpa$network)))
message(sprintf("VPA enrichment: %d GO / %d pathway (padj < 0.05)",
                sum(enrich_vpa$type == "GO"      & enrich_vpa$padj < 0.05, na.rm = TRUE),
                sum(enrich_vpa$type == "pathway" & enrich_vpa$padj < 0.05, na.rm = TRUE)))

# --- Acetaminophen -----------------------------------------------------------

message("=== Acetaminophen ===")
sig_apap <- signature_maison$acetaminophen_all_all

res_apap <- get_route(
  network         = liver_1.8_network,
  closest_dfr     = liver_1.8_rwr_closest_dfr,
  signature_vids  = sig_apap,
  target_type     = target_layers,
  mandatory_nodes = "DB00316",
  do_permutation  = TRUE,
  n_permut        = n_permut,
  seed_permut     = seed_permut
)

enrich_apap <- compute_diffusion_enrichment(
  route_result    = res_apap,
  target_types    = c("GO", "pathway"),
  min_visit_count = 1L,
  padj_method     = "BH"
)

res_report_aceta <- report(res_apap, complete_network = liver_1.8_network)

message(sprintf("APAP: %d nodes, %d edges",
                igraph::vcount(res_apap$network),
                igraph::ecount(res_apap$network)))
message(sprintf("APAP enrichment: %d GO / %d pathway (padj < 0.05)",
                sum(enrich_apap$type == "GO"      & enrich_apap$padj < 0.05, na.rm = TRUE),
                sum(enrich_apap$type == "pathway" & enrich_apap$padj < 0.05, na.rm = TRUE)))

# --- Save --------------------------------------------------------------------

save(
  sig_apap, res_apap, enrich_apap, res_report_aceta,
  sig_vpa,  res_vpa,  enrich_vpa,  res_report_valpro,
  file = path_out
)

message("Saved: ", path_out)
res_apap$network  %>% vertex_attr()  %>% as.data.frame()  %>% group_by(type)  %>% summarise(apap = n()) %>% 
  left_join(res_vpa$network  %>% vertex_attr()  %>% as.data.frame()  %>% group_by(type)  %>% summarise(vpa = n()), by = "type")

#> enrich_apap  %>% left_join(vertex_attr(liver_1.8_network)  %>% as.data.frame(), by = c("term" = "name"))  %>% dplyr::select(type.x, term, display_name, visit_count, mean_perm, sd_perm, z_score, pvalue, padj)  %>% arrange(type.x, padj)

# The xlsx export comparing diffusion / ORA / GSEA now lives at the end of
# data-raw/comparison_ora_gsea.R (it needs the ORA/GSEA results computed there).
