#!/usr/bin/env Rscript
# =============================================================================
# data-raw/liver_1.8/02_rwr_1.8.R
# -----------------------------------------------------------------------------
# Precalcule la diffusion RWR pour liver_1.8 et ecrit
# data/liver_1.8_rwr_closest_dfr.rda (meme structure que liver_1.7_rwr_closest_dfr).
#
# Utilise la fonction interne pipeline_RWR() du package (RWR_build_complete +
# generate_closest_dfr + get_all_routes), restart = 0.7. Le calcul materialise
# une matrice dense n x n (~29 800 noeuds connectes) -> plusieurs Go de RAM.
#
# Executer depuis la racine du package CENetwork.
# =============================================================================

suppressMessages({
  library(igraph)
  pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE,
                    attach_testthat = FALSE)   # expose pipeline_RWR (interne)
})

PKG <- "/Users/antoine/Documents/DEV/CENetwork"
load(file.path(PKG, "data/liver_1.8_network.rda"))  # -> liver_1.8_network

message("pipeline_RWR (restart=0.7) sur liver_1.8 ...")
t0  <- Sys.time()
res <- pipeline_RWR(liver_1.8_network)              # restart = 0.7 par defaut
liver_1.8_rwr_closest_dfr <- res$closest_dfr
message(sprintf("Termine en %.1f min | closest_dfr : %d lignes",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                nrow(liver_1.8_rwr_closest_dfr)))

save(liver_1.8_rwr_closest_dfr,
     file = file.path(PKG, "data/liver_1.8_rwr_closest_dfr.rda"))
message("Ecrit : data/liver_1.8_rwr_closest_dfr.rda")
