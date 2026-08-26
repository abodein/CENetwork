#!/usr/bin/env Rscript
# =============================================================================
# data-raw/liver_1.8/03_rwr_1.8_full.R
# -----------------------------------------------------------------------------
# Precalcule la RWR pour liver_1.8_full (reseau de benchmark, ~49k noeuds) et
# ecrit data/liver_1.8_full_rwr_closest_dfr.rda.
#
# VARIANTE "seed-only" : au lieu de l'inverse dense n x n de pipeline_RWR
# (~49004^2*8o ~ 19 Go, infaisable), on ne calcule les scores RWR QUE depuis les
# noeuds seed pertinents. get_route ne diffuse et ne permute que sur des seeds de
# type "gene" (permutation : sample dans tous les noeuds gene du reseau), donc il
# suffit des lignes des noeuds GENE connectes.
#
# Maths (identiques a RWR_build_complete, mais lignes restreintes aux seeds S) :
#   W = A * diag(1/deg)           (colonne-stochastique)
#   M = I - (1-r) * t(W)
#   P[seed,target] = r * (M^-1)[seed,]              -> on veut seulement seed in S
#   (M^-1)[S,] = t( solve(t(M), E_S) )              (E_S = colonnes indicatrices de S)
# Cout : solve creux (LU) + |S| back-solves ; resultat dense |S| x n (~1.5 Go).
#
# Executer depuis la racine du package CENetwork.
# =============================================================================

suppressMessages({
  library(igraph); library(Matrix)
  pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE,
                    attach_testthat = FALSE)   # generate_closest_dfr, get_all_routes (internes)
})

PKG     <- "/Users/antoine/Documents/DEV/CENetwork"
RESTART <- 0.7
load(file.path(PKG, "data/liver_1.8_full_network.rda"))  # -> liver_1.8_full_network
g <- liver_1.8_full_network

# ---- transition matrix (comme RWR_build_complete) ----
Xi      <- CENetwork:::remove_unconnected_nodes(g)
seed_xi <- V(Xi)$name
n       <- length(seed_xi)
A       <- as_adjacency_matrix(Xi, sparse = TRUE)
cs      <- colSums(A); cs[cs == 0] <- 1
W       <- A %*% Diagonal(x = 1 / cs)
M       <- Diagonal(n) - (1 - RESTART) * t(W)

# ---- seeds = tous les noeuds gene connectes ----
S   <- intersect(V(Xi)$name[V(Xi)$type == "gene"], seed_xi)
idx <- match(S, seed_xi)
E   <- sparseMatrix(i = idx, j = seq_along(S), x = 1, dims = c(n, length(S)))
message(sprintf("liver_1.8_full : %d noeuds connectes | %d seeds gene", n, length(S)))

# ---- solve seed-only : (M^-1)[S,] = t(solve(t(M), E_S)) ----
message("solve creux (seed-only) ...")
t0 <- Sys.time()
Xsol     <- Matrix::solve(Matrix::t(M), E)          # n x |S|
res_seed <- RESTART * t(as.matrix(Xsol))            # |S| x n : row=seed, col=target
rownames(res_seed) <- S
colnames(res_seed) <- seed_xi
res_seed <- as.data.frame(res_seed, check.names = FALSE)  # garder les noms R-HSA/ENSG intacts
message(sprintf("  solve termine en %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- closest par type + routes (fonctions du package) ----
closest <- generate_closest_dfr(network = g, res_tmp_matrix = res_seed)
closest <- get_all_routes(network = g, closest_dfr = closest)
liver_1.8_full_rwr_closest_dfr <- closest

save(liver_1.8_full_rwr_closest_dfr,
     file = file.path(PKG, "data/liver_1.8_full_rwr_closest_dfr.rda"))
message(sprintf("Termine en %.1f min | closest_dfr : %d lignes",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(closest)))
message("Ecrit : data/liver_1.8_full_rwr_closest_dfr.rda")
