#' RWR scores
#'
#' Generate all random walk scores for all seeds simultaneously using the
#' analytical solution: P = r * solve(I - (1-r) * W^T), where W is the
#' column-stochastic transition matrix. Uses sparse LU factorisation (UMFPACK)
#' via the Matrix package so the factorisation step stays sparse even for
#' large networks (~6k nodes).
#' Return the RWR for the connected nodes in the network as a squared data.frame.

#' @param X (igraph) network
#' @param restart (numeric) restart probability, default = 0.7
#' @param verbose (logical) default = TRUE
#'
#' @return a data.frame (seed in rows, target in columns)

#' @importFrom igraph V as_adjacency_matrix
#' @importFrom Matrix Diagonal colSums t solve
#' @import netOmics
#'
#' @examples
#' library(igraph)
#'
#' set.seed(42)
#' n <- 200
#' g <- igraph::sample_gnm(n, n * 3, directed = FALSE)
#' igraph::V(g)$name <- paste0("node", seq_len(n))
#' igraph::V(g)$type <- rep(c("A", "B", "C", "D"), length.out = n)
#'
#' rwr_matrix <- RWR_build_complete(g)
#' # rwr_matrix is a 50x50 data.frame: rwr_matrix[seed, target]
#' dim(rwr_matrix)


RWR_build_complete <- function(X, restart = 0.7, verbose = TRUE)
{
    Xi <- netOmics:::remove_unconnected_nodes(X)
    seed_xi <- igraph::V(Xi)$name
    n <- length(seed_xi)

    # Sparse column-stochastic transition matrix: W[i,j] = A[i,j] / degree(j)
    A <- igraph::as_adjacency_matrix(Xi, sparse = TRUE)
    col_sums <- Matrix::colSums(A)
    col_sums[col_sums == 0] <- 1
    W <- A %*% Matrix::Diagonal(x = 1 / col_sums)

    # Analytical RWR for all seeds simultaneously:
    # p*_s = r * (I - (1-r)*W)^{-1} * e_s
    # Full matrix (row = seed, col = target): P = r * (I - (1-r)*W^T)^{-1}
    # M is sparse — Matrix::solve() uses sparse LU (UMFPACK), factorises once,
    # then does n back-solves. Result is dense (n x n) but factorisation is not.
    M <- Matrix::Diagonal(n) - (1 - restart) * Matrix::t(W)
    res_matrix <- restart * as.matrix(Matrix::solve(M))
    dimnames(res_matrix) <- list(seed_xi, seed_xi)

    order_row_col <- sort(colnames(res_matrix))
    res_matrix <- res_matrix[order_row_col, order_row_col]
    return(as.data.frame(res_matrix))  # target as column, seed as row
}


#' Closest dfr
#'
#' From RWR result matrix generate closest dfr based on "type" attribue
#'
#' @param network igraph
#' @param res_tmp_matrix squared data.frame with seed in row and targets in columns
#'
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom magrittr %>%
#' @importFrom igraph vertex_attr
#' @importFrom purrr set_names imap_dfr
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select
#'
#' @examples
#' library(igraph)
#'
#' set.seed(42)
#' n <- 200
#' g <- igraph::sample_gnm(n, n * 3, directed = FALSE)
#' igraph::V(g)$name <- paste0("node", seq_len(n))
#' igraph::V(g)$type <- rep(c("A", "B", "C", "D"), length.out = n)
#'
#' rwr_matrix <- RWR_build_complete(g)
#' closest <- generate_closest_dfr(network = g, res_tmp_matrix = rwr_matrix)
#' # closest is a data.frame: SeedNode, type.seed, NodeNames, type.target
#' head(closest)
#'
generate_closest_dfr <- function(network, res_tmp_matrix){

    vecteur_group <- colnames(res_tmp_matrix)

    tmp_group <- igraph::vertex_attr(network) %>% # full network
        as.data.frame() %>% dplyr::select(name, type) %>% # or other column
        tibble::column_to_rownames("name") %>% .[colnames(res_tmp_matrix),,drop = FALSE]  %>% pull(type)

    groups <- unique(tmp_group)

    closest_target <- list()
    cli::cli_progress_bar("generate_closest_dfr", total = nrow(res_tmp_matrix))
    for(i in 1:nrow(res_tmp_matrix)){
        closest_target[[i]] <- list()
        for(j in groups) {
            tmp <- res_tmp_matrix[i, tmp_group == j]
            closest_target[[i]][[j]] <- names(tmp)[which(tmp == max(tmp))]
        }
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    names(closest_target) <- rownames(res_tmp_matrix)
    closest.dfr <- purrr::imap_dfr(closest_target, ~lapply(groups, function(gr){
        data.frame("NodeNames" = .x[[gr]],
                   "type" = rep(gr, times = length(.x[[gr]])),
                   "SeedNode" = rep(.y,times = length(.x[[gr]])))
    }))

    va <- vertex_attr(network) %>% as.data.frame()

    closest.dfr <- closest.dfr %>% dplyr::select(NodeNames, SeedNode) %>%
        left_join(va %>% rename_with(~paste0(.x, ".target")),by = c("NodeNames" = "name.target")) %>%
        left_join(va %>% rename_with(~paste0(.x, ".seed")),by = c("SeedNode" = "name.seed")) %>%
        dplyr::select(SeedNode, type.seed, NodeNames, type.target)

    return(closest.dfr)
}


#' Get all routes
#'
#' For each (SeedNode, NodeNames) pair in closest_dfr, compute the shortest
#' path on the network. Grouped by unique seed so BFS runs once per seed,
#' not once per pair.
#'
#' @param network (igraph) network
#' @param closest_dfr (data.frame) output of generate_closest_dfr
#'
#' @return a tibble with an additional list column \code{shortest_path},
#'   each element a character vector of node names along the path.
#'
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom igraph shortest_paths
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#'
#' @examples
#' library(igraph)
#'
#' set.seed(42)
#' n <- 1000
#' g <- igraph::sample_gnm(n, n * 3, directed = FALSE)
#' igraph::V(g)$name <- paste0("node", seq_len(n))
#' igraph::V(g)$type <- rep(c("A", "B", "C", "D"), length.out = n)
#'
#' rwr_matrix <- RWR_build_complete(g)
#' closest    <- generate_closest_dfr(network = g, res_tmp_matrix = rwr_matrix)
#' routes     <- get_all_routes(network = g, closest_dfr = closest)
#' # routes is a tibble with a shortest_path list column
#' routes$shortest_path[[1]]
#'
get_all_routes <- function(network, closest_dfr) {
    unique_seeds <- unique(closest_dfr$SeedNode)
    result_sp <- vector("list", nrow(closest_dfr))

    cli::cli_progress_bar("get_all_routes", total = length(unique_seeds))
    for (s in unique_seeds) {
        idx     <- which(closest_dfr$SeedNode == s)
        targets <- closest_dfr$NodeNames[idx]

        paths <- try(
            igraph::shortest_paths(network, from = s, to = targets)$vpath,
            silent = TRUE
        )

        if (inherits(paths, "try-error")) {
            for (i in idx) result_sp[[i]] <- character(0)
        } else {
            for (j in seq_along(idx)) result_sp[[idx[j]]] <- names(paths[[j]])
        }
        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    tibble::as_tibble(closest_dfr) %>%
        dplyr::mutate(shortest_path = result_sp)
}


#' Pipeline RWR
#'
#' Run RWR_build_complete, generate_closest_dfr and get_all_routes
#'
#' @param X network
#'
pipeline_RWR <- function(X){
    res_tmp_matrix <- RWR_build_complete(X)
    closest_dfr    <- generate_closest_dfr(network = X, res_tmp_matrix = res_tmp_matrix)
    closest_dfr    <- get_all_routes(network = X, closest_dfr = closest_dfr)
    return(list("rwr_matrix"  = res_tmp_matrix,
                "closest_dfr" = closest_dfr))
}
