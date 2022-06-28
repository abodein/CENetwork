#` RWR scores
#
# Generate all random walk scores. Iteratively, each node become the seed.
# Return the RWR for the connected nodes in the network as a squared data.frame.

#' @param X (igraph) network
#' @param verbose (logical) default = TRUE
#'
#' @return a data.frame

#' @import progress
#' @import magrittr
#' @importFrom igraph V
#' @importFrom RandomWalkRestartMH create.multiplex
#' @importFrom RandomWalkRestartMH compute.adjacency.matrix
#' @importFrom RandomWalkRestartMH normalize.multiplex.adjacency
#' @importFrom RandomWalkRestartMH Random.Walk.Restart.Multiplex
#' @importFrom purrr set_names imap_dfr
#' @importFrom tibble column_to_rownames
#' @importFrom netOmics remove_unconnected_nodes


RWR_build_complete <- function (X, verbose = TRUE)
{
    # X <- check_graph(X)

    # comp <- igraph::components(X)

    Xi <- netOmics:::remove_unconnected_nodes(X)
    seed_xi <- igraph::V(Xi)$name


    multiplex <- RandomWalkRestartMH::create.multiplex(LayersList = list(L1 = Xi),
                                                       Layers_Name = "layers_name")
    adj_matrix <- RandomWalkRestartMH::compute.adjacency.matrix(x = multiplex)
    adj_matrix_norm <- RandomWalkRestartMH::normalize.multiplex.adjacency(x = adj_matrix)


    res_tmp <- list()
    p <- progress::progress_bar$new(total = length(seed_xi))
    for (seed_xi_i in seed_xi) {

        rwr_res <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(x = adj_matrix_norm,
                                                                      MultiplexObject = multiplex, Seeds = seed_xi_i)
        res_tmp[[seed_xi_i]] <- rwr_res
        p$tick()
    }


    res_tmp_matrix <- purrr::imap_dfr(res_tmp, ~{.x$RWRM_Results %>%
            tibble::column_to_rownames("NodeNames") %>%
            purrr::set_names(.y) %>% t %>% as.data.frame})

    order_row_col <- sort(colnames(res_tmp_matrix))
    res_tmp_matrix <- res_tmp_matrix[order_row_col,order_row_col]
    return(res_tmp_matrix)  # target as column, seed as row to apply generate_closest_dfr (rowwise)
}


#' Closest dfr
#'
#' From RWR result matrix generate closest dfr based on "type" attribue
#'
#' @param res_tmp_matrix
#'
#' @import progress
#' @import magrittr
#' @importFrom igraph vertex_attr
#' @importFrom purrr set_names imap_dfr
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select
#'
generate_closest_dfr <- function(network, res_tmp_matrix){

    vecteur_group <- colnames(res_tmp_matrix)

    tmp_group <- igraph::vertex_attr(network) %>% # full network
        as.data.frame() %>% dplyr::select(name, type) %>% # or other column
        tibble::column_to_rownames("name") %>% .[colnames(res_tmp_matrix),,drop = FALSE]  %>% pull(type)

    groups <- unique(tmp_group)

    closest_target <- list()
    p <- progress::progress_bar$new(total =nrow(res_tmp_matrix))
    for(i in 1:nrow(res_tmp_matrix)){
        closest_target[[i]] <- list()
        for(j in groups) {
            #closest_target[[i]][[j]] <- names(which.max(res_tmp_matrix[i, tmp_group == j]))
            tmp <- res_tmp_matrix[i, tmp_group == j]
            closest_target[[i]][[j]] <- names(tmp)[which(tmp == max(tmp))]
        }
        p$tick()
    }
    names(closest_target) <- rownames(res_tmp_matrix)
    closest.dfr <- purrr::imap_dfr(closest_target, ~lapply(groups, function(gr){
        data.frame("NodeNames" = .x[[gr]],
                   "type" = rep(gr, times = length(.x[[gr]])),
                   "SeedNode" = rep(.y,times = length(.x[[gr]])))
    }))
    return(closest.dfr)
}


#' Pipeline RWR
#'
#' Run RWR_build_complete and generate_closest_dfr
#' @param X network
#'
pipeline_RWR <- function(X){
    res_tmp_matrix <- RWR_build_complete(X)
    closest_dfr <- generate_closest_dfr(network = X, res_tmp_matrix = res_tmp_matrix)
    return(list("rwr_matrix" = res_tmp_matrix,
                "closest_dfr" = closest_dfr))
}
