#' get_route
#'
#' @importFrom netOmics random_walk_restart
#' @importFrom netOmics rwr_find_closest_type
#' @importFrom igraph induced_subgraph
#' @importFrom purrr map2
#' @export
get_route <- function(network, signature, input_type = c("ENSEMBL"),
                      target_type = c("UNIPROT", "SMPDB", "DrugBank","CHEMBL")){

    # check if signature is not empty
    stopifnot(is_empty(signature))

    #1 -- get RWR
    RWR_res <- netOmics::random_walk_restart(X = network, seed = signature)
    # can be computed
    ## /!\ Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102
    ## works with no PPI ego

    closest <- lapply(target_type,
                      function(x){netOmics::rwr_find_closest_type(X = RWR_res,
                                                                  seed = signature,
                                                                  attribute = "type",
                                                                  value = x) %>%
            map_dfr(~.x %>% as.data.frame)}) %>%
        map_dfr(~.x)

    # get SP
    input <- closest$SeedName
    target <- closest$NodeNames

    all_sp <-  purrr::map2(input, target, ~{shortest_paths(graph=network, from=.x, to=.y) %>% .$vpath %>% lapply(names) %>% unlist()})
    sub.res <- igraph::induced_subgraph(network, unlist(all_sp))

    return(sub.res)
}

