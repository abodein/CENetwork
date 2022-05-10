#' get_route
#'
#' @importFrom netOmics random_walk_restart
#' @importFrom netOmics rwr_find_closest_type
#' @importFrom igraph induced_subgraph shortest_paths
#' @importFrom purrr map2
#' @export

get_route <- function(network, signature, input_type = c("ENSEMBL"),
                      target_type = c("UNIPROT", "SMPDB", "DrugBank","CHEMBL"),
                      num_cores = 4){

    # check if num_cores ok

    # check if signature is not empty
    stopifnot(!purrr::is_empty(signature))


    #1 -- get RWR
    print("get RW scores")
    RWR_res <- netOmics::random_walk_restart(X = network, seed = signature)
    # can be computed
    ## /!\ Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102
    ## works with no PPI ego

    # SP_res <- shortest_paths(graph = network, from = signature, to = test_target)
    print("find closest")
    closest.dfr <- parallel::mclapply(target_type, function(x){
        netOmics::rwr_find_closest_type(X = RWR_res,
                                        seed = signature,
                                        attribute = "type",
                                        value = x) %>%
            map_dfr(~.x %>% as.data.frame)}, mc.cores = num_cores)

    closest <- map_dfr(closest.dfr, ~.x)

    # closest <- netOmics::rwr_find_closest_type(X = RWR_res, seed = signature, attribute = "type", value = target_type[[i]])
    # get SP
    input <- closest$SeedName
    target <- closest$NodeNames

    print("extract shortest path")
    all_sp <-  map2(input, target, ~{igraph::shortest_paths(graph=network, from=.x, to=.y) %>% .$vpath %>% lapply(names) %>% unlist()})
    # SP <- shortest_paths(graph = network, from = closest$SeedName[1], to = closest$NodeNames[1])
    # igraph::induced_subgraph(network, SP$vpath[[1]]) %>% plot
    #
    induced.g <- igraph::induced_subgraph(network, unlist(all_sp))
    vertex_attr(induced.g) <- vertex_attr(induced.g) %>% as.data.frame() %>%
        mutate(input_signature = name %in% signature)


    result <- list(network = induced.g, closest = closest)
    class(result) <- "get_route.res"
    return(result)

}
