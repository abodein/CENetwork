#' Diffusion
#'
#' This function spreads a signal over the network from an entrance point (seeds) on the network via the random walk strategy.
#' A particle moves randomly from neighbour to neighbour based on the available connections.
#' In this implementation, it is possible to force the passage through the different layers of the network.
#'
#' For each seed, the nearest neighbours per layer are selected.
#'
#' Finally, a sub-network is returned by shortest paths including the nearest neighbours per seed.
#'
#' @param network (igraph) the network for diffusion
#' @param closest_dfr (data.frame)  precalculated closest neighbors (the data.frame must contains the columns "SeedNode","type.seed","NodeNames","type.target", see the vignette)
#' @param signature_vids (character) list of seeds, names of vertices in the network
#' @param target_type (character) list of layers to force the diffusion; closest_dfr is filtered based on target types.
#'
#' @return
#' Return a list with 2 items:
#' \item{input_diffusion}{signature_vids inside the network}
#' \item{network}{the resulting subnetork (igraph). In addition, input seeds are tagged in the network (input_diffusion attribute). Also input_gene_signature and input_protein_signature attributes are set to TRUE if corresponding nodes are present in the subnetowrk.}
#' \item{closest}{closest_dfr filtered based on seed and target types}

#'
#'
#' @examples
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#'
#' signature_vids <- signature_maison$acetaminophen_all_all
#'
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, signature_vids)
#'
#' @importFrom netOmics random_walk_restart
#' @importFrom netOmics rwr_find_closest_type
#' @importFrom igraph induced_subgraph shortest_paths vertex_attr
#' @importFrom purrr map2
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @importFrom magrittr %>%
#'
#' @export
get_route <- function(network, closest_dfr, signature_vids,
                      target_type =c("pathway", "drug/compound")){


    # Check argument ----------------------------------------------------------
    #
    ## network
    stopifnot(is(network, "igraph"))

    ## closest dfr
    stopifnot(is(closest_dfr, "data.frame"))
    stopifnot(colnames(closest_dfr) %in% c("SeedNode", "type.seed", "NodeNames", "type.target"))
    # SeedNode and NodeNa,es are present in network # too long to test
    # stopifnot(c(closest_dfr$SeedNode, closest_dfr$NodeNames) %in% V(network)$name)

    # check if signature is not empty
    stopifnot(is(signature_vids, "character"))
    stopifnot(!purrr::is_empty(signature_vids))

    # target_type
    layers <-  unique(c(closest_dfr$type.seed, closest_dfr$type.target))
    stopifnot(is(target_type, "character"))
    stopifnot(all(target_type %in% layers))


    # Process -----------------------------------------------------------------
    #
    # filter RWR results to get only signature
    closest_res <- closest_dfr %>% filter(SeedNode %in% signature_vids) %>%
        filter(type.target %in% target_type)

    input <- closest_res$SeedNode
    target <- closest_res$NodeNames

    # extract shortest path
    all_sp <-  purrr::map2(input, target, ~{igraph::shortest_paths(graph=network, from=.x, to=.y) %>% .$vpath %>% lapply(names) %>% unlist()})

    # update autoroute
    all_sp_df_distance <- purrr::map2(input, all_sp, ~data.frame(start = .x, path = .y, distance_from_start = seq_along(.y))) %>%
        bind_rows

    induced.g <- igraph::induced_subgraph(network, unlist(all_sp))
    igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>% as.data.frame() %>%
        dplyr::mutate(input_diffusion = name %in% signature_vids) %>%
        dplyr::mutate(input_gene_signature = (name %in% signature_vids) & type == "gene") %>%
        as.list()

    # get gene input signature and convert to protein
    input_gene_signature <- igraph::vertex_attr(induced.g) %>% as.data.frame() %>%
        dplyr::filter(input_gene_signature) %>%
        dplyr::pull(name)
    if(length(input_gene_signature)){
        lut <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db,
                                     keys = input_gene_signature,
                                     keytype = "ENSEMBL", columns = "UNIPROT")
        igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>% as.data.frame() %>%
            dplyr::mutate(input_protein_signature = name %in% lut$UNIPROT) %>% as.list
    } else {
        igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>% as.data.frame() %>%
            dplyr::mutate(input_protein_signature = NA) %>% as.list
    }

    # return results
    result <- list(network = induced.g, closest = closest_res, input = input,
                   autoroute = all_sp_df_distance)
    class(result) <- "get_route.res"
    return(result)
}
