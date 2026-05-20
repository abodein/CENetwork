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
#' @param mandatory_nodes (character) optional list of node names that must be included in the route as waypoints.
#'   Shortest paths from seeds and targets to each mandatory node are computed and added to the induced subnetwork.
#' @param do_permutation (logical) whether to run permutation analysis to assess significance
#' @param n_permut (integer) number of permutations (default 100)
#' @param seed_permut (integer) random seed for reproducibility (default 123)
#'
#' @return
#' Return a list with items:
#' \item{input_diffusion}{signature_vids inside the network}
#' \item{network}{the resulting subnetwork (igraph). Input seeds are tagged via the input_diffusion attribute. input_gene_signature and input_protein_signature attributes are set for corresponding nodes.}
#' \item{closest}{closest_dfr filtered based on seed and target types}
#' \item{permutation_results}{(only when do_permutation=TRUE) list of per-permutation summary stats: n_seeds_found, n_targets, n_nodes, visit_counts (named numeric vector of visit counts for each observed node)}
#' \item{observed_visit_counts}{(only when do_permutation=TRUE) named numeric vector of visit counts for the observed network}
#' \item{z_scores}{(only when do_permutation=TRUE) named numeric vector of z-scores per node}
#'
#' @examples
#' data(liver_1.7_network)
#' data(liver_1.7_rwr_closest_dfr)
#' data(signature_maison)
#'
#' signature_vids <- signature_maison$acetaminophen_all_all
#'
#' res.diffusion <- get_route(liver_1.7_network, liver_1.7_rwr_closest_dfr, signature_vids)
#'
#' # With mandatory nodes (e.g. force a specific drug node into the route)
#' res.diffusion <- get_route(liver_1.7_network, liver_1.7_rwr_closest_dfr, signature_vids,
#'                            mandatory_nodes = "DB00316")
#'
#' # With permutation analysis (assess significance of the observed route)
#' res.diffusion <- get_route(liver_1.7_network, liver_1.7_rwr_closest_dfr, signature_vids,
#'                            do_permutation = TRUE, n_permut = 10, seed_permut = 123)
#'
#' # Inspect permutation results
#' permut_stats <- do.call(rbind, lapply(res.diffusion$permutation_results, function(x) {
#'   data.frame(n_seeds_found = x$n_seeds_found,
#'              n_targets     = x$n_targets,
#'              n_nodes       = x$n_nodes)
#' }))
#'
#' # Compare observed vs null distribution (e.g. number of nodes)
#' observed_n_nodes <- igraph::vcount(res.diffusion$network)
#' p_value_global <- mean(permut_stats$n_nodes >= observed_n_nodes)
#'
#' # Z-score per node: (observed_visit_count - mean_permut) / sd_permut
#' # visit_count = number of seed->target shortest paths passing through the node
#' # high z-score => node is visited much more than expected by chance
#' head(sort(res.diffusion$z_scores, decreasing = TRUE))
#'
#' # Both visit_count and z_score are also stored as vertex attributes
#' igraph::vertex_attr(res.diffusion$network, "visit_count")
#' igraph::vertex_attr(res.diffusion$network, "z_score")
#'
#' @importFrom netOmics random_walk_restart
#' @importFrom netOmics rwr_find_closest_type
#' @importFrom igraph induced_subgraph shortest_paths vertex_attr as_data_frame
#' @importFrom purrr map2 keep is_empty
#' @importFrom AnnotationDbi select
#' @importFrom dplyr select bind_rows filter mutate pull
#' @import org.Hs.eg.db
#' @importFrom magrittr %>%
#'
#' @export
get_route <- function(network, closest_dfr, signature_vids,
                      target_type = c("pathway", "drug/compound"),
                      mandatory_nodes = NULL,
                      do_permutation = FALSE, n_permut = 100, seed_permut = 123) {

    # Check arguments ----------------------------------------------------------

    stopifnot(is(network, "igraph"))

    stopifnot(is(closest_dfr, "data.frame"))
    stopifnot(colnames(closest_dfr) %in% c("SeedNode", "type.seed", "NodeNames", "type.target", "distance_from_start", "shortest_path"))

    stopifnot(is(signature_vids, "character"))
    stopifnot(!purrr::is_empty(signature_vids))

    layers <- unique(c(closest_dfr$type.seed, closest_dfr$type.target))
    stopifnot(is(target_type, "character"))
    stopifnot(all(target_type %in% layers))

    if (!is.null(mandatory_nodes)) {
        stopifnot(is(mandatory_nodes, "character"))
        missing_nodes <- mandatory_nodes[!mandatory_nodes %in% closest_dfr$SeedNode]
        if (length(missing_nodes) > 0) {
            stop(paste("The following mandatory_nodes are not found as SeedNode in closest_dfr:",
                       paste(missing_nodes, collapse = ", ")))
        }
    }

    # Process -----------------------------------------------------------------

    if (!"shortest_path" %in% colnames(closest_dfr)) {
        closest_dfr <- get_all_routes(network, closest_dfr)
    }

    effective_seeds <- unique(c(signature_vids, mandatory_nodes))

    closest_res <- closest_dfr %>%
        dplyr::filter(SeedNode %in% effective_seeds) %>%
        dplyr::filter(type.target %in% target_type)

    induced_nodes <- closest_res$shortest_path %>%
        purrr::keep(~!purrr::is_empty(.x)) %>%
        unlist() %>%
        unique()

    induced.g <- igraph::induced_subgraph(network, induced_nodes)
    igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>%
        as.data.frame() %>%
        dplyr::mutate(input_diffusion      = name %in% signature_vids) %>%
        dplyr::mutate(input_gene_signature = (name %in% signature_vids) & type == "gene") %>%
        as.list()

    # Gene → protein conversion
    input_gene_signature <- igraph::vertex_attr(induced.g) %>%
        as.data.frame() %>%
        dplyr::filter(input_gene_signature) %>%
        dplyr::pull(name)

    if (length(input_gene_signature)) {
        lut <- AnnotationDbi::select(x        = org.Hs.eg.db::org.Hs.eg.db,
                                     keys     = input_gene_signature,
                                     keytype  = "ENSEMBL",
                                     columns  = "UNIPROT")
        igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>%
            as.data.frame() %>%
            dplyr::mutate(input_protein_signature = name %in% lut$UNIPROT) %>%
            as.list()
    } else {
        igraph::vertex_attr(induced.g) <- igraph::vertex_attr(induced.g) %>%
            as.data.frame() %>%
            dplyr::mutate(input_protein_signature = NA) %>%
            as.list()
    }

    result <- list(
        network   = induced.g,
        closest   = closest_res,
        input     = unique(closest_res$SeedNode),
        autoroute = closest_res
    )
    class(result) <- "get_route.res"

    # Permutations ------------------------------------------------------------

    if (do_permutation) {
        set.seed(seed_permut)

        all_genes <- igraph::as_data_frame(network, what = "vertices") %>%
            dplyr::filter(type == "gene") %>%
            dplyr::pull(name)

        # Only permute the gene part of the signature; mandatory_nodes are always fixed
        gene_seeds <- signature_vids[signature_vids %in% all_genes]
        n_gene_seeds <- max(length(gene_seeds), 1L)
        pool <- setdiff(all_genes, mandatory_nodes)

        signature_permut <- lapply(seq_len(n_permut), function(i) {
            unique(c(mandatory_nodes, sample(pool, n_gene_seeds)))
        })

        observed_nodes <- igraph::V(result$network)$name

        # Observed visit count: number of seed->target paths that pass through each node
        observed_paths <- closest_res$shortest_path %>% purrr::keep(~!purrr::is_empty(.x))
        observed_visit_counts <- vapply(observed_nodes, function(node) {
            sum(vapply(observed_paths, function(p) node %in% p, logical(1L)))
        }, numeric(1L))

        # Permutation visit counts: same metric on randomised signatures
        permutation_results <- lapply(signature_permut, function(sig) {
            tryCatch({
                cr <- closest_dfr %>%
                    dplyr::filter(SeedNode %in% sig) %>%
                    dplyr::filter(type.target %in% target_type)

                paths <- cr$shortest_path %>% purrr::keep(~!purrr::is_empty(.x))

                visit_counts <- vapply(observed_nodes, function(node) {
                    sum(vapply(paths, function(p) node %in% p, logical(1L)))
                }, numeric(1L))

                list(
                    n_seeds_found = length(unique(cr$SeedNode)),
                    n_targets     = length(unique(cr$NodeNames)),
                    n_nodes       = sum(visit_counts > 0L),
                    visit_counts  = visit_counts
                )
            }, error = function(e) NULL)
        })

        # Z-score per node: (observed - mean_permut) / sd_permut
        permut_matrix <- do.call(rbind, lapply(permutation_results, function(x) {
            if (is.null(x)) rep(NA_real_, length(observed_nodes)) else x$visit_counts
        }))

        permut_means <- colMeans(permut_matrix, na.rm = TRUE)
        permut_sds   <- apply(permut_matrix, 2L, sd, na.rm = TRUE)

        z_scores <- (observed_visit_counts - permut_means) / permut_sds
        z_scores[is.nan(z_scores) | !is.finite(z_scores)] <- NA_real_

        igraph::vertex_attr(result$network, "visit_count") <- observed_visit_counts
        igraph::vertex_attr(result$network, "z_score")     <- z_scores

        result$permutation_results    <- permutation_results
        result$observed_visit_counts  <- observed_visit_counts
        result$z_scores               <- z_scores
    }

    return(result)
}


# get all routes
get_all_routes <- function(network, closest_dfr) {

    input  <- closest_dfr$SeedNode
    target <- closest_dfr$NodeNames

    all_sp <- list()
    for (i in seq_along(input)) {
        print(i)
        all_sp[[i]] <- try(
            igraph::shortest_paths(graph = network, from = input[i], to = target[i]) %>%
                .$vpath %>% lapply(names) %>% unlist(),
            silent = TRUE
        )
    }

    new_closest_dfr <- closest_dfr %>% as_tibble() %>% mutate(shortest_path = all_sp)
    return(new_closest_dfr)
}
