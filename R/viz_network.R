#' Interactive visualisation of a diffusion subnetwork
#'
#' Renders the subnetwork returned by \code{get_route()} as an interactive
#' \pkg{visNetwork} widget: nodes are coloured by layer type, labelled with
#' their human-readable name, and seed nodes are highlighted. The widget is
#' self-contained (drag, zoom, hover for details, click to highlight a node's
#' neighbourhood) and can be saved to a standalone HTML file with
#' \code{htmlwidgets::saveWidget()}.
#'
#' This is a lightweight alternative to \code{export_to_cytoscape()} for
#' users who do not have Cytoscape installed or want a shareable, static-file
#' visualisation.
#'
#' @param route_result Output of \code{get_route()}.
#' @param seed_color Colour used for seed nodes (nodes flagged
#'   \code{input_diffusion == TRUE}). Default \code{"firebrick"}.
#' @param height,width Widget dimensions, passed to
#'   \code{visNetwork::visNetwork()}. Default \code{"800px"} and \code{"100\%"}.
#'
#' @return A \code{visNetwork} htmlwidget.
#'
#' @examples
#' \dontrun{
#' data(liver_1.8_network)
#' data(liver_1.8_rwr_closest_dfr)
#' data(signature_maison)
#'
#' res_apap <- get_route(
#'   network        = liver_1.8_network,
#'   closest_dfr    = liver_1.8_rwr_closest_dfr,
#'   signature_vids = signature_maison$acetaminophen_all_all,
#'   target_type    = c("drug/compound", "pathway", "GO", "side_effect")
#' )
#' viz_network(res_apap)
#' }
#'
#' @export
viz_network <- function(route_result, seed_color = "firebrick",
                         height = "800px", width = "100%") {

    if (!requireNamespace("visNetwork", quietly = TRUE)) {
        stop("Package 'visNetwork' is required for viz_network(). ",
             "Install it with install.packages('visNetwork').")
    }
    stopifnot(is(route_result, "get_route.res"))

    g <- route_result$network

    type_colors <- c(
        "gene"          = "#8dd3c7",
        "protein"       = "#80b1d3",
        "GO"            = "#fdb462",
        "pathway"       = "#b3de69",
        "drug/compound" = "#fccde5",
        "side_effect"   = "#bebada"
    )

    vdf <- igraph::as_data_frame(g, what = "vertices")
    edf <- igraph::as_data_frame(g, what = "edges")

    display <- ifelse(!is.na(vdf$display_name) & nzchar(vdf$display_name),
                       vdf$display_name, vdf$name)

    is_seed <- if (!is.null(vdf$input_diffusion)) {
        vapply(vdf$input_diffusion, isTRUE, logical(1L))
    } else {
        rep(FALSE, nrow(vdf))
    }

    node_color <- unname(type_colors[vdf$type])
    node_color[is.na(node_color)] <- "#d9d9d9"
    node_color[is_seed] <- seed_color

    nodes <- data.frame(
        id    = vdf$name,
        label = display,
        group = vdf$type,
        color = node_color,
        title = paste0("<b>", display, "</b><br>type: ", vdf$type,
                        ifelse(is_seed, "<br><i>seed</i>", "")),
        stringsAsFactors = FALSE
    )

    edges <- data.frame(from = edf$from, to = edf$to, stringsAsFactors = FALSE)

    widget <- visNetwork::visNetwork(nodes, edges, height = height, width = width) %>%
        visNetwork::visNodes(font = list(size = 14)) %>%
        visNetwork::visEdges(color = list(color = "#cccccc", opacity = 0.5)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                                nodesIdSelection = TRUE) %>%
        visNetwork::visLegend(useGroups = TRUE) %>%
        visNetwork::visPhysics(stabilization = FALSE) %>%
        visNetwork::visInteraction(navigationButtons = TRUE)

    widget
}
