#' Combine network layers
#'
#' Merge two layers of a multi-layer network, or connect a layer to new nodes
#' through a table of interactions. Nodes are matched by name; shared vertex
#' attributes are coalesced, the value of \code{graph1} taking precedence.
#'
#' Reimplementation of \code{netOmics::combine_layers()} and of the internal
#' \code{netOmics:::merge_graphs()} (netOmics 1.x, GPL-3, Bodein et al.),
#' restricted to the two cases used to build the liver networks. netOmics is
#' deprecated in Bioconductor (removed in 3.21); internalising these keeps
#' CENetwork installable and its data reproducible in the long run. The original
#' behaviour is preserved as-is, including attribute type coercion, so that
#' networks rebuilt with this function are identical to the published ones.
#'
#' @param graph1 (igraph) layer to extend.
#' @param graph2 (igraph) layer to merge into \code{graph1}. If \code{NULL},
#'   \code{interaction.df} is used instead.
#' @param interaction.df (data.frame) interactions with columns \code{from} and
#'   \code{to}; only rows with at least one endpoint in \code{graph1} are kept.
#'
#' @return an igraph object.
#'
#' @examples
#' library(igraph)
#'
#' g1 <- igraph::make_ring(3)
#' igraph::V(g1)$name <- c("a", "b", "c")
#' igraph::V(g1)$type <- "layer1"
#'
#' g2 <- igraph::make_ring(3)
#' igraph::V(g2)$name <- c("c", "d", "e")
#' igraph::V(g2)$type <- "layer2"
#'
#' g <- combine_layers(g1, g2)
#' igraph::vertex_attr(g, "type")  # "c" keeps the type of graph1
#'
#' @export
combine_layers <- function(graph1, graph2 = NULL, interaction.df = NULL) {
    stopifnot(inherits(graph1, "igraph"))
    if (!is.null(graph2)) {
        stopifnot(inherits(graph2, "igraph"))
        return(merge_graphs(graph1, graph2))
    }
    stopifnot(is.data.frame(interaction.df),
              all(c("from", "to") %in% colnames(interaction.df)))
    keep <- interaction.df$from %in% igraph::V(graph1)$name |
        interaction.df$to %in% igraph::V(graph1)$name
    interaction.graph <- igraph::graph_from_data_frame(interaction.df[keep, , drop = FALSE],
                                                       directed = FALSE)
    merge_graphs(graph1, interaction.graph)
}


#' Merge two graphs by node name
#'
#' Union of two igraph objects, coalescing the vertex attributes they share
#' (value of \code{graph1} first, then \code{graph2}). Internal;
#' see \code{\link{combine_layers}}.
#'
#' @param graph1,graph2 (igraph)
#'
#' @return an igraph object of class \code{merged.igraph}.
#'
#' @noRd
merge_graphs <- function(graph1, graph2) {
    shared_attr <- intersect(names(igraph::vertex_attr(graph1)),
                             names(igraph::vertex_attr(graph2)))
    shared_attr <- shared_attr[shared_attr != "name"]
    merged <- igraph::union(graph1, graph2)
    ma <- igraph::vertex_attr(merged)
    for (sa in shared_attr) {
        # `vector(length = n)` starts as a logical vector: the final attribute
        # type follows from the coercion at the first assignment. netOmics
        # behaviour, kept as-is.
        ma[[sa]] <- vector(length = igraph::vcount(merged))
        for (i in seq_along(ma[[sa]])) {
            ma[[sa]][i] <- ifelse(!is.na(ma[[paste0(sa, "_1")]][i]),
                                  ma[[paste0(sa, "_1")]][i],
                                  ma[[paste0(sa, "_2")]][i])
        }
        merged <- igraph::delete_vertex_attr(merged, paste0(sa, "_1"))
        merged <- igraph::delete_vertex_attr(merged, paste0(sa, "_2"))
        merged <- igraph::set_vertex_attr(merged, sa, value = ma[[sa]])
    }
    class(merged) <- c("merged.igraph", "igraph")
    merged
}


#' Remove isolated nodes
#'
#' Drop the nodes that have no neighbour once multiple edges and self-loops are
#' collapsed. Reimplementation of \code{netOmics:::remove_unconnected_nodes()}
#' (netOmics 1.x, GPL-3); see \code{\link{combine_layers}} for the rationale.
#'
#' @param X (igraph) network
#'
#' @return an igraph object.
#'
#' @noRd
remove_unconnected_nodes <- function(X) {
    isolated <- which(igraph::degree(igraph::simplify(X)) == 0)
    igraph::delete_vertices(X, isolated)
}
