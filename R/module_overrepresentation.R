#' Set-level over-representation of toxic modules among diffusion-reached modules
#'
#' Tests whether a set of \dQuote{toxic} modules (e.g. PPI modules enriched in
#' hepato-toxic GO terms, DILI targets or liver-associated side effects) is
#' over-represented among the modules \emph{reached} by a diffusion subnetwork.
#'
#' This complements the per-module enrichment test (\code{report()}), which is
#' conservative once corrected across all modules and therefore under-powered
#' for small signatures. Pooling the evidence at the set level asks a more
#' powerful question: rather than \dQuote{is any single module significantly
#' over-represented?}, it asks \dQuote{do the modules reached by the diffusion
#' contain more toxic modules than expected by chance?}. A module is considered
#' \emph{reached} (touched) when at least one of the diffusion subnetwork nodes
#' of \code{node_type} belongs to it.
#'
#' @param route_result Output of \code{get_route()} (class \code{get_route.res}).
#' @param toxic_modules Character vector of module identifiers considered toxic
#'   (e.g. \code{radar_plots_data$module_GOTox$module} from \code{report()}).
#' @param modules Data frame with columns \code{module} and \code{molecule}
#'   mapping each module to its member molecules. Defaults to the package data
#'   \code{module_ppi_M1_20231211}.
#' @param node_type Node type of the diffusion subnetwork used to define
#'   \dQuote{reached} molecules (default \code{"protein"}).
#' @param id_attr Vertex attribute holding the molecule identifier that matches
#'   \code{modules$molecule} (default \code{"UNIPROT"}).
#' @param alternative Alternative hypothesis for \code{fisher.test}
#'   (default \code{"greater"}, i.e. over-representation).
#'
#' @return A list with:
#'   \item{n_modules}{total number of modules}
#'   \item{n_toxic}{number of toxic modules}
#'   \item{n_touched}{number of modules reached by the diffusion}
#'   \item{n_toxic_touched}{number of toxic modules reached}
#'   \item{expected}{toxic modules expected among reached under the null}
#'   \item{odds_ratio}{Fisher odds ratio}
#'   \item{p_value}{Fisher exact test p-value}
#'   \item{toxic_touched}{identifiers of the toxic modules reached}
#'
#' @examples
#' \dontrun{
#' res <- get_route(liver_network, closest_dfr, signature_vids)
#' rep <- report(res, complete_network = liver_network)
#' tox <- rep$radar_plots_data$module_GOTox$module
#' module_overrepresentation(res, toxic_modules = tox)
#' }
#'
#' @importFrom igraph as_data_frame
#' @importFrom stats fisher.test na.omit
#' @export
module_overrepresentation <- function(route_result, toxic_modules,
                                      modules = NULL,
                                      node_type = "protein",
                                      id_attr = "UNIPROT",
                                      alternative = "greater") {

    if (!inherits(route_result, "get_route.res")) {
        stop("'route_result' must be the output of get_route().")
    }
    if (!is.character(toxic_modules) || length(toxic_modules) < 1L) {
        stop("'toxic_modules' must be a non-empty character vector of module ids.")
    }

    if (is.null(modules)) {
        env <- new.env()
        utils::data("module_ppi_M1_20231211", package = "CENetwork", envir = env)
        modules <- env[["module_ppi_M1_20231211"]]
    }
    if (!all(c("module", "molecule") %in% colnames(modules))) {
        stop("'modules' must contain the columns 'module' and 'molecule'.")
    }

    # molecules reached by the diffusion (subnetwork nodes of the given type)
    va <- igraph::as_data_frame(route_result$network, what = "vertices")
    if (!id_attr %in% colnames(va)) {
        stop(sprintf("vertex attribute '%s' not found in the subnetwork.", id_attr))
    }
    reached <- unique(stats::na.omit(va[[id_attr]][va$type == node_type]))

    all_modules   <- unique(modules$module)
    touched       <- unique(modules$module[modules$molecule %in% reached])
    toxic         <- intersect(toxic_modules, all_modules)
    toxic_touched <- intersect(touched, toxic)

    n_total       <- length(all_modules)
    n_toxic       <- length(toxic)
    n_touched     <- length(touched)
    n_toxic_touch <- length(toxic_touched)

    # 2 x 2 contingency: reached (yes/no) x toxic (yes/no)
    a <- n_toxic_touch                 # reached & toxic
    b <- n_touched - a                 # reached & non-toxic
    c <- n_toxic  - a                  # not reached & toxic
    d <- n_total  - n_touched - c      # not reached & non-toxic
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE),
                             alternative = alternative)

    list(
        n_modules       = n_total,
        n_toxic         = n_toxic,
        n_touched       = n_touched,
        n_toxic_touched = n_toxic_touch,
        expected        = n_touched * n_toxic / n_total,
        odds_ratio      = unname(ft$estimate),
        p_value         = ft$p.value,
        toxic_touched   = sort(toxic_touched)
    )
}
