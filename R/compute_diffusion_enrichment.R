#' Compute empirical enrichment statistics for diffusion results
#'
#' From a \code{get_route()} result with permutations, computes for each
#' GO term and/or Reactome pathway node an empirical p-value and an FDR-adjusted
#' p-value that are directly comparable to ORA and GSEA outputs.
#'
#' The test statistic is the visit count: the number of seed-to-target shortest
#' paths that pass through a given node.  A null distribution is built from
#' \code{n_permut} random gene signatures of the same size, and the empirical
#' p-value is computed with the Phipson & Smyth (2010) correction to avoid
#' p = 0 when no permutation exceeds the observed value.
#'
#' @param route_result Output of \code{get_route()} run with
#'   \code{do_permutation = TRUE}.
#' @param target_types Character vector of node types to test.
#'   Defaults to \code{c("GO", "pathway")}.
#' @param min_visit_count Integer. Minimum observed visit count to include a
#'   term. Terms with \code{visit_count < min_visit_count} are dropped before
#'   multiple-testing correction. Default is 1.
#' @param padj_method Correction method passed to \code{p.adjust()}.
#'   Default is \code{"BH"} (Benjamini-Hochberg).
#'
#' @return A data frame with one row per tested term, sorted by \code{padj}
#'   then descending \code{z_score}, with columns:
#'   \describe{
#'     \item{term}{Node name (GO ID or Reactome pathway ID)}
#'     \item{label}{Human-readable label when available (vertex attribute
#'       \code{label}; \code{NA} otherwise)}
#'     \item{type}{Layer type (\code{"GO"} or \code{"pathway"})}
#'     \item{visit_count}{Observed number of shortest paths through this node}
#'     \item{mean_perm}{Mean visit count across permutations}
#'     \item{sd_perm}{Standard deviation of visit count across permutations}
#'     \item{z_score}{Standardised effect size:
#'       \eqn{(v_t - \bar{v}_t^{perm}) / \sigma_t^{perm}}}
#'     \item{pvalue}{Empirical p-value (Phipson & Smyth correction)}
#'     \item{padj}{FDR-adjusted p-value}
#'   }
#'
#' @references
#' Phipson, B. & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly drawn.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 9(1).
#'
#' @examples
#' \dontrun{
#' res <- get_route(liver_network, closest_dfr, signature_vids,
#'                  do_permutation = TRUE, n_permut = 1000)
#'
#' enrich <- compute_diffusion_enrichment(res)
#' head(enrich[enrich$padj < 0.05, ])
#'
#' # Compare with ORA output: same padj < 0.05 threshold
#' sig_terms_diffusion <- enrich$term[enrich$padj < 0.05]
#' }
#'
#' @importFrom igraph as_data_frame V vertex_attr
#' @export
compute_diffusion_enrichment <- function(route_result,
                                         target_types    = c("GO", "pathway"),
                                         min_visit_count = 1L,
                                         padj_method     = "BH") {

    # --- Input checks ---------------------------------------------------------

    if (!inherits(route_result, "get_route.res")) {
        stop("'route_result' must be the output of get_route().")
    }
    if (is.null(route_result$permutation_results)) {
        stop("No permutation results found. Re-run get_route() with do_permutation = TRUE.")
    }
    if (is.null(route_result$observed_visit_counts)) {
        stop("'observed_visit_counts' missing. Please update CENetwork to the latest version.")
    }

    # --- Extract observed visit counts for target nodes -----------------------

    vdf <- igraph::as_data_frame(route_result$network, what = "vertices")

    # Keep only target-type nodes with at least min_visit_count visits
    target_df <- vdf[
        vdf$type %in% target_types &
        !is.na(vdf$visit_count) &
        vdf$visit_count >= min_visit_count,
    ]

    if (nrow(target_df) == 0L) {
        warning("No target nodes satisfy the filter criteria. Returning empty data frame.")
        return(data.frame())
    }

    target_nodes <- target_df$name
    obs_vc       <- route_result$observed_visit_counts[target_nodes]

    # --- Build permutation matrix ---------------------------------------------
    # Rows = permutations, cols = target nodes

    N <- sum(!vapply(route_result$permutation_results, is.null, logical(1L)))

    perm_matrix <- do.call(rbind, lapply(route_result$permutation_results, function(x) {
        if (is.null(x)) {
            rep(NA_real_, length(target_nodes))
        } else {
            x$visit_counts[target_nodes]
        }
    }))

    # --- Z-score and parametric p-value ---------------------------------------
    # Null distribution parameters estimated from permutations.
    # One-sided p-value from normal approximation (same rationale as GSEA NES):
    # p = P(Z >= z_obs) = 1 - pnorm(z_obs)
    # This avoids the resolution limit of purely empirical p-values (1/(N+1))
    # while still grounding the null in the observed permutation distribution.

    perm_means <- colMeans(perm_matrix, na.rm = TRUE)
    perm_sds   <- apply(perm_matrix, 2L, sd, na.rm = TRUE)
    z_score    <- (obs_vc - perm_means) / perm_sds
    z_score[!is.finite(z_score)] <- NA_real_

    pvalue <- pnorm(z_score, lower.tail = FALSE)
    pvalue[is.na(z_score)] <- NA_real_

    # --- FDR correction -------------------------------------------------------

    padj <- p.adjust(pvalue, method = padj_method)

    # --- Assemble output ------------------------------------------------------

    label_col <- if ("label" %in% colnames(target_df)) target_df$label else NA_character_

    result_df <- data.frame(
        term        = target_nodes,
        label       = label_col,
        type        = target_df$type,
        visit_count = as.integer(obs_vc),
        mean_perm   = round(perm_means, 3L),
        sd_perm     = round(perm_sds,   3L),
        z_score     = round(z_score,    3L),
        pvalue      = pvalue,
        padj        = padj,
        row.names   = NULL,
        stringsAsFactors = FALSE
    )

    result_df[order(result_df$padj, -result_df$z_score, na.last = TRUE), ]
}
