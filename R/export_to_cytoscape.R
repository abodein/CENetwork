#' @importFrom NetPathMiner plotCytoscapeGML
#' @export
export_to_cytoscape <- function(x, file){

    # check x
    stopifnot(is(x, "igraph"))

    # check file
    stopifnot(checkmate::checkPathForOutput(file))

    # check layout

    NetPathMiner::plotCytoscapeGML(x, file=file)
}
