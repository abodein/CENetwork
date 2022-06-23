#' Check Cytoscpae connection
#'
#' @importFrom RCy3 cytoscapePing
#' @export
check_Cytoscape_connection <- function(){
    RCy3::cytoscapePing()
}
