#' Export to Cytoscape
#'
#' This function export an igraph object to cytoscape..
#'
#' @param x igraph object, the network to export
#' @param collectionName character, name of collection inside Cytoscape
#' @param networkName character, name of network inside Cytoscape

#' @examples
#' \donttest{
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#'
#' signature_vids <- signature_maison$acetaminophen_all_all
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, signature_vids, target_type = c("drug/compound", "pathway", "side_effect", "GO"))
#' x <- res.diffusion
#' export_to_cytoscape(res.diffusion)
#' apply_custom_theme()
#' }
#'
#' @importFrom checkmate assert_string
#' @importFrom RCy3 cytoscapePing
#' @importFrom RCy3 createNetworkFromIgraph
#' @importFrom RCy3 setVisualStyle
#' @export
#'
#'
export_to_cytoscape <- function(x, collectionName = "myNewCollection", networkName = "myNetwork"){

    # check x
    stopifnot(is(x, "get_route.res"))

    checkmate::assert_string(collectionName)
    checkmate::assert_string(networkName)

    network <- x$network

    # RCy3
    RCy3::cytoscapePing()  # stop execution if status != 200

    collectionList.names <- RCy3::getCollectionList()
    if(!(collectionName %in% collectionList.names)){
        message(paste0("The collection ", collectionName, " does not exist and will be created."))
    }
    networkList.names <- RCy3::getNetworkList()
    if(networkName %in% networkList.names){
        paste0("The network ", networkName, " already exists and will be overwritten.")
    }

    # create network
    RCy3::createNetworkFromIgraph(igraph = network, title = networkName,
                                  collection = collectionName)

}

#' apply_custom_theme
#'
#' This function apply custom theme to current network in cytoscape. Here is the description of the theme.
#'
#' \describe{
#'   \item{Fill color}{(based on type) gene = blue, protein = orange, GO = grey, side_effect = green, drug/compound = pink, pathway = red}
#'   \item{Node shape}{input diffusion = ELLIPSE, other = ROUND_RECTANGLE}
#'   \item{Node size}{(based on type) gene = 20, protein = 20, GO = 35, side_effect = 40, drug/compound = 40, pathway = 35.}
#'   \item{Node border}{gene_hepatox = red, other = black}
#'   }
#'
#' @seealso export_to_cytoscape
#'
#' @importFrom tibble tribble
#' @importFrom RCy3 setNodeLabelMapping
#' @importFrom RCy3 setEdgeLineWidthDefault
#' @importFrom RCy3 setEdgeColorDefault
#' @importFrom RCy3 setNodeColorMapping
#' @importFrom RCy3 setNodeShapeMapping
#' @importFrom RCy3 setNodeSizeMapping
#' @importFrom RCy3 setNodeBorderWidthDefault
#' @export
apply_custom_theme <- function(){

    # customized default style ------------------------------------------------
    ## display name
    RCy3::setNodeLabelMapping('display_name')

    RCy3::setEdgeLineWidthDefault(new.width = 0.5)
    RCy3::setEdgeColorDefault(new.color = "#d3d3d3")

    ## color
    fill_color <- tibble::tribble(~type, ~color,
                                  "gene", "#99ffff",
                                  "protein", "#ffcc99",
                                  "GO", "#d9d9d9",
                                  "side_effect", "#b3ff99",
                                  "drug/compound", "#ffccff",
                                  "pathway", "#ff9999")

    RCy3::setNodeColorMapping(table.column = "type",
                              mapping.type = 'd',
                              table.column.values = fill_color$type,
                              colors = fill_color$color)

    ## shape based on input diffusion
    node_shape <- tibble::tribble(~input_diffusion, ~shape,
                                  FALSE, "ROUND_RECTANGLE",
                                  TRUE, "ELLIPSE")

    RCy3::setNodeShapeMapping(table.column = 'input_diffusion',
                              table.column.values = node_shape$input_diffusion,
                              shapes = node_shape$shape,
                              default.shape = "ROUND_RECTANGLE")

    # border color
    border_color <- tibble::tribble(~hepatox, ~color,
                                    "TRUE", "#ff0000",
                                    "FALSE", "#000000",
                                    "NA", "#000000")

    RCy3::setNodeBorderWidthDefault(1)
    RCy3::setNodeBorderColorMapping(table.column = 'gene_hepatox',
                                    table.column.values = border_color$hepatox,
                                    mapping.type = "d",
                                    colors = border_color$color,
                                    default.color = "#000000")

    # size
    node_size <- tibble::tribble(~type, ~size,
                                 "gene", 20,
                                 "protein", 20,
                                 "GO", 35,
                                 "side_effect", 40,
                                 "drug/compound", 40,
                                 "pathway", 35)

    RCy3::setNodeSizeMapping(table.column = "type",
                             mapping.type = 'd',
                             table.column.values = node_size$type,
                             sizes = node_size$size,
                             default.size =10)
}




