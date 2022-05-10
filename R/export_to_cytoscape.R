#' @importFrom NetPathMiner plotCytoscapeGML
#' @export
export_to_cytoscape <- function(x, collectionName = "myNewCollection", networkName = "myNetwork"){

    # check x
    stopifnot(is(x, "get_route.res"))

    checkmate::assert_string(collectionName)
    checkmate::assert_string(networkName)

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
    RCy3::createNetworkFromIgraph(igraph = x$network, title = networkName,
                                  collection = collectionName)
    # from dataframes
    # edge_node_list <- convert_get_route_to_data.frames(x)
    # RCy3::createNetworkFromDataFrames(edge_node_list$node,edge_node_list$edge,
    #                                   title=networkName, collection = collectionName)


    # build and set visual style
    build_visual_style(style.name = "cosmEU_Style", network = networkName)
    RCy3::setVisualStyle("cosmEU_Style", network = networkName)
    # hideEdges()
}

#' build visual style
build_visual_style <- function(style.name = "cosmEU_Style", network){

    defaults <- list(NODE_SHAPE="ROUND_RECTANGLE",
                     NODE_SIZE=30,
                     EDGE_TRANSPARENCY=120,
                     NODE_LABEL_POSITION="C,C,c,0.00,0.00")

    nodeLabels <- RCy3::mapVisualProperty('node label','id','p')

    # fill color
    fill_color <- tribble(
        ~type, ~hex,
      "ENSEMBL", "#96cbff",
      "UNIPROT", "#FFA500",
        "SMPDB", "#f77272",
     "DrugBank", "#fcc2f6",
       "CHEMBL", "#00FF00")

    nodeFills <- RCy3::mapVisualProperty(visual.prop = 'node fill color',
                                   table.column = 'type',
                                   mapping.type = 'd',
                                   table.column.values = fill_color$type,
                                   visual.prop.values = fill_color$hex,
                                   network = network)

    nodeLabels <- RCy3::mapVisualProperty('node label','display_name','p')
    # edgeLineStyle <- mapVisualProperty('EDGE_LINE_TYPE', table.column = "interaction",
    #                                    mapping.type = 'd',
    #                                    table.column.values = c("interacts with","closest"),
    #                                    visual.prop.values = c("SOLID", "DOTS"),
    #                                    network = network)

    # arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',
    #                                  c("activates","inhibits","interacts"),c("Arrow","T","None"))
    # edgeWidth <- mapVisualProperty('edge width','weight','p')

    # and then create the style
    RCy3::createVisualStyle(style.name, defaults, list(nodeFills, nodeLabels))
    # createVisualStyle(style.name, defaults, list(nodeFills, nodeLabels, edgeLineStyle))
}

#' @importFrom NetPathMiner plotCytoscapeGML
export_to_cytoscape_w_netPathMiner <- function(x, file, collectionName = "myNewCollection", networkName = "myNetwork"){

    # check x
    stopifnot(is(x, "igraph"))

    # check file
    stopifnot(checkmate::checkPathForOutput(file))

    # check layout
    NetPathMiner::plotCytoscapeGML(x, file=file)
}

convert_igraph_to_dataframes <- function(x){
    Vx <- vertex_attr(x) %>% as.data.frame()
    Ex <- x %>% as_long_data_frame() %>% dplyr::select(from_name, to_name, from_type, to_type) %>%
        mutate(interaction = "interacts with")

    return(list(V = Vx, E = Ex))
}

convert_get_route_to_data.frames <- function(x){
    tmp <- convert_igraph_to_dataframes(x$network)

    join.df <- as.data.frame(tmp$V) %>% dplyr::select(name, type)

    Ex_closest <- x$closest %>% dplyr::select(SeedName, NodeNames) %>%
        left_join(join.df, c("SeedName" = "name")) %>%
        dplyr::rename("from_type" = "type") %>%
        left_join(join.df, c("NodeNames" = "name")) %>%
        purrr::set_names("from_name", "to_name", "from_type", "to_type") %>%
        mutate(interaction = "closest")

    new_E <- rbind(tmp$E, Ex_closest) %>%
        purrr::set_names("source", "target", "from_type", "to_type", "interaction")
    new_V <- tmp$V %>% mutate(id = name)

    return(list(node = tmp$V, edge = new_E))
}

