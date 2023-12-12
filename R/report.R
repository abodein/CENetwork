#' Report
#'
#' Reporting remarkable nodes within the network. The function produce a list of data.frame containing the valuable informations.
#'
#' @return
#' A list of data.frame:
#' \describe{
#'   \item{`degree_signature`}{degree of input nodes and hepatox infos}
#'   \item{`drugs`}{drug infos with DILI score, IC50}
#'   \item{`drugs.targets`}{drug to protein target}
#'   \item{`drugs.side_effect`}{drug side effects}
#'   \item{`drug.dist.signature`}{closest input node(s) for each drug and their distance}
#'   \item{`pathways`}{pathway infos}
#'   \item{`pathways.dist.signature`}{closest input node(s) for each pathway and their distance}
#'   \item{`GOs`}{GO infos}
#'   \item{`GO.dist.signature`}{closest input node(s) for each GO and their distance}
#'   }
#'
#' @param x diffusion results from get_route()
#'
#' @examples
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#' signature_vids <-  signature_maison$`acetaminophen_all_all`
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, c(signature_vids, "DB00313"), target_type = c("drug/compound", "pathway"))
#' res_report <- report(res.diffusion)
#'
#' # with module enrichment
#'
#'
#' @importFrom igraph vertex_attr make_ego_graph distances degree
#' @importFrom dplyr filter select mutate bind_rows rename_with left_join pull group_by top_n ungroup arrange across
#' @importFrom purrr imap set_names is_logical
#' @importFrom tibble rownames_to_column
#' @importFrom  tidyr pivot_longer replace_na
#' @export
report <- function(x){
#    checkmate::check_access("")

    stopifnot(is(x, "get_route.res"))

    va <- igraph::vertex_attr(x$network) %>% as.data.frame()

    # 1. find drugs and DILI scores
    drugs <- va %>% dplyr::filter(type == "drug/compound") %>%
        dplyr::select(name, drug_name, drugbank_id, chembl_id, DILI_severity_class, vDILIConcern, drug_has_side_effect)

    # 1.1 find protein target
    drugs.ego <- igraph::make_ego_graph(graph = x$network, order = 1, nodes = drugs$name)
    names(drugs.ego) <- drugs$name
    drugs.targets <- purrr::imap(drugs.ego,
                                 ~{igraph::vertex_attr(.x) %>% as.data.frame %>%
                                         dplyr::filter(type == "protein") %>%
                                         dplyr::select(name) %>%
                                         dplyr::mutate(drug_seed = .y)}) %>%
        dplyr::bind_rows() %>%
        purrr::set_names(c("drug_target", "name"))  %>%  # name = drug id
        dplyr::left_join(va %>% dplyr::select(name, display_name, protein_name) %>% dplyr::rename_with(~paste0(.x, ".target")), by = c("drug_target" = "name.target")) %>%
        dplyr::left_join(va %>% dplyr::select(name, display_name) %>% dplyr::rename_with(~paste0(.x, ".drug")), by = c("name" = "name.drug")) %>%
        dplyr::select(display_name.drug, name, display_name.target, protein_name.target) %>%
        purrr::set_names("drug_name", "drug_id", "target_id", "target_name")

    # 1.2 drug distance to signature
    drug.dist.signature.gene <- igraph::distances(graph = x$network, to = drugs$name,
                                                  v = va %>% dplyr::filter(input_diffusion) %>% dplyr::pull(name)) %>%
        as.data.frame() %>%
    # purrr::imap_dfr(drug.dist.signature.gene, ~ data.frame(signature_vids = rownames(drug.dist.signature.gene)[which(.x == min(.x))]) %>%
    #                 mutate(drug = .y)) %>%
        tibble::rownames_to_column("signature_vids") %>%
        tidyr::pivot_longer(names_to = "drug", values_to = "distance", -signature_vids) %>%
        dplyr::group_by(drug) %>% dplyr::top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
        # left_join with va
        dplyr::left_join(va %>% dplyr::select(name, display_name, gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% dplyr::rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
        dplyr::left_join(va %>% dplyr::select(name, display_name) %>%
                         dplyr::rename_with(~paste0(.x, ".drug")), by = c("drug" = "name.drug")) %>%
        dplyr::select(display_name.drug, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, drug, signature_vids) %>%
        purrr::set_names(c("drug_name","signature_name", "distance", "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "drug_node_id", "signature_node_id"))

    # 1.3 drug side effects
    drugs.se <- purrr::imap(drugs.ego,
                                 ~{igraph::vertex_attr(.x) %>% as.data.frame %>%
                                     dplyr::filter(type == "side_effect") %>%
                                         dplyr::select(name) %>%
                                     dplyr::mutate(drug_seed = .y)}) %>%
        dplyr::bind_rows() %>%
        purrr::set_names(c("drug_side_effect", "name"))  %>%  # name = drug id
      dplyr::left_join(va %>% dplyr::select(name, display_name) %>% dplyr::rename_with(~paste0(.x, ".target")), by = c("drug_side_effect" = "name.target")) %>%
      dplyr::left_join(va %>% dplyr::select(name, display_name) %>% dplyr::rename_with(~paste0(.x, ".drug")), by = c("name" = "name.drug")) %>%
        dplyr::select(display_name.drug, name, display_name.target) %>%
        purrr::set_names("drug_name", "drug_id", "side_effect")

     # 2. degree signature
     degree_signature <- data.frame(degree = igraph::degree(graph = x$network, v = va %>% dplyr::filter(input_diffusion) %>% dplyr::pull(name))) %>%
         tibble::rownames_to_column("node_id") %>% dplyr::left_join(va, by = c("node_id" = "name")) %>%
         dplyr::select("node_id", "display_name", "degree", "gene_hepatox_Toxicology2014","gene_hepatox_ToxicologyInVitro2020", "input_protein_signature") %>%
         dplyr::arrange(desc(degree)) %>%
         purrr::set_names("node_id", "display_name", "degree", "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "coding_protein_is_present")


    # 3. pathway
    # 3.1 get pathways
     pathways <- va %>% dplyr::filter(type == "pathway") %>%
         dplyr::select(name, pathway_name, pathway_id, pathway_database)

     # 3.2 get pathways
     pathway.dist.signature <- igraph::distances(graph = x$network, to = pathways$name,
                                                 v = va %>% dplyr::filter(input_diffusion) %>% pull(name)) %>%
         as.data.frame() %>%
         tibble::rownames_to_column("signature_vids") %>%
         tidyr::pivot_longer(names_to = "pathway", values_to = "distance", -signature_vids) %>%
       dplyr::group_by(pathway) %>% dplyr::top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
         # left_join with va
       dplyr::left_join(va %>% dplyr::select(name, display_name,  gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% dplyr::rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
       dplyr::left_join(va %>% dplyr::select(name, display_name) %>% dplyr::rename_with(~paste0(.x, ".pathway")), by = c("pathway" = "name.pathway")) %>%
         dplyr::select(display_name.pathway, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, pathway, signature_vids) %>%
         purrr::set_names(c("pathway_name","signature_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "display_node_id", "signature_node_id"))

     # 4) GO terms
     #
     # 4.1 GO info
     GOs <- va %>% dplyr::filter(type == "GO") %>%
         dplyr::select(name, go_id, go_ontology, go_term_name)

     # 4.2 GO distance to signature
     #signature_in_gene <- va %>% filter(input_gene_signature) %>% pull(name)

     if(!purrr::is_empty(GOs$name)){

     GO.dist.signature <- igraph::distances(graph = x$network, to = GOs$name,
                                                 v = va %>% dplyr::filter(input_diffusion) %>% dplyr::pull(name)) %>%

     # signature.dist.GO <- igraph::distances(graph = x$network, to = signature_in_gene,
     # v = GOs$name) %>%
         as.data.frame()  %>%
         tibble::rownames_to_column("signature_vids") %>%
         tidyr::pivot_longer(names_to = "GO", values_to = "distance", -signature_vids) %>%
       dplyr::group_by(GO) %>%
         # shortest distance
         dplyr::top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
         # left_join with va
         left_join(va %>% dplyr::select(name, display_name,  gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% dplyr::rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
         left_join(va %>% dplyr::select(name, display_name) %>% dplyr::rename_with(~paste0(.x, ".go")), by = c("GO" = "name.go")) %>%
         dplyr::select(display_name.go, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, GO, signature_vids) %>%
         purrr::set_names(c("GO_term","signature_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "display_node_id", "signature_node_id"))
     } else {
       GO.dist.signature <- data.frame()
}

     # autoroute, most visited nodes
     visite_autoroute <- x$autoroute %>% dplyr::group_by(path) %>% dplyr::summarise(val = n()) %>% dplyr::arrange(dplyr::desc(val)) %>%
         purrr::set_names(c("name", "nb.visite")) %>%
         left_join(va, by = "name") %>%
         dplyr::select(name, nb.visite, type, display_name)


     ## update: enrichment in module
     # convert signature to uniprot
     enrich_input <- AnnotationDbi::select(x = org.Hs.eg.db,
                                           keys = c(x$input, va %>% dplyr::pull(name)),
                                           columns = "UNIPROT", keytype = "ENSEMBL") %>%
         na.omit() %>% pull(UNIPROT)


     enrich_res <- module_ppi_M1_20231211 %>% mutate(is_present = molecule %in% enrich_input) %>%
         group_by(module) %>%
         summarise(k = sum(is_present), # nb of molecule in input in the module
                   n = n()) %>% # taille module
         mutate(K = sum(k),  # nb of molecule in all modules
                N = sum(n)) %>%   # size all modules
         nest(enrich_module_count = -module) %>%
         mutate(enrich_module_contingency = imap(enrich_module_count, ~make_contingency(k = .x$k, K = .x$K, n = .x$n, N = .x$N))) %>%
         mutate(enrich_module_p.value = map_dbl(enrich_module_contingency, ~fisher.test(.x, alternative = "greater")$p.value)) %>%
         mutate(enrich_module_p.value_adj = p.adjust(.$enrich_module_p.value, method = "fdr"))
     dplyr::select(c(module, enrich_module_count, enrich_module_p.value, enrich_module_p.value_adj)) %>% unnest(cols = c(enrich_module_count)) %>% arrange(enrich_module_p.value)



     to_return <- list()
     to_return[["drug.dist.signature"]] <- drug.dist.signature.gene
     to_return[["pathway.dist.signature"]] <- pathway.dist.signature
     to_return[["drugs.targets"]] <- drugs.targets
     to_return[["degree_signature"]] <- degree_signature
     to_return[["drugs"]] <- drugs
     to_return[["pathways"]] <- pathways
     to_return[["drugs.side_effect"]] <- drugs.se
     to_return[["GOs"]] <- GOs
     to_return[["GO.dist.signature"]] <- GO.dist.signature
     to_return[["visite_autoraute"]] <- visite_autoroute
     to_return[["module_enrichment"]] <- enrich_res


    # replace NA where logical
     to_return <- lapply(to_return, function(x) {
         x %>% mutate(dplyr::across(where(purrr::is_logical), ~tidyr::replace_na(data = .x, replace = FALSE)))
     })

     return(to_return)
}

#' Produce diffusion report
#'
#' Produce diffusion report based remarkable nodes in the network
#'
#' @param res_report `report()` results
#' @param report_title (character) name of the report
#' @param report_out_filepath (character) filepath where to write the .Rmd report
#' @param overwrite (logical, default = FALSE) if TRUE, overwrite out filepath
#' @param render (logical, default = FALSE) if TRUE, render the html file
#'
#'
#' @seealso
#' report
#'
#' @examples
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#' signature_vids <-  signature_maison$`acetaminophen_all_all`
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, c(signature_vids), target_type = c("drug/compound", "pathway", "side_effect"))
#'
#' res_report <- report(res.diffusion)
#' report_info_df <- produce_diffusion_report(res_report = res_report,  report_title = "acetaminophen", report_out_filepath = "report_acetaminophen_example.Rmd",
#' render = TRUE, overwrite = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom checkmate assert_logical assertPathForOutput
#' @importFrom purrr imap_dfr
#' @importFrom rmarkdown render
#' @importFrom stringr str_remove_all str_replace_all
#' @export
produce_diffusion_report <- function(res_report,
                                     report_title,
                                     report_out_filepath = NA,
                                     overwrite = FALSE,
                                     render = FALSE){


    # Check arguments ---------------------------------------------------------

    ## res_report
    stopifnot(is(res_report, "list"))
    stopifnot(all(unlist(lapply(res_report, function(x)is(x, "data.frame")))))
    stopifnot(all(names(res_report) %in%  c('drug.dist.signature','pathway.dist.signature',
                                            'drugs.targets','degree_signature','drugs','pathways','drugs.side_effect','GOs','GO.dist.signature')))

    #report title
    stopifnot(is(report_title, "character"))

    # render and overwrite
    checkmate::assert_logical(x = render)
    checkmate::assert_logical(x = overwrite)


    # report_out_filepath
    if(!is.na(report_out_filepath)){
        checkmate::assertPathForOutput(x = report_out_filepath, overwrite = overwrite)
    }


    # Process -----------------------------------------------------------------


    filpath_report_res <- file.path(tempdir(),
                                    paste0(as.character(Sys.time()) %>% stringr::str_replace_all(" ", "-") %>% stringr::str_remove_all(":"),
                                           "_res_diffusion.Rds"))
    saveRDS(object = res_report, file = filpath_report_res)

    report_info_df = list()
    counter_report <- 1

intro <- paste0("---
title: '", report_title, "'
date: \"`r Sys.Date()`\"
output:
  BiocStyle::html_document:
      df_print: paged

---

```{r, echo = FALSE}
knitr::opts_chunk$set( echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center')
library(DT)
```

")

    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             # value = paste0(
                                             #     "---\ntitle: '", report_title, "'\ndate: \"`r Sys.Date()`\"\noutput: html_document\n---\n\n",
                                             #     "```{r, echo = FALSE}\n knitr::opts_chunk$set( echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center')\n
                                             #     library(datatable)```"))
                                             value = intro)

    # get data or load data
load_data <- paste0("
```{r}
res_report <- readRDS('", filpath_report_res,"')
```")
    counter_report <- counter_report + 1
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             value =  load_data)

# show data
core <- "
# Most connected node from signature input

The following table returns the list of molecules as input for diffusion, ranked from highest number of connections (degree) to lowest.

Columns description:

* node_id: Node identifier
* display_name: The node name displayed on the network
* degree: The number of node connections
* gene_hepatox_Toxicology2014: True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* coding_protein_is_present: True or False, if a gene node has its coding protein in the subnetwork.

```{r}
datatable(res_report$degree_signature)
```

# Drugs

Focus on drug nodes.
The following section provides information on the *drug* nodes returned in the subnetwork.

## Drug info

The following table returns general information on drugs contained in the subnetwork.

Columns description:


* name: Node identifier
* drug_name: Drug most commom name, if aivailable
* drugbank_id: DrugBank database identifier
* chembl_id: ChEMBL database identifier
* DILI_severity_class: Drug Induced Liver Injury severity class (from LTKB database)
* vDILIConcern: Drug Induced Liver Injury concern (from LTKB database)
* drug_has_side_effect: True or False, if the drug has reported side effects.

```{r}
datatable(res_report$drugs)
```

## Drug target

The following table lists all the drug targets in the subnetwork.

Column description:

* drug_id: drug node identifier
* drug_name: drug most common name
* target_id: drug target Uniprot identifier
* target_name: target function

```{r}
datatable(res_report$drugs.targets)
```

## Drug Side Effects

The following table lists all the drug side effects in the subnetwork.

Column description:

* drug_name: drug most common name
* drug_id: drug node identifier
* side_effect: drug induced side effect present in the subnetwork.

```{r}
datatable(res_report$drugs.side_effect)
```

## Drug shortest distance to signature

In order to refine the subnetwork, the table below returns the closest node in the diffusion input for each drug.

Column description:

* drug_name: drug most common name
* signature_name: displayed input node for diffusion, closest node from the drug node.
* distance: geodesic distance between the drug node and input node.
* gene_hepatox_Toxicology2014: True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* drug_node_id: drug node identifier
* signature_node_id: input node identifier for diffusion.

```{r}
datatable(res_report$drug.dist.signature %>% arrange(distance))
```

# Pathways

The following section provides information on the *pathway* nodes returned in the subnetwork.

Focus on drug nodes.

## Pathway infos

General informations about the pathway nodes in the subnetworks.

Column description:

* name: displayed pathway name
* pathway_name: pathway description
* pathway_id: node identifier
* pathway_database: database of pathway origin

```{r}
datatable(res_report$pathways)
```

## Pathways shortest distance to signature

To refine the subnetwork, the table below returns the closest node in the diffusion input from each pathway node.

Column description:

* pathway_name: pathway description
* signature_name: displayed input node for diffusion, closest node from the pathway node.
* distance: geodesic distance between the drug node and input node.
* gene_hepatox_Toxicology2014: True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* display_node_id: displayed pathway node id.
* signature_node_id: input node identifier for diffusion.

```{r}
datatable(res_report$pathway.dist.signature)
```

# GO terms hepatox

The following section provides information on the *GO_term* nodes returned in the subnetwork.


## GO terms infos

General informations about the GO terms nodes in thesubnetwork.

Columns description:

* name: GO displayed name
* go_id: GO node identifier
* go_ontology: GO term ontology
* go_term_name: GO term description

```{r}
datatable(res_report$GOs)

```
## GO shortest distance to signature

To refine the subnetwork, the table below returns the closest node in the diffusion input from each pathway node.

Column description:

* GO_term: GO term description
* signature_name: displayed input node for diffusion, closest node from the GO node.
* distance: geodesic distance between the drug node and input node.
* gene_hepatox_Toxicology2014: True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* display_node_id: displayed pathway node id.
* signature_node_id: input node identifier for diffusion.

```{r}
datatable(res_report$GO.dist.signature)
```

# Most visited nodes

Column description:

* name :
* nb.visite : number of time
* type : type of node
* display_name : displayed node name in the network.

Each input seed generates a path to the final subnetwork, and some parts of the entire network may be visited multiple times. In this table, we list the most visited nodes during the diffusion process.
```{r}
datatable(res_report$visite_autoraute)
```

# Module enrichment

Enrichment analysis from diffusion input and diffusion results against modules.
Node names are converted to UNIPROT_ID to perform enrichment.

Column description:

* module: module identifier.
* k: number of diffusion targets inside the module
* n: size of the module
* K: sum of diffusion targets inside all modules
* N: size of all module
* enrich_module_p.value: Hypergeometric test pvalue
* enrich_module_p.value_adj: Corrected pvalue (fdr)

```{r}
datatable(res_report$module_enrichment)
```
"

    counter_report <- counter_report + 1
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             value =  core)

    report_info_df <- purrr::imap_dfr(report_info_df, ~.x)

    lines <- produce_report(report_infos = report_info_df,
                            report_filename = report_out_filepath)

    if(render){
        rmarkdown::render(report_out_filepath)
    }

    return(invisible(report_info_df))
}



#' @importFrom purrr map
produce_report <- function(report_infos, report_filename = "report.Rmd") {

    stopifnot(is(report_infos, "data.frame") | is(report_infos, "character"))
    stopifnot(nrow(report_infos) > 0)
    stopifnot(all(c("add", "value") %in% colnames(report_infos)))
    stopifnot(all(report_infos$add %in% c("text", "file", "plot")))

    for (i in 1:nrow(report_infos)) {
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]
        if (current_add == "text") {
            stopifnot(is(current_value, "character"))
        }

        if (!is.null(report_filename)) {
            stopifnot(is(report_filename, "character"))
            if (dirname(report_filename) != ".") {
                stopifnot(dir.exists(dirname(report_filename)))
            }
        }
    }


    # 2. Parse report_infos

    produce_lines <- function(i) {
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]
        if (current_add == "text") {
            lines <- c(current_value, "\n")
        }
        return(lines)
    }
    all_lines <- purrr::map(1:nrow(report_infos), produce_lines) %>% unlist


    if (!is.null(report_filename)) {
        file_conn<-file(report_filename)
        writeLines(all_lines, file_conn)
        close(file_conn)
    }

    return(invisible(lines))

}

make_contingency <- function(k, K, n, N){
    dat <- matrix(data = c(k, K-k,
                           n-k, N-n-k-K), ncol = 2, byrow = TRUE)

    colnames(dat) <- c("gene.in.interest", "gene.not.interest")
    rownames(dat) <- c("in.category", "not.in.category")
    return(dat)
}
