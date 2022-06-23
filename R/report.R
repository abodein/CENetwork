#' Report
#'
#' Reporting remarkable nodes within the network
#'
#' @param x
#'
#' @examples
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#' signature_vids <-  signature_maison$`valproic acid_all_all`
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, c(signature_vids, "DB00313"), target_type = c("drug/compound", "pathway", "side effect"))
#'
#' x <- res.diffusion
#' res_report <- report(res.diffusion)
#'
report <- function(x, pathout = ""){
#    checkmate::check_access("")

    stopifnot(is(x, "get_route.res"))

    va <- vertex_attr(x$network) %>% as.data.frame()

    # 1. find drugs and DILI scores
    drugs <- va %>% dplyr::filter(type == "drug/compound") %>%
        dplyr::select(name, drug_name, drugbank_id, chembl_id, DILI_severity_class, vDILIConcern, drug_has_side_effect)

    # 1.1 find protein target
    drugs.ego <- igraph::make_ego_graph(graph = x$network, order = 1, nodes = drugs$name)
    names(drugs.ego) <- drugs$name
    drugs.targets <- purrr::imap(drugs.ego,
                                 ~{vertex_attr(.x) %>% as.data.frame %>%
                                         filter(type == "protein") %>%
                                         dplyr::select(name) %>%
                                         mutate(drug_seed = .y)}) %>%
        dplyr::bind_rows() %>%
        purrr::set_names(c("drug_target", "name"))  %>%  # name = drug id
        left_join(va %>% dplyr::select(name, display_name, protein_name) %>% rename_with(~paste0(.x, ".target")), by = c("drug_target" = "name.target")) %>%
        left_join(va %>% dplyr::select(name, display_name) %>% rename_with(~paste0(.x, ".drug")), by = c("name" = "name.drug")) %>%
        dplyr::select(display_name.drug, name, display_name.target, protein_name.target) %>%
        purrr::set_names("drug_name", "drug_id", "target_id", "target_name")

    # 1.2 drug distance to signature
    drug.dist.signature.gene <- igraph::distances(graph = x$network, to = drugs$name,
                                                  v = va %>% dplyr::filter(input_diffusion) %>% pull(name)) %>%
        as.data.frame() %>%
    # purrr::imap_dfr(drug.dist.signature.gene, ~ data.frame(signature_vids = rownames(drug.dist.signature.gene)[which(.x == min(.x))]) %>%
    #                 mutate(drug = .y)) %>%
        tibble::rownames_to_column("signature_vids") %>%
        tidyr::pivot_longer(names_to = "drug", values_to = "distance", -signature_vids) %>%
        group_by(drug) %>% top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
        # left_join with va
        left_join(va %>% dplyr::select(name, display_name, gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
        left_join(va %>% dplyr::select(name, display_name) %>%
                      rename_with(~paste0(.x, ".drug")), by = c("drug" = "name.drug")) %>%
        dplyr::select(display_name.drug, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, drug, signature_vids) %>%
        purrr::set_names(c("drug_name","signature_name", "distance", "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "drug_node_id", "signature_node_id"))

    # 1.3 drug side effects
    drugs.se <- purrr::imap(drugs.ego,
                                 ~{vertex_attr(.x) %>% as.data.frame %>%
                                         filter(type == "side effect") %>%
                                         dplyr::select(name) %>%
                                         mutate(drug_seed = .y)}) %>%
        dplyr::bind_rows() %>%
        purrr::set_names(c("drug_side_effect", "name"))  %>%  # name = drug id
        left_join(va %>% dplyr::select(name, display_name) %>% rename_with(~paste0(.x, ".target")), by = c("drug_side_effect" = "name.target")) %>%
        left_join(va %>% dplyr::select(name, display_name) %>% rename_with(~paste0(.x, ".drug")), by = c("name" = "name.drug")) %>%
        dplyr::select(display_name.drug, name, display_name.target) %>%
        purrr::set_names("drug_name", "drug_id", "side_effect")

     # 2. degree signature
     degree_signature <- data.frame(degree = igraph::degree(graph = x$network, v = va %>% dplyr::filter(input_diffusion) %>% pull(name))) %>%
         tibble::rownames_to_column("node_id") %>% left_join(va, by = c("node_id" = "name")) %>%
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
         group_by(pathway) %>% top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
         # left_join with va
         left_join(va %>% dplyr::select(name, display_name,  gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
         left_join(va %>% dplyr::select(name, display_name) %>% rename_with(~paste0(.x, ".pathway")), by = c("pathway" = "name.pathway")) %>%
         dplyr::select(display_name.pathway, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, pathway, signature_vids) %>%
         purrr::set_names(c("pathway_name","signature_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "display_node_id", "signature_node_id"))

     # 4) GO terms
     #
     # 4.1 GO info
     GOs <- va %>% dplyr::filter(type == "GO") %>%
         dplyr::select(name, go_id, go_ontology, go_term_name)

     # 4.2 GO distance to signature
     #signature_in_gene <- va %>% filter(input_gene_signature) %>% pull(name)

     GO.dist.signature <- igraph::distances(graph = x$network, to = GOs$name,
                                                 v = va %>% dplyr::filter(input_diffusion) %>% pull(name)) %>%

     # signature.dist.GO <- igraph::distances(graph = x$network, to = signature_in_gene,
     # v = GOs$name) %>%
         as.data.frame()  %>%
         tibble::rownames_to_column("signature_vids") %>%
         tidyr::pivot_longer(names_to = "GO", values_to = "distance", -signature_vids) %>%
         group_by(GO) %>%
         # shortest distance
         top_n(wt = distance, n = -1) %>% dplyr::ungroup() %>%
         # left_join with va
         left_join(va %>% dplyr::select(name, display_name,  gene_hepatox_Toxicology2014, gene_hepatox_ToxicologyInVitro2020) %>% rename_with(~paste0(.x, ".sig")), by = c("signature_vids" = "name.sig")) %>%
         left_join(va %>% dplyr::select(name, display_name) %>% rename_with(~paste0(.x, ".go")), by = c("GO" = "name.go")) %>%
         dplyr::select(display_name.go, display_name.sig, distance,  gene_hepatox_Toxicology2014.sig, gene_hepatox_ToxicologyInVitro2020.sig, GO, signature_vids) %>%
         purrr::set_names(c("GO_term","signature_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "display_node_id", "signature_node_id"))




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

    # replace NA where logical
     to_return <- lapply(to_return, function(x) {
         x %>% mutate(across(where(is_logical), ~replace_na(data = .x, replace = FALSE)))
     })

     return(to_return)
}

#' Produce diffusion report
#'
#' Produce diffusion report based remarkable nodes in the network
#'
#' @param res_report
#' @param report_title
#' @param report_out_filepath
#' @param render = FALSE
#'
#' @examples
#' data(liver_1.3_network)
#' data(liver_1.3_rwr_closest_dfr)
#' data(signature_maison)
#' signature_vids <-  signature_maison$`acetaminophen_all_all`
#' res.diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, c(signature_vids, "DB00313"), target_type = c("drug/compound", "pathway", "side effect"))
#'
#' res_report <- report(res.diffusion)
#' report_info_df <- produce_diffusion_report(res_report = res_report,  report_title = "acetaminophen", report_out_filepath = "report_acetaminophen_example.Rmd",
#' render = TRUE, overwrite = TRUE)
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
                                    paste0(as.character(Sys.time()) %>% str_replace_all(" ", "-") %>% str_remove_all(":"),
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
```{r}
datatable(res_report$degree_signature)
```

# Drugs
## Drug info
```{r}
datatable(res_report$drugs)
```
## Drug target
```{r}
datatable(res_report$drugs.targets)
```
## Drug Side Effects
```{r}
datatable(res_report$drugs.side_effect)
```
## Drug shortest distance to signature
```{r}
datatable(res_report$drug.dist.signature)
```

# Pathways
## Pathway infos
```{r}
datatable(res_report$pathways)
```
## Pathways shortest distance to signature
```{r}
datatable(res_report$pathway.dist.signature)
```

# GO terms hepatox
## GO terms infos
```{r}
datatable(res_report$pathway.dist.signature)
```
## GO shortest distance to signature
```{r}
datatable(res_report$pathway.dist.signature)
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
#' produce report from rnaseq
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

