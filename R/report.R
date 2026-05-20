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
#' @param complete_network igraph, same used in diffusion() to compute some distance metrics

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
report <- function(x, complete_network = NULL){

    stopifnot(is(x, "get_route.res"))
    if (!is.null(complete_network)) stopifnot(is(complete_network, "igraph"))

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
        purrr::set_names("drug_display_name", "drug_name", "target_display_name", "target_name")

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
        purrr::set_names(c("drug_display_name","input_display_name", "distance", "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "drug_name", "input_name"))

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
         purrr::set_names("name", "display_name", "degree", "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "coding_protein_is_present")


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
         purrr::set_names(c("pathway_display_name","input_display_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "pathway_name", "input_name"))

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
         purrr::set_names(c("GO_display_name","input_display_name", "distance",  "gene_hepatox_Toxicology2014", "gene_hepatox_ToxicologyInVitro2020", "GO_name", "input_name"))
     } else {
       GO.dist.signature <- data.frame()
}

     # most visited nodes — prefer pre-computed visit_count vertex attribute (from permutations),
     # fall back to counting node occurrences across shortest_path lists
     if ("visit_count" %in% colnames(va) && !all(is.na(va$visit_count))) {
         visite_autoroute <- va %>%
             dplyr::filter(!is.na(visit_count), visit_count > 0L) %>%
             dplyr::arrange(dplyr::desc(visit_count)) %>%
             dplyr::select(name, visit_count, type, display_name) %>%
             dplyr::rename(nb.visite = visit_count)
     } else {
         all_nodes <- x$autoroute$shortest_path %>%
             purrr::keep(~!purrr::is_empty(.x)) %>%
             unlist()
         visite_autoroute <- if (length(all_nodes) > 0L) {
             data.frame(name = all_nodes, stringsAsFactors = FALSE) %>%
                 dplyr::count(name, name = "nb.visite") %>%
                 dplyr::arrange(dplyr::desc(nb.visite)) %>%
                 dplyr::left_join(va %>% dplyr::select(name, type, display_name), by = "name")
         } else {
             data.frame()
         }
     }


     # 5. GO / pathway enrichment from permutations
     GO_enrichment <- NULL
     if (!is.null(x$permutation_results) && !is.null(x$observed_visit_counts)) {
         GO_enrichment <- tryCatch(
             compute_diffusion_enrichment(x, target_types = c("GO", "pathway"), min_visit_count = 1L),
             error = function(e) NULL
         )
     }

     ## update: enrichment in module (requires complete_network)
     # convert signature to uniprot
     enrich_input <- ""
     try(
     enrich_input <- AnnotationDbi::select(x = org.Hs.eg.db,
                                           keys = c(x$input, va %>% dplyr::pull(name)),
                                           columns = "UNIPROT", keytype = "ENSEMBL") %>%
         na.omit() %>% pull(UNIPROT), silent = TRUE)


     enrich_res <- module_ppi_M1_20231211 %>% mutate(is_present = molecule %in% enrich_input) %>%
         group_by(module) %>%
         summarise(k = sum(is_present), # nb of molecule in input in the module
                   n = n()) %>% # taille module
         mutate(K = sum(k),  # nb of molecule in all modules
                N = sum(n)) %>%   # size all modules
         nest(enrich_module_count = -module) %>%
         mutate(enrich_module_contingency = imap(enrich_module_count, ~make_contingency(k = .x$k, K = .x$K, n = .x$n, N = .x$N))) %>%
         mutate(enrich_module_p.value = map_dbl(enrich_module_contingency, ~fisher.test(.x, alternative = "greater")$p.value)) %>%
         mutate(enrich_module_p.value_adj = p.adjust(.$enrich_module_p.value, method = "fdr")) %>%
         dplyr::select(c(module, enrich_module_count, enrich_module_p.value, enrich_module_p.value_adj)) %>%
         unnest(cols = c(enrich_module_count)) %>%
         arrange(enrich_module_p.value)

     # add radar / ORA metrics — only when complete_network is provided
    if (is.null(complete_network)) {
        ORA_value        <- NULL
        radar_plots_data <- NULL
        enrich_res       <- data.frame()
    } else {

    va_complete <-  complete_network %>% vertex_attr() %>% as.data.frame()
    complete_edge_type <-  complete_network %>% igraph::as_long_data_frame() %>% dplyr::select(from_name, from_type, from_display_name, to_name, to_type, to_display_name)

    gotox_mapping <- list()
    gotox_mapping[["GO"]] <- va_complete %>% filter(type == "GO") %>% pull(name)
    gotox_mapping[["gene"]] <- complete_edge_type %>% filter(to_type == "GO", from_type == "gene") %>% pull(from_name) %>% unique()
    gotox_mapping[["protein"]] <-  complete_edge_type %>% filter(to_type == "GO", from_type == "protein") %>% pull(from_name) %>% unique()
    gotox_mapping[["pathway"]] <- complete_edge_type %>% filter(from_name %in% gotox_mapping[["protein"]], to_type == "pathway") %>% pull(to_name) %>% unique
    gotox_mapping[["drug"]] <- complete_edge_type %>% filter(from_name %in% gotox_mapping[["protein"]], to_type == "drug/compound") %>% pull(to_name) %>% unique
    gotox_mapping[["side_effect"]] <- complete_edge_type %>% filter(from_name %in% gotox_mapping[["drug"]], to_type == "side_effect") %>% pull(to_name) %>% unique

    liver_SE_mapping <- list()
    # liver SE
    liver_SE_mapping[["side_effect"]] <-  va_complete %>% filter(type == "side_effect") %>% filter(str_detect(name, "(?i)liver|hepa|steat"), !str_detect(name, "(?i)delivery")) %>% pull(name)
    # drugs connected to liver SE
    liver_SE_mapping[["drug"]] <- complete_edge_type %>% filter(to_name %in% liver_SE_mapping[["side_effect"]], from_type == "drug/compound") %>% pull(from_name) %>% unique
    # proteins connected to drugs connected to liver SE
    liver_SE_mapping[["protein"]] <- complete_edge_type %>% filter(to_name %in% liver_SE_mapping[["drug"]], from_type == "protein") %>% pull(from_name) %>% unique
    # genes connected to proteins connected to drugs connected to liver SE
    liver_SE_mapping[["gene"]] <-  complete_edge_type %>% filter(to_name %in% liver_SE_mapping[["protein"]], from_type == "gene") %>% pull(from_name) %>% unique
    # GOtox connected to proteins connected to drugs connected to liver SE
    liver_SE_mapping[["GO"]] <-  complete_edge_type %>% filter(from_name %in% liver_SE_mapping[["protein"]], to_type == "GO") %>% pull(to_name) %>% unique
    # pathway
    liver_SE_mapping[["pathway"]] <-  complete_edge_type %>% filter(from_name %in% liver_SE_mapping[["protein"]], to_type == "pathway") %>% pull(to_name) %>% unique

    # do ORA
    # get stats

    stat_enrich_interest <- tribble(
        ~name, ~ki, ~Ki, ~ni, ~Ni,
        # gene, ki = GOTOX gene in diff, Ki = GOTOX gene in network, ni = taille diff, Ni = taille network. ## ni and Ni could be length gene in diff/network
        "Gene_GOTox", intersect(gotox_mapping[["gene"]], va$name), va_complete %>% filter(type == "gene") %>% pull(name), va$name, va_complete$name,
        # proteins, ki = GOTOX prot in diff, Ki = GOTOX prot in network, ni = taille diff, Ni = taille network
        "Protein_GOTox",  intersect(gotox_mapping[["protein"]], va$name), va_complete %>% filter(type == "protein") %>% pull(name), va$name, va_complete$name,
        # pathway, pathway linked to GOTOX protein, pathway in networks
        "Pathway_GOTox",  intersect(gotox_mapping[["pathway"]], va$name), va_complete %>% filter(type == "pathway") %>% pull(name), va$name, va_complete$name,
        # drug, drug linked to GOTOX protein, drug in networks
        "Drug_GOTox", intersect(gotox_mapping[["drug"]], va$name), va_complete %>% filter(type == "drug/compound") %>% pull(name), va$name, va_complete$name,
        # side effect, SE linked to drugs linked to GOTOX protein, SE in networks
        "Side_effect_GOTox", intersect(gotox_mapping[["side_effect"]], va$name), va_complete %>% filter(type == "side_effect") %>% pull(name), va$name, va_complete$name,
        # GOTOX
        "GO", intersect(gotox_mapping[["GO"]], va$name),va_complete %>% filter(type == "GO") %>% pull(name), va$name, va_complete$name,

        ## SELiver
        # gene, gene connected to ... to SE_liver, gene in network
        "Gene_SELiver", intersect(liver_SE_mapping[["gene"]], va$name), va_complete %>% filter(type == "gene") %>% pull(name),  va$name, va_complete$name,
        # proteins, GOTOX protein in network, protein in network
        "Protein_SELiver",  intersect(liver_SE_mapping[["protein"]], va$name), va_complete %>% filter(type == "protein") %>% pull(name),  va$name, va_complete$name,
        # pathway, pathway linked to GOTOX protein, pathway in networks
        "Pathway_SELiver",  intersect(liver_SE_mapping[["pathway"]], va$name), va_complete %>% filter(type == "pathway") %>% pull(name), va$name, va_complete$name,
        # drug, drug linked to GOTOX protein, drug in networks
        "Drug_SELiver", intersect(liver_SE_mapping[["drug"]], va$name),va_complete %>% filter(type == "drug/compound") %>% pull(name), va$name, va_complete$name,
        # side effect, SE linked to drugs linked to GOTOX protein, SE in networks
        "Side_effect_SELiver", intersect(liver_SE_mapping[["side_effect"]], va$name),va_complete %>% filter(type == "side_effect") %>% pull(name), va$name, va_complete$name,
        # GO
        "GO_SELiver", intersect(liver_SE_mapping[["GO"]], va$name),va_complete %>% filter(type == "GO") %>% pull(name), va$name, va_complete$name,

        # Toxic pathways: pathway that touch dili drugs
        "drugs with Less DILI", filter(va, vDILIConcern == "vLess-DILI-Concern") %>% pull(name), filter(va_complete, vDILIConcern == "vLess-DILI-Concern") %>% pull(name), va$name, va_complete$name,
        "drugs with Ambigous DILI", filter(va, vDILIConcern == "Ambiguous DILI-concern") %>% pull(name), filter(va_complete, vDILIConcern == "Ambiguous DILI-concern") %>% pull(name), va$name, va_complete$name,
        "drugs with Most DILI", filter(va, vDILIConcern == "vMost-DILI-Concern") %>% pull(name), filter(va_complete, vDILIConcern == "vMost-DILI-Concern") %>% pull(name), va$name, va_complete$name,
        "drugs with DILI", filter(va, vDILIConcern %in% c("vMost-DILI-Concern", "Ambiguous DILI-concern", "vMost-DILI-Concern")) %>% pull(name), filter(va_complete, vDILIConcern %in% c("vMost-DILI-Concern", "Ambiguous DILI-concern", "vMost-DILI-Concern")) %>% pull(name), va$name, va_complete$name,
        "Pathways with DILI targets", intersect(va$name, filter(complete_edge_type, from_name %in% unique(filter(complete_edge_type, to_name %in% filter(va, vDILIConcern %in% c("vMost-DILI-Concern", "Ambiguous DILI-concern", "vMost-DILI-Concern"))$name)$from_name), to_type == "pathway") %>% pull(to_name) %>% unique), filter(complete_edge_type, from_name %in% unique(filter(complete_edge_type, to_name %in% filter(va, vDILIConcern %in% c("vMost-DILI-Concern", "Ambiguous DILI-concern", "vMost-DILI-Concern"))$name)$from_name), to_type == "pathway") %>% pull(to_name) %>% unique, va$name, va_complete$name
    )

    ORA_value <- stat_enrich_interest %>%
        #proportion
        mutate(k = map_dbl(ki, ~length(.x)),
               n = map_dbl(ni, ~length(.x)),
               K = map_dbl(Ki, ~length(.x)),
               N = map_dbl(Ni, ~length(.x))) %>%
        group_by(name) %>%
        nest(count = c(k, K, n, N)) %>%
        # ORA
        mutate(contingency = imap(count, ~make_contingency(k = .x$k, K = .x$K, n = .x$n, N = .x$N))) %>%
        mutate(p.value = map_dbl(contingency, ~fisher.test(.x, alternative = "greater")$p.value)) %>%
        mutate(padj = p.adjust(p.value, method = "fdr")) %>%
        mutate(ki_length = imap_dbl(ki, ~length(.x)), ni_length = imap_dbl(ni, ~length(.x)),
               Ki_length = imap_dbl(Ki, ~length(.x)), Ni_length = imap_dbl(Ni, ~length(.x))) %>%
        mutate(display_name = paste0(name, "\n", ki_length, "/", Ki_length))


    radar_plots_data <- list()
    radar_plots_data[["interest"]] <- ORA_value
    radar_plots_data[["module_all"]] <- enrich_res %>%
        mutate(name = module) %>%
        mutate(p.value = enrich_module_p.value) %>%
        mutate(padj = enrich_module_p.value_adj) %>%
        mutate(ki_length = k, ni_length = n, Ki_length = K, Ni_length = N) %>%
        mutate(display_name = paste0(name, "\n", ki_length, "/", Ki_length)) %>%
        arrange(module)

    module_enrich_gotox <- module_ppi_M1_20231211 %>% mutate(is_present = molecule %in% gotox_mapping[["protein"]]) %>%
        group_by(module) %>%
        summarise(k = sum(is_present), # nb of molecule in input in the module
                  n = n()) %>% # taille module
        mutate(K = sum(k),  # nb of molecule in all modules
               N = sum(n)) %>%   # size all modules
        nest(enrich_module_count = -module) %>%
        mutate(enrich_module_contingency = imap(enrich_module_count, ~make_contingency(k = .x$k, K = .x$K, n = .x$n, N = .x$N))) %>%
        mutate(enrich_module_p.value = map_dbl(enrich_module_contingency, ~fisher.test(.x, alternative = "greater")$p.value)) %>%
        mutate(enrich_module_p.value_adj = p.adjust(.$enrich_module_p.value, method = "fdr")) %>%
        dplyr::select(c(module, enrich_module_count, enrich_module_p.value, enrich_module_p.value_adj)) %>%
        unnest(cols = c(enrich_module_count)) %>%
        arrange(enrich_module_p.value) %>%
        filter(enrich_module_p.value_adj < 0.05) %>%
        pull(module)

    radar_plots_data[["module_GOTox"]] <- radar_plots_data[["module_all"]] %>%
        filter(name %in% c(module_enrich_gotox))

    } # end if (!is.null(complete_network))

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
     to_return[["GO_enrichment"]] <- GO_enrichment

    # replace NA where logical
     to_return <- lapply(to_return, function(x) {
         if (is.data.frame(x) && nrow(x) > 0L)
             x %>% mutate(dplyr::across(where(purrr::is_logical), ~tidyr::replace_na(data = .x, replace = FALSE)))
         else x
     })

     to_return[["ORA_value"]] <- ORA_value
     to_return[["radar_plots_data"]] <- radar_plots_data

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
    #stopifnot(all(unlist(lapply(res_report, function(x)is(x, "data.frame")))))
    #stopifnot(all(names(res_report) %in%  c('drug.dist.signature','pathway.dist.signature',
                                            # 'drugs.targets','degree_signature','drugs','pathways','drugs.side_effect','GOs','GO.dist.signature',
                                            # 'visite_autoraute', 'module_enrichment')))

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

```{r setup, echo = FALSE, message = FALSE, warning = FALSE}
knitr::opts_chunk$set( echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center')
library(DT)
library(tidyverse)
library(fmsb)
library(plotly)
```

")

    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             # value = paste0(
                                             #     "---\ntitle: '", report_title, "'\ndate: \"`r Sys.Date()`\"\noutput: html_document\n---\n\n",
                                             #     "```{r, echo = FALSE}\n knitr::opts_chunk$set( echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center')\n
                                             #     library(datatable)\n
                                             #      library(tidyverse)\n
                                             #      library(fmsb)\n
                                             #.     library(plotly)```"))
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

Columns description (columns with \\* are not shown in the Cytoscape file):

* name: Node identifier
* display_name: The node name displayed on the network
* \\* degree: The number of node connections
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

* drug_display_name *(column display_name)*: drug common name
* drug_name *(column name)*: drug node identifier
* target_display_name *(column display_name)*: drug target gene Symbol
* target_name *(column name)*: drug target node identifier

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

To refine the subnetwork, the table below returns the closest node in the diffusion input for each drug.

Columns description (columns with \\* are not shown in the Cytoscape file):

* drug_display_name *(display_name column)*: drug displayed name
* input_display_name *(display_name column)*: displayed input node for diffusion, closest node from the drug node.
* \\* distance: geodesic distance between the drug and input node.
* gene_hepatox_Toxicology2014: True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* drug_name *(name column)*: drug node identifier
* input_name *(name column)*: input node identifier for diffusion.

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

Columns description (column with \\* are not shown in the Cytoscape file):

* pathway_display_name *(display_name column)*: pathway displayed name
* input_display_name *(display_name column)*: displayed input node for diffusion, closest node from the pathway node.
* \\* distance: geodesic distance between the pathway and input node.
* gene_hepatox_Toxicology2014: True or False, if the closest node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the closest node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* GO_name *(name column)*: GO node identifier
* input_name *(name column)*: input node identifier for diffusion.

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

Columns description (columns with \\* are not shown in the Cytoscape file):

* GO_display_name *(display_name column)*: GO term description
* input_display_name *(display_name column)*: displayed input node for diffusion, closest node from the GO node.
* \\* distance: geodesic distance between the go_term node and input node.
* gene_hepatox_Toxicology2014: True or False, if the closest node is connected to a hepatotoxic function (GO terms) from the list given in Toxicology2014
* gene_hepatox_ToxicologyInVitro2020:  True or False, if the closest node is connected to a hepatotoxic function (GO terms) from the list given in ToxicologyInVitro2020
* GO_name *(name column)*: GO node identifier
* input_name *(name column)*: input node identifier for diffusion.

```{r}
datatable(res_report$GO.dist.signature)
```

# Most visited nodes

Columns description (columns with \\* are not shown in the Cytoscape file):

* name : node identifier
* \\* nb.visite : number of time
* type : type of node
* display_name : displayed node name in the network.

Each input seed generates a path to the final subnetwork, and some parts of the entire network may be visited multiple times. In this table, we list the most visited nodes during the diffusion process.
```{r}
datatable(res_report$visite_autoraute)
```

# Module enrichment

Enrichment analysis from diffusion input and diffusion results against modules.
Node names are converted to UNIPROT_ID to perform enrichment.

Columns description (columns with \\* are not shown in the Cytoscape file):

* \\* module: module identifier.
* \\* k: number of diffusion targets inside the module
* \\* n: size of the module
* \\* K: sum of diffusion targets inside all modules
* \\* N: size of all module
* \\* enrich_module_p.value: Hypergeometric test pvalue
* \\* enrich_module_p.value_adj: Corrected pvalue (fdr)

```{r}
datatable(res_report$module_enrichment)
```

```{r}
# if(!is(res_report$module_enrichment, 'try-errr')){
#   x <- res_report$module_enrichment
# # ONE version
# x.module_gotox <- x %>%  arrange(module) %>% mutate(pval = 1-enrich_module_p.value) %>%
#   dplyr::select(module, pval) %>% mutate(max = 1, min = 0) %>%
#   as.data.frame() %>% tibble::column_to_rownames('module') %>% t %>%
#   as.data.frame() %>%
#   .[c('max', 'min', 'pval'), sort(c('mod_M1_98', 'mod_M1_79', 'mod_M1_205', 'mod_M1_96', 'mod_M1_236', 'mod_M1_213', 'mod_M1_23', 'mod_M1_243', 'mod_M1_206', 'mod_M1_129', 'mod_M1_101', 'mod_M1_9', 'mod_M1_215', 'mod_M1_278', 'mod_M1_28', 'mod_M1_195', 'mod_M1_26', 'mod_M1_247'))]
#
# radarchart(x.module_gotox, title= 'Diffusion in GO tox enriched modules')
# }


try(res_report$radar_plots_data$interest %>%
    make_ggradar(max.upper.limit = 30, interactive = TRUE))

try(res_report$radar_plots_data$module_GOTox %>%
    make_ggradar(max.upper.limit = 30, interactive = TRUE))

try(res_report$radar_plots_data$module_all %>%
    make_ggradar(max.upper.limit = 30, interactive = TRUE))
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
                           n-k, N-K-(n-k)), ncol = 2, byrow = TRUE)

    colnames(dat) <- c("gene.in.interest", "gene.not.interest")
    rownames(dat) <- c("in.category", "not.in.category")
    return(dat)
}

#' @importFrom ggradar ggradar
make_ggradar <- function(x, max.upper.limit, mode_pvalue_highlight = c("-10log10pval", "1-pval"), interactive = FALSE){
    # make ggradar df
    max.upper.limit <- 30
    # mode_pvalue_highlight <- "-10log10pval" # c("-10log10pval", "1-pval"),
    mode_pvalue_highlight = match.arg(mode_pvalue_highlight)

    if(mode_pvalue_highlight == "-10log10pval"){
        ggradar_df <- x %>% mutate(signif = 0.05) %>%
            mutate(trans.pval =  -10*log10(p.value)) %>%
            mutate(max = max.upper.limit, min = 0) %>%
            mutate(trans.pval = pmin(trans.pval, max)) %>%
            mutate(trans.signif = -10*log10(signif))
    } else { # 1 - p.value
        ggradar_df <- x %>% mutate(signif = 0.05) %>%
            mutate(trans.pval =  (1 - p.value)) %>%
            mutate(max = 1, min = 0) %>%
            mutate(trans.signif = (1 - signif))
    }

    axis.label <- ggradar_df$display_name
    if(interactive){
        fig <- plot_ly( type = 'scatterpolar', mode = "markers") %>%
            add_trace(fill = 'toself', mode = "markers",
                r = ggradar_df$trans.pval,
                theta = ggradar_df$display_name,
                color = I("#00AFBB"),
                name = paste0("ORA ", ggradar_df$Ki_length[[1]], "/", ggradar_df$Ni_length[[1]])
            ) %>%
            add_trace(r = ggradar_df$trans.signif, mode = "lines",
                      theta = ggradar_df$display_name,
                      fill = "none",
                      line = list(smoothing = 1,
                                  color = "red",
                                  shape = "spline"),
                      hoverinfo = "skip",
                      showlegend = TRUE,
                      name = "Signif.") %>%
            layout(
                polar = list(
                    radialaxis = list(
                        visible = T,
                        range = c(-4,max.upper.limit)
                    )
                ),
                font=list("size" = 7)
            )
        return(fig)
    } else {
        p <- ggradar_df %>% dplyr::select(name, trans.pval) %>% column_to_rownames("name") %>% t %>%
            as.data.frame %>% rownames_to_column("group") %>%
            ggradar(fill = TRUE,
                    group.point.size = 2,
                    group.colours = "#00AFBB",
                    group.line.width = 1, grid.mid = 0.43,
                    values.radar = c(0, 13, 30),

                    gridline.min.linetype = "solid",
                    gridline.max.linetype = "solid",
                    gridline.mid.linetype = "solid",

                    axis.label.size = 3, base.size = 2, grid.label.size = 3,
                    gridline.mid.colour = "red",
                    axis.labels = axis.label,
                    background.circle.colour = "#D7D6D1",
                    gridline.label.offset = 0, )

        return(p)
    }
}

export_report_tables <- function(res_report, filepath){

    to_write = list(
        degree_signature = res_report$degree_signature,
        drugs = res_report$drugs,
        drugs.targets = res_report$drugs.targets,
        drugs.side_effect = res_report$drugs.side_effect,
        drug.dist.signature = res_report$drug.dist.signature,
        pathways = res_report$pathways,
        pathway.dist.signature = res_report$pathway.dist.signature,
        GOs = res_report$GOs,
        GO.dist.signature = res_report$GO.dist.signature,
        visite_autoraute = res_report$visite_autoraute,
        module_enrichment = res_report$module_enrichment)

   # to_write <- imap(to_write, ~.x %>% as.data.frame %>% rownames_to_column("item"))

    openxlsx::write.xlsx(to_write, file = filepath)
}


#' Produce a Quarto dashboard from CENetwork report results
#'
#' Generates a self-contained `.qmd` Quarto dashboard file from the output of
#' \code{report()} and \code{get_route()}.  Optionally renders it to HTML via
#' \code{quarto::quarto_render()}.
#'
#' @param res_report list returned by \code{report()}.
#' @param res_diffusion \code{get_route.res} returned by \code{get_route()}.
#' @param complete_network igraph; the full reference network used during
#'   diffusion (e.g. \code{liver_1.7_network}).
#' @param title character; title shown in the dashboard header.
#' @param output_filepath character; path to write the \code{.qmd} file.
#' @param se_to_severity data.frame mapping side-effect names to severity
#'   levels (columns \code{se} and \code{max_severity}).  Required.
#' @param render logical; if \code{TRUE} call \code{quarto::quarto_render()}
#'   after writing the file (requires the \pkg{quarto} package).
#' @param overwrite logical; if \code{FALSE} (default) an error is raised when
#'   \code{output_filepath} already exists.
#'
#' @return Invisibly returns the output file path.
#'
#' @examples
#' data(liver_1.7_network)
#' data(liver_1.7_rwr_closest_dfr)
#' data(signature_maison)
#' signature_vids <- signature_maison$`acetaminophen_all_all`
#' res_diffusion <- get_route(liver_1.7_network, liver_1.7_rwr_closest_dfr, signature_vids,
#'                            target_type    = c("drug/compound", "pathway", "GO", "side_effect"),
#'                            do_permutation = TRUE, n_permut = 100)
#' # optional: inspect enrichment results before generating the dashboard
#' # (the dashboard recomputes this internally)
#' enrich <- compute_diffusion_enrichment(res_diffusion,
#'                                        target_types = c("GO", "pathway"),
#'                                        min_visit_count = 1L)
#' head(enrich[enrich$padj < 0.05, ])
#' res_report <- report(res_diffusion, complete_network = liver_1.7_network)
#' produce_quarto_dashboard(
#'   res_report       = res_report,
#'   res_diffusion    = res_diffusion,
#'   complete_network = liver_1.7_network,
#'   title            = "Acetaminophen v1.7",
#'   output_filepath  = "dashboard_acetaminophen.qmd",
#'   se_to_severity   = se_to_severity,
#'   render           = TRUE,
#'   overwrite        = TRUE
#' )
#'
#' @importFrom checkmate assert_logical assertPathForOutput
#' @export
produce_quarto_dashboard <- function(res_report,
                                     res_diffusion,
                                     complete_network,
                                     title,
                                     output_filepath,
                                     se_to_severity,
                                     render = FALSE,
                                     overwrite = FALSE) {

    # Validate -----------------------------------------------------------------
    stopifnot(is(res_report, "list"))
    stopifnot(is(res_diffusion, "get_route.res"))
    stopifnot(is(complete_network, "igraph"))
    stopifnot(is(title, "character"), nchar(title) > 0)
    stopifnot(is(output_filepath, "character"))
    checkmate::assert_logical(render)
    checkmate::assert_logical(overwrite)
    checkmate::assertPathForOutput(output_filepath, overwrite = overwrite)

    # Save data objects to temp files ------------------------------------------
    ts <- format(Sys.time(), "%Y%m%d%H%M%S")

    path_res_report       <- file.path(tempdir(), paste0(ts, "_res_report.rds"))
    path_res_diffusion    <- file.path(tempdir(), paste0(ts, "_res_diffusion.rds"))
    path_complete_network <- file.path(tempdir(), paste0(ts, "_complete_network.rds"))
    path_se_to_severity   <- file.path(tempdir(), paste0(ts, "_se_to_severity.rds"))

    saveRDS(res_report,       path_res_report)
    saveRDS(res_diffusion,    path_res_diffusion)
    saveRDS(complete_network, path_complete_network)
    saveRDS(se_to_severity,   path_se_to_severity)

    # Read template and substitute placeholders --------------------------------
    # template_path <- system.file("templates", "dashboard_template.qmd",
    #                              package = "CENetwork")
    # temporairement
    template_path <- "/Users/antoine/Documents/DEV/CENetwork/inst/templates/dashboard_template.qmd"
    if (!nzchar(template_path)) {
        stop("Dashboard template not found. Is CENetwork installed / loaded with devtools::load_all()?")
    }

    qmd <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

    replacements <- list(
        "{{DASHBOARD_TITLE}}"      = title,
        "{{PATH_RES_REPORT}}"      = path_res_report,
        "{{PATH_RES_DIFFUSION}}"   = path_res_diffusion,
        "{{PATH_COMPLETE_NETWORK}}" = path_complete_network,
        "{{PATH_SE_TO_SEVERITY}}"  = path_se_to_severity
    )

    for (token in names(replacements)) {
        qmd <- gsub(token, replacements[[token]], qmd, fixed = TRUE)
    }

    # Write QMD ----------------------------------------------------------------
    writeLines(qmd, output_filepath)
    message("Dashboard written to: ", output_filepath)

    # Optionally render --------------------------------------------------------
    if (render) {
        if (!requireNamespace("quarto", quietly = TRUE)) {
            warning("Package 'quarto' not available; skipping render.")
        } else {
            quarto::quarto_render(output_filepath)
        }
    }

    invisible(output_filepath)
}
