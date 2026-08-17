#' Liver specific multi-layer network v1.3
#'
#' The liver network is composed of 6 layers (gene, protein, drug/compound, pathway, side effect and Hepatox GO terms).
#'
#' First, the protein-protein interaction network layout wos build based on BioGRID (https://thebiogrid.org) and only proteins expressed in liver were kept (https://www.proteinatlas.org/humanproteome/tissue/liver).
#'
#' Proteins were connected to an in-house gene coregulation network (ARACNE https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7) with gene to protein coding information and Transctipted Factor (TF) to Targeted Genes (TG) with TF2DNA (https://www.fiserlab.org/tf2dna_db/), TRRUST https://www.grnpedia.org/trrust/ and Dorothea (https://saezlab.github.io/dorothea/). Only interactions between the present genes and proteins were included.
#'
#' Drugs were extracted from DrugBank and were linked to their protein targets (https://drugbank.ca).
#' We used CHEMBL to add IC50 information when available for the HepG2 cell line and if the compound and its targets were present in the network.
#' Side Effects were extracted from SIDER and were linked to drugs (http://sideeffects.embl.de).
#'
#' Proteins were also linked to Reactome pathways (https://reactome.org).
#'
#' Finaly, GO terms linked to hepato-toxicicity were connected to gene and protein.
#'
#' @format an igraph object
#'
"liver_1.3_network"


#' Closest Targets dfr
#'
#' Closest target from each node (each node is a seed), and for each layer.
#'
#' @format data.frame with 4 columns: Seed, Target, Seed.type, Target.type
#'
"liver_1.3_rwr_closest_dfr"

#' Liver specific multi-layer network v1.7
#'
#' Updated version of the liver network (v1.7), composed of 6 layers
#' (gene, protein, drug/compound, pathway, side effect and hepatotoxicity GO terms).
#' Compared to v1.3, this version extends the drug layer (DrugBank update),
#' refines the GO term list, and includes pre-computed shortest paths
#' across all seed–target pairs (see \code{liver_1.7_rwr_closest_dfr}).
#'
#' @format an igraph object
#'
"liver_1.7_network"


#' Pre-calculated diffusion scores v1.7
#'
#' Pre-calculated closest neighbours from each seed node to every target layer,
#' matching \code{liver_1.7_network}. The \code{shortest_path} column stores
#' the ordered list of node names along each seed-to-target shortest path,
#' enabling fast subnetwork extraction in \code{get_route()} without runtime
#' igraph calls.
#'
#' @format a data.frame with columns: SeedNode, type.seed, NodeNames,
#'   type.target, distance_from_start, shortest_path
#'
"liver_1.7_rwr_closest_dfr"


#' Liver specific multi-layer network v1.8
#'
#' Current reference version of the liver network. Identical to v1.7 except for
#' the pathway layer, pruned from 2018 to 814 Reactome pathways by keeping only
#' those mapping to hepatotoxicity key characteristics (Tsai et al. 2025).
#' Pruning removes the generic root pathways (Metabolism, Signal Transduction, …)
#' that acted as diffusion hubs and dominated the enrichment results.
#'
#' @format an igraph object
#'
"liver_1.8_network"


#' Pre-calculated diffusion scores v1.8
#'
#' Pre-calculated closest neighbours from each seed node to every target layer,
#' matching \code{liver_1.8_network}. Same format as
#' \code{liver_1.7_rwr_closest_dfr}.
#'
#' @format a data.frame with columns: SeedNode, type.seed, NodeNames,
#'   type.target, shortest_path
#'
"liver_1.8_rwr_closest_dfr"


#' Custom signature
#'
#' Custom gene signature from in house study.
#'
#' @format a list of char
#'
"signature_maison"


#' Module
#'
#' Precalulated PPI module using M1 olgo from MONET
#'
#' @format a data.frame
#'
"module_ppi_M1_20231211"


#' ORA and GSEA comparison results for APAP and VPA
#'
#' Pre-computed enrichment results for acetaminophen (APAP) and valproic acid
#' (VPA) from three methods — network diffusion, ORA, and GSEA — all filtered
#' to the hepatotoxicity-specific background of \code{liver_1.8_network}.
#' Intended for use in the vignette comparison section.
#'
#' ORA was run with \code{clusterProfiler::enrichGO()} /
#' \code{ReactomePA::enrichPathway()}, GSEA with \code{gseGO()} /
#' \code{gsePathway()}, both with \code{pAdjustMethod = "none"} and
#' \code{pvalueCutoff = 1}, then re-filtered to the network background and
#' FDR-adjusted (\code{qvalue < 0.05}).
#'
#' @format A named list with 12 elements:
#' \describe{
#'   \item{diffusion_go_apap}{Character vector of significant GO term IDs (APAP)}
#'   \item{diffusion_go_vpa}{Character vector of significant GO term IDs (VPA)}
#'   \item{diffusion_pathway_apap}{Character vector of significant Reactome pathway IDs (APAP)}
#'   \item{diffusion_pathway_vpa}{Character vector of significant Reactome pathway IDs (VPA)}
#'   \item{ora_go_apap}{data.frame of significant ORA GO results (APAP)}
#'   \item{ora_go_vpa}{data.frame of significant ORA GO results (VPA)}
#'   \item{gsea_go_apap}{data.frame of significant GSEA GO results (APAP)}
#'   \item{gsea_go_vpa}{data.frame of significant GSEA GO results (VPA)}
#'   \item{ora_reactome_apap}{data.frame of significant ORA Reactome results (APAP)}
#'   \item{ora_reactome_vpa}{data.frame of significant ORA Reactome results (VPA)}
#'   \item{gsea_reactome_apap}{data.frame of significant GSEA Reactome results (APAP)}
#'   \item{gsea_reactome_vpa}{data.frame of significant GSEA Reactome results (VPA)}
#' }
#'
"comparison_ora_gsea"

