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

#' Custom signature
#'
#' Custom gene signature from in house study.
#'
#' @format a list of char
#'
"signature_maison"


