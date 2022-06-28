## ----setup, echo = FALSE, warning=FALSE, message=FALSE------------------------
knitr::opts_chunk$set(fig.align = "center", message = FALSE, warning = FALSE)
library(tidyverse)
library(igraph)

## ----install, eval = FALSE----------------------------------------------------
#  # install devtools
#  install.packages("devtools")
#  # install the package (last version)
#  devtools::install_github("abodein/CENetwork")

## ----load, eval = TRUE--------------------------------------------------------
# load the package
library(CENetwork)

## ---- eval = FALSE, echo=TRUE-------------------------------------------------
#  # return an error if Cytoscape can not be found
#  RCy3::cytoscapePing()

## ----load_data----------------------------------------------------------------
data(liver_1.3_network)  # the network
data(liver_1.3_rwr_closest_dfr) # pre-calculated diffusion scores
data("signature_maison") # custom gene signature
signature_acetaminophen <- signature_maison$`acetaminophen_all_all`
signature_acetaminophen

## ----getroute-----------------------------------------------------------------
# load data RWR scores
diffusion.res <- get_route(network = liver_1.3_network, # the network
            closest_dfr = liver_1.3_rwr_closest_dfr, # pre-calculated diffusion scores
            signature_vids = signature_acetaminophen, # input seeds
            target_type = c("drug/compound", "pathway", "side_effect")) # layers to reached

## ----cyto, eval = FALSE-------------------------------------------------------
#  # export to cytoscape
#  export_to_cytoscape(diffusion.res)
#  apply_custom_theme()

## ----report, eval = FALSE-----------------------------------------------------
#  # generate report
#  diffusion.report <- report(diffusion.res) # result from get_route()
#  
#  # produce the table
#  produce_diffusion_report(res_report = diffusion.report, # result from report()
#                           report_title = "Acetaminophen",
#                           report_out_filepath = "Acetaminophen_report.Rmd", # Rmd file
#                           overwrite = TRUE, # overwrite report_out_filepath
#                           render = TRUE) # produce html

