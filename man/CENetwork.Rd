% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CENetwork-package.R
\docType{package}
\name{CENetwork}
\alias{CENetwork}
\title{CENetwork: Network exploration}
\description{
This package provides a multi-omics network and allows diffusion analysis from a gene list of nodes.
}
\details{
The main diffusion method is based on the Random Walk with Restart principle on a multi-layer network.
We have modified the algorithm to force the passage through distant layers of the network (user defined parameter).
From each element of a list of nodes (seeds), we calculate the RWR scores, and then identify the closest targets for each of the selected layers. Finally, we extract the shortest paths between the seed and these targets and return a sub-network with the intermediate nodes.

This subnetwork is exportable to Cytoscape.

To help the interpretation, we build automated reports to highlight notable features (high connected nodes with hepatotoxic functions, distance drug - gene signature, distance pathways signature, ...)
}
\section{Diffusion}{

\describe{
  \item{`get_route`}{get diffusion result from a list of seeds}
  }
}

\section{Export to Cytoscape}{

\describe{
  \item{`export_to_cytoscape`}{Export to cytoscape (need to be open)}
  \item{`apply_custom_theme`}{Apply custom theme such as node file color, border color, ... on cytoscape graph}
 }
}

\section{Diffusion reports}{

\describe{
  \item{`report`}{Generate report tables}
  \item{`produce_diffusion_report`}{Produce html report from `report()` results.`}
 }
}

