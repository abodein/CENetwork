# Package: Cosmetic Europe Network exploration

This package provides a multi-omics network and proposes the exploration of the network by propagation from a gene list (signature).

## informations about the network

## Install

```r 
# install devtools
install.packages("devtools")
# install the package (last version)
devtools::install_github("abodein/CENetwork")

# load the package
library(CENetwork)
```
## Install Cytoscape

Install the latest version of Cytoscape: https://cytoscape.org/download.html

Cytoscape must be open during the exportation process.
Please check your connection via

```r
RCy3::cytoscapePing()
```
