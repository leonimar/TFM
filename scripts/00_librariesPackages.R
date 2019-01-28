## Package installation

package_name <- "devtools"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

## First Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")

if (!requireNamespace("hipathia", quietly = TRUE))
  BiocManager::install("hipathia", version = "3.8")

package_name <- "limma"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "GEOquery"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "preprocessCore"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "affy"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "org.Hs.eg.db"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "genefilter"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "biomaRt"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "KEGGgraph"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

package_name <- "magrittr"
if (!requireNamespace(package_name, quietly = TRUE))
  BiocManager::install(package_name)

# otherlibraries
package_name <- "igraph"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "here"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "gplots"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "gdata"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "dplyr"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "Matching"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "xlsx"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "stringr"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "reshape2"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "ggplot2"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "XML"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

package_name <- "knitr"
if (!requireNamespace(package_name, quietly = TRUE))
  install.packages(package_name)

# devtools libraries
package_name <- "drugbankR"
if (!requireNamespace(package_name, quietly = TRUE))
  devtools::install_github("yduan004/drugbankR")

package_name <- "HPO.db"
if (!requireNamespace(package_name, quietly = TRUE))
  devtools::install_github("cran/HPO.db")

package_name <- "HPOSim"
if (!requireNamespace(package_name, quietly = TRUE))
  devtools::install_github("cran/HPOSim")


## Load of all libraries
library(AnnotationDbi)
library(magrittr) # pipe operator %>%
library(igraph)
library(KEGGgraph)
library(dplyr)
library(biomaRt)
library(gdata)
library(xlsx) 
library(here)
library(drugbankR)
library(stringr)
library("org.Hs.eg.db")
library(Matching)
library("GEOquery")
library("affy")
library("limma")
library("genefilter")
library("gplots")
library(preprocessCore)
library(hipathia)
library(here)
library(AnnotationDbi)
library(HPO.db)
library(HPOSim)
library(org.Hs.eg.db)
library(reshape2)
library(ggplot2)
library(XML)
library(knitr)
