
        #### CREATION OF NEW METAGINFO (INCLUDING DISEASE MAP)####

  # Adapted from hpAnnot pre Bioconductor version

library(igraph)

source(here("library", "graphs.R"))
source(here("library", "layout.R"))

data_folder <- here("data")
sif_folder <- file.path(data_folder, "sif")
tmp_folder <- file.path(data_folder, "tmp")
dir.create(tmp_folder, showWarnings = FALSE)

ammend_file <- file.path(data_folder, "sif_amendments.txt")

spe = "hsa"

pathway_names <- unique(gsub(".sif", "", list.files(sif_folder,
                                                    pattern="sif")))
titles <- pathway_names
titles_df <- cbind(gsub(spe, "", titles), titles)

# Load pathways from created SIF files
pgs <- load.graphs(sif_folder, spe, pathway_names)
save(pgs, file = file.path(tmp_folder, "pgs.RData"))

# Ammend pathways
apgs <- amend.kegg.pathways(ammend_file, pgs, spe)
save(apgs, file = file.path(tmp_folder, "apgs.RData"))

# Add final functions to the pathways

# Load annotations
dbannot <- hipathia:::load_annots("uniprot", spe)
entrez2hgnc <- hipathia:::load_entrez_hgnc(spe)

fpgs <- add.functions.to.pathigraphs(apgs, entrez2hgnc, dbannot,
                                     maxiter = 1000)
save(fpgs, file = file.path(tmp_folder, "fpgs.RData"))

# Compute Path Normalization Values
metaginfo <- create.metaginfo.object(fpgs, spe)
save(metaginfo,
     file = file.path(tmp_folder, paste0("meta_graph_info_", spe, ".RData")))
