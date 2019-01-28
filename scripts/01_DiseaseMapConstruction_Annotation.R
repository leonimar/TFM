
source(here("scripts/marina_funct.R"))

                  #### DISEASE MAP CONSTRUCTION ####

## INPUT: KEGG pathway ID of rare disease and extra genes and interactions we wish to add
  # !!!!!enter variable path_id <- "XXXX"
  # for Fanconi Anemia (FA) analysis: path_id <- "03460"  as in TEST_PIPELINE.R
  vertexfile = here("data", paste0(path_id, "_vertex.xls"))
  edgesfile = here("data", paste0(path_id, "_edges.xls"))

## Build up of the graph from KEGG pathway ID and add the extra genes and interactions
  # (extra_genes and extra_interaction = vertexfile and edgesfile respectively)
  # OUTPUT: vertex dataframe and edges dataframe that conform the graph of the Disease Map
  db <- kegg2db(path_id)
  lst_df <- kegg2df(db)

  datos_finales <- adding_ve(lst_df)
  vertex_final <- datos_finales$vertex_final
  edges_final <- datos_finales$edges_final
  vertex_extra <- datos_finales$vertex_extra
  edges_extra <- datos_finales$edges_extra
  vertex_final <- vertex_final[!duplicated(vertex_final$label),]
  rownames(vertex_final) <- vertex_final$kegg_id
  complete_ve <- adding_ve(lst_df)
  vertex_final <- complete_ve$vertex_final
  edges_final <- complete_ve$edges_final
  check_genes(edges_final, vertex_final)

  # Construct the final graph
  # For correct viewing click -> view -> fit to screen
  g <- convert2igraph(vertex_final, edges_final)
  tkplot(g)

  # # the current version needs all .sif and .att files in one place
  #TODO: Curation and automatization of the generation of .sif and .att files


                    ### DISEASE MAP ANNOTATION ####

## DRUGBANK ANNOTATION
  # Convert "hsa" codes of genes in Disease Map into UniProtKB codes to query them with DrugBankDB.
  entrez_ids <- strsplit(gsub("hsa: ","",vertex_final$inside),split = ";")
  entrez_ids <- lapply(entrez_ids, trimws, which = "both")
  names(entrez_ids) <- vertex_final$label
  myuni <- sapply(entrez_ids, conversor)
  names(myuni) <- vertex_final$label
  vertex_final$uniprot <- sapply(myuni, paste, collapse = "; ")

  # Retrieve Drugbank Data Base latest version and filter UniProtKB codes in targets colummn
  xml_fpath <- here("data", "drugbank.xml")
  drugbank_dataframe <- dbxml2df(xmlfile = xml_fpath, version = "5.1.1")
  drugbank_UniProt <- str_match(drugbank_dataframe$targets, "UniProtKB[A-Z]+[0-9]+") %>% gsub("UniProtKB", "", .)
  drugbank_dataframe$targets <- drugbank_UniProt

  # Find matched UniProtKB  codes  of genes in the Disease Map in Drugbank_dataframe$targets
  ignore <- NA
  HITS <- t(as.data.frame(lapply(myuni, function(x) {checking_targets(x, targets = drugbank_dataframe$targets)}))) %>% cbind(.,drugbank_dataframe$name[match(.,drugbank_dataframe$targets,  nomatch = NA_integer_ , incomparables = ignore )]) %>% as.data.frame(.)
      # output: matched items of ids in targets (when hit = TRUE) and store them as a list in HIT with the drug that corresponds to the target
    colnames(HITS)[1] <- "Uniprot_gene"
    colnames(HITS)[2] <- "Drug"
  HITS

## HPO ANNOTATION
  # Download and anotate Protein coding genes with HPO data base
  hsa_genes_fpath <- here("data", "genes.xls")
  Human_genes <- as.data.frame(read.csv(hsa_genes_fpath, sep = "\t"))
  hbo_db_fpath <- here("data", "HPO_OMIM-IVA.xls")
  HPO_DB <- as.data.frame(read.csv(hbo_db_fpath, sep = "\t"))
  Genes_HPO <- aggregate(cbind(as.character(HPO_DB$HPO.ID)) ~ HPO_DB$gene.symbol,
                         data = HPO_DB,
                         FUN = paste)
  names(Genes_HPO)[1] <- "gene_sym"
  names(Genes_HPO)[2] <- "HPO"
    #Annotation of Disease Map (DM) genes with existing HPO codes associated to each gene
  DM_genes <- vertex_final
  DM_genes_HPO <- t(as.data.frame(lapply(vertex_final$label, function(x) {
      checking_targets(x, targets = Genes_HPO$gene_sym)
  }))) %>%
      cbind(., Genes_HPO$HPO[match(., Genes_HPO$gene_sym,  nomatch = NA_integer_)]) %>%
      as.data.frame(.)
