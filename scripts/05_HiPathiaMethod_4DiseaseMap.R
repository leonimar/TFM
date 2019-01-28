
          ##### HIPATHIA METHOD FOR DISEASES MAP (EXTENDED PATHWAY CREATED) ####

 ### 1. NEW METAGINFO CREATION ###

## Creation of new metaginfo including the Disease Map (DM) (New pathway not contemplated in actual metaginfo of HiPathia
## metaginfo containing--> .sif and .att files created after Disease Map construction.
## metaginfo contains the information about the pathwaysthat HiPathia needs to compute the activation values

    source(here("scripts", "metaginfo_create_from_sif.R"))


 ### 2. HIPATHIA METHOD ###

  ## PREPROCESSMENT: Data scaling and normalization
    exp_data_FAmarina <- normalize_data(gExp)
    targets$NiceClasses <- targets$Classes
    targets$NiceClasses[grep("Normal", targets$Classes,ignore.case = T)] <- "Normal"
    targets$NiceClasses[grep("FAaplas", targets$Classes,ignore.case = T)] <- "Fanconi"
    table(targets$NiceClasses)
    boxplot(gExp)
    exp_data_FAmarina <- normalize_data(gExp, by_quantiles = TRUE)
    boxplot(exp_data_FAmarina)

 ## Loading Pathways
    pathways_FAmarina <-  metaginfo

 ## Using Hipathia to compute the signal (subpathways)
    results_FAmarina <- hipathia(exp_data_FAmarina, pathways_FAmarina, decompose = FALSE, verbose = FALSE)
    results_FAmarina

    path_vals_FAmarina <- get_paths_data(results_FAmarina, matrix = TRUE)

    tab <- t(sapply(c("hsa", "mmu", "rno"), function(species){
      p <- suppressMessages(metaginfo)
      effs <- sum(sapply(p$pathigraphs, function(pathi) length(
        pathi$effector.subgraphs)))
      decs <- sum(sapply(p$pathigraphs, function(pathi) length(pathi$subgraphs)))
      n <- length(p$pathigraphs)
      c(n, effs, decs)
    }))
    colnames(tab) <- c("Pathways", "Effector subpathways", "Decomposed subpathways")
    knitr::kable(tab)

 ## Function activation computation with UniPort Terms and GO terms
    uniprot_vals_FAmarina <- quantify_terms(results_FAmarina, pathways_FAmarina, dbannot = "uniprot")
    go_vals_FAmarina <- quantify_terms(results_FAmarina, pathways_FAmarina, dbannot = "GO")
    sample_group_FAmarina <- targets$NiceClasses
    comp_FAmarina <- do_wilcoxon(path_vals_FAmarina, sample_group_FAmarina, g1 = "Fanconi", g2 = "Normal")
    top_genes_FAmarina <- comp_FAmarina[comp_FAmarina$FDRp.value <= 0.05,]

    pathways_summary_FAmarina <- get_pathways_summary(comp_FAmarina, pathways_FAmarina)
    head(pathways_summary_FAmarina, 4)

 ## PCA analysis
    ranked_path_vals_FAmarina <- path_vals_FAmarina[order(comp_FAmarina$p.value, decreasing = FALSE), ]
    pca_model_FAmarina <- do_pca(ranked_path_vals_FAmarina[1:ncol(ranked_path_vals_FAmarina), ])

 ## Heat Map
    heatmap_plot(path_vals_FAmarina, group = sample_group_FAmarina)

    heatmap_plot(uniprot_vals_FAmarina, group = sample_group_FAmarina, colors = "hipathia",
                 variable_clust = TRUE)
    heatmap_plot(go_vals_FAmarina, group = sample_group_FAmarina, colors = "redgreen",
                 variable_clust = TRUE)

 ##PCA visualization
    pca_plot(pca_model_FAmarina, sample_group_FAmarina, legend = TRUE)

    pca_plot(pca_model_FAmarina, group = rep(1:5, 8), main = "Random types",
             legend = TRUE)
    multiple_pca_plot(pca_model_FAmarina, sample_group_FAmarina, cex = 3, plot_variance = TRUE)

 ## Pathway comparision
    pathway_comparison_plot(comp_FAmarina, metaginfo = pathways_FAmarina, pathway = "hsa03460_marina")
    colors_de_FAmarina <- node_color_per_de(results_FAmarina, pathways_FAmarina, sample_group_FAmarina, "Fanconi",
                                            "Normal")
    pathway_comparison_plot(comp_FAmarina, metaginfo = pathways_FAmarina, pathway = "hsa03460_marina",
                            node_colors = colors_de_FAmarina)
    colors_de_hipathia_FAmarina <- node_color_per_de(results_FAmarina, pathways_FAmarina, sample_group_FAmarina,
                                                     "Fanconi", "Normal", colors = "hipathia")
    pathway_comparison_plot(comp_FAmarina, metaginfo = pathways_FAmarina, pathway = "hsa03460_marina",
                            node_colors = colors_de_hipathia_FAmarina, colors = "hipathia")

