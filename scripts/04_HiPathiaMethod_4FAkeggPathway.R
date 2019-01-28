
##### HIPATHIA METHOD FOR FA KEGG PATHWAY ####

 ## PREPROCESSMENT: Data scaling and normalization 
  exp_data <- normalize_data(gExp)
  targets$NiceClasses <- targets$Classes
  targets$NiceClasses[grep("Normal",targets$Classes,ignore.case = T)] <- "Normal"
  targets$NiceClasses[grep("FAaplas",targets$Classes,ignore.case = T)] <- "Fanconi"
  table(targets$NiceClasses)
  boxplot(gExp)
  boxplot(exp_data)
  exp_data <- normalize_data(gExp, by_quantiles = TRUE)
  boxplot(exp_data)
  exp_data <- normalize_data(gExp, percentil = TRUE)
  boxplot(exp_data)
  exp_data <- normalize_data(gExp, truncation_percentil = 0.95)
  boxplot(exp_data)

 ## Loading Pathways we are interested in checking
  pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03460"))

 ## Using Hipathia to compute the signal (subpathways)
  results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)
  results
  
  path_vals <- get_paths_data(results, matrix = TRUE)
  
  
  tab <- t(sapply(c("hsa", "mmu", "rno"), function(species){
    p <- suppressMessages(load_pathways(species))
    effs <- sum(sapply(p$pathigraphs, function(pathi) length(
      pathi$effector.subgraphs)))
    decs <- sum(sapply(p$pathigraphs, function(pathi) length(pathi$subgraphs)))
    n <- length(p$pathigraphs)
    c(n, effs, decs)
  }))
  colnames(tab) <- c("Pathways", "Effector subpathways", "Decomposed subpathways")
  knitr::kable(tab)

 ## Function activation computation UniPort Terms and GO terms
  uniprot_vals <- quantify_terms(results, pathways, dbannot = "uniprot")
  go_vals <- quantify_terms(results, pathways, dbannot = "GO")
  sample_group <- targets$NiceClasses
  comp <- do_wilcoxon(path_vals, sample_group, g1 = "Fanconi", g2 = "Normal")
  top_genes<- comp[comp$FDRp.value<=0.05,]
  hhead(comp)
  pathways_summary <- get_pathways_summary(comp, pathways)
  

 ## PCA analysis
  ranked_path_vals <- path_vals[order(comp$p.value, decreasing = FALSE),]
  if(ncol(ranked_path_vals)>nrow(ranked_path_vals)){ 
    pca_model <- do_pca(ranked_path_vals)
  }else{
    pca_model <- do_pca(ranked_path_vals[1:ncol(ranked_path_vals),])
  }

 ## Heat Map
  heatmap_plot(path_vals, group = sample_group)
  
  heatmap_plot(uniprot_vals, group = sample_group, colors="hipathia", 
               variable_clust = TRUE)
  heatmap_plot(go_vals, group = sample_group, colors="redgreen", 
               variable_clust = TRUE)
  
 ##PCA visualization
  pca_plot(pca_model, sample_group, legend = TRUE)
  
  pca_plot(pca_model, group = rep(1:5, 8), main = "Random types", 
           legend = TRUE)
  multiple_pca_plot(pca_model, sample_group, cex=3, plot_variance = TRUE)

 ## Pathway comparision
  pathway_comparison_plot(comp, metaginfo = pathways, pathway = "hsa03460")
  colors_de <- node_color_per_de(results, pathways, sample_group, "Fanconi", 
                                 "Normal")
  pathway_comparison_plot(comp, metaginfo = pathways, pathway = "hsa03460", 
                          node_colors = colors_de)
  colors_de_hipathia <- node_color_per_de(results, pathways, sample_group, 
                                          "Fanconi", "Normal", colors = "hipathia")
  pathway_comparison_plot(comp, metaginfo = pathways, pathway = "hsa03460", 
                          node_colors = colors_de_hipathia, colors = "hipathia")