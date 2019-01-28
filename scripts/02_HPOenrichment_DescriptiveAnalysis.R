
                     #### HPO ENRICHMENT ANALYSES ####

  # Enrichment of HPO terms in Disease Map list of genes
FA_HPO_enrich <- HPOGeneEnrichment(vertex_final$entrzIDs, filter = 1, cutoff = 0.05, background = getGeneDefaultBackground())
FA_HPO_enrich50 <- FA_HPO_enrich[FA_HPO_enrich$qvalue < 1e-13, ]
HPO_Enriched_FA <- getTerm(as.character(FA_HPO_enrich50$HPOID))
x = HPO_Enriched_FA[[1]]
FA_HPO_enrich50$Term <- as.vector(sapply(HPO_Enriched_FA, function(x) x@Term))
out_fpath <- here("OUTPUT", "FA_HPO_enrich50.csv")
write.csv(FA_HPO_enrich50, file = out_fpath)


  # Enrichment of HPO terms in a Disease Map interaction-network of genes
edges1 <- gsub("hsa:", "", edges_final$V1)
edges2 <- gsub("hsa:", "", edges_final$V2)
gene_network <- cbind(edges1, edges2)
gene_network_fpath <- here("OUTPUT", "gene_network.csv")
write.csv(gene_network, file = gene_network_fpath)
FA_Network_HPOenrich <- HPOGeneNOAWholeNetEnrichment(gene_network_fpath, cutoff = 0.005, filter = 5)
Terms_Network <- getTerm(as.character(FA_Network_HPOenrich$HPOID))
y = Terms_Network[[1]]
FA_Network_HPOenrich_terms <- as.data.frame(sapply(Terms_Network, function(y) y@Term))
fpath <- here("OUTPUT", "FA_Network_HPOenrich")
write.csv(FA_Network_HPOenrich, file = fpath)
fpath <- here("OUTPUT", "FA_Network_HPOenrich_Terms")
write.csv(FA_Network_HPOenrich_terms, file = fpath)


                          #### DESCRIPTIVE ANALYSIS ####

 # From table of all genes in Rare Disease retrieved from ORPHANET and HiPathia genes
fpath <- here("data", "RD_ORPHA-pathways.xlsx")
RD_ORPHA <- read.xlsx(file = fpath, 1)
RD_ORPHA <- subset(RD_ORPHA, disease_genes_in_pathways > 3)
rownames(RD_ORPHA) <- RD_ORPHA$ORPHA_disease_term_name
RD_ORPHA <- RD_ORPHA[, -2]

 # Stacked graph comparing total RD genes versus total genes in HiPathia
 # For correct viewing click --> Zoom
graph_DF <- data.frame(t(RD_ORPHA[, -c(1,4)]), check.names = F)
graph_DF <- graph_DF[, order(apply(graph_DF,2,sum), decreasing = F)]
graph_DF$genetype <- rownames(graph_DF)
graph_DF <- melt(graph_DF, id.vars = "genetype")
p <- ggplot(graph_DF, aes(x = variable, y = value, fill = genetype,  cex.axis = 15, cex.lab = 15)) +
  geom_bar(stat = "identity") +
  xlab("\nDiseases") +
  ylab("Number of genes") +
  guides(fill = FALSE)
theme_bw()

p <- p +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = rel(2)),
          axis.text.y = element_text(size = rel(2))) +
    coord_flip() +
    guides(fill = guide_legend(), size = rel(4))
plot(p)
