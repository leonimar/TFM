
#### DATA PREPARATION FOR HIPATHIA VALIDATION ####

## 1.NORMALIZATION OF AFFIMETRIX EXPRESSION SET FROM GEO

  # Download gene expression set files (GSEA16334 for Fanconi Anemia gene expression set)
  #df = getGEOSuppFiles('GSEA16334'). NOT RUNNING.
  # We have downloaded the files from GEO web interface and stored in /data/raw

  # Import targets.txt file
fa_targets_fpath <- here("data", "raw", "FA_targets.txt")
targets <- readTargets(fa_targets_fpath, row.names = "FileName")

  # Import .CEL files
data <- ReadAffy(filenames = targets$FileName,
                 celfile.path=fa_targets_fpath <- here("data", "raw"))
# data is an object of AffyBatch class

  # Normalize with RMA
  #generates object eset (class ExprSet),
  #expresso function provides intensities in log scale
eset <- expresso(data,
                 bg.correct = TRUE,
                 bgcorrect.method = "rma",
                 normalize = TRUE,
                 normalize.method = "quantiles",
                 pmcorrect.method = "pmonly",
                 summary.method = "medianpolish",
                 verbose = TRUE
)

  # Generate BOXPLOTS before and after normalization

  #boxplot for raw data
boxplot(data,
        main = "Boxplot Before Normalization",
        col = "lightgrey")


  #boxplot for normalized data
exprseset <- as.data.frame(exprs(eset))
boxplot(data.frame(exprseset),
        main = "Boxplot After Normalization (log scale)",
        col = "white")

  # Data filtering using IQR.
esetIQR <- varFilter(eset, var.func = IQR, var.cutoff = 0.5, filterByQuantile = TRUE)
exprsesetIQR <- as.data.frame(exprs(esetIQR))

## 2. PROBES OBTAINED AFTER NORMALIZATION OF GEO DATASETS TO GENE TRANFORMATION

  # Disease genes extraction and translation of probes to genes for HiPathia method
fanconigenes_fpath <- here("data", "fanconigenes.txt")
write.table(unique(gsub("hsa:", "", vertex_final$kegg_id)),
            file = fanconigenes_fpath,
            quote = F,
            row.names = F)

affyNorm <- as.matrix(exprsesetIQR)
export_fpath <- here("data", "raw", "mart_export.txt")
biomart <- read.delim(export_fpath,
                      sep = "\t",
                      stringsAsFactors = F)
head(biomart)
biomart <- biomart[!is.na(biomart$NCBI.gene.ID), ]
biomart <- biomart[which(biomart$AFFY.HG.U133A.probe %in% rownames(affyNorm)), ]
NCBIgenes <- unique(biomart$NCBI.gene.ID)

gExp <- mat.or.vec(nr=length(NCBIgenes),nc = ncol(affyNorm))
gExp[, ] <- -99999
multipleprobe <- names(which(table(biomart$AFFY.HG.U133A.probe) > 1))

for (n in 1:length(NCBIgenes)) {
  ng <- NCBIgenes[n]
  probes <- unique(biomart$AFFY.HG.U133A.probe[biomart$NCBI.gene.ID == ng])

  # if there is 1 probe for 1 gene, the probe expression is assigned to gene
  if (length(probes) == 1) {
    myexp <- affyNorm[match(probes,rownames(affyNorm)), ]
    gExp[n, ] <- myexp
  }

  # discard probes which are shared between different genes (disProbe).
  # when you discard them if there is no any probe for a gene keep
  # all probes and use median of their expression as gene expression.
  if (length(probes) > 1) {
    idx_mprobe <- match(probes, multipleprobe)
    numberofprobes <- length(idx_mprobe[is.na(idx_mprobe)])

  # After discarding disProbes if 1 specific probe remains use its expression for gene
  # if many probe remains use percentile 90 of probe expressions as gene expression
    if (numberofprobes >= 1) {
      probes <- probes[is.na(idx_mprobe)]
      if (length(probes) == 1) {
        myexp <- affyNorm[match(probes,rownames(affyNorm)), ]
      } else {
        myexp <- affyNorm[match(probes,rownames(affyNorm)), ]
        myexp <- apply(myexp, 2, function(x) quantile(x, probs = 0.9))
      }
    }else{
      # if all probes of a gene shared with other genes then median of probes' expressions are assigned to gene
      myexp <- affyNorm[match(probes, rownames(affyNorm)), ]
      myexp <- apply(myexp, 2, function(x) median(x))
    }
    gExp[n, ] <- myexp
  }
  if ((n %% 1000) == 0)
    print(paste0(n, "... genes calculated"))
}
dimnames(gExp) <- list(NCBIgenes, colnames(affyNorm))

par(bg = 'white', col = "black")
boxplot(gExp)

  ## gene espression dataset prepared for HiPathia
gExp <- normalize.quantiles(gExp)
dimnames(gExp) <- list(NCBIgenes, colnames(affyNorm))
boxplot(gExp)

