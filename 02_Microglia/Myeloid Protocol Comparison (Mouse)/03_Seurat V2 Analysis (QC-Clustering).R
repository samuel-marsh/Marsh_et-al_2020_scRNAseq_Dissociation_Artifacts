# Two cells have extremly high Mito % (>40%) while everything else is <20.  Filter now so that plots scale appropriately for pre-filter plots 
exp17_micro_ALL <- FilterCells(object = exp17_micro_ALL, subset.names = c("nGene", "percent_mito", "nUMI"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(6000, 0.3, 12000))

# Filter based on data
exp17_micro_ALL <- FilterCells(object = exp17_micro_ALL, subset.names = c("nGene", "percent_mito", "nUMI"), low.thresholds = c(600, -Inf, -Inf), high.thresholds = c(2000, 0.1, 3500))

# Normalize 
exp17_micro_ALL <- NormalizeData(object = exp17_micro_ALL, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)

# Find Variable genes in the data
exp17_micro_ALL <- FindVariableGenes(object = exp17_micro_ALL, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)

# Check how many variable genes were identified
length(x = exp17_micro_ALL@var.genes)

# Scale & Regress
exp17_micro_ALL <- ScaleData(object = exp17_micro_ALL, vars.to.regress = c("nUMI", "percent_mito"))

# PCA & Evaluation of PCs
exp17_micro_ALL <- RunPCA(object = exp17_micro_ALL, pc.genes = exp17_micro_ALL@var.genes, pcs.compute = 40)

exp17_micro_ALL <- ProjectPCA(object = exp17_micro_ALL, do.print = FALSE)

PCElbowPlot(object = exp17_micro_ALL, num.pc = 40)

# Clustering
exp17_micro_ALL <- FindClusters(object = exp17_micro_ALL, reduction.type = "pca", dims.use = 1:15, resolution = .2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

# Run TSNE 
exp17_micro_ALL <- RunTSNE(object = exp17_micro_ALL, dims.use = 1:15)

TSNEPlot(object = exp17_micro_ALL, do.label = TRUE)

write_rds(exp17_micro_ALL, "OBJECT_NAME.RDS")