# =============================================================================
# Weighted Gene Co-expression Network Analysis (WGCNA)
# Staphylococcus aureus RNA-seq data (RPKM-normalized expression)
# =============================================================================

library(WGCNA)
library(flashClust)
library(doParallel)

# Enable multi-threaded calculations for WGCNA
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# Set working directory (update to your local project folder)
setwd('../Data')

# -----------------------------------------------------------------------------
# 1. Load expression data
# -----------------------------------------------------------------------------
# RPKM expression matrix: rows = genes (LocusTag), columns = samples
df <- read.csv('WGCNARPKMSaureus.txt', header = TRUE, sep = '\t')
df[is.na(df)] <- 0
head(df)

EGData <- df

# Transpose so that rows = samples and columns = genes, as required by WGCNA
datExpr0 <- as.data.frame(t(EGData[, -c(1:1)]))
names(datExpr0) <- EGData$LocusTag
rownames(datExpr0) <- names(EGData)[-c(1:1)]

# -----------------------------------------------------------------------------
# 2. Quality control: genes and samples
# -----------------------------------------------------------------------------
# Check for genes/samples with too many missing values or zero variance
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# Hierarchical clustering of samples to visually detect outliers
sampleTree <- hclust(dist(datExpr0), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Height cutoff used to flag outlying samples
abline(h = 1000000, col = "red")

# Keep samples that fall into the main cluster below the cutoff
clust <- cutreeStatic(sampleTree, cutHeight = 1000000, minSize = 2)
table(clust)
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Robust PCA (rrcov) to perform a second, more sensitive round of outlier
# detection on the retained samples
library(rrcov)
datcomp <- datExpr
pcaHub <- PcaHubert(datcomp)
pcaHub$flag
outliers <- which(pcaHub@flag == 'FALSE')

# Keep only samples flagged as non-outliers by the robust PCA
datExpr <- datcomp[pcaHub$flag, ]

# -----------------------------------------------------------------------------
# 3. Soft-thresholding power selection
# -----------------------------------------------------------------------------
# Candidate powers to evaluate for approximate scale-free topology
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Biweight midcorrelation (bicor) is more robust to outliers than Pearson
sft <- pickSoftThreshold(datExpr, networkType = "signed", powerVector = powers,
                          verbose = 5, corFnc = "bicor")

sizeGrWindow(10, 5)
par(mar = c(6, 6, 6, 2), mfrow = c(1, 2))
cex1 <- 1.5

# Scale-free topology fit (R^2) as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence",
     cex.lab = 2, cex.axis = 2, cex.main = 2)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# Reference line for an R^2 cutoff of 0.8
abline(h = 0.8, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity",
     cex.lab = 2, cex.axis = 2, cex.main = 2)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# Selected soft-thresholding power (yields R^2 ~ 0.8; see manuscript Methods)
softPower <- 18

# -----------------------------------------------------------------------------
# 4. Network construction (adjacency and TOM)
# -----------------------------------------------------------------------------
# Signed adjacency matrix using bicor correlation
adjacency <- adjacency(datExpr, type = "signed", power = softPower, corFnc = "bicor")

# Topological Overlap Matrix (TOM), which reduces noise/spurious connections
TOM <- TOMsimilarityFromExpr(adjacency, networkType = "signed", power = softPower,
                              corType = "bicor")
dissTOM <- 1 - TOM

# -----------------------------------------------------------------------------
# 5. Gene clustering and module detection
# -----------------------------------------------------------------------------
# Hierarchical clustering of genes based on TOM dissimilarity
geneTree <- flashClust(as.dist(dissTOM), method = "average")

sizeGrWindow(12, 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# Minimum number of genes required to form a module
minModuleSize <- 20

# Dynamic tree cut to identify modules; label 0 = unassigned genes
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                              pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric module labels into colors (grey = unassigned genes)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     main = "Gene dendrogram and module colors")

# -----------------------------------------------------------------------------
# 6. Merge modules with highly similar eigengenes
# -----------------------------------------------------------------------------
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Dissimilarity between module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Modules whose eigengenes are closer than this threshold are merged
MEDissThres <- 0.10
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 5)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

sizeGrWindow(14, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                     c("Dynamic Tree Cut", "Merged dynamic"),
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     cex.axis = 1.5, cex.lab = 1.5, cex.main = 2,
                     cex.colorLabels = 1.5, marAll = c(2, 10, 3, 0))

# Final module colors and eigengenes after merging
moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

# Save results for downstream analyses
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Staphy-RPKM-Signed-Bicor-newman.RData")

# -----------------------------------------------------------------------------
# 7. Export the full network for Cytoscape
# -----------------------------------------------------------------------------
# Exclude the grey (unassigned) module
modules <- unique(moduleColors)
modules <- modules[modules != 'grey']

probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]

modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = "Signed-edges-Staphy-RPKM-bicor00Newman.txt",
  nodeFile = "Signed-nodes-Staphy-RPKM-bicor00Newman.txt",
  weighted = TRUE,
  threshold = 0.0,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

# Gene-to-module color table
dfrm <- data.frame(NameGen = modProbes, Mcolors = moduleColors[inModule])
write.table(dfrm, file = "ModulesStaphy-RPKM-Signed-bicorNewman.csv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------
# 8. Identify hub genes (highest intramodular connectivity per module)
# -----------------------------------------------------------------------------
hubs <- chooseTopHubInEachModule(datExpr, colorh = moduleColors, power = 18,
                                  type = "signed", corFnc = "bicor")
hubs

# -----------------------------------------------------------------------------
# 9. Export a sub-network for specific modules of interest (Cytoscape)
# -----------------------------------------------------------------------------
modules <- c("royalblue")
probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]

modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste0("StaphyCytoscapeInput-edges-100", paste(modules, collapse = "-"), ".txt"),
  nodeFile = paste0("StaphyCytoscapeInput-nodes-100", paste(modules, collapse = "-"), ".txt"),
  weighted = TRUE,
  threshold = 0.031,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

# -----------------------------------------------------------------------------
# 10. Visualize expression profiles of genes within selected modules
# -----------------------------------------------------------------------------
library(magrittr)
library(tidyverse)

module_df <- data.frame(
  gene_id = modProbes,
  colors = moduleColors[inModule]
)

# Modules to plot
modules_of_interest <- c("black", "magenta", "lightgreen", "royalblue")

# Genes belonging to the modules of interest
submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

# Log-transformed expression matrix (genes x samples)
dataexp <- log(t(datExpr))

write.table(datExpr, file = "datExpr.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Expression of the selected genes, reshaped to long format for ggplot
subexpr <- dataexp[submod$gene_id, ]
submod_df <- data.frame(subexpr) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id) %>%
  mutate(module = module_df[gene_id, ]$colors)

# One expression-profile panel per module, genes as individual lines
submod_df %>%
  ggplot(., aes(x = name, y = value, group = gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "none"
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Treatment", y = "log(expression level)") +
  scale_color_manual(values = c("#000000", "#ff00ff", "#90ee90", "#4169e1"))
