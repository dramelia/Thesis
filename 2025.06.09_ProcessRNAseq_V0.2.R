# See https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html and
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog
# for help.
# NOTE: will need to manually alter sample.txt after creation by prepareRsemForTximport.py to add contrasts that need to be explored
# Will then need to add those to code below, at DESeqDataSetFromTximport etc

# NOTE2: Have removed rlog from this script as often too slow
# IMPORTANT -- The Group names need to only contain "_" and ".", definitely not "-"!

pacman::p_load(tidyverse)
pacman::p_load(RColorBrewer)
pacman::p_load(DESeq2)
pacman::p_load(pheatmap)
pacman::p_load(DEGreport)
pacman::p_load(ggplot2)
pacman::p_load(tximport)
pacman::p_load(ggrepel)
pacman::p_load(umap)
pacman::p_load(sva)
pacman::p_load(apeglm)
pacman::p_load(factoextra)

#--------------------------- Setup paths and parameters ------------------------#
# Set working directory to where your nf-core results are located
dir <- file.path("my", "file", "path", "Isoforms_Sun") # change file path
setwd(dir)

# Input file names (edit as needed)
metaFile <- "samples_test.txt" # Sample metadata file
rsemFolder <- "RSEM_Quant" # Folder with RSEM quantification results
tx2geneFile <- "tx2gene_nrGencodev48.csv" # Transcript-to-gene mapping file

# Output folder names
baseResultsDir <- "DESeq2_Results"
dataFolder <- "Data_All"
diffFolder <- "Diff_All"

# Output file naming conventions and filtering parameters
prefix <- "RNAseq_CLL_"
suffix <- ".txt"
MIN_COUNT_PER_ROW <- 1 # Minimum count per gene to keep
MIN_MEAN_PER_ROW <- 1 # Minimum mean count per gene to keep
FILTER_BY_MIN_COUNT <- TRUE # If FALSE, filter by mean instead

# Subsetting
getGroups <- c(1:3) # CHANGE ME

#  EXAMPLE MANUAL COLOURS/SHAPES
# For good colours see https://www.thinkingondata.com/something-about-viridis-library/
# colours <- c("#000000", "#440D54", "#440D54", "#39568C", "#39568C", "#1F968B", "#1F968B", "#73D055", "#73D055")
colours <- c("#000000", "#7FFFD4", "#440D54", "#B4DE2C", "#829c30", "#596e19")
shapes <- c(9, 19, 19, 19, 19, 19)
MANUAL_COLOURS_SHAPES <- FALSE # Set to TRUE to use manual colors/shapes

# PCA, MDS, UMAP settings
GROUP_BY <- "1" # Grouping variable for plots ("1" = no grouping)
UMAP_STATE <- 1234 # Random seed for UMAP
PCA_GENE_NUM <- 2000 # Number of genes for PCA
UMAP_PCA_MAX <- 10 # Number of PCs for UMAP

# Plot dimensions
CORR_HEATMAP_WIDTH <- 60
CORR_HEATMAP_HEIGHT <- 60
PCA_MDS_WIDHT <- 10
PCA_MDS_HEIGHT <- 10
UMAP_WIDTH <- 10
UMAP_HEIGHT <- 10
PCA_POINT_SIZE <- 4
UMAP_POINT_SIZE <- 4
MA_Y_LIMITS <- c(-3, 3) # Y-axis limits for MA plots

#--------------------------- Output folder setup -------------------------------#
# Create output directories if they do not exist
if (MIN_COUNT_PER_ROW) {
  dataFolder <- paste(dataFolder, "_MCPR", MIN_COUNT_PER_ROW, sep = "")
  diffFolder <- paste(diffFolder, "_MCPR", MIN_COUNT_PER_ROW, sep = "")
} else {
  dataFolder <- paste(dataFolder, "_MMPR", MIN_MEAN_PER_ROW, sep = "")
  diffFolder <- paste(diffFolder, "_MMPR", MIN_MEAN_PER_ROW, sep = "")
}
baseDir <- file.path(dir, baseResultsDir)
dir.create(baseDir)
diffFolder <- file.path(baseDir, diffFolder)
dir.create(diffFolder)
dataFolder <- file.path(baseDir, dataFolder)
dir.create(dataFolder)

################################################################################
#                               Functions                                      #
################################################################################

#' Plot PCA of samples using top variable genes
#'
#' @param object SummarizedExperiment or DESeq2 object
#' @param intgroup Column name in colData to group samples by
#' @param ntop Number of top variable genes to use
#' @param pc1,pc2 Principal components to plot
#' @param returnData If TRUE, return data frame instead of plotting
#' @param fileName Output PDF file name
#' @param width,height Plot dimensions
#' @param textLabels Labels for points
#' @param annotatePoints If TRUE, add text labels to points
#' @param combat_exp Optional batch-corrected expression matrix
#' @param colours,shapes Colors and shapes for groups
#' @param pointSize Size of points
#' @param manualColours,manualShapes Use manual colors/shapes
#' @param forceSquarePlot If TRUE, force square aspect ratio
#' @return PCA object (from prcomp)
plotPCA <- function(object, intgroup = "condition", ntop = 500, pc1 = 1, pc2 = 2,
                    returnData = FALSE, fileName = "temp.pdf", width = 8, height = 8,
                    textLabels = colnames(object), annotatePoints = FALSE, combat_exp = NULL,
                    colours = c("#ff7f00", "#e377c2", "#17becf"), shapes = c(25),
                    pointSize = 3, manualColours = FALSE, manualShapes = FALSE, forceSquarePlot = TRUE) {
  # Select top variable genes for PCA
  if (!is.null(combat_exp)) {
    rv <- rowVars(combat_exp)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(combat_exp[select, ]))
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  } else {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  }

  # Handle grouping (if intgroup is "1", treat all samples as one group)
  if (identical(intgroup, "1")) {
    group <- factor(rep("all", ncol(object)))
  } else {
    group <- factor(colData(object)[[intgroup]], levels = unique(colData(object)[[intgroup]]))
  }

  d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc1:pc2]
    return(d)
  }

  # Plot PCA using ggplot2
  aes_mapping <- if (manualShapes) {
    aes_string(x = "PC1", y = "PC2", color = "group", shape = "group")
  } else {
    aes_string(x = "PC1", y = "PC2", color = "group")
  }

  p <- ggplot(data = d, mapping = aes_mapping) +
    geom_point(size = pointSize) +
    xlab(paste0("PC", pc1, ": ", round(percentVar[pc1] * 100), "% variance")) +
    ylab(paste0("PC", pc2, ": ", round(percentVar[pc2] * 100), "% variance")) +
    coord_fixed() +
    theme_classic() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL)

  if (manualShapes) {
    p <- p + scale_shape_manual(values = shapes)
  }
  if (annotatePoints) {
    p <- p + geom_text_repel(aes(label = textLabels), max.overlaps = Inf)
  }
  if (manualColours) {
    p <- p + scale_color_manual(values = colours)
  }
  if (forceSquarePlot) {
    p <- p + theme(aspect.ratio = 1)
  }

  ggsave(fileName, width = width, height = height, useDingbats = FALSE)
  return(pca)
}

#' Plot MDS (Multi-dimensional Scaling) of samples
#'
#' @param object SummarizedExperiment or DESeq2 object
#' @param intgroup Column name in colData to group samples by
#' @param fileName Output PDF file name
#' @param width,height Plot dimensions
#' @param textLabels Labels for points
#' @param annotatePoints If TRUE, add text labels to points
#' @param combat_exp Optional batch-corrected expression matrix
#' @param colours,shapes Colors and shapes for groups
#' @param pointSize Size of points
#' @param manualColours,manualShapes Use manual colors/shapes
#' @param forceSquarePlot If TRUE, force square aspect ratio
plotMDS <- function(object, intgroup = "Condition", fileName = "temp.pdf", width = 8, height = 8,
                    textLabels = colnames(object), annotatePoints = FALSE, combat_exp = NULL,
                    colours = c("#ff7f00", "#e377c2", "#17becf"), shapes = c(25),
                    pointSize = 3, manualColours = FALSE, manualShapes = FALSE, forceSquarePlot = TRUE) {
  # Compute sample distances and perform MDS
  if (!is.null(combat_exp)) {
    sampleDists <- dist(t(combat_exp))
    sampleDistMatrix <- as.matrix(sampleDists)
    mdsData <- data.frame(cmdscale(sampleDistMatrix))
  } else {
    sampleDists <- dist(t(assay(object)))
    sampleDistMatrix <- as.matrix(sampleDists)
    mdsData <- data.frame(cmdscale(sampleDistMatrix))
  }

  # Handle grouping
  if (identical(intgroup, "1")) {
    group <- factor(rep("all", ncol(object)))
  } else {
    group <- factor(colData(object)[[intgroup]], levels = unique(colData(object)[[intgroup]]))
  }

  d <- data.frame(MDS1 = mdsData$X1, MDS2 = mdsData$X2, group = group, name = colnames(object))

  # Plot MDS using ggplot2
  aes_mapping <- if (manualShapes) {
    aes_string(x = "MDS1", y = "MDS2", color = "group", shape = "group")
  } else {
    aes_string(x = "MDS1", y = "MDS2", color = "group")
  }

  p <- ggplot(data = d, mapping = aes_mapping) +
    geom_point(size = pointSize) +
    xlab("MDS1") +
    ylab("MDS2") +
    coord_fixed() +
    theme_classic() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL)

  if (manualShapes) {
    p <- p + scale_shape_manual(values = shapes)
  }

  if (annotatePoints) {
    p <- p + geom_text_repel(aes(label = textLabels), max.overlaps = Inf)
  }

  if (manualColours) {
    p <- p + scale_color_manual(values = colours)
  }

  if (forceSquarePlot) {
    p <- p + theme(aspect.ratio = 1)
  }

  ggsave(fileName, width = width, height = height, useDingbats = FALSE)
}

#' Plot UMAP of samples
#'
#' @param object SummarizedExperiment or DESeq2 object
#' @param labels Sample labels
#' @param intgroup Column name in colData to group samples by
#' @param colours,shapes Colors and shapes for groups
#' @param returnData If TRUE, return data frame instead of plotting
#' @param fileName Output PDF file name
#' @param width,height Plot dimensions
#' @param textLabels Labels for points
#' @param annotatePoints If TRUE, add text labels to points
#' @param combat_exp Optional batch-corrected expression matrix
#' @param umap.config UMAP configuration
#' @param pointSize Size of points
#' @param manualColours,manualShapes Use manual colors/shapes
#' @param forceSquarePlot If TRUE, force square aspect ratio
plotUMAP <- function(object, labels, intgroup = "condition",
                     colours = c("#ff7f00", "#e377c2", "#17becf"), shapes = c(25),
                     returnData = FALSE, fileName = "temp.pdf", width = 8, height = 8,
                     textLabels = colnames(object), annotatePoints = FALSE, combat_exp = NULL, umap.config = NULL,
                     pointSize = 3, manualColours = FALSE, manualShapes = FALSE, forceSquarePlot = TRUE) {
  # Set UMAP configuration
  if (is.null(umap.config)) {
    custom.config <- umap.defaults
    custom.config$random_state <- 1234
  } else {
    custom.config <- umap.config
  }

  # Run UMAP on expression matrix
  if (!is.null(combat_exp)) {
    umap.obj <- umap(t(combat_exp), custom.config)
  } else {
    umap.obj <- umap(t(assay(object)), custom.config)
  }

  layout <- umap.obj$layout

  # Handle grouping
  if (identical(intgroup, "1")) {
    group <- factor(rep("all", ncol(object)))
  } else {
    group <- factor(colData(object)[[intgroup]], levels = unique(colData(object)[[intgroup]]))
  }

  d <- data.frame(umapX = layout[, 1], umapY = layout[, 2], group = group, name = colnames(object))
  if (returnData) {
    return(d)
  }

  # Plot UMAP using ggplot2
  aes_mapping <- if (manualShapes) {
    aes_string(x = "umapX", y = "umapY", color = "group", shape = "group")
  } else {
    aes_string(x = "umapX", y = "umapY", color = "group")
  }

  p <- ggplot(data = d, mapping = aes_mapping) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    geom_point(size = pointSize) +
    coord_fixed() +
    theme_classic() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL)

  if (manualShapes) {
    p <- p + scale_shape_manual(values = shapes)
  }
  if (annotatePoints) {
    p <- p + geom_text_repel(aes(label = textLabels), max.overlaps = Inf)
  }
  if (manualColours) {
    p <- p + scale_color_manual(values = colours)
  }
  if (forceSquarePlot) {
    p <- p + theme(aspect.ratio = 1)
  }

  ggsave(fileName, width = width, height = height, useDingbats = FALSE)
}

#' Plot UMAP using principal components as input
#'
#' @param object SummarizedExperiment or DESeq2 object
#' @param labels Sample labels
#' @param intgroup Column name in colData to group samples by
#' @param colours,shapes Colors and shapes for groups
#' @param returnData If TRUE, return data frame instead of plotting
#' @param fileName Output PDF file name
#' @param width,height Plot dimensions
#' @param textLabels Labels for points
#' @param annotatePoints If TRUE, add text labels to points
#' @param combat_exp Optional batch-corrected expression matrix
#' @param umap.config UMAP configuration
#' @param pointSize Size of points
#' @param manualColours,manualShapes Use manual colors/shapes
#' @param ntop Number of top variable genes for PCA
#' @param pcaMax Number of principal components for UMAP
#' @param forceSquarePlot If TRUE, force square aspect ratio
plotUMAP_PCA <- function(object, labels, intgroup = "condition",
                         colours = c("#ff7f00", "#e377c2", "#17becf"), shapes = c(25),
                         returnData = FALSE, fileName = "temp.pdf", width = 8, height = 8,
                         textLabels = colnames(object), annotatePoints = FALSE, combat_exp = NULL, umap.config = NULL,
                         pointSize = 3, manualColours = FALSE, manualShapes = FALSE, ntop = 500, pcaMax = 10,
                         forceSquarePlot = TRUE) {
  # Run PCA on top variable genes
  if (!is.null(combat_exp)) {
    rv <- rowVars(combat_exp)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(combat_exp[select, ]))
  } else {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
  }

  if (pcaMax > length(pca$sdev)) {
    pcaMax <- length(pca$sdev)
  }
  print(paste("Using the top", pcaMax, "PC for UMAP"))

  # Set UMAP configuration
  if (is.null(umap.config)) {
    custom.config <- umap.defaults
    custom.config$random_state <- 1234
  } else {
    custom.config <- umap.config
  }

  # Run UMAP on principal components
  umap.obj <- umap(pca$x[, 1:pcaMax], custom.config)
  layout <- umap.obj$layout

  # Handle grouping
  if (identical(intgroup, "1")) {
    group <- factor(rep("all", ncol(object)))
  } else {
    group <- factor(colData(object)[[intgroup]], levels = unique(colData(object)[[intgroup]]))
  }

  d <- data.frame(umapX = layout[, 1], umapY = layout[, 2], group = group, name = colnames(object))
  if (returnData) {
    return(d)
  }

  # Plot UMAP using ggplot2
  aes_mapping <- if (manualShapes) {
    aes_string(x = "umapX", y = "umapY", color = "group", shape = "group")
  } else {
    aes_string(x = "umapX", y = "umapY", color = "group")
  }

  p <- ggplot(data = d, mapping = aes_mapping) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    geom_point(size = pointSize) +
    coord_fixed() +
    theme_classic() +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL)

  if (manualShapes) {
    p <- p + scale_shape_manual(values = shapes)
  }
  if (annotatePoints) {
    p <- p + geom_text_repel(aes(label = textLabels), max.overlaps = Inf)
  }
  if (manualColours) {
    p <- p + scale_color_manual(values = colours)
  }
  if (forceSquarePlot) {
    p <- p + theme(aspect.ratio = 1)
  }

  ggsave(fileName, width = width, height = height, useDingbats = FALSE)
}

#' Plot total counts per sample from tximport object
#'
#' @param txiO tximport object
#' @param fileName Output file name (without extension)
#' @param width,height Plot dimensions
plot_TXI_countSum <- function(txiO, fileName = "temp.pdf", width = 8, height = 8) {
  # Sum counts for each sample
  cc <- txiO$counts
  dim(cc)
  countSumPerSample <- apply(cc, 2, sum)
  countSumPerSampleDF <- as_data_frame(countSumPerSample, rownames = "Sample")

  # Barplot of total counts per sample
  p <- ggplot(data = countSumPerSampleDF, aes(x = Sample, y = value)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("TXI Count Sum")
  ggsave(paste(fileName, ".pdf", sep = ""), width = width, height = height)
  write.table(countSumPerSampleDF, paste(fileName, ".txt", sep = ""), sep = "\t", row.names = FALSE)
}

################################################################################
#                        Data Import and Processing                            #
################################################################################

# List files in working directory (for debugging)
list.files(dir)

# Read sample metadata
samples <- read.table(file.path(dir, metaFile), header = TRUE, sep = "\t")
rownames(samples) <- samples$Sample

# Generate list of RSEM quantification files for each sample
files <- file.path(dir, rsemFolder, paste(samples$run, ".isoforms.results", sep = ""))
names(files) <- samples$run

# Check that all files exist
files
all(file.exists(files))

# Read transcript-to-gene mapping
tx2gene <- read_csv(file.path(dir, tx2geneFile))
head(tx2gene)

# Select all samples for analysis (edit getGroupsA as needed)
names(files)
getGroupsA <- c(0:length(files))
selectF <- files[getGroupsA]
selectSamplesA <- samples[getGroupsA, ]
names(selectF)

################################################################################
# Import quantification data using tximport
txi <- tximport(selectF, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)
names(txi)
head(txi$counts)
head(txi$abundance)

# Write out raw counts and TPMs for inspection
write.table(as.data.frame("Gene_ID" = rownames(txi$counts), txi$counts), paste(dataFolder, "/", prefix, "RSEM_Gene_count", suffix, sep = ""), sep = "\t", row.names = FALSE)
write.table(as.data.frame("Gene_ID" = rownames(txi$abundance), txi$abundance), paste(dataFolder, "/", prefix, "RSEM_Gene_TPM", suffix, sep = ""), sep = "\t", row.names = FALSE)

# Baseline QC: plot total counts per sample
plot_TXI_countSum(txi, paste(dataFolder, "/", prefix, "RSEM_Gene_TXIcountSum", sep = ""), width = 20, height = 8)

################################################################################
# Import data into DESeq2 for normalization and variance stabilization

# Create DESeq2 dataset (design formula can be changed as needed)
dds <- DESeqDataSetFromTximport(txi, colData = selectSamplesA, design = formula(paste("~", GROUP_BY))) # CHANGE ME

# Filter out lowly-expressed genes
nrow(dds)
if (FILTER_BY_MIN_COUNT) {
  keep <- rowSums(counts(dds)) >= MIN_COUNT_PER_ROW
} else {
  keep <- rowMeans(counts(dds)) >= MIN_MEAN_PER_ROW
}
dds <- dds[keep, ]
nrow(dds)

# If avgTxLength is present, use it for normalization
nm <- assays(dds)[["avgTxLength"]]
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds) / nm)
sizeFactors(dds)

# Retrieve normalized counts
normalized_counts <- counts(dds, normalize = TRUE)

# Save normalized counts to file
write.table(normalized_counts, file = paste(dataFolder, "/", prefix, "RSEM_Gene_normalisedCounts", suffix, sep = ""), sep = "\t", quote = F, col.names = NA)

# Variance stabilizing transformation (VST) for downstream analysis
vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)

# Save VST-transformed data
write.table(data.frame("Gene_ID" = rownames(assay(vsd)), assay(vsd)), paste(dataFolder, "/", prefix, "RSEM_Gene_count_VST", suffix, sep = ""), sep = "\t", row.names = FALSE)

# Select all samples for visualization
names(files)
selectF <- files[getGroupsA]
selectSamples <- samples[getGroupsA, ]
vsd_S <- vsd[, getGroupsA]
names(selectF)
selectSamples

################################################################################
#                        Visualization and QC Plots                            #
################################################################################

#--------------------------- UMAP configuration -------------------------------#
umap.config <- umap.defaults
umap.config$random_state <- UMAP_STATE
umap.config$metric <- "euclidean"

#--------------------------- PCA, MDS, UMAP plots -----------------------------#
if (MANUAL_COLOURS_SHAPES) {
  # PCA, MDS, UMAP with manual colors and shapes
  res.pca <- plotPCA(vsd_S,
    intgroup = GROUP_BY, pc1 = 1, pc2 = 2, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_PCA-ntop", PCA_GENE_NUM, ".pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, ntop = PCA_GENE_NUM, pointSize = PCA_POINT_SIZE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, colours = colours
  )
  res.pca.2 <- plotPCA(vsd_S,
    intgroup = GROUP_BY, pc1 = 1, pc2 = 2, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_PCA-ntop", PCA_GENE_NUM, "WithLabels.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, ntop = PCA_GENE_NUM, pointSize = PCA_POINT_SIZE, annotatePoints = TRUE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, colours = colours
  )
  plotMDS(vsd_S,
    intgroup = GROUP_BY, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_MDS.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, pointSize = PCA_POINT_SIZE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, colours = colours
  )
  plotMDS(vsd_S,
    intgroup = GROUP_BY, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_MDS_WithLabels.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, pointSize = PCA_POINT_SIZE, annotatePoints = TRUE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, colours = colours
  )
  plotUMAP(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP.pdf", sep = ""), width = UMAP_WIDTH, height = UMAP_HEIGHT, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes
  )
  plotUMAP(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP_WithLables.pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, annotatePoints = TRUE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes
  )
  plotUMAP_PCA(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP-PCA-ntop", PCA_GENE_NUM, "-pc", UMAP_PCA_MAX, ".pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, pcaMax = UMAP_PCA_MAX, ntop = PCA_GENE_NUM
  )
  plotUMAP_PCA(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP-PCA-ntop", PCA_GENE_NUM, "-pc", UMAP_PCA_MAX, "_WithLables.pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, annotatePoints = TRUE, manualColours = TRUE,
    manualShapes = TRUE, shapes = shapes, pcaMax = UMAP_PCA_MAX, ntop = PCA_GENE_NUM
  )
} else {
  # PCA, MDS, UMAP with default colors and shapes
  res.pca <- plotPCA(vsd_S,
    intgroup = GROUP_BY, pc1 = 1, pc2 = 2, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_PCA-ntop", PCA_GENE_NUM, ".pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, ntop = PCA_GENE_NUM, pointSize = PCA_POINT_SIZE, manualColours = FALSE
  )
  res.pca.2 <- plotPCA(vsd_S,
    intgroup = GROUP_BY, pc1 = 1, pc2 = 2, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_PCA-ntop", PCA_GENE_NUM, "WithLabels.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, ntop = PCA_GENE_NUM, pointSize = PCA_POINT_SIZE, annotatePoints = TRUE, manualColours = FALSE
  )
  plotMDS(vsd_S,
    intgroup = GROUP_BY, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_MDS.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, pointSize = PCA_POINT_SIZE
  )
  plotMDS(vsd_S,
    intgroup = GROUP_BY, fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_MDS_WithLabels.pdf", sep = ""),
    width = PCA_MDS_WIDHT, height = PCA_MDS_HEIGHT, pointSize = PCA_POINT_SIZE, annotatePoints = TRUE
  )
  plotUMAP(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP.pdf", sep = ""), width = UMAP_WIDTH, height = UMAP_HEIGHT
  )
  plotUMAP(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP_WithLables.pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, annotatePoints = TRUE
  )
  plotUMAP_PCA(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP-PCA-ntop", PCA_GENE_NUM, "-pc", UMAP_PCA_MAX, ".pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, pcaMax = UMAP_PCA_MAX, ntop = PCA_GENE_NUM
  )
  plotUMAP_PCA(vsd_S,
    colours = colours, intgroup = GROUP_BY, umap.config = umap.config, pointSize = UMAP_POINT_SIZE,
    fileName = paste(dataFolder, "/", prefix, "RSEM_count_VST_UMAP-PCA-ntop", PCA_GENE_NUM, "-pc", UMAP_PCA_MAX, "_WithLables.pdf", sep = ""),
    width = UMAP_WIDTH, height = UMAP_HEIGHT, annotatePoints = TRUE, pcaMax = UMAP_PCA_MAX, ntop = PCA_GENE_NUM
  )
}

#--------------------------- Sample distance heatmap ---------------------------#
# Compute sample distances using VST data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$SampleName, vsd$Type, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot heatmap of sample distances
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors, filename = paste(dataFolder, "/", prefix, "RSEM_count_VST_Heatmap_distances.pdf", sep = ""),
  width = CORR_HEATMAP_WIDTH,
  height = CORR_HEATMAP_HEIGHT
)

#--------------------------- Sample correlation heatmap ------------------------#
# Extract VST matrix and compute sample correlations
vsd_mat <- assay(vsd) # VST-transformed expression matrix
vsd_cor <- cor(vsd_mat) # Compute pairwise sample correlations
head(vsd_cor) # Print first few rows for inspection

# Plot heatmap of sample correlations
pheatmap(vsd_cor,
  filename = paste(dataFolder, "/", prefix, "RSEM_count_VST_Heatmap_corr.pdf", sep = ""),
  width = CORR_HEATMAP_WIDTH,
  height = CORR_HEATMAP_HEIGHT
)

################################################################################
#                        PCA Component Exploration                             #
################################################################################

# Use factoextra to extract PCA variable information
var <- get_pca_var(res.pca)

# Uncomment below to inspect PCA variable coordinates, cos2, and contributions
# head(var$coord)   # Coordinates of variables on principal components
# head(var$cos2)    # Quality of representation (cos2) on factor map
# head(var$contrib) # Contributions to principal components

# Write out cos2 (quality of representation) for each gene
write.table(var$cos2, paste(dataFolder, "/", prefix, "RSEM_count_VST_PCA-ntop", PCA_GENE_NUM, "_Cos2.txt", sep = ""), sep = "\t")
#####################################

