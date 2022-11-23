## Installation
devtools::install_github("navinlabcode/copykit")
packageVersion("copykit")

library(copykit)
library(BiocParallel)

## set up parallelization
register(MulticoreParam(workers = 16, progressbar = T), default = T)
bpparam()

## run CBS segmentation from bam files 
### (This part takes ~10 mins thus will be skipped in the workshop)
# tumor <- runVarbin(
#           dir = "~/path/to/bam/files/",
#           genome = "hg38",
#           resolution = "220kb", # c("220kb", "55kb", "110kb", "195kb", "280kb", "500kb", "1Mb", "2.8Mb")
#           remove_Y = FALSE,
#           is_paired_end = FALSE,
#           seed = 17,
#           BPPARAM = bpparam()
#         )

## example of basic segmentation and calculate met
copykit_obj <- mock_bincounts(ncells = 5,
                              ncells_diploid = 1,
                              position_gain = c(1,5,1000))
copykit_obj <- runMetrics(copykit_obj)
colData(copykit_obj)
plotMetrics(copykit_obj, metric = "overdispersion", label = "ground_truth")

## Load in sample data
tumor <- readRDS("sample_obj.rds")
head(tumor)
rowData(tumor)
segment_ratios(tumor)[1:10,1:5]

## QC
### Find aneuploid cells
tumor <- findAneuploidCells(tumor, seed = 17)
### Find outlier cells
tumor <- findOutliers(tumor, k = 5, resolution = 0.9)
colData(tumor)
### visual the results
plotHeatmap(tumor, 
            order_cells = 'hclust',
            label = c('is_aneuploid','outlier'),
            row_split = 'outlier', 
            n_threads = 40)
## remove the unwanted cells
tumor <- tumor[,colData(tumor)$is_aneuploid == TRUE]
tumor <- tumor[,colData(tumor)$outlier == FALSE]

## Clustering
### Dimension reduction
tumor <- runUmap(tumor, seed = 17)
reducedDim(tumor)
plotUmap(tumor)
### determine k
tumor <- findSuggestedK(tumor, k_range = 5:9)
plotSuggestedK(tumor)
plotSuggestedK(tumor, geom = 'tile')
### hdbscan clustering 
tumor <- findClusters(tumor, k_subclones = 6)
### Visual
tumor <- calcConsensus(tumor);tumor <- runConsensusPhylo(tumor)
plotUmap(tumor, label = 'subclones')
plotHeatmap(tumor, 
            label = 'subclones', 
            n_threads = 40)
### remove cluster outliers
tumor_c0 <- tumor[,colData(tumor)$subclones == 'c0']
tumor <- tumor[,colData(tumor)$subclones != 'c0']
colData(tumor)$subclones <- droplevels(colData(tumor)$subclones)
### merge dataset
merged <- cbind(tumor, tumor_c0)

## Ploidy and Integer estimation
### calculate from fixed ploidy
tumor <- calcInteger(tumor, 
                     method = 'fixed', 
                     ploidy_value = 4.3)
### visual
plotHeatmap(tumor, 
            assay = 'integer',
            label = 'subclones', 
            n_threads = 40)
### plot rounding errors
plotHeatmap(tumor, 
            assay = 'integer',
            rounding_error = TRUE,
            label = 'subclones', 
            n_threads = 40)
### profiles of a cell from c1
plotRatio(tumor, sample_name = 'PMTC6LiverC281AL6L7S5_S1433_L004_R1_001')

## Phylogeny
### single cell tree
tumor <- runPhylo(tumor,
                  assay = 'segment_ratios',
                  metric = 'manhattan',
                  n_threads = 40)
phylo(tumor)
### visual
plotPhylo(tumor, label = 'subclones')

### consensus tree
tumor <- calcConsensus(tumor, assay = 'segment_ratios')
consensus(tumor)[1:10,]
### visual
tumor <- runConsensusPhylo(tumor, root = 'neutral')
consensusPhylo(tumor)
plotPhylo(tumor, label = 'subclones', consensus = TRUE)

### heatmap with gene annotation
genes = c("CDKN2A",
          "FGFR1",
          "TP53",
          "PTEN",
          "MYC",
          "CDKN1A",
          "MDM2",
          "AURKA",
          "PIK3CA",
          "CCND1",
          "KRAS")

tumor <- calcConsensus(tumor, assay = 'integer')
plotHeatmap(tumor, 
            consensus = TRUE,
            label = 'subclones',
            genes = genes,
            n_threads = 40)

### check copy number across genes
plotGeneCopy(tumor, 
             genes = genes,
             label = 'subclones',
             dodge.width = 0.8)

### customize ggplot2 object
plotGeneCopy(tumor, 
             genes = genes,
             label = 'subclones',
             dodge.width = 0.8) +
  ggplot2::ylim(c(0,4)) +
  ggplot2::ggtitle("Copy number status across genes") +
  ggplot2::theme_linedraw()

## session info
sessionInfo()









