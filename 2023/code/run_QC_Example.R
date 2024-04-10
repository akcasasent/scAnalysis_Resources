
library(Seurat) ## read 10X
library(ggplot2) ## ggplot
library(cowplot) ## plot_grid
library(matrixStats) ## rowMedians 
library(dplyr) ## %>% 

# general variables 
# colors and such
barfill <- "#4271AE"
barlines <- "#1F3552"
gene_no_filter=2000


##### Basic QC ##### 
## sample path -- Please EDIT
Sample1.path <- "Desktop/Class_Prep/Workshop/datasets/GSE205013_RAW/P01/"
Sample1.name <- "P01"

## read in a sample 
Sample1.mtx <- Read10X(data.dir = paste0(Sample1.path)) # sets up matrix 

# get the count of genes per cell
Sample1.gene <- apply(Sample1.mtx, 2, function(x) sum(x>0)) # number of genes per cell

# just provide text to be added to graph
Sample1_num <-paste0(sum(Sample1.gene>gene_no_filter)," / ",length(Sample1.gene)," cells >2000 genes")

# Graph of Detected genes per cell
Sample1_gene_plot <- ggplot() + 
  geom_histogram(aes(x = Sample1.gene),bins = 200, colour = barlines, fill = barfill) + 
  scale_x_continuous(name = "detected genes per cell",breaks=seq(0,10000,2000),limits=c(0,10000)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(paste("Detected genes per cell in", Sample1.name))
print(Sample1_gene_plot)

# Graph of UMI expression
Sample1_exp <- ggplot() + 
  geom_histogram(aes(x = colSums(as.matrix(Sample1.mtx))), bins = 800,colour = barlines, fill = barfill) + 
  scale_x_continuous(name = "UMI expression",breaks=seq(10000,80000,10000), limits=c(0,80000)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(paste("UMI expression in", Sample1.name))
print(Sample1_exp)

# sets up as data frame of the sample 
# you might get an warning about the vector size.
Sample1_df <- data.frame(nUMI=colSums(as.matrix(Sample1.mtx)), nGene=Sample1.gene) 

# Graph for dot plot comparing UMI and nGenes
# Question what do you expect? 
Sample1_UMI_Gene <- ggplot(Sample1_df, aes(nUMI, nGene)) + 
  geom_point(size=0.5) + scale_x_continuous(breaks=seq(0,200000,100000), limits=c(0,200000)) + 
  theme_bw()+ 
  ggtitle(paste("number of UMI vs number of genes in", Sample1.name))
print(Sample1_UMI_Gene)

# plot in one grid 
Sample1_QC_raw1 <- plot_grid(Sample1_gene_plot, Sample1_exp, Sample1_UMI_Gene, nrow=1)
print(Sample1_QC_raw1)

# create a seurat object 
# this creates the object without restriction on cells or features (genes)
# Include features detected in at least this many cells
min_cells_required = 0
# Include cells where at least this many features are detected
min_features_required = 0
Sample1_seurat <- CreateSeuratObject(Sample1.mtx, min.cells = min_cells_required, 
                                     min.features = min_features_required, 
                                     project = Sample1.name)

# genes per same cell (nFeature_RNA).
# number of UMI reads detected per cell (nCount_RNA)

# Seurat object info
Sample1_seurat

# calculate the mitochondrial DNA 
Sample1_seurat[["percent.mt"]] <- PercentageFeatureSet(object = Sample1_seurat, pattern = "^MT-")

# calculate the ribosomal DNA
Sample1_seurat[["percent.ribo"]] <- PercentageFeatureSet(object = Sample1_seurat, pattern = "^RP[SL]")

# Plot of nCount_RNA and  percent mitochondrial DNA
plot_nRNA_by_mitoP <- FeatureScatter(Sample1_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Plot of nCount_RNA and  percent ribosomal DNA
plot_nRNA_by_riboP <- FeatureScatter(Sample1_seurat, feature1 = "nCount_RNA", feature2 = "percent.ribo")

# Plot of percent mitochondrial DNA and  percent ribosomal DNA
plot_mitoP_by_riboP <- FeatureScatter(Sample1_seurat, feature1 = "percent.mt", feature2 = "percent.ribo")

# Plot of RNA count by feature 
plot_nRNA_by_Feature <- FeatureScatter(Sample1_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#print the three plots

# plot in one grid 
Sample1_QC_process <- plot_grid(plot_nRNA_by_mitoP, plot_nRNA_by_riboP, 
                                plot_mitoP_by_riboP, plot_nRNA_by_Feature, nrow=2)
print(Sample1_QC_process)

# vilion plots of the number of features, count of RNA per cells, mito per cells, brecent prib
# just to show
# VlnPlot(Sample1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

# 
Sample1_mito_plot =VlnPlot(object = Sample1_seurat, features = c("percent.mt"), ncol = 1, pt.size=0.2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Sample1_ribo_plot =VlnPlot(object = Sample1_seurat, features = c("percent.ribo"), ncol = 1, pt.size=0.2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# this is a set cut off of 20% -- note that this cut off is standard for our QC 
# however, it might depend on your organ time. 
mito_cutoff = 20
num_cut=paste0(sum(Sample1_seurat$percent.mt>mito_cutoff), paste0(" cells > ",mito_cutoff,"%"))

# This gives an idea of the cell viabilities based on computational estimates 
#Sample1_cellViability2 <-ggplot(Sample1_seurat@meta.data,aes(x=percent.mt)) + 
#  stat_ecdf(aes(colour=orig.ident)) + 
#  scale_x_continuous(name = "percent.Mitochondrial",breaks=seq(0,100,10),limits=c(0,100)) + 
#  theme_bw() + ylab("The percentage of cells") + 
#  theme(legend.position="none")

# shows the cut of of 20 for the mito in yellow
# shows the number of cells in the percentate that would be dropped... 
Sample1_cellViability1 <- Sample1_mito_plot + geom_hline(yintercept = mito_cutoff, linewidth = 1, colour="yellow", linetype = "dashed")+
  geom_text(aes(y=30), label=num_cut, colour="blue",size =5,x=0.9)+
  ylab("nCell")+theme_bw()+ theme(legend.position="none")+xlab("")

# shows the ribosomal count no line for cut off

Sample1_cellViability3 <- Sample1_ribo_plot + ylab("nCell") + theme_bw() + 
  theme(legend.position="none")+xlab("")

# shows the violin plots of these 
Sample1_QC_raw2 <- plot_grid(Sample1_cellViability1, NULL, Sample1_cellViability3, 
                             ncol=1,rel_heights = c(4,0.5,4),labels = c("A","","B"))
 
Sample1_QC_raw2

#Genewise topGenes QC
getExprProp <- function(expr.data){
  expr.frac <- t(t(expr.data) / Matrix::colSums(expr.data))
  return(expr.frac)
}

#this tests for the genes that you might want to remove 
Sample1_topExpressed <- rowMedians(as.matrix(GetAssayData(object = Sample1_seurat,slot = "counts")))
names(Sample1_topExpressed)<-rownames(GetAssayData(object = Sample1_seurat,slot = "counts"))
mito.genes <- grep('^MT-', names(Sample1_topExpressed), value = TRUE)
ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', names(Sample1_topExpressed), value = TRUE)
rbc.genes<-grep('^HBB|^HBB-BS|^HBB-BT', names(Sample1_topExpressed), value = TRUE) # the paper also included these
dis.genes <- grep('^HBA1|^HBA2|^HBM|^ALAS2)', names(Sample1_topExpressed), value = TRUE) # genes from the paper see methods
Sample1_gene_show <- Sample1_topExpressed %>% sort(.,decreasing=TRUE) %>% names %>% head(n=50)

# note you might get an error or a warning related to size...
expr.frac<-getExprProp(as.matrix(Sample1_seurat[['RNA']]@counts))

# wait for it to finish
Sample1_rate.df.plot <- reshape2::melt(as.data.frame(as.matrix(t(expr.frac[Sample1_gene_show, ]))), id.vars = NULL)
annotation <- rep("other", dim(Sample1_rate.df.plot)[1])
annotation[Sample1_rate.df.plot$variable %in% mito.genes] <- "mitochondrial"
annotation[Sample1_rate.df.plot$variable %in% ribo.genes] <- "ribosome"
annotation[Sample1_rate.df.plot$variable %in% rbc.genes] <- "RBC"
annotation[Sample1_rate.df.plot$variable %in% dis.genes] <- "erythroid"

Sample1_rate.df.plot$annotation<-annotation
Sample1_rate.df.plot$annotation <- factor(Sample1_rate.df.plot$annotation,
                                          levels = c("mitochondrial", "ribosome","RBC","erythroid","other"), ordered = T)
Sample1_rate.df.plot$variable <- factor(Sample1_rate.df.plot$variable,
                                        levels = rev(Sample1_gene_show), ordered = TRUE)

# create a plot of what these genes look like. 
gene.colors <- c(other = "gray70", ribosome = "#FF7F00", erythroid = "#65a55d", mitochondrial = "#3778bf",RBC="#E31A1C")
plot_topgeneOrder <- ggplot() + coord_flip() +
  scale_color_manual(values = "red", name = "") +
  scale_fill_manual(values = gene.colors, name = "Gene:") +
  scale_y_continuous(breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                     trans = 'log10')+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face="bold"),
        axis.text.y =element_text(size=3.5)
  )
plot_topgeneOrder1 <- plot_topgeneOrder + geom_boxplot(data = Sample1_rate.df.plot,
                       mapping = aes(x = variable, y = value, fill = annotation),
                       outlier.size = 0.1,
                       alpha = 0.5) +
  xlab("") + ylab(paste0("Gene proportion of total UMI (1-50)"))
# this plot the top 50 genes and colors them based on status
# what are some items that authors might want to consider?
plot_topgeneOrder1
