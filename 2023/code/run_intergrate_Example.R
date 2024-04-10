
library(Seurat) ## read 10X
library(ggplot2) ## ggplot
library(cowplot) ## plot_grid
library(matrixStats) ## rowMedians 
library(dplyr) ## %>% 

##### Basic intergration ##### 

### read in each sample 
# the best set up for this is a loop. 
# you will need to set for your sytem
base_path <- "Desktop/Class_Prep/Workshop/datasets/"
general_path <- file.path(base_path, "GSE205013_RAW")
sample_info <- read.csv(file.path(base_path,"Sample_Info.csv"))
# Include features detected in at least this many cells
min_cells_required = 0
# Include cells where at least this many features are detected
# this should match what the paper did -- but it might not
min_features_required = 500

# We are going to just select 2 patients 
# You can try other is you want. ;-) 
sam_names_list <- c("P03","P04") 

# set up loop to go through the list of items you want to work with
# if I was doing this for a general setting you can also set saves for RSD files
# Usually this lets you take your samples that you are intrested in add meta data
# assign them and save them 
# I have comment out a number of "images" you might want to look at if you were doing QC
# Usually that would be done before this loop
assign_name_list <- NULL
for(sample_name in sam_names_list)
{
  #path for each sample
  sample_path <- file.path(general_path,sample_name)
  sample_mtx <- Read10X(data.dir = sample_path) 
  sample_seurat <- CreateSeuratObject(sample_mtx, min.cells = min_cells_required, 
                                      min.features = min_features_required, 
                                      project = sample_name)

  meta <- sample_seurat@meta.data
  dim(meta)
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, pattern = "^MT-")
  sample_seurat[["percent.rb"]] <- PercentageFeatureSet(sample_seurat, pattern = "^RP[SL]")
  
 #Plot before filtering 
#VlnPlot(sample_seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
#    theme(plot.title = element_text(size=10))
  #FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
  #FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

  # only keep those with less than 15% of transcripts coming from mitochondrial genes
  sample_seurat_sub = subset(x = sample_seurat, cells = which(sample_seurat@meta.data$percent.mt < 15))
  
  # remove Cells with >1% of transcripts representing erythroid genes (HBA1, HBA2, HBB, HBM, and ALAS2)
  sample_seurat_sub[["percent.HB"]] <- PercentageFeatureSet(sample_seurat_sub, pattern = "^HBA1|^HBA2|^HBB|^HBM|^ALAS2")
  range(sample_seurat_sub@meta.data$percent.HB)
  sample_seurat_sub = subset(x = sample_seurat_sub, cells = which(sample_seurat_sub@meta.data$percent.HB < 1))
  
  # add meta data 
  sample_seurat_sub@meta.data[["Treatment"]] <- sample_info$Treatment[which(sample_info$Sample==sample_name)]
  sample_seurat_sub@meta.data[["Procedure"]] <- sample_info$Procedure[which(sample_info$Sample==sample_name)]
  sample_seurat_sub@meta.data[["Stage"]] <- sample_info$Stage[which(sample_info$Sample==sample_name)]
  sample_seurat_sub@meta.data[["Moffit"]] <- sample_info$Moffit[which(sample_info$Sample==sample_name)]
  
  # assigns to a new alias 
  # this good way to make objects in loops
  print(paste0(sample_name,"_seurat")) 
  assign_name <- paste0(sample_name,"_seurat")
  assign(x=assign_name, value = sample_seurat_sub)
  
  #compiles the list of all names
  if(!is.null(assign_name_list)){
       assign_name_list <- c(assign_name_list,assign_name)
  }
  if(is.null(assign_name_list)){
    assign_name_list <- assign_name
  }

}
# clean up loop info 
rm(sample_path, sample_mtx, sample_seurat, sample_name, assign_name)

# I am using this to combine all the samples 
# add object to mergeed object set 
# note merge DOES not integrate 
pdac.merged <- merge(get(assign_name_list[1]), y = sapply(assign_name_list[-1],function(x) get(x)),
                       add.cell.ids = sam_names_list, 
                       treatment = sample_info$Treatment,procedure = sample_info$Procedure, 
                       moffit = sample_info$Moffit, stage = sample_info$Stage,
                       project = "PDAC")

# check size of object 
pdac.merged
#check first sample
head(colnames(pdac.merged))
#check last sample 
tail(colnames(pdac.merged))
# tells you what samples are all there 
table(pdac.merged$orig.ident)

# clean up extra 
# you are having trouble you can clean up the extra items.
#rm(P03_seurat, P04_seurat)

# I am showing you RPCA because it is recommended for when a 
#  * A substantial fraction of cells in one dataset have no matching type in the other 
#  * Datasets originate from the same platform (i.e. multiple lanes of 10x genomics) 
#  * There are a large number of datasets or cells to integrate (see INSERT LINK for more tips on integrating large datasets)

# Note for this we required PCA individually 
# While the list of commands is nearly identical, this workflow requires users to run principal components analysis (PCA) individually on each dataset prior to integration. 

### Try this with just 2 samples 
# by patient 
pdac.list <- SplitObject(pdac.merged, split.by = "orig.ident")

# you could also go by something else like a meta data treatment 
# by treatment 
# pdac.list <- SplitObject(pdac.merged, split.by = "Treatment")

# this is used to try and find what we call anchors 
# normalize and identify variable features for each dataset independently
pdac.list <- lapply(X = pdac.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pdac.list)
pdac.list <- lapply(X = pdac.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Perform integration
# this requires that we find achor points / genes
pdac.anchors <- FindIntegrationAnchors(object.list = pdac.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
pdac.combined <- IntegrateData(anchorset = pdac.anchors)

# you need to make sure you know what you are working with so you need to tell it
# Note that original unmodified data still resides in the 'RNA' assay -- if you use that you are using the raw data
# they won't work the same 
# so we tell the defaults what we want to work with.
DefaultAssay(pdac.combined) <- "integrated"

# all of these steps will increase with number of samples 
# I suggest you only try and intregrate 2-4 samples on your computers
# anything above that you need to consider looking a working on a sever
# if you use different samples your results will be different 
# can you match them up from what you have? 

# Run the standard workflow for visualization and clustering
# Note that you start will scaling the data 
pdac.combined <- ScaleData(pdac.combined, verbose = FALSE)
# then you run PCA -- in this case it used 30 - you can change that 
pdac.combined <- RunPCA(pdac.combined, npcs = 30, verbose = FALSE)
# this set migth take a while but it should be printing to your screen
pdac.combined <- RunUMAP(pdac.combined, reduction = "pca", dims = 1:30)
pdac.combined <- FindNeighbors(pdac.combined, reduction = "pca", dims = 1:30)

# Visualization look at the interations 
# using umaps
# these are the new standards you will noticed older data used tsne - 
plot_by_treatment <- DimPlot(pdac.combined, reduction = "umap", group.by = "Treatment")
# if you are looking at it by treatment
plot_by_treatment

plot_by_sample <- DimPlot(pdac.combined, reduction = "umap", group.by = "orig.ident")
# Here you can see how well the integration did as it was run by sample (orig.ident)
plot_by_sample


### beside intergration another major important item is 
# what are my clusters 

# resolution is something that changes a lot based on what you are doing 
# default training res is usually 0.5 
# for some items we like increasing it 
# test and look at what the different resolutions look like for 
# O.3 0.5 and 0.8 
# why would you use one instead of another? 
pdac.combined <- FindClusters(pdac.combined, resolution = 0.2) 


# you often want to look at the "markers" to determine the names of the clusters
# however it is good to know some of what the 
# find all markers of cluster 0
cluster0.markers <- FindMarkers(pdac.combined, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)

# you can look at these using VlnPlot
# this will look at is for a specific marker list 
VlnPlot(pdac.combined, features = c("MMP7", "LYZ"))

# or the "raw" data as count 
# make sure to set it to count 
VlnPlot(pdac.combined, features = c("MMP7", "LYZ"), slot = "counts", log = TRUE)

# we are going to run this with s lower resolution set 
plot_by_cluster_0.2 <- DimPlot(pdac.combined, reduction = "umap", label = TRUE, repel = TRUE)
# pdac.combined
plot_by_cluster_0.2

#pdac.combined <- FindClusters(pdac.combined, resolution = 0.4) 
#plot_by_cluster_0.4 <- DimPlot(pdac.combined, reduction = "umap", label = TRUE, repel = TRUE)
# pdac.combined
#plot_by_cluster_0.4

#note I have overwritten it 
#pdac.combined <- FindClusters(pdac.combined, resolution = 0.8) 
#plot_by_cluster_0.8 <- DimPlot(pdac.combined, reduction = "umap", label = TRUE, repel = TRUE)
# pdac.combined
#plot_by_cluster_0.8

# here we are trying to find markers that are conserved 
# Identify conserved cell type markers
# so for this we switch back to the original data
DefaultAssay(pdac.combined) <- "RNA"

# T/NK "clusrter ID from paper
paper_tnk_markers <-  c("CCL3", "GNLY", "NKG7", "CCL4", "IL7R")
FeaturePlot(pdac.combined, features =paper_tnk_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_tnk_markers)

# Epithelial "clusrter ID from paper
paper_epi_markers <- c("TFF1", "LCN2", "KRT19", "SPINK1", "TFF2")
FeaturePlot(pdac.combined, features =paper_epi_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_epi_markers)

# Endothelial "clusrter ID from paper
paper_endo_markers <- c("SPARCL1", "ACKR1", "CCL14", "VWF", "COL4A1")
FeaturePlot(pdac.combined, features = paper_endo_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_endo_markers)

# Myeliod "clusrter ID from paper 
paper_myeloid_markers <- c("SPP1", "APOE", "S100A9", "IFI30", "HLA-DRA")
FeaturePlot(pdac.combined, features = paper_myeloid_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_myeloid_markers)

# Prolif Epithelial "clusrter ID from paper
paper_pEpi_markers <- c("UBE2C", "AKR1C2", "UBE2S", "UBE2S", "CCNB1")
FeaturePlot(pdac.combined, features = paper_pEpi_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_pEpi_markers)

# Prolif Lymphoid "clusrter ID from paper
paper_pLym_markers <- c("HIST1H4C", "HMGB2", "STMN1", "MKI67", "CENPF")
FeaturePlot(pdac.combined, features = paper_pLym_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_pLym_markers)

# B/Plasma "clusrter ID from paper
paper_B_markers <- c("IGKC", "IGLC2", "IGHG1", "IGLC3", "IGLC1")
FeaturePlot(pdac.combined, features = paper_B_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_B_markers)

# Myeloid "clusrter ID from paper
paper_pMye_markers <-c("MMP9", "CTSK", "ACP5", "CST3", "RGS10")
FeaturePlot(pdac.combined, features =paper_pMye_markers , min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_pMye_markers)

# CAFs "clusrter ID from paper 
paper_cafs_markers <- c("COL3A1", "COL1A2", "COL1A1", "SFRP2", "DCN")
FeaturePlot(pdac.combined, features = paper_cafs_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_cafs_markers)

# Mast "clusrter ID from paper 
paper_mast_markers <- c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT")
FeaturePlot(pdac.combined, features = paper_mast_markers, min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_mast_markers)


# Make the assignments 
# sometimes this part is confusing 
# here we are "guided" by the paper 

# using the clust
genes_types <- read.csv(file.path(base_path,"genes_types.csv"))

genes_types
# if you were doing this for your own data 
# you would want to look at the the clusters items like 

# however other times you will probably want to be guided by comparing the different clusters using 
# cluster markers 

# this compares against ALL clusters 
# looks like the cluster 0
for(clusterID in levels(pdac.combined@meta.data$integrated_snn_res.0.2))
{
  print(clusterID)
  clusterX.markers <- FindMarkers(pdac.combined, ident.1 = as.numeric(clusterID), min.pct = 0.25)
  
  assign_name <- paste0("cluster",clusterID ,".markers")
  print(assign_name)
  assign(x=assign_name, value = clusterX.markers)
  print(head(clusterX.markers))
  
  # match 
  rownames(clusterX.markers)
  intGenes <- intersect(rownames(clusterX.markers),genes_types$Gene)
  
  genes_types$cluster <- NA
  for(gene in intGenes)
  {
    #print(gene)
    # Get the rank 
    # you could do this differently but we are going for simple
    genes_types$cluster[genes_types$Gene==gene] <- which(rownames(clusterX.markers) == gene)
  }
  print(table(genes_types$Celltype[!is.na(genes_types$cluster)]))
  colnames(genes_types)[colnames(genes_types)=="cluster"] <- assign_name
  
}
# clean up loop info 
rm(clusterX.markers, clusterID, assign_name)

# now it printed out but sometimes you want to look at it more 

# looking for the most that are matched with the genes of interst
# looking for those with the best rank 
# some clusters should have similaritlies 

# cluster0.markers
table(genes_types$Celltype[!is.na(genes_types$cluster0.markers)])
cluster0Tab <- genes_types[!is.na(genes_types$cluster0.markers),]
cluster0Tab <- cluster0Tab[order(cluster0Tab$cluster0.markers),]
cluster0Tab
# calling cluster
# 0 -- T/NK 

VlnPlot(pdac.combined, features = paper_tnk_markers)
FeaturePlot(pdac.combined, features =paper_tnk_markers , min.cutoff = "q9")
plot_by_cluster_0.2

# cluster1.markers
table(genes_types$Celltype[!is.na(genes_types$cluster1.markers)])
cluster1Tab <- genes_types[!is.na(genes_types$cluster1.markers),]
cluster1Tab <- cluster1Tab[order(cluster1Tab$cluster1.markers),]
cluster1Tab
# unclear but either myeloid or Epithelial
# therefore look at graphs 
plot_by_cluster_0.2
VlnPlot(pdac.combined, features = paper_myeloid_markers)
FeaturePlot(pdac.combined, features =paper_myeloid_markers , min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_epi_markers)
FeaturePlot(pdac.combined, features =paper_epi_markers , min.cutoff = "q9")
# 1 -- Epithelial 

#cluster2.markers
table(genes_types$Celltype[!is.na(genes_types$cluster2.markers)])
cluster2Tab <- genes_types[!is.na(genes_types$cluster2.markers),]
cluster2Tab <- cluster2Tab[order(cluster2Tab$cluster2.markers),]
cluster2Tab
# either Myeliod or pMyeloid
plot_by_cluster_0.2
VlnPlot(pdac.combined, features = paper_myeloid_markers)
FeaturePlot(pdac.combined, features =paper_myeloid_markers , min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_pMye_markers)
FeaturePlot(pdac.combined, features =paper_pMye_markers , min.cutoff = "q9")
# 2 -- pMyeloid

#cluster3.markers
table(genes_types$Celltype[!is.na(genes_types$cluster3.markers)])
cluster3Tab <- genes_types[!is.na(genes_types$cluster3.markers),]
cluster3Tab <- cluster3Tab[order(cluster3Tab$cluster3.markers),]
cluster3Tab
# 3 -- CAFS
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features =paper_cafs_markers , min.cutoff = "q9")
VlnPlot(pdac.combined, features = paper_cafs_markers)

#cluster4.markers
table(genes_types$Celltype[!is.na(genes_types$cluster4.markers)])
cluster4Tab <- genes_types[!is.na(genes_types$cluster4.markers),]
cluster4Tab <- cluster4Tab[order(cluster4Tab$cluster4.markers),]
cluster4Tab
# 4 -- Epithelial
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features =paper_epi_markers , min.cutoff = "q9")

#cluster5.markers
table(genes_types$Celltype[!is.na(genes_types$cluster5.markers)])
cluster5Tab <- genes_types[!is.na(genes_types$cluster5.markers),]
cluster5Tab <- cluster5Tab[order(cluster5Tab$cluster5.markers),]
cluster5Tab
# 5 -- Epithelial 
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features =paper_epi_markers , min.cutoff = "q9")


#cluster6.markers
table(genes_types$Celltype[!is.na(genes_types$cluster6.markers)])
# Myeliod 
cluster6Tab <- genes_types[!is.na(genes_types$cluster6.markers),]
cluster6Tab <- cluster6Tab[order(cluster6Tab$cluster6.markers),]
cluster6Tab
# 6 -- Myeliod 
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features =paper_myeloid_markers , min.cutoff = "q9")
FeaturePlot(pdac.combined, features =paper_pMye_markers , min.cutoff = "q9")

#cluster7.markers
table(genes_types$Celltype[!is.na(genes_types$cluster7.markers)])
# Prolif cluster -- 
cluster7Tab <- genes_types[!is.na(genes_types$cluster7.markers),]
cluster7Tab <- cluster7Tab[order(cluster7Tab$cluster7.markers),]
cluster7Tab
# 7 -- Prolif Epithelial  
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features = paper_pLym_markers , min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_pEpi_markers , min.cutoff = "q9")

#cluster8.markers
table(genes_types$Celltype[!is.na(genes_types$cluster8.markers)])
# Endothelial
cluster8Tab <- genes_types[!is.na(genes_types$cluster8.markers),]
cluster8Tab <- cluster8Tab[order(cluster8Tab$cluster8.markers),]
cluster8Tab
# 8 -- Endothelial
plot_by_cluster_0.2
FeaturePlot(pdac.combined, features = paper_endo_markers , min.cutoff = "q9")

#cluster9.markers
table(genes_types$Celltype[!is.na(genes_types$cluster9.markers)])
# unclear 
cluster9Tab <- genes_types[!is.na(genes_types$cluster9.markers),]
cluster9Tab <- cluster9Tab[order(cluster9Tab$cluster9.markers),]
cluster9Tab
# 9 -- Epithelial
plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_epi_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_epi_markers[1], min.cutoff = "q9")

#cluster10.markers
table(genes_types$Celltype[!is.na(genes_types$cluster10.markers)])
# Mast  
cluster10Tab <- genes_types[!is.na(genes_types$cluster10.markers),]
cluster10Tab <- cluster10Tab[order(cluster10Tab$cluster10.markers),]
cluster10Tab
# 10 -- Mast
plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_mast_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_mast_markers[1], min.cutoff = "q9")


#cluster11.markers
table(genes_types$Celltype[!is.na(genes_types$cluster11.markers)])
# unclear -- Endothelial or CAFs
cluster11Tab <- genes_types[!is.na(genes_types$cluster11.markers),]
cluster11Tab <- cluster11Tab[order(cluster11Tab$cluster11.markers),]
cluster11Tab
# 11 -- CAFs

plot_by_cluster_0.2
FeaturePlot(pdac.combined, features = paper_endo_markers, min.cutoff = "q9")
#FeaturePlot(pdac.combined, features = paper_endo_markers[2], min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_cafs_markers, min.cutoff = "q9")
#FeaturePlot(pdac.combined, features = paper_cafs_markers[1], min.cutoff = "q9")


#cluster12.markers
table(genes_types$Celltype[!is.na(genes_types$cluster12.markers)])
# Prolif Lymphoid
cluster12Tab <- genes_types[!is.na(genes_types$cluster12.markers),]
cluster12Tab <- cluster12Tab[order(cluster12Tab$cluster12.markers),]
cluster12Tab
# 12 -- Prolif Lymphoid 
plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_pLym_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_pLym_markers[1], min.cutoff = "q9")

plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_endo_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_endo_markers[1], min.cutoff = "q9")


#cluster13.markers
table(genes_types$Celltype[!is.na(genes_types$cluster13.markers)])
# B/Plasma
cluster13Tab <- genes_types[!is.na(genes_types$cluster13.markers),]
cluster13Tab <- cluster13Tab[order(cluster13Tab$cluster13.markers),]
cluster13Tab
# 13 -- B/Plasma
plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_B_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_B_markers[1], min.cutoff = "q9")


#cluster14.markers
table(genes_types$Celltype[!is.na(genes_types$cluster14.markers)])
# unclear
cluster14Tab <- genes_types[!is.na(genes_types$cluster14.markers),]
cluster14Tab <- cluster14Tab[order(cluster14Tab$cluster14.markers),]
cluster14Tab
# 14 -- CAFs
plot_by_cluster_0.2
#FeaturePlot(pdac.combined, features = paper_B_markers, min.cutoff = "q9")
FeaturePlot(pdac.combined, features = paper_cafs_markers[1], min.cutoff = "q9")

# Based on this Estimate are 

# 0 -- T/NK 
# 1 -- Epithelial  -- Yes 
# 2 -- Myeloid -- Yes -- Prolif maybe
# 3 -- CAFS -- Yes
# 4 -- Epithelial  -- Yes
# 5 -- Epithelial  -- Yes
# 6 -- Myeliod -- Yes 
# 7 -- Prolif Epithelial -- Yes
# 8 -- Endothelial -- Yes
# 9 -- Epithelial  -- Yes
# 10 -- Mast -- Yes
# 11 -- CAFs -- Yes
# 12 -- Prolif Lymphoid -- Yes 
# 13 -- B/Plasma -- Yes
# 14 -- CAFs -- Yes

## there can be some agruement on these assignments for example 
## I would suggest that 2 pMyeloid might be just Myeloid 
## but for this we will say it is its own group

#Not we combined some at this point
# note this will change a number of items so make sure you want to do
# I general like to save before I do the names or combines
pdac.combined <- RenameIdents(pdac.combined, 
                              `0` = "T/NK", 
                              `1` = "Epithelial", 
                              `2` = "Myeliod",
                              `3` = "CAFs", 
                              `4` = "Epithelial", 
                              `5` = "Epithelial", 
                              `6` = "Myeliod", 
                              `7` = "Prolif Epithelial", 
                              `8` = "Endothelial", 
                              `9` = "Epithelial",
                              `10` = "Mast", 
                              `11` = "CAFs", 
                              `12` = "Prolif Lymphoid", 
                              `13` = "B/Plasma", 
                              `14` = "CAFs")

# we can now do "guided" clustering and assignments 
# it should be noted that this works "mostly"
# it also should be noted this is why it is "good" to have some genes in mind

# Now these are only 2 of the samples used in the sample 
# and I randomly grabbed these samples 
# and we are using a different integration 
# but we are still going to try and name them 
# based on the "Graphs" 
# not if you want to discuss alternative options, 
# or if you decided to use different starter samples -- you will not get the same results. 



# look at new plot of this 
DimPlot(pdac.combined, label = TRUE)

### Examine cluster differences using these 

FeaturePlot(pdac.combined, features = paper_pMye_markers[1], min.cutoff = "q9")