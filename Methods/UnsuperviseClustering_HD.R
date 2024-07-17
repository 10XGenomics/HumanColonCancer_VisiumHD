# Unsupervised Clustering

library(tidyverse,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(dplyr,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(SeuratObject,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(Seurat,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(BPCells,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(Matrix,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(scattermore,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(ggplot2,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")

options(future.globals.maxSize = 1000 * 1024^2)

source("/mnt/home/juanpablo.romerorioj/Builds/HumanColonCancer_VisiumHD/Methods/AuxFunctions.R")


SampleInfo<-data.frame(Sample=c("P1CRC","P2CRC","P5CRC"),
                       H5=c("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P2CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P5CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"),
                       Deconv=c("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/DeconvolutionResults_P1CRC.csv.gz",
                                "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P2CRC/DeconvolutionResults_P2CRC.csv.gz",
                                "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P5CRC/DeconvolutionResults_P5CRC.csv.gz"))

MD<-readRDS("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/MetaDataCombined.rds")
MetaData<-vector("list",length = nrow(SampleInfo))
Matrices<-vector("list",length = nrow(SampleInfo))

for(jj in 1:nrow(SampleInfo))
{
  PathX<-paste0("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/",SampleInfo$Sample[jj])
  
  #MatData <- open_matrix_10x_hdf5(path = SampleInfo$H5[jj])
  #write_matrix_dir(mat = MatData,dir = PathX,overwrite = TRUE)
  
  mat <- open_matrix_dir(dir = PathX)
  colnames(mat)<-paste0(colnames(mat),"_",SampleInfo$Sample[jj])
  
  Genes<-read.delim(file="/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header = F)
  rownames(mat)<-make.unique(Genes$V2[match(rownames(mat),Genes$V1)])
  
  MetaData[[jj]] <- MD[MD$Patient==SampleInfo$Sample[jj],]
  Matrices[[jj]] <- mat
  
  
}

names(Matrices)<-c("P1CRC","P2CRC","P5CRC")
names(MetaData)<-c("P1CRC","P2CRC","P5CRC")

MetaData <- Reduce(rbind, MetaData)

merged.object <- CreateSeuratObject(counts = Matrices, meta.data = MetaData)
merged.object<-JoinLayers(merged.object)

AllDeconvolution<-vector("list",length=nrow(SampleInfo))
for(jj in 1:nrow(SampleInfo))
{
  DecTmp<-read.delim(SampleInfo$Deconv[jj],sep=",") %>% na.omit()
  DecTmp$barcode<-paste0(DecTmp$barcode,"_",SampleInfo$Sample[jj])
  
  AllDeconvolution[[jj]]<-DecTmp
  
}

AllDeconvolution <- Reduce(rbind, AllDeconvolution)


merged.object$DeconvolutionClass<-AllDeconvolution$DeconvolutionClass[match(colnames(merged.object),AllDeconvolution$barcode)]
merged.object$DeconvolutionLabel1<-AllDeconvolution$DeconvolutionLabel1[match(colnames(merged.object),AllDeconvolution$barcode)]
merged.object$DeconvolutionLabel2<-AllDeconvolution$DeconvolutionLabel2[match(colnames(merged.object),AllDeconvolution$barcode)]


saveRDS(object = merged.object,file = "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/CombinedObject_BPCells.rds")

merged.object<-readRDS("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/CombinedObject_BPCells.rds")

#

merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object)
merged.object <- SketchData(object = merged.object,ncells = 240000,method = "LeverageScore",sketched.assay = "sketch")

DefaultAssay(merged.object) <- "sketch"

merged.object <- FindVariableFeatures(merged.object)
merged.object <- ScaleData(merged.object)
merged.object <- RunPCA(merged.object)
merged.object <- FindNeighbors(merged.object, dims = 1:20)
merged.object <- FindClusters(merged.object, resolution = 0.8)
merged.object <- RunUMAP(merged.object, dims = 1:20,return.model = T)

IdentsLvl1<-c("Tumor","Intestinal Epithelial","Endothelial","Smooth Muscle","Tumor","Tumor","T cells",
              "Fibroblast","B cells","Myeloid","Fibroblast","Tumor","Unknown","Intestinal Epithelial","Fibroblast",
              "Unknown","Myeloid","B cells","Tumor","Tumor","Neuronal","Tumor","Unknown")

merged.object$Level1<-IdentsLvl1[as.numeric(as.vector(merged.object$seurat_clusters))+1]
merged.object$Level1<-factor(merged.object$Level1,levels=sort(unique(merged.object$Level1)))

merged.object<-SetIdent(merged.object,value = "Level1")

Mks<-FindAllMarkers(merged.object,logfc.threshold = 0.2,min.diff.pct = 0.2,only.pos = T)
Mks<-Mks[Mks$p_val_adj<0.05,]

# Iterative sub-clustering for level2
MarkersSubcluster<-vector("list",length=length(IdentsLvl1))
names(MarkersSubcluster)<-IdentsLvl1

MetaDataSubClusters<-vector("list",length=length(IdentsLvl1))
names(MetaDataSubClusters)<-IdentsLvl1

for(ClusterID in levels(merged.object$Level1))
{
  message(ClusterID)
  Subset<-subset(merged.object,idents = ClusterID)
  Subset <- FindVariableFeatures(Subset)
  Subset <- ScaleData(Subset)
  Subset <- RunPCA(Subset)
  Subset <- FindNeighbors(Subset, dims = 1:25)
  Subset <- FindClusters(Subset, resolution = 0.1)
  Subset$Level2 <- paste0(gsub(" ","",ClusterID),"_",Subset$seurat_clusters)
  Subset<-SetIdent(Subset,value="Level2")
  
  MetaDataSubClusters[[ClusterID]]<-Subset@meta.data
  
  # Get Markers
  SubMks<-FindAllMarkers(Subset,min.diff.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
  SubMks<-SubMks[SubMks$p_val_adj<0.05,]
  
  MarkersSubcluster[[ClusterID]]<-SubMks
  
}

# Convert Markers and MD lists to data.frames
MarkersSubcluster<-do.call(rbind,MarkersSubcluster)
MetaDataSubClusters<-do.call(rbind,MetaDataSubClusters)

MetaDataSubClusters$Barcode<-sapply(strsplit(rownames(MetaDataSubClusters),"[.]"),function(X){return(X[2])})

IdentsLevel2<-sort(unique(MetaDataSubClusters$Level2))

MetaDataSubClusters$Level2<-factor(MetaDataSubClusters$Level2,levels = IdentsLevel2)
MarkersSubcluster$cluster<-factor(MarkersSubcluster$cluster,levels = IdentsLevel2)

merged.object$Level2<-MetaDataSubClusters$Level2[match(colnames(merged.object),MetaDataSubClusters$Barcode)]
#merged.object$Level2<-factor(merged.object$Level2,levels=sort(unique(merged.object$Level2)))

# Create color palettes
ColsL1<-paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(merged.object$Level1))]
names(ColsL1)<-levels(MetaData$Level1)

ColsL2<-paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(merged.object$Level2))]
names(ColsL2)<-levels(merged.object$Level2)


P1<-DimPlot(merged.object,group.by = "Level1",label = T,label.size = 4,cols=ColsL1)
P2<-DimPlot(merged.object,group.by = "Level2",label = T,label.size = 4,cols=ColsL2)

P1+P2

# Project to Full dataset

merged.object <- Seurat:::ProjectData(
  object = merged.object,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(L1 = "Level1",L2 = "Level2")
)
