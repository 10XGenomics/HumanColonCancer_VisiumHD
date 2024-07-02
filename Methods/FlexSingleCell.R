#### Visium HD Manuscript
### 

## METHOD: Single Cell Flex Analysis

## Load required Packages
library(Seurat)
library(Azimuth)
library(BPCells)

# Load Auxiliary Functions
source("~/Sections/AuxFunctions.R")

# Read Aggr'd dataset as a Seurat V5 object
FlexOutPath <- "~/AggrOutput/outs" # Path to cellranger aggr output folder
ColonFlex.data <- open_matrix_10x_hdf5(path = paste0(FlexOutPath,"/count/filtered_feature_bc_matrix.h5")) 
write_matrix_dir(mat = ColonFlex.data,dir = '~/Outputs/FlexSeurat/') 
Flex.mat <- open_matrix_dir(dir = "~/Outputs/FlexSeurat/")
Flex.mat <- ConvertEnsembleToSymbol(mat = Flex.mat, species = "human") 

# Read aggregation_csv.csv file to be used as MetaData ( Patient, etc)
MetaData<-read.csv(paste0(FlexOutPath,"aggregation_csv.csv"))
MetaData$Patient<-sapply(strsplit(MetaData$sample_id,"_"),function(X){return(X[6])})
MetaData$BC<-sapply(strsplit(MetaData$sample_id,"_"),function(X){return(X[7])})
MetaData$Condition<-gsub("P[0-9]","",MetaData$Patient)
Index<-as.numeric(sapply(strsplit(colnames(Flex.mat),"-"),function(X){return(X[2])}))
MetaData<-MetaData[Index,3:4]
MetaData$Barcode<-colnames(Flex.mat)
rownames(MetaData)<-MetaData$Barcode

# Cell Ranger Results
UMAP<-read.csv(paste0(FlexOutPath,"/analysis_csv/umap/gene_expression_2_components/projection.csv"))
MetaData<-cbind(MetaData,UMAP[match(MetaData$Barcode,UMAP$Barcode),2:3])

# Create Seurat Object
ColonCancer_Flex<-CreateSeuratObject(Flex.mat,meta.data = MetaData)

# Add % MT
ColonCancer_Flex[["MT.percent"]] <- PercentageFeatureSet(ColonCancer_Flex, pattern = "^MT-")

# Add UMAP projection
ColonCancer_Flex[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(MetaData[,c("UMAP.1","UMAP.2")]), key = "UMAP_", assay = DefaultAssay(ColonCancer_Flex))

# UMI and Gene Threshold (by removing 5%)
UMI_TH<-quantile(ColonCancer_Flex$nCount_RNA,c(0.025,0.975))
Gene_TH<-quantile(ColonCancer_Flex$nFeature_RNA,c(0.025,0.975))

# Add variable with QC filter status
ColonCancer_Flex$QCFilter<-ifelse(ColonCancer_Flex$MT.percent < 25 & 
                                    ColonCancer_Flex$nCount_RNA > UMI_TH[1] & ColonCancer_Flex$nCount_RNA < UMI_TH[2] & 
                                    ColonCancer_Flex$nFeature_RNA > Gene_TH[1] & ColonCancer_Flex$nFeature_RNA < Gene_TH[2],"Keep","Remove")

# Remove bcs that failed QC
ColonCancer_Flex<-subset(ColonCancer_Flex,cells=colnames(ColonCancer_Flex)[ColonCancer_Flex$QCFilter=="Keep"])

# Seurat V5 Processing (We will sketch and then project) at ~ 15% of the full data
ColonCancer_Flex <- NormalizeData(ColonCancer_Flex)
ColonCancer_Flex <- FindVariableFeatures(ColonCancer_Flex)
ColonCancer_Flex <- SketchData(object = ColonCancer_Flex,ncells = 37000,
                               method = "LeverageScore",sketched.assay = "sketch")

DefaultAssay(ColonCancer_Flex) <- "sketch"
ColonCancer_Flex <- FindVariableFeatures(ColonCancer_Flex)
ColonCancer_Flex <- ScaleData(ColonCancer_Flex)
ColonCancer_Flex <- RunPCA(ColonCancer_Flex)
ElbowPlot(ColonCancer_Flex,ndims=40)
ColonCancer_Flex <- FindNeighbors(ColonCancer_Flex, dims = 1:25)
ColonCancer_Flex <- FindClusters(ColonCancer_Flex, resolution = 0.6)

# Find Markers
Mks<-FindAllMarkers(ColonCancer_Flex,min.diff.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
Mks<-Mks[Mks$p_val_adj<0.05,]

# Rename Clusters (Manual Annotation)
NewClusterIDs <- c("Tumor","Smooth Muscle","Tumor","Myeloid","Fibroblast","Fibroblast",
                   "Intestinal Epithelial","Endothelial","Smooth Muscle","T cells","Intestinal Epithelial",
                   "B cells","T cells","B cells","Fibroblast","Neuronal","Fibroblast","Myeloid",
                   "Myeloid","Neuronal","Myeloid","Endothelial","Neuronal","Tumor",
                   "Intestinal Epithelial","Intestinal Epithelial")

names(NewClusterIDs) <- levels(ColonCancer_Flex)

ColonCancer_Flex <- RenameIdents(ColonCancer_Flex, NewClusterIDs)
ColonCancer_Flex@meta.data$Level1<-NewClusterIDs[match(ColonCancer_Flex$seurat_clusters,names(NewClusterIDs))]

# Iterative sub-clustering for level2
MarkersSubcluster<-vector("list",length=length(NewClusterIDs))
names(MarkersSubcluster)<-NewClusterIDs

MetaDataSubClusters<-vector("list",length=length(NewClusterIDs))
names(MetaDataSubClusters)<-NewClusterIDs

for(ClusterID in unique(NewClusterIDs))
{
  message(ClusterID)
  Subset<-subset(ColonCancer_Flex,idents = ClusterID)
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

# Rename Clusters
Level2Clusters<-data.frame(ID=sort(unique(MetaDataSubClusters$Level2)),
Level2Clusters<-data.frame(ID=sort(unique(MetaDataSubClusters$Level2)),
                           Label=c("Plasma","Mature B","Plasma","Proliferating Immune II","Endothelial","Endothelial","Lymphatic Endothelial","CAF",
                                   "Fibroblast","Pericytes","Fibroblast","Myofibroblast","CAF","Myofibroblast","Fibroblast","CAF","Goblet","Enterocyte",
                                   "Goblet","Tuft","Macrophage","Neutrophil","Mast","Macrophage","Macrophage","mRegDC","pDC","Enteric Glial","Enteric Glial",
                                   "Neuroendocrine","Enteric Glial","Adipocyte","Epithelial","SM Stress Response","Smooth Muscle","vSM","Unknown III (SM)",
                                   "SM Stress Response","Adipocyte","CD4 T cell","CD8 Cytotoxic T cell","CD4 T cell","CD8 Cytotoxic T cell",
                                   "CD4 T cell","Tumor I","Tumor III","Tumor II","Tumor V","Tumor I","Proliferating Immune II"))

MetaDataSubClusters$Level2<-Level2Clusters$Label[match(MetaDataSubClusters$Level2,Level2Clusters$ID)]

ColonCancer_Flex$Level2<-MetaDataSubClusters$Level2[match(colnames(ColonCancer_Flex),MetaDataSubClusters$Barcode)]
ColonCancer_Flex<-SetIdent(ColonCancer_Flex,value = "Level2")

# Project to Full dataset

ColonCancer_Flex <- Seurat:::ProjectData(
  object = ColonCancer_Flex,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(ProjectAll_L1 = "Level1",ProjectAll_L2 = "Level2")
)

ColonCancer_Flex$Level1<-factor(ColonCancer_Flex$Level1,levels = sort(unique(ColonCancer_Flex$Level1)))
ColonCancer_Flex$Level2<-factor(ColonCancer_Flex$Level2,levels = sort(unique(ColonCancer_Flex$Level2)))

# Save processed Seurat Object & MetaData
saveRDS(ColonCancer_Flex,file='~/Outputs/Flex/FlexSeuratV5.rds')
saveRDS(ColonCancer_Flex@meta.data,file='~/Outputs/Flex/FlexSeuratV5_MetaData.rds')