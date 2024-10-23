#### Visium HD Manuscript
## Supplemental Figure 6

# load packages
library(Seurat)
library(scattermore)
library(tidyverse)
library(data.table)
library(wesanderson)
library(patchwork)
library(RColorBrewer)
library(furrr)
library(paletteer)
library(arrow)
library(pheatmap)
library(RColorBrewer)
library(distances)
library(enrichR)

# load aux functions
source("~/Methods/AuxFunctions.R")


options(future.globals.maxSize = 1000 * 1024^2)

SampleInfo<-data.frame(Sample=c("P1CRC","P2CRC","P5CRC"),
                       H5=c("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "~/VisiumHD/PatientCRC2/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "~/VisiumHD/PatientCRC5/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"),
                       Deconv=c("~/MetaData/DeconvolutionResults_P1CRC.csv.gz",
                                "~/MetaData/DeconvolutionResults_P2CRC.csv.gz",
                                "~/MetaData/DeconvolutionResults_P5CRC.csv.gz"))

MetaData<-vector("list",length = nrow(SampleInfo))
Matrices<-vector("list",length = nrow(SampleInfo))

for(jj in 1:nrow(SampleInfo))
{

  MatData <- open_matrix_10x_hdf5(path = SampleInfo$H5[jj])
  write_matrix_dir(mat = MatData,dir = getwd(),overwrite = TRUE)
  
  mat <- open_matrix_dir(dir = PathX)
  colnames(mat)<-paste0(colnames(mat),"_",SampleInfo$Sample[jj])
  
  Genes<-read.delim(file="~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header = F)
  rownames(mat)<-make.unique(Genes$V2[match(rownames(mat),Genes$V1)])
  
  Matrices[[jj]] <- mat
  
}

names(Matrices)<-c("P1CRC","P2CRC","P5CRC")

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

# Sketch Full Seurat Object
merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object)
merged.object <- SketchData(object = merged.object,ncells = 240000,method = "LeverageScore",sketched.assay = "sketch")

# Analysis on the Sketch subset
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
  SubMks<-FindAllMarkers(Subset,min.diff.pct = 0.1,logfc.threshold = 0.1,only.pos = T)
  SubMks<-SubMks[SubMks$p_val_adj<0.05,]
  
  MarkersSubcluster[[ClusterID]]<-SubMks
  
}

# Convert Markers and MD lists to data.frames
MarkersSubcluster<-do.call(rbind,MarkersSubcluster)
MetaDataSubClusters<-do.call(rbind,MetaDataSubClusters)

MetaDataSubClusters$Barcode<-sapply(strsplit(rownames(MetaDataSubClusters),"[.]"),function(X){return(X[2])})

IdentsLevel2<-c("Bcells_0","Bcells_1","Bcells_2","Bcells_3","Bcells_4",
                "Endothelial_0","Endothelial_1","Endothelial_2","Endothelial_3",
                "Fibroblast_0","Fibroblast_1","Fibroblast_2",
                "IntestinalEpithelial_0","IntestinalEpithelial_1","IntestinalEpithelial_2","IntestinalEpithelial_3","IntestinalEpithelial_4","IntestinalEpithelial_5",
                "Myeloid_0","Myeloid_1","Myeloid_2",
                "Neuronal_0","Neuronal_1",
                "SmoothMuscle_0","SmoothMuscle_1","SmoothMuscle_2","SmoothMuscle_3","SmoothMuscle_4",
                "Tcells_0","Tcells_1","Tcells_2","Tcells_3",
                "Tumor_0","Tumor_1","Tumor_2","Tumor_3",
                "Unknown_0","Unknown_1","Unknown_2","Unknown_3","Unknown_4","Unknown_5","Unknown_6")

MetaDataSubClusters$Level2<-factor(MetaDataSubClusters$Level2,levels = IdentsLevel2)
MarkersSubcluster$cluster<-factor(MarkersSubcluster$cluster,levels = IdentsLevel2)

merged.object$Level2<-MetaDataSubClusters$Level2[match(colnames(merged.object),MetaDataSubClusters$Barcode)]

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


MetaData.merged$L2<-factor(MetaData.merged$L2,levels=sort(unique(MetaData.merged$L2)))
LabelsL2<-MetaData.merged %>% group_by(L2) %>% summarise(X=median(fullumap_1),Y=median(fullumap_2))
LabelsL1<-MetaData.merged %>% group_by(L1) %>% summarise(X=median(fullumap_1),Y=median(fullumap_2))


# Create color palettes
MetaData.merged$L1<-factor(MetaData.merged$L1,levels=sort(unique(MetaData.merged$L1)))
ColsL1<-paletteer::paletteer_d("ggsci::category10_d3")[1:length(levels(MetaData.merged$L1))]
names(ColsL1)<-levels(MetaData.merged$L1)
ColsL1["T cells"]<-"#8C564BFF"
ColsL1["Neuronal"]<-"#7F7F7FFF"

MetaData.merged$L2<-factor(MetaData.merged$L2,levels=sort(unique(MetaData.merged$L2)))
ColsL2<-paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(MetaData.merged$L2))]
names(ColsL2)<-levels(MetaData.merged$L2)


# Proportion Heatmap using only Deconvolution singlet bins
ix_Singlet<-MetaData.merged$DeconvolutionClass=="singlet"
ConfusionMatrix<-prop.table(table(MetaData.merged$DeconvolutionLabel1[ix_Singlet],MetaData.merged$L2[ix_Singlet]),margin = 1)*100

ConfusionMatrix[which(is.na(ConfusionMatrix))]<-0

AnnotRow<-data.frame(CellType=rownames(ConfusionMatrix))
rownames(AnnotRow)<-rownames(ConfusionMatrix)

AnnotCol<-data.frame(Level1=sapply(strsplit(colnames(ConfusionMatrix),"-"),function(X){return(X[1])}),
                     Level2=colnames(ConfusionMatrix))
rownames(AnnotCol)<-colnames(ConfusionMatrix)

ColorsAnnotation<-list(CellType=ColorPalette(),Level1=ColsL1,Level2=ColorsExtra())
names(ColorsAnnotation$Level1)<-gsub(" ","",names(ColorsAnnotation$Level1))

ColorHM<-colorRampPalette(c("white","red"))(25)

pheatmap(ConfusionMatrix,annotation_row = AnnotRow,annotation_col = AnnotCol,color=ColorHM,
         annotation_legend = FALSE,border_color=NA,annotation_colors = ColorsAnnotation,cluster_cols = F)




# SpatialPlots (P1CRC as an example)
Patient<-"P1CRC"

BarcodeData<-split(MetaData.merged,MetaData.merged$Patient)
BarcodeData<-BarcodeData[[Patient]]

DF<-GenerateSampleData("~/VisiumHD/PatientCRC1/outs/")$bcs

DF$Level1<-BarcodeData$L1[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]
DF$Level1<-factor(DF$Level1,levels = sort(unique(DF$Level1)))

DF$Level2<-BarcodeData$L2[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]
DF$Level2<-factor(DF$Level2,levels = sort(unique(DF$Level2)))

DF$DeconvolutionClass<-BarcodeData$DeconvolutionClass[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]
DF$DeconvolutionL1<-BarcodeData$DeconvolutionLabel1[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]
DF$DeconvolutionL1<-factor(DF$DeconvolutionL1,levels = sort(unique(DF$DeconvolutionL1)))

DF$nCount_RNA<-BarcodeData$nCount_RNA[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]
DF$nFeature_RNA<-BarcodeData$nFeature_RNA[match(paste0(DF$barcode,paste0("_",Patient)),rownames(BarcodeData))]


L1_A<-DF  %>% filter(tissue==1) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Level1)) +
  geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=ColsL1)+
  labs(color="Level 1")+ggtitle("")+NoLegend()


L2_A<-DF  %>% filter(tissue==1) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Level2)) +
  geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=ColorsExtra())+
  labs(color="Level 2")+ggtitle("")+NoLegend()


L2_D<-DF  %>% filter(tissue==1 & DeconvolutionClass == "singlet") %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=DeconvolutionL1)) +
  geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=ColorPalette())+
  labs(color="Deconvolution")+ggtitle("")+NoLegend()

L1_A+L2_A+L2_D

