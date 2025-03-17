#### Visium HD Manuscript
## Supplemental Figure 9

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

# load aux functions
source("~/Methods/AuxFunctions.R")



SampleData<-data.frame(Patient = c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/VisiumHD/PatientCRC1/outs/","~/VisiumHD/PatientCRC2/outs/","~/VisiumHD/PatientCRC5/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Periphery/P1CRC_TME_Barcodes.rds","~/Outputs/Periphery/P2CRC_TME_Barcodes.rds","~/Outputs/Periphery/P5CRC_TME_Barcodes.rds"),
                       Tumor = c("Tumor II","Tumor III","Tumor IV"))

#Colocalization Macrophage and T cells
# We will use PCRC1 as an example
index<-1

TumorCluster<-SampleData$Tumor[index]
Sample_path <- SampleData$PathSR[index]

# Create data.frame and add Deconvolution Info
bcDF<-GenerateSampleData(Sample_path)$bcs
DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
bcDF<-AddDeconvolutionInfo(bcDF,DeconvolutionHD,AddWeights=FALSE)

# Detect Tumor microenviroment (barcodes within 50 microns of tumor)
# Either load the results from Figure 4 or compute again
if(file.exists(SampleData$PathPeriphery[index]))
{
  PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
  
}else{
  
  PeripheryBCs<-SelectPeripheryDiscrete(bcs = bcDF,CellType = TumorCluster,distance=50,PATH = SampleData$PathSR[index])
}

bcDF$Periphery<-NA
bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==TumorCluster]<-"Tumor"

# Keep only singlet bins to run Macrophage analysis
bcDF <- bcDF %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()

# Create Seurat Object
SeuratObj <- Load10X_Spatial(Sample_path, bin.size = 8)

# Subset only to Macrophages in the TME
SeuratObj<-subset(SeuratObj,cells=bcDF$barcode[bcDF$DeconvolutionLabel1=="Macrophage" & bcDF$Periphery=="50 micron"])
SeuratObj$CellType<-bcDF$DeconvolutionLabel1[match(colnames(SeuratObj),bcDF$barcode)]

# Standard Seurat Processing
SeuratObj<-NormalizeData(SeuratObj)
SeuratObj<-FindVariableFeatures(SeuratObj)
SeuratObj<-ScaleData(SeuratObj)
SeuratObj<-RunPCA(SeuratObj)
SeuratObj<-FindNeighbors(SeuratObj,dims=1:10)
SeuratObj<-FindClusters(SeuratObj,resolution=0.2)

# Get markers
CtMarkers<-FindAllMarkers(SeuratObj,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T)
CtMarkers<-CtMarkers[CtMarkers$p_val_adj<0.05,]

# Rename Clusters
SPP1_Cluster<-as.vector(CtMarkers$cluster[match("SPP1",CtMarkers$gene)])
SeuratObj$Group<-ifelse(SeuratObj$seurat_clusters==SPP1_Cluster,"SPP1+","SELENOP+")

# Split barcodes by group
TypeInfo<-split(rownames(SeuratObj@meta.data),SeuratObj$Group)

# Re compute 
bcDF<-GenerateSampleData(Sample_path)$bcs
bcDF<-AddDeconvolutionInfo(bcDF,DeconvolutionHD,AddWeights=FALSE)

# Detect Tumor microenviroment (barcodes within 50 microns of tumor)
# Either load the results from Figure 4 or compute again
if(file.exists(SampleData$PathPeriphery[index]))
{
  PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
  
}else{
  
  PeripheryBCs<-SelectPeripheryDiscrete(bcs = bcDF,CellType = TumorCluster,distance=50,PATH = SampleData$PathSR[index])
}

bcDF$Periphery<-NA
bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==TumorCluster]<-"Tumor"

# Only keep Singlet and Doublet Certain barcodes
bcDF <- bcDF %>% filter(tissue==1 & DeconvolutionClass %in% c("singlet","doublet_certain")) %>% na.omit()

#Rename Macrophages
bcDF$DeconvolutionLabel1[match(TypeInfo[[1]],bcDF$barcode)]<-paste0("Macrophage ",names(TypeInfo)[1])
bcDF$DeconvolutionLabel1[match(TypeInfo[[2]],bcDF$barcode)]<-paste0("Macrophage ",names(TypeInfo)[2])

# Custom function to plot colocalization of macrophages and T cells

Plot1<-PlotColocalization(bcDF,"CD8 T cell",c("Macrophage SELENOP+","Macrophage SPP1+"),
                          Tumor=TumorCluster,BarcodeSet=bcDF$barcode[bcDF$Periphery=="50 micron"],option = "CD8")

Plot2<-PlotColocalization(bcDF,"CD4 T cell",c("Macrophage SELENOP+","Macrophage SPP1+"),
                          Tumor=TumorCluster,BarcodeSet=bcDF$barcode[bcDF$Periphery=="50 micron"],option = "CD4")

Plot1/Plot2

### Here we list BC and gene relationship for zoom-ins

#               Barcode               
#P1CRC  s_008um_00176_00486-1 
#P2CRC  s_008um_00460_00073-1
#P5CRC  s_008um_00365_00224-1


# Generate sample data.frame
Sample_path <- SampleData$PathSR

DataHD<-GenerateSampleData(Sample_path)
bcsHD<-DataHD$bcs %>% filter(tissue==1) %>% na.omit()


# Get Zoom in Region
Subset<-GetSquare(Spot = "s_008um_00176_00486-1",SizeMicrons = 200,BarcodeDF = bcsHD)

# Import Zoom Image (TIF)
image <- readTIFF("img_Stardist.tif")
image[image>1]<-1
image[image<0]<-0
grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
images_tibble <- tibble(Path = factor("img_Stardist.tif"), 
                        grob = list(grobs), height = nrow(image), width = ncol(image))

# Transform to 2um bins
Transformed2um<-TransformBarcodes(Subset,SizeOriginal = 8,SizeNew = 2)

# Create 2um data.frame
bcDF_2um<-GenerateSampleData(SampleData$PathSR,size="002um")
bcDF_2um<-bcDF_2um$bcs %>% filter(tissue == "1")

# Add equivalent 8um bin
bcDF_2um$BC8um<-Transformed2um$Original[match(bcDF_2um$barcode,Transformed2um$Transformed)]

# Read nuclei segmentation results
NucleiBins<-read.csv(file=paste0("~/Outputs/Deconvolution/",SampleData$Patient,"/Nuclei_Barcode_Map.csv"),header = T)

# Remove bins that overlap (i.e. a pixel that is located in two different Nuclei)
NucleiBins<- NucleiBins %>% filter(is_overlap=="False")

# Add info to 2um data.frame
bcDF_2um$NucleiID<-NucleiBins$id[match(bcDF_2um$barcode,NucleiBins$barcode)]
bcDF_2um$isNuclei<-!is.na(bcDF_2um$NucleiID)

# Read H5 (2um) and subset to BCs in the data.frame
SeuratObj<-Read10X_h5(paste0(SampleData$PathSR,"/binned_outputs/square_002um/filtered_feature_bc_matrix.h5"))[,bcDF_2um$barcode]

# Select only bins that are within Nuclei masks
SeuratObj<-SeuratObj[,bcDF_2um$barcode[bcDF_2um$isNuclei]]

# Transform UMI matrix from Gene x Bins to Gene x Nuclei
Group<-bcDF_2um$NucleiID[match(colnames(SeuratObj),bcDF_2um$barcode)]
Group<-factor(Group,levels = unique(Group))
group_mat <- sparse.model.matrix(~ 0 + Group) 
colnames(group_mat) <- colnames(group_mat) %>% str_extract("(?<=^Group).+")
SeuratObj <- SeuratObj %*% group_mat

# Create Seurat Object and Normalize
SeuratObj<-CreateSeuratObject(SeuratObj)
SeuratObj<-NormalizeData(SeuratObj)

# Add Expression to data.frame
Genes<-c("TRAC","CD3E","PECAM1","IGKC","COL1A1","SPP1","SELENOP","CEACAM5")
Exp<-FetchData(SeuratObj,Genes)
bcDF_2um<-cbind(bcDF_2um,Exp[match(bcDF_2um$NucleiID,rownames(Exp)),])

MaxExpValue<-ceiling(max(Exp))

# Create Plots
PlotList<-vector("list",length=length(Genes))
names(PlotList)<-Genes
for(Gx in Genes)
{
  
  bcDF_2um$Expression<-bcDF_2um[,Gx]

  PlotList[[Gx]]<-bcDF_2um %>% filter(isNuclei & !is.na(BC8um)) %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Expression)) +
    geom_spatial(data = images_tibble, aes(grob = grob), x = 0.5, y = 0.5) +
    geom_point(shape=15,size=1.5)+
    coord_cartesian(expand = FALSE) +
    scale_color_viridis(limits=c(0,MaxExpValue))+
    xlab("") +
    ylab("") + theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + ggtitle(Gx) 
  
  
}


wrap_plots(PlotList,ncol=4,nrow=2 ,guides = "collect")
