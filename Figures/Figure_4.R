#### Visium HD Manuscript
## Figure 4

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


# Create data.frame with Patient information (We use PCRC1 as an example)
SampleData<-data.frame(Patient = "PatientCRC1",
                       PathSR="~/VisiumHD/PatientCRC1/outs/",
                       PathDeconvolution="~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds")


# Generate sample data.frame
Sample_path <- SampleData$PathSR
  
DataHD<-GenerateSampleData(Sample_path)
bcsHD<-DataHD$bcs
  
# Read Deconvolution results and add to DF
DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution)
bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  

### Here we list BCsfor zoom-ins

#P1CRC  s_008um_00556_00532-1 (Tumor III)
#P2CRC  s_008um_00587_00568-1 (Tumor II)
#P5CRC  s_008um_00348_00165-1 (Tumor IV)

TumorCluster<-"Tumor III"
  
bcsHD <- bcsHD %>% dplyr::filter(tissue==1) %>% na.omit()
  
# Subset Seurat object to barcodes in data.frame and Normalize Data
SrtObj<-Load10X_Spatial(Sample_path, bin.size = 8)
SrtObj<-subset(SrtObj,cells=bcsHD$barcode)
SrtObj<-NormalizeData(SrtObj)
  
# Add Expression of genes to data.frame
bcsHD<-AddExpression(bcsHD,SrtObj,c("CSF3R","C1QC","COL1A1"))
  
# Function to obtain all bins that are within a given distance to a specific Cluster/Cell Type.
# Paraneters for SelectPeripheryDiscrete :
#                                   - bcs = data.frame with Section info (Generated via GenerateSampleData)
#                                   - CellType = cell type for which to calculate distance from 
#                                   - distance = distance in microns to select bins
#                                   - PATH = path to SR outs folder
  
PeripheryBCs<-SelectPeripheryDiscrete(bcs = bcDF,CellType = TumorCluster,distance=50,PATH = SampleData$PathSR)

# Save results to avoid computing again
saveRDS(PeripheryBCs,file=paste0("~/Outputs/Periphery/",SampleData$Patient,"_TME_Barcodes.rds"))

# Add results to bcsHD
bcsHD$Periphery<-NA
bcsHD$Periphery[bcsHD$barcode%in%PeripheryBCs]<-"50 micron"
bcsHD$Periphery[is.na(bcsHD$Periphery)]<-"Tissue"
bcsHD$Periphery[bcsHD$Periphery=="Tissue" & bcsHD$DeconvolutionLabel1==TumorCluster]<-"Tumor"
  
# Get square cutouts for zoom in plots (500 micron side length)
Slice<-GetSquare("s_008um_00587_00568-1",SizeMicrons = 500,BarcodeDF = bcsHD)
  
# Get square limits to plot
DFRect<-data.frame(Xmin=min(bcsHD$imagecol_scaled[match(Slice,bcsHD$barcode)]),
                     Xmax=max(bcsHD$imagecol_scaled[match(Slice,bcsHD$barcode)]),
                     Ymin=min(-bcsHD$imagerow_scaled[match(Slice,bcsHD$barcode)]),
                     Ymax=max(-bcsHD$imagerow_scaled[match(Slice,bcsHD$barcode)]),
                     imagecol_scaled=NA,
                     image_row_scaled=NA)
  
# Full section with periphery results
Plot1<-bcsHD %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Periphery)) +
    geom_scattermore(pointsize = 2,pixels = rep(2000,2))+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_color_manual(values=c("dodgerblue","lightgray","firebrick1"))+
    labs(color="")+NoLegend()+geom_rect(aes(xmin=DFRect$Xmin, xmax=DFRect$Xmax, ymin=DFRect$Ymax, ymax=DFRect$Ymin),fill=NA,color="black")
  
# Zoom in with periphery results
Plot2<-bcsHD %>% filter(barcode %in% Slice) %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=Periphery)) +
    geom_point(shape=22,size=3,color="black",stroke=0)+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_fill_manual(values=c("dodgerblue","lightgray","firebrick1"))+
    labs(color="")+NoLegend()
  
# Expression of C1QC in Zoom in region
A1<-bcsHD[bcsHD$barcode%in%Slice,] %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=C1QC)) +
    geom_point(shape=22,size=3,stroke=0)+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_fill_gradient(low = "lightgray",high = "firebrick1")
  
# Expression of COL1A1 in Zoom in region
A2<-bcsHD[bcsHD$barcode%in%Slice,] %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=COL1A1)) +
    geom_point(shape=22,size=3,stroke=0)+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_fill_gradient(low = "lightgray",high = "firebrick1")
  
# Periphery plots (full section + zoom-ins)
Plot1+Plot2+A1+A2+plot_layout(ncol=4)
  
  
# Get Proportions for Periphery (only for singlet bins)
bcsHD<-bcsHD %>% filter(DeconvolutionClass == "singlet")
PeripheryProp<-as.data.frame(prop.table(table(bcsHD$Periphery,bcsHD$DeconvolutionLabel1),margin = 1)*100)
PeripheryProp$Patient<-SampleData$Patient

ggplot(PeripheryProp,aes(x=Var2,y=Freq,color=Var1,shape=Patient))+
  geom_point(size=3)+theme_classic()+scale_color_manual(values=c("dodgerblue","lightgray","firebrick1"))+
  xlab("")+ylab("Proportion (%)")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+coord_flip()

