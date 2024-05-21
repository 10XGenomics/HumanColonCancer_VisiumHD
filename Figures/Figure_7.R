#### Visium HD Manuscript
## Figure 7

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


#HD Zoom-in on same region as Xenium

# Import Zoom in TIF (Output of the NucleiSegmentation.py Script)
imagePath<-"~/Outputs/NucleiSegmentation/img_Stardist_P5CRC_Xenium.tif"
image <- readTIFF(imagePath)

# Black and White image
image[image>1]<-1
image[image<0]<-0

# Create tibble with image data for geom_spatial
grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
images_tibble <- tibble(Path = factor(imagePath), grob = list(grobs), height = nrow(image), width = ncol(image))

# Create data.frame and add Deconvolution Info
bcDF<-GenerateSampleData("~/VisiumHD/PatientCRC5/outs/")$bcs
bcDF<- bcDF %>% filter(tissue==1) %>% na.omit()

# Get Zoom in Region
Subset<-GetSquare(Spot = "s_008um_00367_00233-1",SizeMicrons = 200,BarcodeDF = bcDF)

# Transform to 2um bins
Transformed2um<-TransformBarcodes(Subset,SizeOriginal = 8,SizeNew = 2)

# Zoom in square limits (To add yellow box to plots)
DFRect<-data.frame(Xmin=min(bcDF$imagecol_scaled[match(unique(Transformed2um$Original),bcDF$barcode)]),
                   Xmax=max(bcDF$imagecol_scaled[match(unique(Transformed2um$Original),bcDF$barcode)]),
                   Ymin=min(-bcDF$imagerow_scaled[match(unique(Transformed2um$Original),bcDF$barcode)]),
                   Ymax=max(-bcDF$imagerow_scaled[match(unique(Transformed2um$Original),bcDF$barcode)]),
                   imagecol_scaled=NA,
                   image_row_scaled=NA)

# Function used to display the parameters to be used with the Python (NucleiSegmentation.py)
NucleiSegmentationScript(bcDF,Transformed2um)

# Create 2um DF
bcDF_2um<-GenerateSampleData("~/VisiumHD/PatientCRC5/outs/",size="002um")$bcs
bcDF_2um<-bcDF_2um %>% filter(tissue==1)

# Add corresponding 8um bins. Only for the zoon in region
bcDF_2um$BC8um<-Transformed2um$Original[match(bcDF_2um$barcode,Transformed2um$Transformed)]

# Read Segmentation Results (Output from NucleiSegmentation.py)
NucleiBins<-read.csv(file="/Outputs/NucleiSegmentation/Nuclei_Barcode_Map_P5CRC.csv",header = T)

# Remove overlaps
NucleiBins<-NucleiBins %>% filter(is_overlap=="False")

# Add Nuclei Segmentation Info to the 2um data.frame
bcDF_2um$NucleiID<-NucleiBins$id[match(bcDF_2um$barcode,NucleiBins$barcode)]
bcDF_2um$isNuclei<-!is.na(bcDF_2um$NucleiID)

# Read 2um .h5 UMI matrix
SeuratObj<-Read10X_h5("~/VisiumHD/PatientCRC5/outs/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")[,bcDF_2um$barcode]

# Transform the UMI matrix. From Gene x Bins to Gene x Nuclei
SeuratObj<-SeuratObj[,bcDF_2um$barcode[bcDF_2um$isNuclei]]
Group<-bcDF_2um$NucleiID[match(colnames(SeuratObj),bcDF_2um$barcode)]
Group<-factor(Group,levels = unique(Group))

group_mat <- sparse.model.matrix(~ 0 + Group) 
colnames(group_mat) <- colnames(group_mat) %>% str_extract("(?<=^Group).+")

SeuratObj <- SeuratObj %*% group_mat

# Create Seurat Object and Nornalize Data with the transformed matrix
SeuratObj<-CreateSeuratObject(SeuratObj)
SeuratObj<-NormalizeData(SeuratObj)

# Get normalized expression of selected genes
Genes<-c("CEACAM5","SELENOP","C1QC","JCHAIN","TRAC","CXCL9") 
Exp<-FetchData(SeuratObj,Genes) 
bcDF_2um<-cbind(bcDF_2um,Exp[match(bcDF_2um$NucleiID,rownames(Exp)),])

Plts<-vector("list",length=length(Genes))
names(Plts)<-Genes
for(Gx in Genes)
{
  Tmp<-bcDF_2um
  Tmp$Expression<-Tmp[,Gx]
  
  Plts[[Gx]]<-Tmp %>% filter(tissue == "1" & isNuclei & !is.na(BC8um)) %>%
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Expression)) +
    geom_spatial(data = images_tibble, aes(grob = grob), x = 0.5, y = 0.5) +
    geom_point(shape=15,size=2.5)+
    coord_cartesian(expand = FALSE) +
    scale_color_viridis(limits=c(0,5)) +
    xlab("") +
    ylab("") + theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank()) +
    ggtitle(Gx) 

}

wrap_plots(Plts,ncol=3,nrow =2 ,guides = "collect")