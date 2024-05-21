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

# load aux functions
source("~/Methods/AuxFunctions.R")


# Create data.frame with Patient information (We use PCRC1 as an example)
SampleData<-data.frame(Patient = "PatientCRC1",
                       PathSR="~/VisiumHD/PatientCRC1/outs/",
                       PathDeconvolution="~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds")

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
image <- readTIFF("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Results/NucleiSeg/F00150029/img_Stardist.tif")
image[image>1]<-1
image[image<0]<-0
grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
images_tibble <- tibble(Path = factor("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Results/NucleiSeg/F00150029/img_Stardist.tif"), 
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
