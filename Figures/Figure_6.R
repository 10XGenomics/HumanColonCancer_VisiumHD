#### Visium HD Manuscript
## Figure 6

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
library(Matrix)
library(tiff)

# load aux functions
source("~/Methods/AuxFunctions.R")

# Create data.frame with Patient information 
SampleData<-data.frame(Patient = c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/VisiumHD/PatientCRC1/outs/","~/VisiumHD/PatientCRC2/outs/","~/VisiumHD/PatientCRC5/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Periphery/P1CRC_TME_Barcodes.rds","~/Outputs/Periphery/P2CRC_TME_Barcodes.rds","~/Outputs/Periphery/P5CRC_TME_Barcodes.rds"),
                       Tumor = c("Tumor II","Tumor III","Tumor IV"))


# T cell deconvolution class
Result<-vector("list",length=nrow(SampleData))
names(Result)<-SampleData$Patient
for(index in 1:nrow(SampleData))
{
  bcsHD<-GenerateSampleData(SampleData$PathSR[index])$bcs
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  # Detect Tumor microenviroment (barcodes within 50 microns of tumor)
  # Either load the results from Figure 4 or compute again
  if(file.exists(SampleData$PathPeriphery[index]))
  {
    PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
    
  }else{
    
    PeripheryBCs<-SelectPeripheryDiscrete(bcs = bcDF,CellType = TumorCluster,distance=50,PATH = SampleData$PathSR[index])
  }
  
  bcsHD$Periphery<-NA
  bcsHD$Periphery[bcsHD$barcode%in%PeripheryBCs]<-"50 micron"
  bcsHD$Periphery[is.na(bcsHD$Periphery)]<-"Tissue"
  bcsHD$Periphery[bcsHD$Periphery=="Tissue" & bcsHD$DeconvolutionLabel1==SampleData$Tumor[index]]<-"Tumor"
  
  
  Ix<-as.data.frame(prop.table(table(bcsHD$Periphery,bcsHD$DeconvolutionClass),margin=1)*100)
  Ix$Patient<-names(Result)[index]
  
  Result<-rbind(Result,Ix)
  
}

# Calculate mean and SD for each class to make the barplot with errorbars
# Calculate mean and SD for each class to make the barplot with errorbars
Result2<-rbind(Result %>% filter(Var2=="reject") %>% dplyr::group_by(Var1) %>% dplyr::summarise(MD=mean(Freq),SD=sd(Freq),Group="reject"),
              Result %>% filter(Var2=="doublet_uncertain") %>% dplyr::group_by(Var1) %>% dplyr::summarise(MD=mean(Freq),SD=sd(Freq),Group="doublet_uncertain"),
              Result %>% filter(Var2=="doublet_certain") %>% dplyr::group_by(Var1) %>% dplyr::summarise(MD=mean(Freq),SD=sd(Freq),Group="doublet_certain"),
              Result %>% filter(Var2=="singlet") %>% dplyr::group_by(Var1) %>% dplyr::summarise(MD=mean(Freq),SD=sd(Freq),Group="singlet"))


Result2$Group<-factor(Result2$Group,levels = c("singlet","doublet_certain","doublet_uncertain","reject"))

ggplot(Result2, aes(x=Group, y=MD, fill=Var1)) + geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=MD-SD, ymax=MD+SD), width=.2,position=position_dodge(.9)) +scale_fill_manual(values=c("dodgerblue","lightgray","firebrick1"))+
  theme_classic()+xlab("") + ylab("Proportion (%)")+labs(fill="Group")+geom_point(data=Result,aes(x=Var2,y=Freq),position = position_dodge(width = .9))





# C) 8 um and 2 um Nuclei segnemntation integration results

### Here we list BC used for zoom-ins

#                  BC 
#P1CRC  s_008um_00176_00486-1
#P2CRC  s_008um_00460_00073-1 
#P5CRC  s_008um_00365_00224-1 

index<-1

# Import Zoom in TIF (Output of the NucleiSegmentation.py Script)
imagePath<-"~/Outputs/NucleiSegmentation/img_Stardist_P1CRC.tif"
image <- readTIFF(imagePath)

# Black and White image
image[image>1]<-1
image[image<0]<-0

# Create tibble with image data for geom_spatial
grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
images_tibble <- tibble(Path = factor(imagePath), grob = list(grobs), height = nrow(image), width = ncol(image))

# Create data.frame and add Deconvolution Info
bcDF<-GenerateSampleData(SampleData$PathSR[index])$bcs
DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
bcDF<-AddDeconvolutionInfo(bcDF,DeconvolutionHD,AddWeights=FALSE)
bcDF<- bcDF %>% filter(tissue==1) %>% na.omit()

# Get Zoom in Region
Subset<-GetSquare(Spot = "s_008um_00176_00486-1",SizeMicrons = 200,BarcodeDF = bcDF)

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


# Plot to show full deconvolution Labels with yellow box for zoom in
PlotA<-bcDF %>% filter(DeconvolutionClass=="singlet") %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=DeconvolutionLabel1)) +
  geom_scattermore(pointsize = 2,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  geom_rect(aes(xmin=DFRect$Xmin, xmax=DFRect$Xmax, ymin=DFRect$Ymax, ymax=DFRect$Ymin),fill=NA,color="yellow",linewidth=1)+
  scale_color_manual(values=ColorPalette())+NoLegend()

bcDF$DeconvolutionLabel1[which(is.na(bcDF$DeconvolutionClass) | bcDF$DeconvolutionClass =="reject")]<-NA

# Plot to show deconvolved singlet 8um bins in the zoom in Region
PlotB<-bcDF  %>% filter(barcode %in% Subset & DeconvolutionClass=="singlet") %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=DeconvolutionLabel1)) +
  geom_spatial(data = images_tibble, aes(grob = grob), x = 0.5, y = 0.5) +
  geom_point(shape=15,size=6)+
  coord_cartesian(expand = FALSE) +
  scale_color_manual(values = ColorPalette(),na.value = NA) +
  xlab("") +
  ylab("") + theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank())+NoLegend()

# Draw Plot
PlotA+PlotB

# Create 2um DF
bcDF_2um<-GenerateSampleData(SampleData$PathSR[index],size="002um")
bcDF_2um<-bcDF_2um$bcs %>% filter(tissue==1)

# Add corresponding 8um bins. Only for the zoon in region
bcDF_2um$BC8um<-Transformed2um$Original[match(bcDF_2um$barcode,Transformed2um$Transformed)]

# Read Segmentation Results (Output from NucleiSegmentation.py)
NucleiBins<-read.csv(file="/Outputs/NucleiSegmentation/Nuclei_Barcode_Map_P1CRC.csv",header = T)

# Remove overlaps
NucleiBins<-NucleiBins %>% filter(is_overlap=="False")

# Add Nuclei Segmentation Info to the 2um data.frame
bcDF_2um$NucleiID<-NucleiBins$id[match(bcDF_2um$barcode,NucleiBins$barcode)]
bcDF_2um$isNuclei<-!is.na(bcDF_2um$NucleiID)

# Read 2um .h5 UMI matrix
SeuratObj<-Read10X_h5(paste0(SampleData$PathSR[index],"/binned_outputs/square_002um/filtered_feature_bc_matrix.h5"))[,bcDF_2um$barcode]

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

# Get normalized expression of CD4 and CD8
Exp<-FetchData(SeuratObj,c("CD4","CD8A"))
bcDF_2um<-cbind(bcDF_2um,Exp[match(bcDF_2um$NucleiID,rownames(Exp)),])

# Create plots
PlotExp1<-bcDF_2um %>% filter(isNuclei & !is.na(BC8um)) %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=CD4)) +
  geom_spatial(data = images_tibble, aes(grob = grob), x = 0.5, y = 0.5) +
  geom_point(shape=15,size=1.5)+
  coord_cartesian(expand = FALSE) +
  scale_color_viridis(limits=c(0,7)) +
  xlab("") +
  ylab("") + theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

PlotExp2<-bcDF_2um %>% filter(isNuclei & !is.na(BC8um)) %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=CD8A)) +
  geom_spatial(data = images_tibble, aes(grob = grob), x = 0.5, y = 0.5) +
  geom_point(shape=15,size=1.5)+
  coord_cartesian(expand = FALSE) +
  scale_color_viridis(limits=c(0,7)) +
  xlab("") +
  ylab("") + theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#Draw plots
PlotExp1+PlotExp2
