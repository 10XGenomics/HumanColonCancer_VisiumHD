#### Visium HD Manuscript
## Figure 3

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

# Plot unsupervised clustering results
PlotCluster<-bcsHD  %>% filter(!is.na(Cluster)) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=as.factor(Cluster))) +
  geom_scattermore(pointsize = 2,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=ColorsClusters()[14:39])+
  labs(color="Cluster")+NoLegend()

# Plot Deconvolution Results
PlotDeconvolution<-bcsHD  %>% na.omit() %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=DeconvolutionLabel1)) +
  geom_scattermore(pointsize = 2.5,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=ColorPalette())+
  labs(color="Cell Type")+NoLegend()


# A and B
PlotCluster/PlotDeconvolution

# C
# Get clustering - deconvolution confusion matrix
Mat<-prop.table(table(bcsHD$Cluster,bcsHD$DeconvolutionLabel1),margin=1)*100
Mat[is.nan(Mat)]<-0

HMCol<-colorRampPalette(c("white","firebrick1"))(50)

AnnotColors<- list(Deconvolution = c('Tumor III'='#5050FFFF','Plasma'='#CE3D32FF','Macrophage'='#749B58FF','CD4 T cell'='#F0E685FF','CAF'='#466983FF','vSM'='#BA6338FF',
                                     'Mature B'='#5DB1DDFF','Endothelial'='#802268FF','Tumor I'='#6BD76BFF','CD8 Cytotoxic T cell'='#D595A7FF','Enterocyte'='#924822FF',
                                     'Neutrophil'='#837B8DFF','Proliferating Immune II'='#C75127FF','Pericytes'='#D58F5CFF','Smooth Muscle'='#7A65A5FF','Myofibroblast'='#E4AF69FF',
                                     'Tumor II'='#3B1B53FF','Fibroblast'='#CDDEB7FF','Goblet'='#612A79FF','Lymphatic Endothelial'='#AE1F63FF','Tumor V'='#E7C76FFF',
                                     'Proliferating Macrophages'='#5A655EFF','SM Stress Response'='#CC9900FF','NK'='#99CC00FF','cDC I'='#A9A9A9FF','Tumor IV'='#CC9900FF',
                                     'Proliferating Fibroblast'='#99CC00FF','Epithelial'='#33CC00FF','Tuft'='#00CC33FF','Mast'='#00CC99FF',
                                     'Unknown III (SM)'='#0099CCFF','Adipocyte'='#0A47FFFF','mRegDC'='#4775FFFF','Enteric Glial'='#FFC20AFF',
                                     'pDC'='#FFD147FF','Vascular Fibroblast'='#990033FF','Neuroendocrine'='#991A00FF','Memory B'='#996600FF',
                                     'Unknwon I (Immune)'='#809900FF'),
                   Cluster = c('1'='#7F1786','2'='#EA3323','3'='#5DCBCF','4'='#F09235','5'='#FEFF54','6'='#B72D82','7'='#0C00C5','8'='#75FB4C','9'='#75FB8D','10'='#CA3142',
                               '11'='#56BCF9','12'='#E8A76C','13'='#932CE7','14'='#BEFD5B','15'='#EA33F7','16'='#458EF7','17'='#CD7693','18'='#EEE697','19'='#EA8677','20'='#D4A2D9',
                               '21'='#B6D7E4','22'='#7869E6','23'='#E087E8','24'='#AFF9A2','25'='#A0FBD6'))


AnnotRow<-data.frame(Cluster=rownames(Mat))
AnnotCol<-data.frame(Deconvolution=colnames(Mat))
rownames(AnnotCol)<-AnnotCol$Deconvolution
pheatmap(Mat,cluster_rows = T,color = HMCol,border_color = "black",annotation_colors = AnnotColors,annotation_col = AnnotCol,annotation_legend = FALSE,show_rownames = T,show_colnames = T,annotation_names_row = TRUE,annotation_names_col = FALSE,legend = FALSE)



# D) Expression Zoom Ins

### Here we list BC and gene relationship for zoom-ins

#                  PIGR                CEACAM6               COL1A1
#P1CRC  s_008um_00314_00245-1 s_008um_00194_00640-1 s_008um_00287_00712-1
#P2CRC  s_008um_00228_00549-1 s_008um_00530_00716-1 s_008um_00347_00358-1
#P5CRC  s_008um_00181_00558-1 s_008um_00146_00104-1 s_008um_00383_00204-1


# Define center barcode for zoom ins
Centers<-c("s_008um_00228_00549-1","s_008um_00530_00716-1","s_008um_00347_00358-1")

bcsHD<-bcsHD %>% filter(tissue==1) %>% na.omit()
  
# Subset to barcodes in the bcDF data.frame
SrtObj<-Load10X_Spatial(Sample_path, bin.size = 8)
SrtObj<-subset(SrtObj,cells=bcsHD$barcode)
SrtObj<-NormalizeData(SrtObj)
  
# Add Expression values as columns to data.frame
bcsHD<-AddExpression(bcsHD,SrtObj,c("PIGR","CEACAM6","COL1A1"))
  
# Use GetSquare to make a square cutout of the section. Each side will be 250 microns in length
SliceA<-GetSquare(Centers[1],250,bcsHD)
SliceB<-GetSquare(Centers[2],250,bcsHD)
SliceC<-GetSquare(Centers[3],250,bcsHD)
  
# Define Square limits to be plotted as well
DFRectA<-data.frame(Xmin=min(bcsHD$imagecol_scaled[match(SliceA,bcsHD$barcode)]),
                      Xmax=max(bcsHD$imagecol_scaled[match(SliceA,bcsHD$barcode)]),
                      Ymin=min(-bcsHD$imagerow_scaled[match(SliceA,bcsHD$barcode)]),
                      Ymax=max(-bcsHD$imagerow_scaled[match(SliceA,bcsHD$barcode)]),
                      imagecol_scaled=NA,
                      image_row_scaled=NA)
  
DFRectB<-data.frame(Xmin=min(bcsHD$imagecol_scaled[match(SliceB,bcsHD$barcode)]),
                      Xmax=max(bcsHD$imagecol_scaled[match(SliceB,bcsHD$barcode)]),
                      Ymin=min(-bcsHD$imagerow_scaled[match(SliceB,bcsHD$barcode)]),
                      Ymax=max(-bcsHD$imagerow_scaled[match(SliceB,bcsHD$barcode)]),
                      imagecol_scaled=NA,
                      image_row_scaled=NA)
  
DFRectC<-data.frame(Xmin=min(bcsHD$imagecol_scaled[match(SliceC,bcsHD$barcode)]),
                      Xmax=max(bcsHD$imagecol_scaled[match(SliceC,bcsHD$barcode)]),
                      Ymin=min(-bcsHD$imagerow_scaled[match(SliceC,bcsHD$barcode)]),
                      Ymax=max(-bcsHD$imagerow_scaled[match(SliceC,bcsHD$barcode)]),
                      imagecol_scaled=NA,
                      image_row_scaled=NA)
  
  
# Create All the Plots
(PlotExpression(bcsHD,"PIGR",ptsize = 3)+geom_rect(aes(xmin=DFRectA$Xmin, xmax=DFRectA$Xmax, ymin=DFRectA$Ymax, ymax=DFRectA$Ymin),fill=NA,color="black")+NoLegend())+PlotExpression(bcsHD[bcsHD$barcode%in%SliceA,],"PIGR",ptsize = 4,shape="square")+
    (PlotExpression(bcsHD,"CEACAM6",ptsize = 3)+geom_rect(aes(xmin=DFRectB$Xmin, xmax=DFRectB$Xmax, ymin=DFRectB$Ymax, ymax=DFRectB$Ymin),fill=NA,color="black")+NoLegend())+PlotExpression(bcsHD[bcsHD$barcode%in%SliceB,],"CEACAM6",ptsize = 4,shape="square")+
    (PlotExpression(bcsHD,"COL1A1",ptsize = 3)+geom_rect(aes(xmin=DFRectC$Xmin, xmax=DFRectC$Xmax, ymin=DFRectC$Ymax, ymax=DFRectC$Ymin),fill=NA,color="black")+NoLegend())+PlotExpression(bcsHD[bcsHD$barcode%in%SliceC,],"COL1A1",ptsize = 4,shape="square")+plot_layout(ncol=2)
  
