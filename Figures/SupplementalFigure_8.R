#### Visium HD Manuscript
## Supplemental Figure 8

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

# Create data.frame with Patient information 
SampleData<-data.frame(Patient = c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/VisiumHD/PatientCRC1/outs/","~/VisiumHD/PatientCRC2/outs/","~/VisiumHD/PatientCRC5/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Periphery/P1CRC_TME_Barcodes.rds","~/Outputs/Periphery/P2CRC_TME_Barcodes.rds","~/Outputs/Periphery/P5CRC_TME_Barcodes.rds"),
                       Tumor = c("Tumor III","Tumor II","Tumor IV"))


# T cell deconvolution class
Result<-vector("list",length=nrow(SampleData))
names(Result)<-SampleData$Patient
for(index in 1:nrow(SampleData))
{
  TumorCluster<-SampleData$Tumor[index]
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
  bcsHD$Periphery[bcsHD$Periphery=="Tissue" & bcsHD$DeconvolutionLabel1==TumorCluster]<-"Tumor"
  
  bcsHD <- bcsHD %>% filter(tissue=="1" & DeconvolutionLabel1%in%c("CD4 T cell","CD8 T cell") & Periphery %in% c("50 micron","Tissue"))
  
  bcsHD<-split(bcsHD,bcsHD$Periphery)
  
  Ix<-lapply(bcsHD,function(X){Prop<-as.data.frame(prop.table(table(X$DeconvolutionLabel1,X$DeconvolutionClass),margin = 1)*100)})
  Ix<-do.call(rbind,Ix)
  Ix$Patient<-SampleData$Patient[index]
  
  Result<-rbind(Result,Ix)
  
}

Result$Periphery<-sapply(strsplit(rownames(Ix),"[.]"),function(X){return(X[1])})
# Draw Plot
ggplot(Result,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values=c("orange","purple2","firebrick1","dodgerblue"))+
  facet_grid(rows=vars(Periphery),cols=vars(Patient))+theme_classic()+labs(fill="Group")+ylab("Proportion (%)")+xlab("Cell Type")
