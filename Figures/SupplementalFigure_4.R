#### Visium HD Manuscript
## Supplemental Figure 3

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
                       Tumor = c("Tumor II","Tumor III","Tumor IV"))


# Deconolution Labels Proportions
Result<-vector("list",length=nrow(SampleData))
names(Result)<-SampleData$Patient
for(index in 1:nrow(SampleData))
{
  bcsHD<-GenerateSampleData(SampleData$PathSR[index])$bcs
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  bcsHD <- bcsHD %>% filter(tissue=="1" & DeconvolutionClass=="singlet")
  
  Ix<-as.data.frame(prop.table(table(bcsHD$DeconvolutionLabel1))*100)
  Ix$Patient<-SampleData$Patient[index]
  
  Result<-rbind(Result,Ix)
  
}

# Order factors
Levels<-Result %>% group_by(Var1) %>% summarise(Value=mean(Freq))
Levels<-Levels[order(Levels$Value,decreasing = T),]
Result$Var1<-factor(Result$Var1,levels = Levels$Var1)

# Draw Plot
ggplot(Result,aes(x=Var1,y=Freq,fill=Patient))+geom_bar(stat="identity",position = "dodge")+scale_fill_manual(values=c("#FB6222","#17AB6F","#572D86"))+
  theme_classic()+xlab("Cell Type")+ylab("Proportion (%)")+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

