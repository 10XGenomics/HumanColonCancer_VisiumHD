#### Visium HD Manuscript
## Supplemental Figure 7

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
library(igraph)
library(ggnetwork


# load aux functions
source("~/Methods/AuxFunctions.R")

## This analysis requires to process all samples and combine in order to generate
## the corresponding results. We wrap the full piepeline in a loop.

SampleData<-data.frame(Patient = c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/VisiumHD/PatientCRC1/outs/","~/VisiumHD/PatientCRC2/outs/","~/VisiumHD/PatientCRC5/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Periphery/P1CRC_TME_Barcodes.rds","~/Outputs/Periphery/P2CRC_TME_Barcodes.rds","~/Outputs/Periphery/P5CRC_TME_Barcodes.rds"),
                       Tumor = c("Tumor II","Tumor III","Tumor IV"))


CellType<-"Macrophage"
BCMacro<-vector("list",length=3)
Matrices<-vector("list",length=3)

for(index in 1:nrow(SampleData))
{
  # Generate sample data.frame
  Sample_path <- SampleData$PathSR[index]
  
  TumorCluster<-SampleData$Tumor[index]
  
  DataHD<-GenerateSampleData(Sample_path)
  bcsHD<-DataHD$bcs
  
  # Read Deconvolution results and add to DF
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  # Filter data.frame to keep bins within tissues and are singlets (only 1 cell type)
  bcDF <- bcsHD %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()
  
  PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
  
  # Add periphery results to data.frame
  bcDF$Periphery<-NA
  bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
  bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
  bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==TumorCluster]<-"Tumor"
  
  bcDF$Patient<-SampleData$Patient[index]
  
  Rex<-bcDF %>% filter(DeconvolutionLabel1==CellType) %>% select(c("barcode","Periphery","Patient"))
  Rex$BCMat<-paste0(Rex$barcode,"_",Rex$Patient)
  
  
  mat <- open_matrix_dir(dir = getwd()) # DIR where the expression is saved on disk
  colnames(mat)<-paste0(colnames(mat),"_",SampleData$Patient[index])
  
  Genes<-read.delim(file="~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header = F)
  rownames(mat)<-make.unique(Genes$V2[match(rownames(mat),Genes$V1)])
  
  Matrices[[index]] <- mat[,Rex$BCMat]
  BCMacro[[index]]<-Rex
  
  
}


BCMacro<-do.call(rbind,BCMacro)
BCMacro<-as.data.frame(BCMacro)
rownames(BCMacro)<-BCMacro$BCMat

merged.object <- CreateSeuratObject(counts = Matrices)
merged.object<-JoinLayers(merged.object)

merged.object<-AddMetaData(merged.object,BCMacro)

# Full seurat processing
merged.object<-NormalizeData(merged.object)
merged.object<-FindVariableFeatures(merged.object)
merged.object<-ScaleData(merged.object)
merged.object<-RunPCA(merged.object)
ElbowPlot(merged.object,ndims = 30)
merged.object<-FindNeighbors(merged.object,dims=1:15)
merged.object<-FindClusters(merged.object,resolution=0.2)

Mks<-FindAllMarkers(merged.object,logfc.threshold = 0.05,min.diff.pct = 0.05,only.pos = T) %>% filter(p_val_adj<0.05)


GenesPlot<-c(Mks %>% group_by(cluster) %>% slice_max(order_by = -p_val_adj,n=5) %>% pull(gene),"SELENOP","FOLR2")
GenesPlot<-GenesPlot[c(1:10,12:25,27,26,11)]

DotPlot(merged.object,features = GenesPlot,dot.min = 0.1)+coord_flip()+ylab("Cluster")+xlab("")+
  scale_color_gradientn(colours = colorspace::diverge_hcl(5))


PropsData<-as.data.frame(prop.table(table(merged.object$seurat_clusters,merged.object$Periphery),margin=1)*100)
levels(PropsData$Var2)<-c("<= 50 microns","> 50 microns")
PropsData$Var2[PropsData$Var2=="Tissue"]<-"> 50 microns"

PlotProp1<-ggplot(PropsData,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=c("gray28","lightgray"))+
  xlab("Macrophage Cluster")+ylab("Proportion (%)")+coord_flip()+labs(fill="Distance")

PropsData2<-as.data.frame(prop.table(table(merged.object$seurat_clusters,merged.object$Patient),margin=1)*100)

PlotProp2<-ggplot(PropsData2,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+theme_classic()+
  xlab("Macrophage Cluster")+ylab("Proportion (%)")+coord_flip()+labs(fill="Patient")+scale_fill_manual(values=c("#3b9ab2","#ebcb2a","#f21800"))

PlotProp1+PlotProp2


MacroMeta<-merged.object@meta.data

### Spatial Plots
index<-which(SampleData$Patient=="P1CRC")
MacroIx<-MacroMeta[MacroMeta$Patient==SampleData$Patient[index],]

DataHD<-GenerateSampleData(SampleData$PathSR[index])
bcsHD<-DataHD$bcs

DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)

# Filter data.frame to keep bins within tissues and are singlets (only 1 cell type)
bcDF <- bcsHD %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()

# Create used variables (Patient ID and Tumor Cluster)
Patient<-SampleData$Patient[index]
TumorCluster<-SampleData$Tumor[index]

# Detect Tumor microenviroment (barcodes within 50 microns of tumor)
# Either load the results from Figure 4 or compute again
if(file.exists(SampleData$PathPeriphery[index]))
{
  PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
  
}else{
  
  PeripheryBCs<-SelectPeripheryDiscrete(bcs = bcDF,CellType = TumorCluster,distance=50,PATH = SampleData$PathSR[index])
}

# Add periphery results to data.frame
bcDF$Periphery<-NA
bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==TumorCluster]<-"Tumor"

iix<-match(rownames(bcDF),MacroIx$barcode)

bcDF$MacroNew<-as.vector(MacroIx$seurat_clusters[match(bcDF$barcode,MacroIx$barcode)])
bcDF$MacroNew[bcDF$DeconvolutionLabel1==SampleData$Tumor[index]]<-"Tumor"
bcDF$MacroNew<-factor(bcDF$MacroNew,levels = sort(unique(bcDF$MacroNew)))

bcDF  %>%  
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=MacroNew)) +
  geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
  coord_cartesian(expand = FALSE) +
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10))+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("gold","sienna1","green4","purple","azure4"),na.value = "lightgray")+
  labs(color="Cluster")+ggtitle("")+NoLegend()




