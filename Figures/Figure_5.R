#### Visium HD Manuscript
## Figure 5

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

## This analysis requires to process all samples and combine in order to generate
## the corresponding results. We wrap the full piepeline in a loop.

# Create data.frame with all Patient information 
SampleData<-data.frame(Patient = c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/VisiumHD/PatientCRC1/outs/","~/VisiumHD/PatientCRC2/outs/","~/VisiumHD/PatientCRC5/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds","~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Periphery/P1CRC_TME_Barcodes.rds","~/Outputs/Periphery/P2CRC_TME_Barcodes.rds","~/Outputs/Periphery/P5CRC_TME_Barcodes.rds"),
                       Tumor = c("Tumor III","Tumor II","Tumor IV"))


# List to save plots
ResultPlot<-vector("list",length=length(BarcodeDF))

# Variables to store Combined Results
CombRes<-c()
CombResTumor<-c()


# Specify celltype to be used. 
# Warning: If this is something different than "Macrophage" the function will mostly fail as it assumes there are 2 clusters (SPP1 and SELENOP)
CellType<-"Macrophage"
SeuratObjects<-vector("list",length=nrow(SampleData))
names(SeuratObjects)<-SampleData$Patient

for(index in 1:nrow(SampleData))
{
  # Generate sample data.frame
  Sample_path <- SampleData$PathSR[index]
  
  DataHD<-GenerateSampleData(Sample_path)
  bcsHD<-DataHD$bcs
  
  # Read Deconvolution results and add to DF
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[index])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  # Filter data.frame to keep bins within tissues and are singlets (only 1 cell type)
  bcDF <- bcsHD %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()
  
  # Create used variables (Patient ID and Tumor Cluster)
  Patient<-names(BarcodeDF)[index]
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
  
  # Create Seurat Object (8um)
  SeuratObj <- Load10X_Spatial(Sample_path, bin.size = 8)
  
  # Save object as is needed later
  SeuratObjects[[index]]<-SeuratObj
  
  # Keep only barcodes labeled as Macrophage and within the TME
  SeuratObj<-subset(SeuratObj,cells=bcDF$barcode[bcDF$DeconvolutionLabel1==CellType & bcDF$Periphery=="50 micron"])
  SeuratObj$CellType<-bcDF$DeconvolutionLabel1[match(colnames(SeuratObj),bcDF$barcode)]
  
  # Full seurat processing
  SeuratObj<-NormalizeData(SeuratObj)
  SeuratObj<-FindVariableFeatures(SeuratObj)
  SeuratObj<-ScaleData(SeuratObj)
  SeuratObj<-RunPCA(SeuratObj)
  SeuratObj<-FindNeighbors(SeuratObj,dims=1:10)
  SeuratObj<-FindClusters(SeuratObj,resolution=0.2)
  
  # Conditional based if two (SELENOP and SPP1) clusters are found
  if(length(unique(SeuratObj$seurat_clusters))>1)
  {
    # Identify cluster assignment by finding markers
    CtMarkers<-FindAllMarkers(SeuratObj,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T)
    CtMarkers<-CtMarkers[CtMarkers$p_val_adj<0.05,]
    
    # Label SPP1 and SELENOP macrophages
    SPP1_Cluster<-as.vector(CtMarkers$cluster[match("SPP1",CtMarkers$gene)])
    CtMarkers$cluster<-as.factor(ifelse(CtMarkers$cluster==SPP1_Cluster,"SPP1+","SELENOP+"))
    FinalMarkers<-CtMarkers
    
    # Relabel Seurat Object
    SeuratObj$Group<-ifelse(SeuratObj$seurat_clusters==SPP1_Cluster,"SPP1+","SELENOP+")
    SeuratObj<-SetIdent(SeuratObj,value = "Group")
    
    # Update data.frame labek as well
    iix<-match(colnames(SeuratObj),bcDF$barcode)
    bcDF$DeconvolutionLabel1[iix]<-paste0(bcDF$DeconvolutionLabel1[iix],"-",SeuratObj$Group)
    
    #Result data.frame with barcodes, macrophage type and Patient. To be used when combining the objects
    CombRes<-rbind(CombRes,data.frame(BC=bcDF$barcode[iix],CT=bcDF$DeconvolutionLabel1[iix],Patient=SampleData$Patient[index]))
    
    # Density plots to identify Macrophage enriched regions in the TME
    DensityPlot1<-PlotDensity(bcDF,"Macrophage-SELENOP+",nBins = 4,Tumor = TumorCluster,
                              BarcodeSet = bcDF$barcode[bcDF$Periphery=="50 micron"])
    
    DensityPlot2<-PlotDensity(bcDF,"Macrophage-SPP1+",nBins = 4,Tumor = TumorCluster,
                              BarcodeSet = bcDF$barcode[bcDF$Periphery=="50 micron"])
    
    # Save Plots
    ResultPlot[[index]]<-DensityPlot1+DensityPlot2+plot_layout(ncol=1)
    
    # Obtain the barcodes within the enriched regions 
    RegionsSPP1<-SelectEnrichedRegion("Macrophage-SPP1+",bcDF[bcDF$barcode%in%PeripheryBCs,],SampleData$PathSR[index],Area = 350,N = 1)[[1]]
    RegionsSELENOP<-SelectEnrichedRegion("Macrophage-SELENOP+",bcDF[bcDF$barcode%in%PeripheryBCs,],SampleData$PathSR[index],Area = 350,N = 1)[[1]]
    
    # Label the identified Regions
    bcDF$Region<-"Tissue"
    bcDF$Region[match(RegionsSPP1,bcDF$barcode)]<-"Region SPP1"
    bcDF$Region[match(RegionsSELENOP,bcDF$barcode)]<-"Region SELENOP"
    
    # Select closest tumor spots
    TumorSPP1<-SelectTumor(bcDF,TumorCluster,bcDF %>% filter(Region=="Region SPP1") %>% pull(barcode),SampleData$PathSR[index],50)
    TumorSELENOP<-SelectTumor(bcDF,TumorCluster,bcDF %>% filter(Region=="Region SELENOP") %>% pull(barcode),SampleData$PathSR[index],50)
    
    # Create data.frame that includes closest tumor information to be used in combined analysis
    Gpx<-data.frame(Barcode=c(TumorSPP1,TumorSELENOP),Group=c(rep("SPP1",length(TumorSPP1)),rep("SELENOP",length(TumorSELENOP))),Patient=SampleData$Patient[index])
    CombResTumor<-rbind(CombResTumor,Gpx)
    
    
    
  }else{
    # This is the case where no SPP1+ macrophages aren't found. Should only happen in PCRC5.
    
    # Rename all TME macrophages as SELENOP+
    iix<-match(colnames(SeuratObj),bcDF$barcode)
    bcDF$DeconvolutionLabel1[iix]<-paste0(bcDF$DeconvolutionLabel1[iix],"-SELENOP+")
    
    #Result data.frame with barcodes, macrophage type and Patient. To be used when combining the objects
    CombRes<-rbind(CombRes,data.frame(BC=bcDF$barcode[iix],CT=bcDF$DeconvolutionLabel1[iix],Patient=SampleData$Patient[index]))
    
    # Density plots to identify Macrophage enriched regions in the TME
    DensityPlot1<-PlotDensity(bcDF,"Macrophage-SELENOP+",nBins = 5,Tumor = TumorCluster,BarcodeSet = bcDF$barcode[bcDF$Periphery=="50 micron"])
    
    # Save Plots
    ResultPlot[[index]]<-DensityPlot1+plot_layout(ncol=1)
    
    # Obtain the barcodes within the enriched regions 
    RegionsSELENOP<-SelectEnrichedRegion("Macrophage-SELENOP+",bcDF[bcDF$barcode%in%PeripheryBCs,],PATH=SampleData$PathSR[index],Area = 350,N = 1)[[1]]
    
    # Label the identified Regions
    bcDF$Region<-"Tissue"
    bcDF$Region[match(RegionsSELENOP,bcDF$barcode)]<-"Region SELENOP"
    
    # Select closest tumor spots
    TumorSELENOP<-SelectTumor(bcDF,TumorCluster,bcDF %>% filter(Region=="Region SELENOP") %>% pull(barcode),PATH=SampleData$PathSR[index],50)
    
    # Create data.frame that includes closest tumor information to be used in combined analysis
    Gpx<-data.frame(Barcode=TumorSELENOP,Group=rep("SELENOP",length(TumorSELENOP)),Patient=SampleData$Patient[index])
    CombResTumor<-rbind(CombResTumor,Gpx)
    
    
  }
  
}

# Plot Macrophage Density Estimations for each patient
ResultPlot[[1]]
ResultPlot[[2]]
ResultPlot[[3]]

# Using the combined data.frames we will re-create Seurat Objects and merge
CombRes<-split(CombRes,CombRes$Patient)
CombinedObject<-vector("list",length=2)
i<-1
for(ix in names(CombRes))
{
  
  SeuratObj <- SeuratObjects[[ix]]
  SeuratObj <- subset(SeuratObj,cells=CombRes[[ix]]$BC)
  SeuratObj$CellType<-CombRes[[ix]]$CT[match(colnames(SeuratObj),CombRes[[ix]]$BC)]
  SeuratObj$Patient<-CombRes[[ix]]$Patient[match(colnames(SeuratObj),CombRes[[ix]]$BC)]
  CombinedObject[[i]]<-SeuratObj
  i<-i+1
}

# Merge objects and find DEG for each macrophage subpopulation
Combined<-merge(CombinedObject[[1]],c(CombinedObject[[2]],CombinedObject[[3]]))
Combined<-JoinLayers(Combined)
Combined<-NormalizeData(Combined)
Combined<-SetIdent(Combined,value="CellType")

CtMarkers<-FindAllMarkers(Combined,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T)
CtMarkers<-CtMarkers[CtMarkers$p_val_adj<0.05,]

# Get top 10 DEGs for each macrophage subpopulation
Genes<-CtMarkers %>% group_by(cluster) %>% arrange(-avg_log2FC) %>% slice(1:10) %>% pull(gene)

# Seurat DotPlot, we will use the data to make a customized version
Plot<-DotPlot(Combined,features = unique(Genes),dot.min = 0.15,idents=levels(CtMarkers$cluster),split.by = "Patient",cols = c("blue","red","green"))
Plot<-Plot$data
Plot<- Plot %>% na.omit()
Plot$LogExp<-log1p(Plot$avg.exp)

# Dot Plot Gene Expression
ggplot(Plot,aes(x=features.plot,y=id,size=pct.exp,color=LogExp))+geom_point()+coord_flip()+
  theme_classic()+scale_colour_viridis(limits=c(0,7))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Subpopulation")+guides(size = guide_legend(title = "Percentage Expressed"),color = guide_colorbar(title = "log Mean Expression"))

# Barplot with Enrichment based on the found Markers
EnrichRBarPlot(CtMarkers,"MSigDB_Hallmark_2020",TermsX = 10)

# Same as before, but this time with the Tumor 
CombResTumor<-split(CombResTumor,CombResTumor$Patient)
CombinedTumor<-vector("list",length=2)
i<-1
for(ix in names(CombResTumor))
{
  
  SeuratObj <- SeuratObjects[[ix]]
  SeuratObj <- subset(SeuratObj,cells=CombResTumor[[ix]]$Barcode)
  SeuratObj$CellType<-CombResTumor[[ix]]$Group[match(colnames(SeuratObj),CombResTumor[[ix]]$Barcode)]
  SeuratObj$Patient<-CombResTumor[[ix]]$Patient[match(colnames(SeuratObj),CombResTumor[[ix]]$Barcode)]
  CombinedTumor[[i]]<-SeuratObj
  i<-i+1
}


Combined<-merge(CombinedTumor[[1]],c(CombinedTumor[[2]],CombinedTumor[[3]]))
Combined<-JoinLayers(Combined)
Combined<-NormalizeData(Combined)
Combined<-SetIdent(Combined,value="CellType")

CtMarkers<-FindAllMarkers(Combined,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T)
CtMarkers<-CtMarkers[CtMarkers$p_val_adj<0.05,]

# Create Plots
PlotsExp<-vector("list",length=3)
for(jj in 1:length(BarcodeDF))
{
  
  bcDF<-BarcodeDF[[jj]]
  bcDF <- bcDF %>% filter(tissue==1 ) %>% na.omit()
  
  SeuratObj <- SeuratObjects[[jj]]
  SeuratObj<-subset(SeuratObj,cells=bcDF$barcode)
  SeuratObj<-NormalizeData(SeuratObj)
  
  bcDF<-AddExpression(bcDF,SeuratObj,c("REG1B","REG3A","REG1A","TGFBI"))
  Pa<-PlotExpression(bcDF,"REG1A",ptsize = 2.25)+scale_color_gradient(low = "lightgray",high = "firebrick1")
  Pb<-PlotExpression(bcDF,"TGFBI",ptsize = 2.25)+scale_color_gradient(low = "lightgray",high = "firebrick1")
  
  PlotsExp[[jj]]<-Pa+Pb+plot_layout(ncol=1)
}

# Draw Plots
PlotsExp[[1]]
PlotsExp[[2]]
PlotsExp[[3]]
