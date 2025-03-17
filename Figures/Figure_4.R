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

#P1CRC  s_008um_00556_00532-1 (Tumor II)
#P2CRC  s_008um_00587_00568-1 (Tumor III)
#P5CRC  s_008um_00348_00165-1 (Tumor IV)

TumorCluster<-"Tumor II"
  
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
Slice<-GetSquare("s_008um_00556_00532-1",SizeMicrons = 500,BarcodeDF = bcsHD)
  
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

  # Macrophage Analysis

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

