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
                       Tumor = c("Tumor II","Tumor III","Tumor IV"))

# List to save plots
ResultPlot<-vector("list",length=length(SampleData))

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

# Barplot with Enrichment based on the found Markers
# A)
EnrichRBarPlot(CtMarkers,"MSigDB_Hallmark_2020",TermsX = 15,colsT=c("#01702E","#D35D05"))

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
GenesCluster<-CtMarkers %>% group_by(cluster) %>% arrange(p_val_adj) %>% slice_head(n=5) %>% pull(gene)

ExpData<-FetchData(Combined,GenesCluster)

TumorPlotData<-cbind(Combined@meta.data[,4:5],ExpData)

TumorPlotData<-melt(TumorPlotData)

# B)
TumorPlotData %>% ggplot(aes(x=CellType,y=value,fill=Patient))+geom_violin(scale = "width")+facet_wrap(~variable,ncol = 5) + theme_classic() +
  scale_fill_manual(values=c("darkorange2","darkcyan","darkslateblue"))



## Goblet Cell Analysis

MatZ<-vector("list",length=nrow(SampleData))

for(jj in 1:nrow(SampleData))
{
  # Generate sample data.frame
  Sample_path <- SampleData$PathSR[jj]
  TumorCluster<-SampleData$Tumor[jj]
  DataHD<-GenerateSampleData(Sample_path)
  bcsHD<-DataHD$bcs
  
  # Read Deconvolution results and add to DF
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[jj])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  # Filter data.frame to keep bins within tissues and are singlets (only 1 cell type)
  bcDF <- bcsHD %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()
  BcsX <- bcDF %>% filter(DeconvolutionLabel1=="Goblet") %>% pull(barcode)
  
  
  Cnts<-Read10X_h5(paste0(SampleData$PathSR[jj],"/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"))[,BcsX]
  colnames(Cnts)<-paste0(colnames(Cnts),"-",SampleData$Patient[jj])
  
  MatZ[[jj]]<-Cnts
  
}

Counts<-do.call(cbind,MatZ)

GobletObj<-CreateSeuratObject(Counts)
GobletObj<-NormalizeData(GobletObj)
GobletObj<-FindVariableFeatures(GobletObj)
GobletObj<-ScaleData(GobletObj)
GobletObj<-RunPCA(GobletObj)
GobletObj<-FindNeighbors(GobletObj,dims=1:12)
GobletObj<-FindClusters(GobletObj,resolution=0.2)

Mks<-FindAllMarkers(GobletObj,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T) %>% filter(p_val_adj<0.05)

GobletMD<-GobletObj@meta.data
GobletMD$Barcode<-paste0(sapply(strsplit(rownames(GobletMD),"-"),function(X){return(X[1])}),"-1")
GobletMD$Patient<-sapply(strsplit(rownames(GobletMD),"-"),function(X){return(X[3])})
GobletMD$FinalCluster<-paste0("Goblet-",GobletMD$seurat_clusters)
GobletMD$FinalCluster<-factor(GobletMD$FinalCluster,levels = sort(unique(GobletMD$FinalCluster)))

# C)
DotPlot(GobletObj,c("MUC2","FCGBP","TFF3","CLCA1","OLFM4","DUOX2","DMBT1","REG1A","REG1B"),dot.min = 0.05)+coord_flip()+
  scale_color_gradient2(low="dodgerblue",mid = "white",high = "firebrick1",limits=c(-2.5,2.5))

ColorsX<-c("#CA3142","#F09235","#FEFF54","#0C00C5","#60B177","#EEE697","#74140C")

AllP1<-vector("list",length=3)
for(jj in 1:nrow(SampleData))
{
  
  MDD<-GobletMD %>% filter(Patient==SampleData$Patient[jj])
  
  # Generate sample data.frame
  Sample_path <- SampleData$PathSR[jj]
  TumorCluster<-SampleData$Tumor[jj]
  DataHD<-GenerateSampleData(Sample_path)
  bcsHD<-DataHD$bcs
  
  # Read Deconvolution results and add to DF
  DeconvolutionHD<-readRCTD(SampleData$PathDeconvolution[jj])
  bcsHD<-AddDeconvolutionInfo(bcsHD,DeconvolutionHD,AddWeights=FALSE)
  
  # Filter data.frame to keep bins within tissues and are singlets (only 1 cell type)
  bcDF <- bcsHD %>% filter(tissue==1 & DeconvolutionClass=="singlet") %>% na.omit()
  
  bcDF$isGoblet<-bcDF$DeconvolutionLabel1=="Goblet"
  
  bcDF$NewCT<-NA
  bcDF$NewCT[match(MDD$Barcode,bcDF$barcode)]<-as.vector(MDD$FinalCluster)
  bcDF$NewCT<-factor(bcDF$NewCT,levels = levels(GobletMD$FinalCluster))
  
  AllP1[[jj]]<-bcDF  %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=NewCT)) +
    geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_color_manual(values=ColorsX,na.value = "lightgray")+
    labs(color="Goblet Cluster")+ggtitle("")+NoLegend()
  
  
}

# D)
AllP1[[1]]+AllP1[[2]]+AllP1[[3]]

# Cell Cell Communication
# Load packages
library(liana)
library(circlize)
library(igraph)
library(distances)

AllResults<-vector("list",length=nrow(SampleData))
PlotA<-vector("list",length=nrow(SampleData))
PlotB<-vector("list",length=nrow(SampleData))
names(AllResults)<-names(PlotA)<-names(PlotB)<-SampleData$Patient

for(Patient in SampleData$Patient)
{
  #Generate MetaData
  index<-which(SampleData$Patient==Patient)
  bcDF<-GenerateSampleData(SampleData$PathSR[index])$bcs
  Deconv<-readRCTD(SampleData$PathDeconvolution[index])
  bcDF<-AddDeconvolutionInfo(bcDF,Deconv)
  bcDF<-bcDF %>% filter(tissue==1)
  
  # Add Periphery Results
  PeripheryBCs<-readRDS(SampleData$PathPeriphery[index])
  
  # Add periphery results to data.frame
  bcDF$Periphery<-NA
  bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
  bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
  bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==SampleData$Tumor[index]]<-"Tumor"
  
  # Create Seurat Object (8um)
  object<-Read10X_h5(paste0(SampleData$PathSR[index],"binned_outputs/square_008um/filtered_feature_bc_matrix.h5"))
  CommonBC<-intersect(bcDF$barcode,colnames(object))
  object<-CreateSeuratObject(object[,CommonBC],meta.data=bcDF[match(CommonBC,bcDF$barcode),])
  
  # Analysis To identify Macrophage subtypes
  
  # Keep only barcodes labeled as Macrophage and within the TME
  SeuratObj<-subset(object,cells=bcDF$barcode[bcDF$DeconvolutionLabel1=="Macrophage" & bcDF$Periphery=="50 micron" & bcDF$DeconvolutionClass=="singlet"])
  
  # Full seurat processing
  SeuratObj<-NormalizeData(SeuratObj)
  SeuratObj<-FindVariableFeatures(SeuratObj)
  SeuratObj<-ScaleData(SeuratObj)
  SeuratObj<-RunPCA(SeuratObj)
  SeuratObj<-FindNeighbors(SeuratObj,dims=1:10)
  SeuratObj<-FindClusters(SeuratObj,resolution=0.2)
  
  if(length(unique(SeuratObj$seurat_clusters))>1)
  {
    MkMac<-FindAllMarkers(SeuratObj,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = T)
    SPP1_Cluster<-as.vector(MkMac$cluster[match("SPP1",MkMac$gene)])
    SeuratObj$Group<-ifelse(SeuratObj$seurat_clusters==SPP1_Cluster,"SPP1+","SELENOP+")
    SeuratObj$ID<-paste0(SeuratObj$DeconvolutionLabel1,"_",SeuratObj$Group)
  }else{
    
    SeuratObj$Group<-"SELENOP+"
    SeuratObj$ID<-paste0(SeuratObj$DeconvolutionLabel1,"_",SeuratObj$Group)
    
  }
  
  
  # Add to bcDF
  bcDF$CellType<-bcDF$DeconvolutionLabel1
  bcDF$CellType[match(colnames(SeuratObj),bcDF$barcode)]<-SeuratObj$ID
  
  # Select Closest Tumor to cells in Periphery
  RegionI<-bcDF %>% filter(DeconvolutionClass=="singlet" & Periphery%in%c("50 micron","Tumor"))
  dist_matrix <- distances::distances(as.matrix(RegionI[,c("imagecol","imagerow")]))
  
  ## Chunk to make it efficient
  Values<-which(RegionI$Periphery=="Tumor")
  PeriIndex<-which(RegionI$Periphery=="50 micron")
  Chnks<-split(Values, ceiling(seq_along(Values)/1000))
  ClosestTumorSpot<-vector("list",length=length(Chnks))
  scales <- rjson::fromJSON(file = paste0(SampleData$PathSR[index],"/binned_outputs/square_008um/spatial/scalefactors_json.json"))
  
  for(Idx in 1:length(Chnks))
  {
    print(Idx)
    CHUNK<-Chnks[[Idx]]
    BcChnk<-RegionI$barcode[CHUNK]
    TmpDist<-dist_matrix[CHUNK,PeriIndex]
    TmpDist<-(TmpDist*8)/scales$spot_diameter_fullres
    
    TmpDist<-TmpDist<=50
    iix<-which(rowSums(TmpDist)>=5)
    
    #iix<-apply(TmpDist,1,function(X){any(X<=50)})
    if(length(iix)>0)
    {
      ClosestTumorSpot[[Idx]]<-BcChnk[iix]
    }
    
  }
  
  ClosestTumorSpot<-unique(unlist(ClosestTumorSpot))
  
  bcDF$Selection<-NA
  bcDF$Selection[bcDF$Periphery=="50 micron"]<-"50 micron"
  bcDF$Selection[bcDF$barcode%in%ClosestTumorSpot]<-"Tumor"
  
  # Recreate Seurat Object 
  object<-Read10X_h5(paste0(SampleData$PathSR[index],"binned_outputs/square_008um/filtered_feature_bc_matrix.h5"))
  CommonBC<-intersect(bcDF$barcode,colnames(object))
  object<-CreateSeuratObject(object[,CommonBC],meta.data=bcDF[match(CommonBC,bcDF$barcode),])
  object<-subset(object, subset = Selection %in% c("Tumor","50 micron") & DeconvolutionClass=="singlet")
  object <- NormalizeData(object)
  object<-SetIdent(object,value = "CellType")
  
  # Run C-C using Liana
  liana_result <- liana_wrap(object)%>% liana_aggregate()%>% filter(aggregate_rank <= 0.01)
  liana_result$Patient<-Patient
  AllResults[[Patient]]<-liana_result
  
  CTxs<-c("Macrophage_SPP1+","Macrophage_SELENOP+","CD4 T cell","CD8 T cell",SampleInfo$Tumor[SampleInfo$Sample==Patient])
  
  # Interaction Graph
  PlotB[[Patient]]<-PlotInteractionGraph(liana_result,CellTypes = CTxs,
                                         Colors = c(ColorPalette(),"Macrophage_SPP1+"="cyan2","Macrophage_SELENOP+"="orchid3"))
}

# Plot Networks
# E)
PlotB[[1]]


CC_Result<-do.call(rbind,AllResults)

iix<-(CC_Result$Patient==SampleInfo$Sample[1] & CC_Result$target==SampleInfo$Tumor[1]) |
  (CC_Result$Patient==SampleInfo$Sample[2] & CC_Result$target==SampleInfo$Tumor[2]) |
  (CC_Result$Patient==SampleInfo$Sample[3] & CC_Result$target==SampleInfo$Tumor[3])

CC_Result$target[iix]<-"Tumor"

# F)

CC_Result  %>%
  liana_dotplot(source_groups = c("Macrophage_SPP1+"),
                target_groups = c("CD8 T cell","CD4 T cell","Tumor"),ntop = 50)+
  facet_grid(~Patient)+ggtitle("Macrophage SPP1+")


CC_Result  %>%
  liana_dotplot(source_groups = c("Macrophage_SELENOP+"),
                target_groups = c("CD8 T cell","CD4 T cell","Tumor"),ntop = 50)+
  facet_grid(~Patient)+ggtitle("Macrophage SELENOP+")


