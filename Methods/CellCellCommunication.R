### Visium HD Manuscript

## METHOD:  Cell Cell Communication Analysis

# Load packages
library(cli)
library(dplyr)
library(liana)
library(SeuratObject)
library(Seurat)
library(BPCells)
library(ggplot2)
library(liana)
library(circlize)
library(igraph)
library(distances)

source("~/Methods/AuxFunctions.R")

# Create Sample Info DF
SampleInfo<-data.frame(Patient=c("P1CRC","P2CRC","P5CRC"),
                       PathSR=c("~/P1CRC/outs/",
                            "~/P2CRC/outs/",
                            "~/P5CRC/outs/"),
                       PathDeconvolution=c("~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds",
                                "~/Outputs/Deconvolution/PatientCRC2_Deconvolution_HD.rds",
                                "~/Outputs/Deconvolution/PatientCRC5_Deconvolution_HD.rds"),
                       PathPeriphery=c("~/Outputs/Deconvolution/PatientCRC1_PeripheryBCs.rds",
                                   "~/Outputs/Deconvolution/PatientCRC2_PeripheryBCs.rds",
                                   "~/Outputs/Deconvolution/PatientCRC5_PeripheryBCs.rds"),
                       Tumor=c("Tumor II","Tumor III","Tumor IV"))


AllResults<-vector("list",length=nrow(SampleInfo))
PlotA<-vector("list",length=nrow(SampleInfo))
PlotB<-vector("list",length=nrow(SampleInfo))
names(AllResults)<-names(PlotA)<-names(PlotB)<-SampleInfo$Patient

for(Patient in SampleInfo$Sample)
{
  #Generate MetaData
  index<-which(SampleInfo$Patient==Patient)
  bcDF<-GenerateSampleData(SampleInfo$PathSR[index])$bcs
  Deconv<-readRCTD(SampleInfo$PathDeconvolution[index])
  bcDF<-AddDeconvolutionInfo(bcDF,Deconv)
  bcDF<-bcDF %>% filter(tissue==1)
  
  # Add Periphery Results
  PeripheryBCs<-readRDS(SampleInfo$PathPeriphery[index])
  
  # Add periphery results to data.frame
  bcDF$Periphery<-NA
  bcDF$Periphery[bcDF$barcode%in%PeripheryBCs]<-"50 micron"
  bcDF$Periphery[is.na(bcDF$Periphery)]<-"Tissue"
  bcDF$Periphery[bcDF$Periphery=="Tissue" & bcDF$DeconvolutionLabel1==SampleInfo$Tumor[index]]<-"Tumor"
  
  # Create Seurat Object (8um)
  object<-Read10X_h5(paste0(SampleInfo$PathSR,"binned_outputs/square_008um/filtered_feature_bc_matrix.h5"}))
  CommonBC<-intersect(bcDF$barcode,colnames(object))
  object<-CreateSeuratObject(object[,CommonBC],meta.data=bcDF[match(CommonBC,bcDF$barcode),])
  
  # Analysis To identify Macrophage subtypes [TODO: Replace with saved results]
  
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
  scales <- rjson::fromJSON(file = paste0(PathSR,"/binned_outputs/square_008um/spatial/scalefactors_json.json"))
  
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
  
  # Plot
  PlotA[[Patient]]<-bcDF  %>% filter(DeconvolutionClass=="singlet") %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Selection)) +
    geom_scattermore(pointsize = 3,pixels = rep(2000,2))+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    scale_color_manual(values=c("dodgerblue","firebrick2"),na.value = "lightgray")+
    labs(color="Group")+ggtitle("")
  
  
  # Create Seurat Object (Again?)
  object<-Read10X_h5(SampleInfo$H5[index])
  CommonBC<-intersect(bcDF$barcode,colnames(object))
  object<-CreateSeuratObject(object[,CommonBC],meta.data=bcDF[match(CommonBC,bcDF$barcode),])
  object<-subset(object, subset = Selection %in% c("Tumor","50 micron") & DeconvolutionClass=="singlet")
  object <- NormalizeData(object)
  object<-SetIdent(object,value = "CellType")
  
  # Run C-C using Liana
  liana_result <- liana_wrap(object)%>% liana_aggregate()%>% filter(aggregate_rank <= 0.01)
  liana_result$Patient<-Patient
  AllResults[[Patient]]<-liana_result
  
  # Top Interactions
  #ColX<-colorRampPalette(c("white","firebrick1"))(20)
  #pheatmap(table(liana_result$source,liana_result$target),color = ColX,border_color = "gray60")
  
  CTxs<-c("Macrophage_SPP1+","Macrophage_SELENOP+","CD4 T cell","CD8 Cytotoxic T cell",SampleInfo$Tumor[SampleInfo$Sample==Patient])
  
  # Interaction Graph
  PlotB[[Patient]]<-PlotInteractionGraph(liana_result,CellTypes = CTxs,
                               Colors = c(ColorPalette(),"Macrophage_SPP1+"="cyan2","Macrophage_SELENOP+"="orchid3"))
}



CC_Result<-do.call(rbind,AllResults)

iix<-(CC_Result$Patient==SampleInfo$Sample[1] & CC_Result$target==SampleInfo$Tumor[1]) |
  (CC_Result$Patient==SampleInfo$Sample[2] & CC_Result$target==SampleInfo$Tumor[2]) |
  (CC_Result$Patient==SampleInfo$Sample[3] & CC_Result$target==SampleInfo$Tumor[3])

CC_Result$target[iix]<-"Tumor"

CC_Result  %>%
  liana_dotplot(source_groups = c("Macrophage_SPP1+"),
                target_groups = c("CD8 Cytotoxic T cell","CD4 T cell","Tumor"),ntop = 50)+
  facet_grid(~Patient)+ggtitle("Macrophage SPP1+")


CC_Result  %>%
  liana_dotplot(source_groups = c("Macrophage_SELENOP+"),
                target_groups = c("CD8 Cytotoxic T cell","CD4 T cell","Tumor"),ntop = 50)+
  facet_grid(~Patient)+ggtitle("Macrophage SELENOP+")


# Access Plots
PlotA[[1]]
PlotB[[1]]
