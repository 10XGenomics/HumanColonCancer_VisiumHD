#### Visium HD Manuscript
## Supplemental Figure 4

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


# Read previously created objects (see Sections/FlexAnalysis.R)
ColonCancer_Flex<-readRDS('~/Outputs/Flex/FlexSeuratV5.rds')
TumorObject<-subset(ColonCancer_Flex,subset=Level2%in%c("Tumor I","Tumor II","Tumor III","Tumor IV","Tumor V") & Condition=="Tumor")

# Standard processing
TumorObject<-NormalizeData(TumorObject)
TumorObject<-FindVariableFeatures(TumorObject)
TumorObject<-ScaleData(TumorObject)

TumorObject<-RunPCA(TumorObject)
ElbowPlot(TumorObject,ndims = 40)
TumorObject<-RunUMAP(TumorObject,dims=1:20)

TumorObject<-FindNeighbors(TumorObject,dims = 1:20)
TumorObject<-FindClusters(TumorObject,resolution=0.4)


# UMAP plots
P2<-DimPlot(TumorObject,label=T,label.size = 4,group.by = "Level2")+ggtitle("Level 2 Clustering")
P3<-DimPlot(TumorObject,label=T,label.size = 4,group.by = "Sample")+ggtitle("Sample")
P2+P3

# Get Markers
TumorObject$Level2<-factor(TumorObject$Level2,levels = c("Tumor I","Tumor II","Tumor III","Tumor IV","Tumor V"))
TumorObject<-SetIdent(TumorObject,value = "Level2")
Markers<-FindAllMarkers(TumorObject,logfc.threshold = 0.2,min.diff.pct = 0.2,only.pos = T)
Markers<-Markers[Markers$p_val_adj<0.05,]

Genes<-Markers %>% group_by(cluster) %>% arrange(-avg_log2FC) %>% slice(1:12) %>% pull(gene)

Plot<-DotPlot(TumorObject,features = unique(Genes),dot.min = 0.15,col.min = -2,col.max = 2,idents=levels(Markers$cluster))
Plot<-Plot$data
Plot<- Plot %>% na.omit()

palGrad <- colorRampPalette(c(hcl(0,100, c(20,100)), hcl(240,100,c(100,20))))

# Modified DotPlot
ggplot(Plot,aes(x=features.plot,y=id,size=pct.exp,color=avg.exp.scaled))+geom_point()+coord_flip()+
  theme_classic()+scale_colour_gradientn(colors=palGrad(50),limits=c(-2,2))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("")+guides(size = guide_legend(title = "Percentage Expressed"),color = guide_colorbar(title = "Scaled Expression"))+
  ggtitle("Tumor Subtypes Comparison")

# Feature Plots
#FeaturePlot(TumorObject,order = T,ncol = 3,cols=c("lightgray","firebrick1"))

Genes<-c("LGR5","FGGY","REG1A","MUC17","DAB2","MME","LCN2","MUC5B","MMP7","UBD","CEACAM5","CEACAM6")
MD<-TumorObject@meta.data
MD<-cbind(MD,TumorObject@reductions$umap@cell.embeddings)
Exp<-FetchData(TumorObject,Genes)
MD<-cbind(MD,Exp)

AllPlots<-vector("list",length=length(Genes))
names(AllPlots)<-Genes
for(GeneX in Genes)
{
  MD$Expression<-Exp[,GeneX]
  
  AllPlots[[GeneX]]<-ggplot(MD[order(MD$Expression,decreasing = F),],aes(x=umap_1,y=umap_2,color=Expression))+geom_point(size=0.7)+
    scale_color_gradient(low = "lightgray",high="red",limits=c(0,6))+theme(axis.text = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                                                           axis.line = element_blank(),axis.ticks = element_blank(),panel.border = element_blank(),
                                                                           panel.background = element_blank())+xlab("")+ylab("")+ggtitle(GeneX)
}

wrap_plots(AllPlots,ncol = 3,nrow = 4,byrow = T,guides = "collect")


