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


# Read single cell MetaData 
MetaData<-read.csv('SingleCell_MetaData.csv.gz') %>% na.omit()
MetaData$Level1<-factor(MetaData$Level1,levels = sort(unique(MetaData$Level1)))
MetaData$Level2<-factor(MetaData$Level2,levels = sort(unique(MetaData$Level2)))
MetaData$Condition<-gsub("P[0-9]","",MetaData$Patient)
MetaData$PatientID<-substr(MetaData$Patient,1,2)

# Create color palettes
ColsL1<-paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(MetaData$Level1))]
names(ColsL1)<-levels(MetaData$Level1)

ColsL2<-ColorPalette()

# Level 1 Proportion Plot
ProportionsL1<-as.data.frame(table(MetaData$Level1))
ProportionsL1<-ProportionsL1[order(ProportionsL1$Freq),]
ProportionsL1$Var1<-factor(ProportionsL1$Var1,levels = ProportionsL1$Var1)
ProportionsL1_Plot<-ggplot(ProportionsL1,aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat="identity")+coord_flip()+theme_classic()+
  geom_text(aes(label=Freq), hjust = 1)+scale_fill_manual(values=ColsL1)+NoLegend()+xlab("")+ylab("Frequency")


# Label Positions for UMAP Aggr'd dataset
LabelCoords_L1<-MetaData %>% group_by(Level1) %>% summarise(UM1 = median(UMAP1, na.rm = TRUE),
                                                            UM2 = median(UMAP2, na.rm = TRUE))

LabelCoords_L2<-MetaData %>% group_by(Level2) %>% summarise(UM1 = median(UMAP1, na.rm = TRUE),
                                                            UM2 = median(UMAP2, na.rm = TRUE))

# Level 1 UMAP plot
Level1UMAP<-MetaData %>% filter(QCFilter=="Keep") %>% ggplot(aes(x=UMAP1,y=UMAP2,color=Level1))+
  geom_scattermore(pointsize = 2,pixels = rep(2000,2))+theme_classic()+
  scale_color_manual(values=ColsL1)+xlab("")+ylab("")+NoLegend()+
  geom_text_repel(data=LabelCoords_L1,aes(x=UM1,y=UM2,label=Level1),fontface='bold',col="black")+
  theme(axis.text = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_blank(),axis.ticks = element_blank())+
  ggtitle("Flex Single cell",subtitle = "Level 1 Annotations")

ProportionsL1_Plot + Level1UMAP

# Level 2 UMAP plot
Level2UMAP<-MetaData %>% filter(QCFilter=="Keep") %>% ggplot(aes(x=UMAP1,y=UMAP2,color=Level2))+
  geom_scattermore(pointsize = 2,pixels = rep(2000,2))+theme_classic()+
  scale_color_manual(values=ColsL2)+xlab("")+ylab("")+NoLegend()+
  geom_text_repel(data=LabelCoords_L2,aes(x=UM1,y=UM2,label=Level2),fontface="bold",col="black")+
  theme(axis.text = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_blank(),axis.ticks = element_blank())+
  ggtitle("Flex Single cell",subtitle = "Level 2 Annotations")

# UMAP plots by Condition (rows) and Patient (columns)
MetaData %>% ggplot(aes(x=UMAP1,y=UMAP2,color=Level2))+
  geom_scattermore(pointsize = 4,pixels = rep(2000,2))+theme_classic()+
  scale_color_manual(values=ColsL2)+xlab("")+ylab("")+NoLegend()+
  theme(axis.text = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_blank(),axis.ticks = element_blank())+
  facet_grid(cols = vars(PatientID),rows = vars(Condition))
