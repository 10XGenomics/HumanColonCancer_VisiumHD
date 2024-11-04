#### Visium HD Manuscript

## Run SpaceXR (Deconvolution)
# This script can be used as a skeleton to run deconvolution for any sample.

## load libraries
library(arrow)
library(Seurat)
library(spacexr) # Should be the modified version [see PR: ADD LINK TO PR]

FlexOutPath <- "~/AggrOutput/outs" # Path to cellranger aggr output folder
ColonFlex.data <- Read10X_h5(paste0(FlexOutPath,"filtered_feature_bc_matrix.h5")) 

## Load Reference Data
FlexRef<-Read10X_h5(paste0(FlexOutPath,"filtered_feature_bc_matrix.h5"))
MetaData<-readRDS('~/Outputs/Flex/FlexSeuratV5_MetaData.rds') # See FlexSingleCell.R if not generated.


# spacexr restriction, Clusters with > 25 cells
KpIdents<-names(which(table(MetaData$Level2)>25))
MetaData<-MetaData[MetaData$Level2%in%KpIdents,]
FlexRef<-FlexRef[,MetaData$Barcode]

## Fix cell type labels as spacexr doesn't allow special characters (i.e. spaces)
CTRef<-MetaData$Level2
CTRef<-gsub("/","_",CTRef)
CTRef<-as.factor(CTRef)
names(CTRef)<-MetaData$Barcode

## Build reference object
reference <- Reference(FlexRef[,names(CTRef)], CTRef , colSums(FlexRef))

# Deconvolve HD Data
counts<-Read10X_h5("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
coords<-read_parquet("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet",as_data_frame = TRUE)
rownames(coords)<-coords$barcode
coords<-coords[colnames(counts),]
coords<-coords[,3:4]
nUMI <- colSums(counts)

puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts)

myRCTD <- create.RCTD(puck, reference, max_cores = 12)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

saveRDS(myRCTD,file="~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds")