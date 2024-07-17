## Cell Cell Communication Analysis
library(cli,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(dplyr,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(liana,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(SeuratObject,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(Seurat,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")
library(BPCells,lib.loc = "/mnt/home/juanpablo.romerorioj/yard/Seurat5/")

#### One Sample Only (P1CRC)

SampleInfo<-data.frame(Sample=c("P1CRC","P2CRC","P5CRC"),
                       H5=c("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P2CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
                            "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P5CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"),
                       Deconv=c("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/DeconvolutionResults_P1CRC.csv.gz",
                                "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P2CRC/DeconvolutionResults_P2CRC.csv.gz",
                                "/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P5CRC/DeconvolutionResults_P5CRC.csv.gz"))

MD<-readRDS("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/MetaDataCombined.rds")
MetaData<-vector("list",length = nrow(SampleInfo))
Matrices<-vector("list",length = nrow(SampleInfo))

PathX<-paste0("/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Revisions/Clustering/",SampleInfo$Sample[1])

mat <- open_matrix_dir(dir = PathX)
colnames(mat)<-paste0(colnames(mat),"_",SampleInfo$Sample[1])

#mat <- ConvertEnsembleToSymbol(mat = mat, species = "human")
Genes<-read.delim(file="/mnt/home/juanpablo.romerorioj/deck/VisiumHD/Online/P1CRC/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header = F)
rownames(mat)<-make.unique(Genes$V2[match(rownames(mat),Genes$V1)])

MD <- MD[MD$Patient==SampleInfo$Sample[1],]

object <- CreateSeuratObject(counts = mat, meta.data = MD)

Periphery<-readRDS("/mnt/home/juanpablo.romerorioj/yard/Projects/VisiumHD/Results/HD/Periphery_50_100_F00142484_All.rds")
object$Periphery<-gsub("_P[0-9]CRC","",colnames(object))%in%Periphery


DecTmp<-read.delim(SampleInfo$Deconv[1],sep=",") %>% na.omit()
DecTmp$barcode<-paste0(DecTmp$barcode,"_",SampleInfo$Sample[1])

object$DeconvolutionClass<-DecTmp$DeconvolutionClass[match(colnames(object),DecTmp$barcode)]
object$DeconvolutionLabel1<-DecTmp$DeconvolutionLabel1[match(colnames(object),DecTmp$barcode)]
object$DeconvolutionLabel2<-DecTmp$DeconvolutionLabel2[match(colnames(object),DecTmp$barcode)]

object <- NormalizeData(object)

TME<-subset(object,cells = colnames(object)[object$Periphery & !is.na(object$DeconvolutionLabel1)])
TME<-SetIdent(TME,value = "DeconvolutionLabel2")

liana_result <- liana_wrap(TME)

liana_result <- liana_result %>% liana_aggregate()

liana_result %>% liana_dotplot(source_groups = c("Tumor III"),target_groups = c("CAF","CD4 T cell","CD8 Cytotoxic T cell","Mature B","Memory V","Plasma","Macrophage",
                                                                                "Proliferating Immune II","Proliferating Macrophages"),ntop = 40)

liana_result %>% liana_dotplot(source_groups = c("Tumor III"),target_groups = c("Macrophage","Proliferating Macrophages"),ntop = 40)


liana_trunc <- liana_result %>% filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
pheatmap(table(liana_trunc$source,liana_trunc$target))


Idents<-c("Macrophage", "CD 4 T cell", "CD8 Cytotoxic T cell","Tumor III")
p <- chord_freq(liana_trunc,
                source_groups = Idents,
                target_groups = Idents)
