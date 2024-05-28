#### Visium HD Manuscript

# Auxiliary Functions

# Function to read a Visium HD output directory and generate a data.frame with all bins and a images_tibble for plotting.
GenerateSampleData<-function(PATH,size="008um")
{
  images_tibble <- make_images_tibble(PATH) %>% bind_rows()
  
  tissue_positions_path <- list.files(paste0(PATH, "/binned_outputs/square_",size,"/spatial"), pattern = "tissue_positions*", full.names = TRUE)
  tissue_positions_file <- basename(tissue_positions_path)
  tissue_positions_df <- arrow::read_parquet(tissue_positions_path)%>%
    dplyr::rename_with(~ c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
  path_scales <- file.path(PATH, paste0("/binned_outputs/square_",size,"/spatial/scalefactors_json.json"))
  scales <- rjson::fromJSON(file = path_scales)
  path_clusters <- file.path(PATH, paste0("/binned_outputs/square_",size,"/analysis_csv/clustering/gene_expression_graphclust/clusters.csv"))
  
  images_tibble_mod <- images_tibble %>% filter(Path == PATH)
  
  
  if(file.exists(path_clusters))
  {
    clusters <- read.csv(path_clusters)
    path_umap <- file.path(PATH, paste0("/binned_outputs/square_",size,"/analysis_csv/umap/gene_expression_2_components/projection.csv"))
    umap <- read.csv(path_umap)
    
    bcs <- tissue_positions_df %>% mutate(imagerow_scaled = imagerow * 
                                            scales$tissue_lowres_scalef, imagecol_scaled = imagecol * 
                                            scales$tissue_lowres_scalef, imagerow_scaled_round = round(imagerow * 
                                                                                                         scales$tissue_lowres_scalef), imagecol_scaled_round = round(imagecol * 
                                                                                                                                                                       scales$tissue_lowres_scalef), tissue = as.factor(tissue)) %>% 
      left_join(clusters, by = c(barcode = "Barcode")) %>% 
      left_join(umap, by = c(barcode = "Barcode")) %>% 
      mutate(height = images_tibble_mod$height, 
             width = images_tibble_mod$width)
    
    
    
  }else{
    clusters<-NA
    umap<-NA
    
    bcs <- tissue_positions_df %>% mutate(imagerow_scaled = imagerow * 
                                            scales$tissue_lowres_scalef, imagecol_scaled = imagecol * 
                                            scales$tissue_lowres_scalef, imagerow_scaled_round = round(imagerow * 
                                                                                                         scales$tissue_lowres_scalef), imagecol_scaled_round = round(imagecol * 
                                                                                                                                                                       scales$tissue_lowres_scalef), tissue = as.factor(tissue)) %>% 
      mutate(height = images_tibble_mod$height, 
             width = images_tibble_mod$width)
  }
  
  
  return(list(images_tibble=images_tibble,bcs=bcs))
}

# Function used to generate images_tibble 
make_images_tibble<-function (PATH) 
{
  image <- get_spatial_files(PATH, "tissue_lowres_image")
  grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  images_tibble <- tibble(Path = factor(PATH), grob = list(grobs), height = nrow(image), width = ncol(image))
  return(images_tibble)
}


# Create paths and read different spatial files found in output directory
get_spatial_files<-function (PATH, type) 
{
  if (type == "tissue_lowres_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_lowres_image.png", sep = ""))
  }
  if (type == "tissue_hires_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_hires_image.png", sep = ""))
  }
  if (type == "tissue_positions_list") {
    x <- read.csv(paste(PATH, "/spatial/tissue_positions_list.csv", sep = ""), col.names = c("barcode", "tissue", "row",  "col", "imagerow", "imagecol"), header = F)
  }
  if (type == "aligned_fiducials") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/aligned_fiducials.jpg", sep = ""))
  }
  if (type == "detected_tissue_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/detected_tissue_image.jpg", sep = ""))
  }
  if (type == "scales") {
    path_scales <- paste(PATH, "/spatial/scalefactors_json.json", sep = "")
    x <- rjson::fromJSON(file = path_scales)
  }
  return(x)
}

# Define color palettes used throughout the manuscript
ColorPalette<-function()
{
  Colors<-c(paletteer::paletteer_d("ggsci::default_igv")[1:39],"black","azure4")
  names(Colors)<-c('Tumor III','Plasma','Macrophage','CD4 T cell','CAF','vSM','Mature B','Endothelial','Tumor I','CD8 Cytotoxic T cell',
                   'Enterocyte','Neutrophil','Proliferating Immune II','Pericytes','Smooth Muscle','Myofibroblast',
                   'Tumor II','Fibroblast','Goblet','Lymphatic Endothelial','Tumor V','Proliferating Macrophages','SM Stress Response',
                   'NK','cDC I','Tumor IV','Proliferating Fibroblast','Epithelial','Tuft','Mast','Unknown III (SM)',
                   'Adipocyte','mRegDC','Enteric Glial','pDC','Vascular Fibroblast','Neuroendocrine','Memory B','Unknwon I (Immune)',"Undetermined","cell")
  
  return(Colors)
  
}

# Read deconvolution results. Wrapped in a function as large RCTD objects slow down the R session
readRCTD<-function(PATH)
{
  Object<-readRDS(PATH)
  Info <- Object@results$results_df
  weights <- get_doublet_weights_modified(Info,Object@results$weights_doublet, Object@cell_type_info$info[[2]])
  Results <- list(DF = Info, Weights = t(weights))
  
  return(Results)
}

# From spaceXR outputs create weights matrix
get_doublet_weights_modified <- function(ResultDF,Weights,CellTypes) {
  
  barcodes <- rownames(ResultDF)
  my_beta <- matrix(0, nrow = length(barcodes), ncol = length(CellTypes))
  rownames(my_beta) <- barcodes
  colnames(my_beta) <- CellTypes
  
  indexRow_Certain<-which(ResultDF$spot_class %in% c('singlet', 'doublet_certain'))
  indexCol_Certain<-match(ResultDF[indexRow_Certain,'first_type'],colnames(my_beta))
  my_beta[cbind(indexRow_Certain,indexCol_Certain)] <- Weights[indexRow_Certain,1]
  
  indexRow_Doublet<-which(ResultDF$spot_class == "doublet_certain")
  indexCol_Doublet<-match(ResultDF[indexRow_Doublet,'second_type'],colnames(my_beta))
  my_beta[cbind(indexCol_Doublet)] <- Weights[indexRow_Doublet,2]
  
  return(my_beta)
  
  
}

# Add deconvolution results to data.frame created with GenerateSampleData
AddDeconvolutionInfo<-function(BCS,Results,AddWeights=FALSE)
{
  ResultsDF<-Results$DF
  
  index<- match(rownames(ResultsDF),BCS$barcode)
  
  BCS$DeconvolutionClass<-NA
  BCS$DeconvolutionClass[index]<-as.vector(ResultsDF$spot_class)
  
  BCS$DeconvolutionLabel1<-NA
  BCS$DeconvolutionLabel1[index]<-ResultsDF$first_type
  
  BCS$DeconvolutionLabel2<-NA
  BCS$DeconvolutionLabel2[index]<-ResultsDF$scond_type
  
  if(AddWeights)
  {
    Weights<-Results$Weights
    
    index<- match(colnames(Weights),BCS$barcode)
    
    Names<-gsub(" ","",rownames(Weights))
    
    for(jj in 1:nrow(Weights))
    {
      BCS[,Names[jj]]<-NA
      BCS[index,Names[jj]]<-Weights[jj,]
    }
  }
  
  
  return(BCS)
  
}

# Add expression of genes to a data.frame generated by GenerateSampleData from a Seurat Object
AddExpression<-function(Barcodes,Seurat,Genes)
{
  
  Exp<-FetchData(Seurat,Genes)
  
  for(Gx in Genes)
  {
    Barcodes[,Gx]<-NA
    Barcodes[match(rownames(Exp),Barcodes$barcode),Gx]<-Exp[,Gx]
  }
  
  return(Barcodes)
  
}

# Create a square of a given size (in microns) whose center is a given barcode
GetSquare<-function(Spot,SizeMicrons,BarcodeDF,binsize=8)
{
  Xcenter<-BarcodeDF$col[match(Spot,BarcodeDF$barcode)]
  Ycenter<-BarcodeDF$row[match(Spot,BarcodeDF$barcode)]
  
  AddFactor<-round(SizeMicrons/(2*binsize))
  
  Xmin<-Xcenter-AddFactor
  Xmax<-Xcenter+AddFactor
  
  Ymin<-Ycenter-AddFactor
  Ymax<-Ycenter+AddFactor
  
  SquareSection<-BarcodeDF %>% filter(col >= Xmin & col <= Xmax & row >= Ymin & row <= Ymax) %>% pull(barcode)
  
  return(SquareSection)
  
}

# Plot gene expression from a data.frame 
PlotExpression<-function(barcodes,Gene,ptsize=2,shape="circle")
{
  barcodes$Expression<-barcodes %>% pull(Gene)
  
  if(shape=="circle")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Expression)) +
      geom_scattermore(pointsize = ptsize,pixels = rep(2000,2))+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_color_gradient(low="lightgray",high = "red")+
      labs(color=paste0(Gene))
    
  }else if(shape=="square")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=Expression)) +
      geom_point(shape=22,size=ptsize,color=alpha("black",0),stroke=0.25)+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_fill_gradient(low="lightgray",high = "red")+
      labs(fill=paste0(Gene))
    
  }else{
    stop("Wrong Shape")
  }
}

# Create an additional color palette to be used
ColorsClusters<-function()
{
  x <- c("128 128 128",
         "89 106 55",
         "150 86 53",
         "116 20 12",
         "42 98 24",
         "128 128 38",
         "69 61 134",
         "96 177 119",
         "55 126 127",
         "85 128 176",
         "4 0 123",
         "165 204 79",
         "210 167 65",
         "127 23 134",
         "234 51 35",
         "93 203 207",
         "240 146 53",
         "254 255 84",
         "183 45 130",
         "12 0 197",
         "117 251 76",
         "117 251 141",
         "202 49 66",
         "86 188 249",
         "232 167 108",
         "147 44 231",
         "190 253 91",
         "234 51 247",
         "69 142 247",
         "205 118 147",
         "238 230 151",
         "234 134 119",
         "212 162 217",
         "182 215 228",
         "120 105 230",
         "224 135 232",
         "175 249 162",
         "160 251 214",
         "250 228 200")
  
  
  
  Cols<-sapply(strsplit(x, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))
  
  return(Cols)
  
}
  
# Function to select the barcodes that are within a given distance from a given cell type (cluster)
SelectPeripheryDiscrete<-function(bcs,CellType,distance=50,PATH)
{
  
  SelectedBCs<-bcs %>% filter(DeconvolutionLabel1==CellType)
  
  Result<-GetSlice(SelectedBCs$barcode,distance,bcs,PATH,CellT=CellType)
  
  if(length(distance)>1)
  {
    for(jj in 1:length(Result))
    {
      Result[[jj]]<-Result[[jj]][Result[[jj]]%!in%SelectedBCs$barcode]
    }
    
    return(Result)
    
  }else{
    Result<-Result[Result%!in%SelectedBCs$barcode]
    
    return(Result)
  }
  
  
}

# Equivalent to GetSquare but for a circle instead.
GetSlice<-function(Spot,SizeMicrons,BarcodeDF,PATH,CellT=NA,size="008um")
{
  
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  Scale<-(SizeMicrons*scales$spot_diameter_fullres)/as.numeric(unlist(strsplit(size,"um"))[1])

  Index<-match(Spot,BarcodeDF$barcode)
  Result<-vector("list",length=length(Index))
  
  for(jj in 1:length(Index))
  {
    Distance<-sqrt(((BarcodeDF$imagecol-BarcodeDF$imagecol[Index[jj]])^2) + ((BarcodeDF$imagerow-BarcodeDF$imagerow[Index[jj]])^2))
    BarcodeDF$Distance<-Distance
    
    if(!is.na(CellT))
    {
      ValTh <- sum(BarcodeDF$DeconvolutionLabel1[BarcodeDF$Distance<min(Scale)]==CellT,na.rm = T)
      if(ValTh < 25)
      {
        next
      }
    }
    
    if(length(Scale)>1)
    {
      Result[[jj]]<-lapply(Scale,function(X){return(BarcodeDF$barcode[BarcodeDF$Distance < X])})
    }else{
      
      Result[[jj]]<-BarcodeDF$barcode[BarcodeDF$Distance < Scale]
    }
    
  }
  
  if(length(Scale)>1)
  {
    Rxx<-vector("list",length=length(Scale))
    names(Rxx)<-as.character(SizeMicrons)
    
    for(ii in 1:length(Scale))
    {
      Rxx[[ii]]<-lapply(Result,function(X){return(X[[ii]])})
      Rxx[[ii]]<-unique(unlist(Rxx[[ii]]))
    }
    
    return(Rxx)
    
  }else{
    Result<-unique(unlist(Result))
    return(Result)
  }
  
}

# Function to plot enrichR results as a barplot
EnrichRBarPlot<-function(Markers,DataBase,TermsX=10,PTh=1e-3,GO=F,colsT=c("firebrick1","dodgerblue"))
{
  Cluss<-as.vector(unique(Markers$cluster))
  NClusters<-length(unique(Markers$cluster))
  
  Result<-vector("list",length = NClusters)
  
  for(jj in 1:length(Cluss))
  {
    
    Res<-enrichr(Markers[Markers$cluster==Cluss[jj],"gene"],DataBase)
    
    if(nrow(Res[[1]])>0)
    {
      Result[[jj]]<-cbind(Res[[1]],Cluss[jj])
    }
    
  }
  
  Result<-do.call(rbind,Result)
  
  if(GO)
  {
    Result$Term<-trimws(sapply(strsplit(Result$Term,"[(]"),function(X){return(X[1])}))
  }
  
  
  colnames(Result)[10]<-"Cluster"
  
  Result<-Result[Result$P.value<PTh,]
  
  Result<-Result %>% group_by(Cluster) %>% slice_max(order_by = P.value, n = TermsX)
  Result$FDR<--log(Result$Adjusted.P.value)
  Result$FDR<-Result$FDR*ifelse(as.numeric(as.factor(Result$Cluster))==1,1,-1)
  Result<-Result[order(Result$FDR),]
  Result$Term<-factor(Result$Term,levels=unique(Result$Term))
  
  Plot<-ggplot(Result,aes(x=Term,y=FDR,fill=Cluster))+geom_bar(stat="identity")+
    theme_classic()+scale_fill_manual(values=colsT)+
    xlab("")+ylab("log(Adj. pvalue)")+theme(axis.title.y=element_blank(),
                                            axis.text.y=element_blank(),
                                            axis.ticks.y=element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_text(face="bold"))+
    geom_text(aes(label = Term,hjust = ifelse(FDR < 0, 0, 1),vjust = 0.5),y=0,size=3)+coord_flip()
  
  return(Plot)
  
}

# Function to create density map and identifies enriched region of a given cell type.
# BarcodeSet is used to restrict the area to a given collection of barcodes, if missing then
# the whole section is used.
PlotDensity<-function(DF,CellType,nBins=3,ptsize=3,Tumor=NA,BarcodeSet=NA)
{
  require(wesanderson)
  
  if(length(CellType)>1)
  {
    stop("Pass only 1 cell type to CellType argument")
  }
  
  # Create DF for density2d
  if(all(!is.na(BarcodeSet)))
  {
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1%in%CellType & DF$barcode %in% BarcodeSet,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }else{
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1%in%CellType,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }
  
  
  DF3 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==Tumor,]
  DF3 <- DF3 %>% na.omit()
  DF3$XX <- DF3$imagecol_scaled
  DF3$YY <- DF3$imagerow_scaled
  
  DF4 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1%in%CellType,]
  DF4 <- DF4 %>% na.omit()
  DF4$XX <- DF4$imagecol_scaled
  DF4$YY <- DF4$imagerow_scaled
  DF4$Grouping<-DF4$DeconvolutionLabel1
  
  PlotX<-DF %>% filter(tissue == "1") %>% na.omit() %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled)) +
    geom_scattermore(pointsize = ptsize,pixels = rep(2000,2),col="lightgray")+
    geom_scattermore(data=DF3,pointsize=ptsize,pixels = rep(2000,2),col="gray65")+
    geom_scattermore(data=DF4,pointsize=ptsize,pixels = rep(2000,2),col="red")+
    geom_density_2d(data=DF2,bins=nBins,linewidth=1.2,linetype=1,aes(colour=after_stat(level)),
                    contour_var = "ndensity")+
    scale_colour_gradientn(colours=wes_palette("Zissou1", 20, type = "continuous"),limit=c(0,1))+
    #scale_colour_gradient2(low="dodgerblue",mid = "dodgerblue2" ,high="dodgerblue4",limit=c(0,1))+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    ggtitle(CellType)+
    labs(color="Scaled Density")
  
  return(PlotX)
  
}

# Collect the barcodes that are within an enriched for a given celltype/cluster
SelectEnrichedRegion<-function(CellType,bcs,PATH,N=5,Area=200)
{
  bcs <- bcs %>% filter(tissue==1)
  
  bcsCT<-bcs %>% filter (DeconvolutionLabel1%in%CellType)
  
  
  Kernel<-MASS::kde2d(bcsCT %>% pull(imagecol),
                      bcsCT %>% pull(imagerow),
                      n = 200)
  
  Peaks<-findPeaks3D(Kernel$z,N)
  Peaks<-data.frame(X=sapply(Peaks,function(X){return(X$x)}),
                    Y=sapply(Peaks,function(X){return(X$y)}),
                    Z=sapply(Peaks,function(X){return(X$z)}))
  
  Peaks$X<-Kernel$x[Peaks$X]
  Peaks$Y<-Kernel$y[Peaks$Y]
  
  Result<-vector("list",length=nrow(Peaks))
  
  for(jj in 1:nrow(Peaks))
  {
    BCTmp<-bcs
    BCTmp<-GetDistance(BCTmp,PATH,XY=c(Peaks$X[jj],Peaks$Y[jj]))
    BCTmp<-BCTmp[order(BCTmp$DistanceMicrons),]
    BarcodeX<-BCTmp$barcode[order(BCTmp$DistanceMicrons)[1]]
    Result[[jj]]<-GetSlice(BarcodeX,Area,BCTmp,PATH)
    
  }
  
  return(Result)
  
}

# Find peaks of a 3d function. Used by SelectEnrichedRegion
findPeaks3D <- function(matrix, N) {
  # Function to check if a point is a local maximum
  isLocalMaximum <- function(mat, x, y) {
    neighbors <- c(
      mat[x-1, y], mat[x+1, y], mat[x, y-1], mat[x, y+1],
      mat[x-1, y-1], mat[x-1, y+1], mat[x+1, y-1], mat[x+1, y+1]
    )
    # Remove NA values (edges)
    neighbors <- neighbors[!is.na(neighbors)]
    return(all(mat[x, y] > neighbors))
  }
  
  numRows <- nrow(matrix)
  numCols <- ncol(matrix)
  
  # List to store peaks
  peaks <- list()
  
  for (x in 2:(numRows - 1)) {
    for (y in 2:(numCols - 1)) {
      if (isLocalMaximum(matrix, x, y)) {
        peaks <- c(peaks, list(list(x = x, y = y, z = matrix[x, y])))
      }
    }
  }
  
  # Sort peaks by height and select top N
  peaks <- peaks[order(sapply(peaks, function(peak) -peak$z))]
  if (length(peaks) > N) {
    peaks <- peaks[1:N]
  }
  
  return(peaks)
}

# Function to get the distance to a given barcode or XY coordinate
GetDistance<-function(BCS,PATH,barcode=NA,XY=NA,size="008um")
{
  
  # Transform Scale
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  
  # Select Either BC or XY coordinate
  
  if(all(!is.na(barcode),all(!is.na(XY))))
  {
    stop("Either barcode or XY coordinate")
  }else if(!is.na(barcode) & all(is.na(XY)))
  {
    BC_Center<-barcode
    Index<-match(BC_Center,BCS$barcode)
    Distance<-sqrt(((BCS$imagecol-BCS$imagecol[Index])^2) + ((BCS$imagerow-BCS$imagerow[Index])^2))
    DistanceMicrons<-(Distance*8)/scales$spot_diameter_fullres
    BCS$DistanceMicrons<-DistanceMicrons
    
  }else if(is.na(barcode) & all(!is.na(XY)))
  {
    BC_Center<-barcode
    Index<-match(BC_Center,BCS$barcode)
    Distance<-sqrt(((BCS$imagecol-XY[1])^2) + ((BCS$imagerow-XY[2])^2))
    DistanceMicrons<-(Distance*8)/scales$spot_diameter_fullres
    BCS$DistanceMicrons<-DistanceMicrons
  }
  
  return(BCS)
  
}

# Function to find the Tumor (cluster) bins that are within a given distance from a selection of barcodes.
SelectTumor<-function(bcDF,Tumor,barcodes,PATH,DistVal=50,size="008um")
{
  # Get scale factors
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  
  # Get the center spot from the given barcodes
  RegionDF <- bcDF %>% filter(barcode%in%barcodes) %>% dplyr::select(barcode,imagerow,imagecol,DeconvolutionLabel1)
  
  MeanSpots<-RegionDF %>% summarise(Xval=(max(imagecol)+min(imagecol))/2,
                                    Yval=(max(imagerow)+min(imagerow))/2)
  
  Distance<-sqrt(((bcDF$imagecol-MeanSpots$Xval)^2) + ((bcDF$imagerow-MeanSpots$Yval)^2))
  CenterSpot<-bcDF$barcode[which.min(Distance)]
  
  # Select Slice 200 microns and Subset the data.frame
  SliceRegion<-GetSlice(CenterSpot,350,bcDF,PATH)
  
  RegionDF <- bcDF %>% filter(barcode%in%SliceRegion | barcode %in%barcodes) %>% dplyr::select(barcode,imagerow,imagecol,DeconvolutionLabel1)
  
  #Get Distance between given barcodes and Tumor Spots in the RegionDF
  RegionDF<- RegionDF %>% filter(DeconvolutionLabel1==Tumor | barcode %in% barcodes)
  TumBC<-RegionDF %>% filter(DeconvolutionLabel1==Tumor) %>% pull(barcode)
  
  XY_Data<-RegionDF[,c("imagecol","imagerow")]
  DistMat<-distances::distances(as.matrix(XY_Data))
  DistMat<-DistMat[match(barcodes,RegionDF$barcode),match(TumBC,RegionDF$barcode)]
  rownames(DistMat)<-barcodes
  colnames(DistMat)<-TumBC
  
  DistMat<-(DistMat*8)/scales$spot_diameter_fullres
  
  # Select for each Region spot the closest tumor Spot
  iix<-apply(DistMat,2,function(X){any(X<DistVal)})
  ClosestTumorSpot<-colnames(DistMat)[iix]
  
  return(ClosestTumorSpot)
  
}


# Plot colocalization of two cell types in  a sample (used for CD4 and CD8 t cells)
PlotColocalization<-function(DF,CellType1,CelltypesID,nBins=3,ptsize=3,Tumor=NA,BarcodeSet=NA,colX=c("darkorchid","gold"),option="CD8")
{
  require(wesanderson)
  require(ggnewscale)
  
  # Create DF for density2d
  if(all(!is.na(BarcodeSet)))
  {
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1==CellType1 & DF$barcode %in% BarcodeSet,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
    
  }else{
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1==CellType1,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }
  
  
  DF3 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==Tumor,]
  DF3 <- DF3 %>% na.omit()
  DF3$XX <- DF3$imagecol_scaled
  DF3$YY <- DF3$imagerow_scaled
  
  DF4A <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==CelltypesID[1],]
  DF4A <- DF4A %>% na.omit()
  DF4A$XX <- DF4A$imagecol_scaled
  DF4A$YY <- DF4A$imagerow_scaled
  DF4A$Grouping<-DF4A$DeconvolutionLabel1
  
  DF4B <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==CelltypesID[2],]
  DF4B <- DF4B %>% na.omit()
  DF4B$XX <- DF4B$imagecol_scaled
  DF4B$YY <- DF4B$imagerow_scaled
  DF4B$Grouping<-DF4B$DeconvolutionLabel1
  
  if(option=="CD8")
  {
    ColsGrad<-brewer.pal(9,"Blues")
  }else{
    ColsGrad<-brewer.pal(9,"Greens")
  }
  
  PlotX<-DF %>% filter(tissue == "1") %>% na.omit() %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled)) +
    geom_scattermore(pointsize = ptsize,pixels = rep(2000,2),col="lightgray")+
    geom_scattermore(data=DF3,pointsize=ptsize,pixels = rep(2000,2),col="gray65")+
    geom_scattermore(data=DF4A,pointsize=ptsize+2,pixels = rep(2000,2),col=colX[1])+
    geom_scattermore(data=DF4B,pointsize=ptsize+2,pixels = rep(2000,2),col=colX[2])+
    geom_scattermore(data=DF2,pointsize=ptsize+2,pixels = rep(2000,2),col=ifelse(option=="CD8","dodgerblue","forestgreen"))+
    geom_density_2d(data=DF2,bins=nBins,linewidth=1.2,linetype=1,aes(colour=after_stat(level)),contour_var = "ndensity")+
    scale_colour_gradientn(colours=ColsGrad,limit=c(0,1))+
    labs(color="Scaled Density")+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  
  return(PlotX)
  
}

# Transform barcode names between sizes
TransformBarcodes<-function(Barcodes,SizeOriginal,SizeNew)
{
  SizeO<-str_pad(paste0(SizeOriginal,"um"),5,pad="0")
  SizeN<-str_pad(paste0(SizeNew,"um"),5,pad="0")
  
  Barcodes<-strsplit(Barcodes,"[_-]")
  Barcodes<-data.frame(s="s",size=SizeO,
                       X=as.numeric(sapply(Barcodes,function(X){return(X[3])})),
                       Y=as.numeric(sapply(Barcodes,function(X){return(X[4])})),
                       end="-1")
  
  
  nBinsSide <- SizeOriginal / SizeNew
  
  Xmins <- Barcodes$X * nBinsSide
  Xmaxs <- (Barcodes$X * nBinsSide) + (nBinsSide - 1)
  Ymins <- Barcodes$Y * nBinsSide
  Ymaxs <- (Barcodes$Y * nBinsSide) + (nBinsSide - 1)
  
  Result<-lapply(1:nrow(Barcodes),function(jj)
  {
    Rx <- expand.grid(x = Xmins[jj]:Xmaxs[jj], y = Ymins[jj]:Ymaxs[jj])
    Rx$X <- Barcodes$X[jj]
    Rx$Y <- Barcodes$Y[jj]
    return(Rx)
  })
  
  Result <- do.call(rbind, Result)
  
  OldBC<-paste0("s_",SizeO,"_",str_pad(Result$X,5,pad="0"),"_",str_pad(Result$Y,5,pad="0"),"-1")
  NewBC<-paste0("s_",SizeN,"_",str_pad(Result$x,5,pad="0"),"_",str_pad(Result$y,5,pad="0"),"-1")
  
  Result<-data.frame(Original=OldBC,Transformed=NewBC)
  #Result<-split(Result$Transformed,Result$Original)
  
  return(Result)
  
}

# Used to provide parameters for the segmentation script
NucleiSegmentationScript<-function(BarcodeDF,TransformedDF)
{
  # Get Row and Column image limits for Nuclei Segmentation
  Rx<-data.frame(C1=round(min(BarcodeDF$imagecol[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
             C2=round(max(BarcodeDF$imagecol[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
             R1=round(min(BarcodeDF$imagerow[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
             R2=round(max(BarcodeDF$imagerow[match(unique(TransformedDF$Original),BarcodeDF$barcode)])))
  
  Script<-paste0(" ./NucleiSegmentation.py -i PATH/image.btf -r1 ",Rx$R1," -r2 ",Rx$R2," -c1 ",Rx$C1," -c2 ",Rx$C2," -x /PATH/TO/square_00Xum/ -o PATH/TO/OUTPUT")
  
  message(Script)
}

# Plotting function to at image
geom_spatial<-function (mapping = NULL, data = NULL, stat = "identity", position = "identity", 
          na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, ...) 
{
  GeomCustom <- ggproto("GeomCustom", Geom, setup_data = function(self, 
                                                                  data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data
  }, draw_group = function(data, panel_scales, coord) {
    vp <- grid::viewport(x = data$x, y = data$y)
    g <- grid::editGrob(data$grob[[1]], vp = vp)
    ggplot2:::ggname("geom_spatial", g)
  }, required_aes = c("grob", "x", "y"))
  layer(geom = GeomCustom, mapping = mapping, data = data, 
        stat = stat, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
                                                 ...))
}

# Negate %in% operator
'%!in%' <- function(x,y)!('%in%'(x,y))



  