# ![HD Image](ext/Figure_CRC.png)

# Characterization of immune cell populations in the tumor microenvironment of colorectal cancer using high definition spatial profiling

Michelli F. Oliveira<sup>†</sup>, Juan P. Romero<sup>†</sup>, Meii Chung<sup>†</sup>, Stephen Williams, Andrew D. Gottscho, Anushka Gupta, Sue Pilipauskas, Syrus Mohabbat, Nandhini Raman, David Sukovich, David Patterson, [Visium HD Development Team](Team), Sarah E. B. Taylor<sup>‡</sup>

<sup>†</sup> These authors contributed equally to this work

<sup>‡</sup> Corresponding author

## Abstract

Colorectal cancer (CRC) is the second-deadliest cancer in the world, yet a deeper understanding of spatial patterns of gene expression in the tumor microenvironment (TME) remains elusive. Here, we introduce the Visium HD platform (10x Genomics) and use it to investigate human CRC and normal adjacent mucosal tissues from formalin fixed paraffin embedded (FFPE) samples. The first assay available on Visium HD is a probe-based spatial transcriptomics workflow that was developed to enable whole transcriptome single cell scale analysis. We demonstrate highly refined unsupervised spatial clustering in Visium HD data that aligns with the hallmarks of colon tissue morphology and is notably improved over earlier Visium assays. Using serial sections from the same FFPE blocks we generate a single cell atlas of our samples, then we integrate the data to comprehensively characterize the immune cell types present in the TME, specifically at the tumor periphery. We observed enrichment of two pro-tumor macrophage subpopulations with differential gene expression profiles that were localized within distinct tumor regions. Further characterization of the T cells present in one of the samples revealed a clonal expansion that we were able to localize in the tissue using in situ gene expression analysis. In situ analysis also allowed us to perform in-depth characterization of the microenvironment of the clonally expanded T cell population and we identified a third macrophage subpopulation with gene expression profiles consistent with an anti-tumor response. Our study provides a comprehensive map of the cellular composition of the CRC TME and identifies phenotypically and spatially distinct immune cell populations within it. We show that the single cell-scale resolution afforded by Visium HD and the whole transcriptome nature of the assay allows investigations into cellular function and interaction at the tumor periphery in FFPE tissues, which has not been previously possible.

## Data

The full dataset used in this repository and in the manuscript can be downloaded from the following link [Dataset](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression/dataset-human-crc)

## Repository

This repository contains the scripts to replicate the findings displayed in the manuscript. It is organized into two folders *Figures* and *Methods*. The *Figures* folder has the scripts to replicate the figures in the manuscript and the files are named accordingly. The *methods* folder contains the different custom methods developed for the manuscript.

The *Figures* files require specific outputs generated with the *Methods* scripts.

## Methods

In this section with provide a description and a start guide for the different methods used and developed for the manuscript.

### AuxFunctions.R

R script with multiple custom R functions used in the manuscript. To load all the functions, we use the source function:

```R
  source("~/HumanColonCancer_VisiumHD/Methods/AuxFunctions.R")
```

### FlexSingleCell.R
R script used to process the FLEX single cell data. It takes the outputs from [*cellranger aggr*](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-aggr).

Given the dataset's large size, we adopted the sketch-based analysis approach in Seurat<sup>1</sup> v5 [sketched-based analysis](https://satijalab.org/seurat/articles/seurat5_sketch_analysis), sampling 15% of the entire dataset (~37,000 cells) for downstream analysis. After completing the analysis on the subsampled data, we extended it to the entire single cell dataset.

The script saves the full processed Seurat object and the Metadata for plotting purposes in the Figure scripts.

```R
saveRDS(ColonCancer_Flex,file='~/Outputs/Flex/FlexSeuratV5.rds') # Full Seurat Object
saveRDS(ColonCancer_Flex@meta.data,file='~/Outputs/Flex/FlexSeuratV5_MetaData.rds') #Meta Data
```

### Deconvolution.R

R script used to run [spaceXR](https://github.com/dmcable/spacexr)<sup>2</sup> for deconvolution. It requires the UMI count matrix from [*cellranger aggr*](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-aggr) and the **MetaData** generated with the `FlexSingleCell.R`  script to generate the reference. For Visium HD the [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest) outs are also required.

Due to the number of barcodes in Visium HD, we modified the source code of [spaceXR](https://github.com/dmcable/spacexr) to improve runtime. The modified version can be found in the following [Pull Request](https://github.com/dmcable/spacexr/pull/206). However, the original version can also be used to deconvolve the Visium HD data.

In the script we use sample P1CRC as a template to run the algorithm, but it can also be used for any other sample.

### NucleiSegmentation.py

Python script used to run nuclei segmentation on H&E images used for the tissue sections processed with Visium HD.
To create the the conda environemnt please see [yml section](#ymlfile)

The script takes an HE image as an imput and performs nuclei segmentation on the full section using the [stardist](https://github.com/stardist/stardist)<sup>3</sup> package. The user can provide a set of coordinates to generate a crop of the image along with the corresponding masks located within that region. The user also provides the path to the outputs directory for a given bin size (i.e. 2µm) and will output a .csv file that assigns all the barcodes located within the all segmentation masks.

The script can be called as follows:

```Python
python ./HumanColonCancer_VisiumHD/Methods/NucleiSegmentation.py -i ./PATH_TO_HE_image -r1 rowmin -r2 rowmax -c1 colmin -c2 colmax -s ./PATH_TO_SR_outs/binned_outputs/square_002um/ -o Output_directory
```

More details on the required inputs:
```Python
    parser.add_argument('-i','--image', type=str, help='Path to HE image')
    parser.add_argument('-r1','--rmin', type=int, help='row min for zoom in')
    parser.add_argument('-r2','--rmax', type=int, help='row max for zoom in')
    parser.add_argument('-c1','--cmin', type=int, help='column min for zoom in')
    parser.add_argument('-c2','--cmax', type=int, help='column max for zoom in')
    parser.add_argument('-s','--srdir', type=str, help='Path to spaceranger outs at a given bin size')
    parser.add_argument('-o',"--out", type=str,help="Directory where to save outputs")
```

Outputs:

- `Nuclei_Barcode_Map.csv` csv file with the Nuclei and barcode relationship for the full section
- `labels_FullSection.pckl` Labels of the identified nuclei for the full section
- `polys_FullSection.pckl` Coordinates of the identified polygons for the full section.
- `img_rois_Stardist_Subset.zip` if coordinates are given, segmented nuclei for the selected region. Can be visualized with [QuPath](https://qupath.github.io/)<sup>4</sup>
- `img_Stardist.tif` if coordinates are given, tif file with the selected zoom in region. Can be visualized with [QuPath](https://qupath.github.io/)<sup>4</sup>


### environment_nucleisegmentation.yml <a name="ymlfile"></a>
 `yml` file to create a conda environment with all the required dependencies for the `NucleiSegmentation.py` script. 
 To create the evironment using the provided file:

```Python
conda env create --name NucleiSeg --file=./HumanColonCancer_VisiumHD/Methods/environment_nucleisegmentation.yml
```

To activate the environment:

```Python
conda activate NucleiSeg 
```
## MetaData

The MetaData folder contains files with the associated metadata used in the manuscript.

#### Single Cell
The `SingleCell_MetaData.csv.gz` contains the following columns:
1. **Barcode** : cell barcode
2. **Patient** : Patient of origin
3. **BC** : Probe barcode to identify sample of origin
4. **QCFilter** : Binary column denoting if a cell was kept or removed during QC
5. **Level1** : Level 1 cell type annotation
6. **Level2** : Level 2 cell type annotation
7. **UMAP1** : UMAP dimension 1 coordinates
8. **UMAP2** : UMAP dimension 2 coordinates

#### Visium HD
The parquet files (i.e `P1CRC_Metadata.parquet`) can be opened in R using the following code:

```R
  library(arrow)
  Data<-read_parquet("~/HumanColonCancer_VisiumHD/MetaData/P1CRC_Metadata.parquet")
```

These parquet files contain the following columns:
1. **barcode** : 8um bin barcode
2. **tissue** : Binary column denoting if the bin is under tissue or not
3. **X** : Spatial X coordinate
4. **Y** : Spatial Y coordinate
5. **DeconvolutionClass** : Deconvolution class for the bin (singlet, doublet, doublet_certain,doublet_uncertain or reject)
6. **DeconvolutionLabel1** : Gives the first cell type predicted on the bin
7. **DeconvolutionLabel2** : Gives the second cell type predicted on the bin (Not valid for reject or doublet_uncertain)
8. **Periphery** : Indicates if the bin is in the 50 micron tumor periphery, in the tumor or rest of the tissue
9. **UnsupervisedL1** : Merged unsupervised clusering annotation (Level 1)
10. **UnsupervisedL2** : Merged unsupervised clusering annotation (Level 2)
11. **MacrophageSubtype** : Subtype of macrophage (SELENOP+ or SPP1+) in the tumor periphery
12. **GobletSubcluster** : Goblet subcluster used in Figure 5

## Figures

The Figures folder contains all the scripts to create the figures used in the manuscript.
Most of the scripts within this folder require outputs generated from the *Methods* section.

The required R packages are common across the files:

```R
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
library(rhdf5)
library(glue)
library(Matrix)
library(ggpubr)
library(ggeasy)
library(arrow)
```

The beginning of each file starts with a data.frame that can be used as a template to generate the output for the different sections. We use P1CRC as a template, but can be replaced with any other section.

```R
SampleData<-data.frame(Patient = "PatientCRC1", # Name of the Sample
                       PathSR="~/VisiumHD/PatientCRC1/outs/", # Path to space ranger outs folder
                       PathDeconvolution="~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds") #spaceXR Deconvolution results
```

## References

1. Hao, Yuhan, et al. "Dictionary learning for integrative, multimodal and scalable single-cell analysis." Nature biotechnology 42.2 (2024): 293-304. [Paper](https://www.nature.com/articles/s41587-023-01767-y)

2. Cable, Dylan M., et al. "Robust decomposition of cell type mixtures in spatial transcriptomics." Nature biotechnology 40.4 (2022): 517-526. [Paper](https://www.nature.com/articles/s41587-021-00830-w)

3. Weigert, Martin, and Uwe Schmidt. "Nuclei instance segmentation and classification in histopathology images with stardist." 2022 IEEE International Symposium on Biomedical Imaging Challenges (ISBIC). IEEE, 2022. [Paper](https://arxiv.org/abs/2203.02284)

4. Bankhead, Peter, et al. "QuPath: Open source software for digital pathology image analysis." Scientific reports 7.1 (2017): 1-7. [Paper](https://www.nature.com/articles/s41598-017-17204-5)


## bioRxiv version
Here we list the git tag and link to the preprint originall submitted to bioRxiv
[Link to Preprint](https://www.biorxiv.org/content/10.1101/2024.06.04.597233v1)
[Tag]()