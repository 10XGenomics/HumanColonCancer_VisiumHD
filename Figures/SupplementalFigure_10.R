
library(tidyverse)
library(Seurat)
library(glue)
library(arrow)
library(data.table)



sample_names <- c("Xenium_P1", "Xenium_P2", "Xenium_P5","Visium_HD_P1","Visium_HD_P2", "Visium_HD_P5")
sample_base_paths <- list(xenium_p1_base_path = "/path/to/Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE/",
                          xenium_p2_base_path = "/path/to/Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE/",
                          xenium_p5_base_path = "/path/to/Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE/",
                          hd_p1_base_path = "/path/to/Visium_HD_Human_Colon_Cancer_P1/",
                          hd_p2_base_path = "/path/to/Visium_HD_Human_Colon_Cancer_P2/",
                          hd_p5_base_path = "/path/to/Visium_HD_Human_Colon_Cancer_P5/")



# Visium Xenium Alignment

hd_paths <- sample_base_paths[str_detect(names(sample_base_paths), "hd")]
get_tissue_positions_8um <- function(base_path) {
  tp_path <- glue("{base_path}/binned_outputs/square_008um/spatial/tissue_positions.parquet")
  return(read_parquet(tp_path))
}

tp <- map(.x = hd_paths, .f = get_tissue_positions_8um)
names(tp) <- hd_paths
tp



get_microscope_image_path <- function(base_path) {
  image_path <- glue("{base_path}/tissue_image.btf")
  return(image_path)
}
#Note: must have image magick installed locally (https://imagemagick.org/)
image_paths <- map(hd_paths, get_microscope_image_path)
image_dims <- map(image_paths, ~ system(glue("magick identify {.x}"), intern = TRUE)[[1]])

# Get image widths and heights
widths <- numeric(length(image_dims))
heights <- numeric(length(image_dims))

# Loop through each item in the list
for (i in seq_along(image_dims)) {
  # Extract the string containing dimensions
  dimension_string <- regmatches(image_dims[[i]], regexpr("\\d+x\\d+", image_dims[[i]]))[1]
  
  # Split the dimension string into width and height
  dimensions <- strsplit(dimension_string, "x")[[1]]
  
  # Store width and height
  widths[i] <- as.numeric(dimensions[1])
  heights[i] <- as.numeric(dimensions[2])
}

# Create a data frame
image_dimensions <- data.frame(samples_hd = sample_names[str_detect(sample_names, "HD")], Width = widths, Height = heights)

# get HD image dimensions
tp_row_max_full_res <- list()
for (i in seq_along(tp)) {
  im_dims <- image_dimensions[i, ]
  tp_row_max_full_res[[i]] <- list("max" = tp[[i]] %>% 
                                     pull(pxl_row_in_fullres) %>% 
                                     max(),
                                   "min" = tp[[i]] %>% 
                                     pull(pxl_row_in_fullres) %>% 
                                     min()
  )
  if(tp_row_max_full_res[[i]][["min"]] < 0) {
    tp_row_max_full_res[[i]][["min"]] = 0
  }
  
  if(tp_row_max_full_res[[i]][["max"]] > im_dims$Height) {
    tp_row_max_full_res[[i]][["max"]] = im_dims
  }
}

names(tp_row_max_full_res) <- sample_names[str_detect(sample_names, "HD")]



tp_col_max_full_res <- list()
for (i in seq_along(tp)) {
  tp_col_max_full_res[[i]] <- list("max" = tp[[i]] %>% 
                                     pull(pxl_col_in_fullres) %>% 
                                     max(),
                                   "min" = tp[[i]] %>% 
                                     pull(pxl_col_in_fullres) %>% 
                                     min()
  )
  if(tp_col_max_full_res[[i]][["min"]] < 0) {
    tp_col_max_full_res[[i]][["min"]] = 0
  }
  
  if(tp_col_max_full_res[[i]][["max"]] > im_dims$Width) {
    tp_col_max_full_res[[i]][["max"]] = im_dims
  }
}

names(tp_col_max_full_res) <- sample_names[str_detect(sample_names, "HD")]


# Get the xenium and HD aligned image files
## This will allow us to only assess xenium cells that are within the visium hd capture area
## cell_info.csv available in the manuscript git repo
xenium_paths <- sample_base_paths[str_detect(sample_base_paths, "Xenium")]
get_xenium_cell_cords <- function(sample_name) {
  xenium_cell_info <- glue("/path/to/Figures/Xenium_Visium_Alignment/{sample_name}_cell_info.csv.gz")
  return(fread(xenium_cell_info))
}

xenium_cell_cords <- map(.x = sample_names,  ~ get_xenium_cell_cords(.x))
names(xenium_cell_cords) <- sample_names[str_detect(sample_names, "Xenium")]


# Select only cells and barcodes within the visium HD capture area
for (i in seq_along(xenium_cell_cords)) {
  xenium_cell_cords[[i]] <- xenium_cell_cords[[i]] %>% 
  filter(x_centroid_visium_scale <= tp_col_max_full_res[[i]]$max &
         x_centroid_visium_scale >= tp_col_max_full_res[[i]]$min) %>% 
  filter(y_centroid_visium_scale <= tp_row_max_full_res[[i]]$max &
        tp_row_max_full_res[[i]]$min >= 0)
}


# Get totalu UMI counts for each feature
get_per_feature_umis <- function(base_path){
  if(str_detect(base_path, "Xenium")){
    print("Reading Xenium Data")
    mat <-  Read10X_h5(glue("{base_path}/cell_feature_matrix.h5"))
    mat <- mat$`Gene Expression`
  } else {
    print("Reading HD Data")
    mat <- Read10X_h5(glue("{base_path}/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")) 
  }
  per_feature_umis <- rowSums(mat) %>% 
        as.data.frame() %>%
        rownames_to_column("feature_name") %>% 
        dplyr::rename("num_umis" = 2 )%>% 
        distinct(feature_name, .keep_all = TRUE)
    return(per_feature_umis)
}

pfumi <- map(sample_base_paths,  ~ get_per_feature_umis(.x))
names(pfumi) <- sample_names



pfumi_df <- pfumi %>% list_flatten() %>% bind_rows(.id = "sample")


# Make Plots
## P1
  plot_data_p1 <- pfumi_df %>% 
    filter(sample %in% c("Xenium_P1","Visium_HD_P1")) %>% 
    spread(sample, num_umis) %>% 
    filter(!if_any(everything(), is.na))
    
    
  
  plot_p1 <- plot_data_p1 %>% 
    ggplot(aes(x = log10(Xenium_P1 + 1), y = log10(Visium_HD_P1+1))) +
    geom_point(color = "#990F20") +
    geom_abline(color = "blue")+
    ggpubr::stat_cor(data = plot_data_p1, color = "black",
                               aes(x = log10(Xenium_P1 + 1),
                                   y = log10(Visium_HD_P1 + 1),
                                   label = after_stat(rr.label)),
                               inherit.aes = FALSE,
                               size = 5, method = "spearman")+
    xlab(glue("log10(UMIs +1): Xenium")) +
    ylab(glue("log10(UMIs +1): Visium_HD")) +
    ggtitle(glue::glue("P1 CRC"))+
    xlim(0,7) +
    ylim(0,7) +
    theme_classic() + 
    ggeasy::easy_all_text_size(size = 20)

## P2
  plot_data_p2 <- pfumi_df %>% 
    filter(sample %in% c("Xenium_P2","Visium_HD_P2")) %>% 
    spread(sample, num_umis) %>% 
    filter(!if_any(everything(), is.na))
    
    
  
  plot_p2 <- plot_data_p2 %>% 
    ggplot(aes(x = log10(Xenium_P2 + 1), y = log10(Visium_HD_P2+1))) +
    geom_point(color = "#990F20") +
    geom_abline(color = "blue")+
    ggpubr::stat_cor(data = plot_data_p2, color = "black",
                               aes(x = log10(Xenium_P2 + 1),
                                   y = log10(Visium_HD_P2 + 1),
                                   label = after_stat(rr.label)),
                               inherit.aes = FALSE,
                               size = 5, method = "spearman")+
    xlab(glue("log10(UMIs +1): Xenium")) +
    ylab(glue("log10(UMIs +1): Visium_HD")) +
    ggtitle(glue::glue("p2 CRC"))+
    xlim(0,7) +
    ylim(0,7) +
    theme_classic() + 
    ggeasy::easy_all_text_size(size = 20)


## P5
  plot_data_p5 <- pfumi_df %>% 
    filter(sample %in% c("Xenium_P5","Visium_HD_P5")) %>% 
    spread(sample, num_umis) %>% 
    filter(!if_any(everything(), is.na))
    
    
  
  plot_p5 <- plot_data_p5 %>% 
    ggplot(aes(x = log10(Xenium_P5 + 1), y = log10(Visium_HD_P5+1))) +
    geom_point(color = "#990F20") +
    geom_abline(color = "blue")+
    ggpubr::stat_cor(data = plot_data_p5, color = "black",
                               aes(x = log10(Xenium_P5 + 1),
                                   y = log10(Visium_HD_P5 + 1),
                                   label = after_stat(rr.label)),
                               inherit.aes = FALSE,
                               size = 5, method = "spearman")+
    xlab(glue("log10(UMIs +1): Xenium")) +
    ylab(glue("log10(UMIs +1): Visium_HD")) +
    ggtitle(glue::glue("P5 CRC"))+
    xlim(0,7) +
    ylim(0,7) +
    theme_classic() + 
    ggeasy::easy_all_text_size(size = 20)

