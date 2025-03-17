library(tidyverse)
library(data.table)
library(patchwork)
library(Seurat)
library(rhdf5)
library(glue)
library(Matrix)
library(ggpubr)
library(ggeasy)

# Read in the data

## Set samples names
sample_names <- c("Visium_HD_p2", "Visium_v2_p2","Visium_HD_p5", "Visium_v2_p5")

sample_base_paths <- list(hd_base_path_p2 = "/path/to/Visium_HD_Human_Colon_Cancer_P2/binned_outputs/square_002um/",
                          v2_base_path_p2 = "/path/to/Visium_V2_Human_Colon_Cancer_P2/",
                          hd_base_path_p5 = "/path/to/Visium_HD_Human_Colon_Normal_P5/binned_outputs/square_002um/",
                          v2_base_path_p5 = "/path/to/Visium_V2_Human_Colon_Normal_P5/")

# Identify the paths for each raw_probe_bc_matrix.h5
per_probe_paths <- list()
for (i in seq_along(sample_base_paths)) {
  per_probe_paths[[i]] <- paste0(sample_base_paths[[i]], "raw_probe_bc_matrix.h5")
}
names(per_probe_paths) <- sample_names

# Funtion to read any given 10x .h5 matrix
read_matrix <- function (base_path, type = c("filtered", "raw", "per_probe")) {
  if (type == "filtered" | type == "raw") {
    file_name <- glue(base_path,
                      "{type}_feature_bc_matrix.h5")
    t(Read10X_h5(file_name))
  } else if (type == "per_probe") {
    file_name <- glue(base_path,
                      "raw_probe_bc_matrix.h5")
    t(Read10X_h5(file_name))
  } else {
    stop("Not a vaild matrix type")
  }
}

# Read raw_probe_bc_matrix.h5 for each sample
per_probe_matrix <- map(.x = sample_base_paths, ~read_matrix(.x, "per_probe"))

# Get the probes filtered by the gDNA algorithm
filter_probes <- map(per_probe_paths, ~ h5read(.x, "matrix/features/filtered_probes"))

per_probe_matrix_filtered <- list()
for (i in seq_along(per_probe_matrix)) {
  per_probe_matrix_filtered[[i]] <- per_probe_matrix[[i]][,as.logical(filter_probes[[i]])]
  
}

# Get per probe UMI counts
per_probe_sums <- function(per_probe_matrix_filtered) {
  counts <- colSums(per_probe_matrix_filtered)
  probe_names <- colnames(per_probe_matrix_filtered)
  per_probe_df <- data.frame(probe = probe_names, probe_sum = counts)
}

# all probes
per_probe_counts <- map(.x = per_probe_matrix_filtered, ~per_probe_sums(.x)) %>% 
  set_names(sample_names) %>%
  bind_rows(.id = "sample_name")

# only probes that pass the gDNA filter
per_probe_counts_unfiltered <- map(.x = per_probe_matrix, ~per_probe_sums(.x)) %>% 
  set_names(sample_names) %>%
  bind_rows(.id = "sample_name")

# Only spliced probes
## Only need to pull in one file since they all use the same one
target_panel_path <- "https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
target_panel <- fread(target_panel_path)

spliced_probes <- target_panel %>% 
  filter(region == "spliced") %>% 
  mutate(probe_split = str_match(probe_id, "\\|(.*)")[, 2]) %>% 
  pull(probe_split)
spliced_probes

# Supplemental Figure 1A
# Plot the data
plot_data_1a <-  per_probe_counts_unfiltered %>%
  spread(sample_name, probe_sum)
plot_unfiltered_a <- plot_data_1a %>% 
  ggplot(aes(x = log10(Visium_v2_p2 + 1), y = log10(Visium_HD_p2+1))) +
  geom_point(color = "#990F20") +
  geom_smooth(method='lm', se = FALSE) +
  stat_cor(data = plot_data_1a, color = "black",
           aes(x = log10(Visium_v2_p2 + 1),
               y = log10(Visium_HD_p2 + 1),
               label = after_stat(rr.label)),
           inherit.aes = FALSE,
           size = 5)+
  xlab("log10(UMIs per probe +1) Visium v2") +
  ylab("log10(UMIs per probe +1) Visium HD") +
  ggtitle(glue("P2 NAT Unfiltered"))+
  xlim(0,7) +
  ylim(0,7) +
  theme_classic()+
  easy_all_text_size(size = 20)


plot_data_2a <-  per_probe_counts_unfiltered %>%
  spread(sample_name, probe_sum) %>% 
  filter(probe %in% spliced_probes)
plot_spliced_a <- plot_data_2a %>% 
  ggplot(aes(x = log10(Visium_v2_p2 + 1), y = log10(Visium_HD_p2+1))) +
  geom_point(color = "#990F20") +
  geom_smooth(method='lm', se = FALSE) +
  stat_cor(data = plot_data_2a, color = "black",
           aes(x = log10(Visium_v2_p2 + 1),
               y = log10(Visium_HD_p2 + 1),
               label = after_stat(rr.label)),
           inherit.aes = FALSE,
           size = 5)+
  xlab("log10(UMIs per probe +1) Visium v2") +
  ylab("log10(UMIs per probe +1) Visium HD") +
  ggtitle(glue("P2 NAT Spliced Probes"))+
  xlim(0,7) +
  ylim(0,7) +
  theme_classic() +
  easy_all_text_size(size = 20)

highres_sup_fig1a <- plot_unfiltered_a +  plot_spliced_a

# Supplemental Figure 1B
# Plot the data
plot_data_1b <-  per_probe_counts_unfiltered %>%
  spread(sample_name, probe_sum)
plot_unfiltered_b <- plot_data_1b %>% 
  ggplot(aes(x = log10(Visium_v2_p5 + 1), y = log10(Visium_HD_p5+1))) +
  geom_point(color = "#990F20") +
  geom_smooth(method='lm', se = FALSE) +
  stat_cor(data = plot_data_1b, color = "black",
           aes(x = log10(Visium_v2_p5 + 1),
               y = log10(Visium_HD_p5 + 1),
               label = after_stat(rr.label)),
           inherit.aes = FALSE,
           size = 5)+
  xlab("log10(UMIs per probe +1) Visium v2") +
  ylab("log10(UMIs per probe +1) Visium HD") +
  ggtitle(glue("P5 NAT Unfiltered"))+
  xlim(0,7) +
  ylim(0,7) +
  theme_classic()+
  easy_all_text_size(size = 20)


plot_data_2b <-  per_probe_counts_unfiltered %>%
  spread(sample_name, probe_sum) %>% 
  filter(probe %in% spliced_probes)
plot_spliced_b <- plot_data_2b %>% 
  ggplot(aes(x = log10(Visium_v2_p5 + 1), y = log10(Visium_HD_p5+1))) +
  geom_point(color = "#990F20") +
  geom_smooth(method='lm', se = FALSE) +
  stat_cor(data = plot_data_2b, color = "black",
           aes(x = log10(Visium_v2_p5 + 1),
               y = log10(Visium_HD_p5 + 1),
               label = after_stat(rr.label)),
           inherit.aes = FALSE,
           size = 5)+
  xlab("log10(UMIs per probe +1) Visium v2") +
  ylab("log10(UMIs per probe +1) Visium HD") +
  ggtitle(glue("P5 NAT Spliced Probes"))+
  xlim(0,7) +
  ylim(0,7) +
  theme_classic() +
  easy_all_text_size(size = 20)

highres_sup_fig1b <- plot_unfiltered_b +  plot_spliced_b

