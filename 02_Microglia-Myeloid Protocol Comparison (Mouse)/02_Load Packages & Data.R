library(tidyverse)
library(Seurat)
library(ggrepel)
library(viridis)
library(beepr)


# Read in data
raw_data <- Read10X("~/Dropbox (BCH)/untitled folder/RAW_DATA_EXP_17/Micro_for_Seurat_2/")
    # Modified Cell Ranger 3 output to work with Seurat V2 (change features.tsv to genes.tsv)

# Create Seurat object 
exp17_micro_ALL <- CreateSeuratObject(raw.data = raw_data, min.cells = 10, min.genes = 200, project = "Exp17_MICRO", names.delim = "-", names.field = 2); beep(sound = 2)

# Modify metadata 
meta_data <- data.frame(exp17_micro_ALL@meta.data)
meta_data_tb <- meta_data %>% 
  rownames_to_column(var = "cell_name") %>% 
  as_tibble()

# Create string of which samples 
INHIB <- c("3", "4", "7", "8", "11", "12")
DOUNCE <- c("1", "3", "5", "7", "9", "11")
BATCH <- c("1", "2", "3", "4", "5", "6")

# Create new column based on string
meta_data_tb <- mutate(meta_data_tb, Transcription = ifelse(orig.ident %in% INHIB, "INHIB", "NO_INHIB"))
meta_data_tb <- mutate(meta_data_tb, Method = ifelse(orig.ident %in% DOUNCE, "DOUNCE", "ENZYMATIC"))
meta_data_tb <- mutate(meta_data_tb, Batch = ifelse(orig.ident %in% BATCH, "BATCH_01", "BATCH_02"))

# Create additional column based on combined method and Transcription
meta_data_tb$Transcription_Method[meta_data_tb$orig.ident == "3" | meta_data_tb$orig.ident == "7"| meta_data_tb$orig.ident == "11"] <- "DOUNCE_INHIB"
meta_data_tb$Transcription_Method[meta_data_tb$orig.ident == "1" | meta_data_tb$orig.ident == "5"| meta_data_tb$orig.ident == "9"] <- "DOUNCE_NONE"
meta_data_tb$Transcription_Method[meta_data_tb$orig.ident == "4" | meta_data_tb$orig.ident == "8"| meta_data_tb$orig.ident == "12"] <- "ENZYMATIC_INHIB"
meta_data_tb$Transcription_Method[meta_data_tb$orig.ident == "2" | meta_data_tb$orig.ident == "6"| meta_data_tb$orig.ident == "10"] <- "ENZYMATIC_NONE"

# Convert back to data frame
meta_data_mod <- meta_data_tb %>% 
  select("cell_name", "Transcription", "Method", "Transcription_Method", "Batch") %>% 
  data.frame() %>% 
  column_to_rownames(var = "cell_name")

# Add new metadata to seurat object
exp17_micro_ALL <- AddMetaData(exp17_micro_ALL,
                               metadata = meta_data_mod,
                               col.name = Transcription)

View(exp17_micro_ALL@meta.data)

# Add mito percentage
mito_genes <- grep(pattern = "^mt-", x = rownames(x = raw_data), value = TRUE)
percent_mito <- Matrix::colSums(exp17_micro_ALL@raw.data[mito_genes, ]) / Matrix::colSums(exp17_micro_ALL@raw.data)

exp17_micro_ALL <- AddMetaData(object = exp17_micro_ALL, metadata = percent_mito, col.name = "percent_mito")

# Reorder sample identies to proper numberical order
# Define order of samples
sample_order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Relevel object@orig.ident
exp17_micro_ALL@meta.data$orig.ident <- factor(x = exp17_micro_ALL@meta.data$orig.ident, levels = sample_order, ordered = TRUE)
exp17_micro_ALL@ident <- factor(x = exp17_micro_ALL@ident, levels = sample_order, ordered = TRUE)
