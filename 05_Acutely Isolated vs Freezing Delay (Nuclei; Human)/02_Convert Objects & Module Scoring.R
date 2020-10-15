# 1.0 Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) #v2.3.4
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Convert to Seurat Objects -------------------------------------------------------
# read liger objects
mg_fresh <- read_rds("~/Desktop/exp17_Dissociation Paper/Tushar/Fresh Tissue/smarsh_mgs/MG/a.mgcleaned.rds")
astro_fresh <- read_rds("~/Desktop/exp17_Dissociation Paper/Tushar/Fresh Tissue/smarsh_mgs/astros/a.astrocleaned.rds")

# Microglia convert and metadata
mg_status <- mg_fresh@cell.data$status
mg_seurat <- ligerToSeurat(object = mg_fresh)
mg_seurat@meta.data$status <- mg_status
write_rds(mg_seurat, "RDS_Seurat Objects/mg_seurat.RDS")

# Astro convert and metadata
astro_status <- astro_fresh@cell.data$status
astro_seurat <- ligerToSeurat(object = astro_fresh)
astro_seurat@meta.data$status <- astro_status
write_rds(astro_seurat, "RDS_Seurat Objects/astro_seurat.RDS")

beep(sound = 2)

# Restart R and install Seurat V3.1.5 -------------------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Update to V3
mg_seurat <- read_rds("RDS_Seurat Objects/mg_seurat.RDS")
astro_seurat <- read_rds("RDS_Seurat Objects/astro_seurat.RDS")

mg_seurat <- UpdateSeuratObject(mg_seurat)
astro_seurat <- UpdateSeuratObject(astro_seurat)

write_rds(mg_seurat, "RDS_Seurat Objects/mg_seurat3.RDS")
write_rds(astro_seurat, "RDS_Seurat Objects/astro_seurat3.RDS")

beep(sound = 2)

# Load Post-Mortem & Mouse Module Score Lists -----------------------------
all_module_score_names <- read_rds("~/Desktop/exp17_Dissociation Paper/02_Online_Liger_PM_Local/subcluster_factor_gene_list/all_module_score_names.RDS")
all_module_score_lists <- read_rds("~/Desktop/exp17_Dissociation Paper/02_Online_Liger_PM_Local/subcluster_factor_gene_list/all_module_score_gene_lists.RDS")

# Load Fresh Objects ------------------------------------------------------
mg_fresh <- read_rds("~/Desktop/exp17_Dissociation Paper/Tushar/Fresh Tissue/smarsh_mgs/MG/mg_annot.rds")
astro_fresh <- read_rds("~/Desktop/exp17_Dissociation Paper/Tushar/Fresh Tissue/smarsh_mgs/astros/astro_annot.rds")

# Load Seurat Versions
mg_seurat3 <- read_rds("RDS_Seurat Objects/mg_seurat3.RDS")
astro_seurat3 <- read_rds("RDS_Seurat Objects/astro_seurat3.RDS")

mg_stats <- Cluster_Stats_All_Samples(mg_seurat3)
astro_stats <- Cluster_Stats_All_Samples(astro_seurat3)

# Extract Lists -----------------------------------------------------------
mouse_myeloid_list <- all_module_score_lists[["mouse_myeloid_list_HUMAN"]]
mouse_all_cns_list <- all_module_score_lists[["mouse_all_cns_list_HUMAN"]]
mouse_comb_list <- all_module_score_lists[["mouse_combined_list_HUMAN"]]

liger_mg_list <- all_module_score_lists[["micro_liger_top25"]]
liger_astro_list <- all_module_score_lists[["astro_liger_top25"]]

# Module Scores Microglia -----------------------------------------------------------
# Create Module Score
mg_seurat3 <- AddModuleScore(object = mg_seurat3, features = list(liger_mg_list), name = "liger_mg_factor")
Idents(mg_seurat3) <- "status"
hr0 <- unlist(CellsByIdentities(mg_seurat3, idents = "0"))
hr6 <- unlist(CellsByIdentities(mg_seurat3, idents = "6"))

# Module Astros -----------------------------------------------------------
astro_seurat3 <- AddModuleScore(object = astro_seurat3, features = list(liger_astro_list), name = "liger_astro_factor")
Idents(astro_seurat3) <- "status"
hr0 <- unlist(CellsByIdentities(astro_seurat3, idents = "0"))
hr6 <- unlist(CellsByIdentities(astro_seurat3, idents = "6"))