# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Load Objects ------------------------------------------------------------
opc_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/opc_rd2_seurat_V3.RDS")

oligo_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/oligo_rd2_seurat_V3_B.RDS")

astro_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/astro_rd2_seurat_V3_B.RDS")

micro_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/micro_rd2_seurat_V3.RDS")

excit_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/excit_rd3_seurat_V3.RDS")

inhib_seurat_v3 <- read_rds("RDS_Subcluster_Seurat/V3_Objects/inhib_rd4_seurat_V3.RDS")


# Load Factor Lists -------------------------------------------------------
opc_factor_09_list <- read_rds("subcluster_factor_gene_list/liger_opc_rd2_gene_list.RDS")

oligo_factor_08_list <- read_rds("subcluster_factor_gene_list/liger_oligo_rd2_gene_list_b.RDS")

astro_factor_20_list <- read_rds("subcluster_factor_gene_list/liger_astro_rd2_gene_list_b.RDS")

micro_factor_06_list <- read_rds("subcluster_factor_gene_list/liger_micro_rd2_gene_list.RDS")

excit_factor_11_list <- read_rds("subcluster_factor_gene_list/liger_excit_rd3_gene_list_02.RDS")

inhib_factor_18_list <- read_rds("subcluster_factor_gene_list/liger_inhib_rd4_gene_list.RDS")


all_factor_list <- bind_cols(opc_factor_09_list, oligo_factor_08_list, astro_factor_20_list, micro_factor_06_list, excit_factor_05_list, excit_factor_11_list, excit_factor_12_list, inhib_factor_18_list)

# save all factor list
write_rds(all_factor_list, "subcluster_factor_gene_list/all_factor_list_b.RDS")
write.csv(all_factor_list, "subcluster_factor_gene_list/all_factor_list_b.csv")
all_factor_list <- read_rds("subcluster_factor_gene_list/all_factor_list_b.RDS")