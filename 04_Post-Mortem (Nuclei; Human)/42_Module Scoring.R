# 1.0 Load packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# 2.0 Load V3 All Datasets ----------------------------------------------------
marsh_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/marsh_seuratv3_RENAMED.RDS")
zhou_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/zhou_seuratv3_RENAMED.RDS")
morabito_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/morabito_seuratv3_RENAMED.RDS")
leng_ec_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/leng_ec_seuratv3_RENAMED.RDS")
leng_sfg_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/leng_sfg_seuratv3_RENAMED_b.RDS")

beep(sound = 2)

# 3.0 Load Colors -------------------------------------------------------------
marsh_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid", "orange", "gold", "gray")
zhou_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
morabito_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
leng_ec_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
leng_sfg_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")

# 4.0 Load Gene Lists ---------------------------------------------------------
# 4.1 Load Mouse DEG Lists ------------------------------------------------
# Add Mouse list from SI Table XX
mouse_myeloid_list <- read.csv("Mouse_Data_Import/mouse_DEG_list_myeloid.csv", stringsAsFactors = FALSE, header = TRUE)

# Pull myeloid list
mouse_myeloid_list <- pull(mouse_myeloid_list, myeloid)

# Human homolog of Hist2h2aa1 is HIST2H2AA4
mouse_myeloid_list <- gsub("Hist2h2aa1", "HIST2H2AA4", mouse_myeloid_list)
print(mouse_myeloid_list)

# Convert to uppercase for human genes
mouse_myeloid_list_HUMAN <- mouse_myeloid_list %>% 
  str_to_upper()

# Pull all CNS list 
mouse_all_cns_list <- read.csv("Mouse_Data_Import/mouse_DEG_list_all_cns.csv", stringsAsFactors = FALSE, header = TRUE)

mouse_all_cns_list <- pull(mouse_all_cns_list, all_cns)
# Human homolog of 1500015O10Rik is C2orf40
mouse_all_cns_list <- gsub("1500015O10Rik", "C2orf40", mouse_all_cns_list)
print(mouse_all_cns_list)

mouse_all_cns_list_HUMAN <- mouse_all_cns_list %>% 
  str_to_upper() %>% 
  gsub("C2ORF40", "C2orf40", .)
print(mouse_all_cns_list_HUMAN)

# Create Combined Mouse Myeloid and All CNS list (No duplicates)
mouse_combined_list_HUMAN <- union(mouse_myeloid_list_HUMAN, mouse_all_cns_list_HUMAN)

# Remove mouse lists
rm(mouse_myeloid_list)
rm(mouse_all_cns_list)

# 4.2 Load LIGER Factor List --------------------------------------------------
# Load LIGER Factor Lists
all_factor_list <- read_rds("subcluster_factor_gene_list/all_factor_list.RDS")

all_factor_list <- data.frame(lapply(all_factor_list, as.character), stringsAsFactors=FALSE)

micro_liger_full <- all_factor_list$micro_rd2_factor6
micro_liger_top25 <- micro_liger_full[1:25]

astro_liger_full <- all_factor_list$astro_rd4_factor15
astro_liger_top25 <- astro_liger_full[1:25]

oligo_liger_full <- all_factor_list$oligo_rd2_factor18_1000
oligo_liger_top25 <- oligo_liger_full[1:25]

opc_liger_full <- all_factor_list$opc_rd2_factor9
opc_liger_top25 <- opc_liger_full[1:25]

excit_liger_full <- all_factor_list$neuron_rd3_factor11
excit_liger_top25 <- excit_liger_full[1:25]

inhib_liger_full <- all_factor_list$inhib_rd4_factor18
inhib_liger_top25 <- inhib_liger_full[1:25]

# 5.0 Add Module Scores ----------------------------------------------------
# 5.1 Marsh Dataset -------------------------------------------------------
# Mouse Lists
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# LIGER Lists
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(micro_liger_top25), name = "micro_liger_top25", search = TRUE)

marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(astro_liger_top25), name = "astro_liger_top25", search = TRUE)

marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(oligo_liger_top25), name = "oligo_liger_top25", search = TRUE)

marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(opc_liger_top25), name = "opc_liger_top25", search = TRUE)

marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(excit_liger_top25), name = "excit_liger_top25", search = TRUE)

marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(inhib_liger_top25), name = "inhib_liger_top25", search = TRUE)

beep(sound = 2)
View(marsh_seurat3@meta.data)

# 5.2 zhou Dataset -------------------------------------------------------
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# LIGER Lists
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(micro_liger_top25), name = "micro_liger_top25", search = TRUE)

zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(astro_liger_top25), name = "astro_liger_top25", search = TRUE)

zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(oligo_liger_top25), name = "oligo_liger_top25", search = TRUE)

zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(opc_liger_top25), name = "opc_liger_top25", search = TRUE)

zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(excit_liger_top25), name = "excit_liger_top25", search = TRUE)

zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(inhib_liger_top25), name = "inhib_liger_top25", search = TRUE)

# 5.3 Morabito Dataset -------------------------------------------------------
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# LIGER Lists
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(micro_liger_top25), name = "micro_liger_top25", search = TRUE)

morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(astro_liger_top25), name = "astro_liger_top25", search = TRUE)

morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(oligo_liger_top25), name = "oligo_liger_top25", search = TRUE)

morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(opc_liger_top25), name = "opc_liger_top25", search = TRUE)

morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(excit_liger_top25), name = "excit_liger_top25", search = TRUE)

morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(inhib_liger_top25), name = "inhib_liger_top25", search = TRUE)

# 5.4 leng ec Dataset -------------------------------------------------------
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# LIGER Lists
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(micro_liger_top25), name = "micro_liger_top25", search = TRUE)

leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(astro_liger_top25), name = "astro_liger_top25", search = TRUE)

leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(oligo_liger_top25), name = "oligo_liger_top25", search = TRUE)

leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(opc_liger_top25), name = "opc_liger_top25", search = TRUE)

leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(excit_liger_top25), name = "excit_liger_top25", search = TRUE)

leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(inhib_liger_top25), name = "inhib_liger_top25", search = TRUE)

# 5.5 leng sfg Dataset -------------------------------------------------------
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# LIGER Lists
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(micro_liger_top25), name = "micro_liger_top25", search = TRUE)

leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(astro_liger_top25), name = "astro_liger_top25", search = TRUE)

leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(oligo_liger_top25), name = "oligo_liger_top25", search = TRUE)

leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(opc_liger_top25), name = "opc_liger_top25", search = TRUE)

leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(excit_liger_top25), name = "excit_liger_top25", search = TRUE)

leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(inhib_liger_top25), name = "inhib_liger_top25", search = TRUE)

beep(sound = 2)

# 6.0 Create Score Lists ----------------------------------------------
all_module_score_names <- c("mouse_myeloid1", "mouse_all_cns1", "mouse_combined1", "micro_liger_top251", "astro_liger_top251", "oligo_liger_top251", "opc_liger_top251", "excit_liger_top251", "inhib_liger_top251")

all_module_score_gene_lists <- list("mouse_myeloid_list_HUMAN" = mouse_myeloid_list_HUMAN, "mouse_all_cns_list_HUMAN" = mouse_all_cns_list_HUMAN, "mouse_combined_list_HUMAN" = mouse_combined_list_HUMAN, "micro_liger_top25" = micro_liger_top25, "astro_liger_top25" = astro_liger_top25, "oligo_liger_top25" = oligo_liger_top25, "opc_liger_top25" = opc_liger_top25, "excit_liger_top25" = excit_liger_top25, "inhib_liger_top25" = inhib_liger_top25)

# 7.0 Save Objects ----------------------------------------------------
# Save Seurat Objects
write_rds(marsh_seurat3, "RDS_SeuratV3/rename_scored/marsh_seuratv3_RENAMED_SCORED.RDS")
write_rds(zhou_seurat3, "RDS_SeuratV3/rename_scored/zhou_seuratv3_RENAMED_SCORED.RDS")
write_rds(morabito_seurat3, "RDS_SeuratV3/rename_scored/morabito_seuratv3_RENAMED_SCORED.RDS")
write_rds(leng_ec_seurat3, "RDS_SeuratV3/rename_scored/leng_ec_seuratv3_RENAMED_SCORED.RDS")
write_rds(leng_sfg_seurat3, "RDS_SeuratV3/rename_scored/leng_sfg_seuratv3_RENAMED_SCORED.RDS")

# Save gene list and score names
write_rds(all_module_score_names, "subcluster_factor_gene_list/all_module_score_names.RDS")
write_rds(all_module_score_gene_lists, "subcluster_factor_gene_list/all_module_score_gene_lists.RDS")
beep(sound = 5)
