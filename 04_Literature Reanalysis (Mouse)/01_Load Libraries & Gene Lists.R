# All literature reanalysis was performed using:
    # R 3.6.1
    # Seurat 3.1.5
    # tidyverse 1.3.0
    # loomR 0.2.1.9000

# 1.0 Load Packages & Scripts ---------------------------------------------
library(tidyverse)
library(readxl)
library(Seurat)
library(viridis)
library(beepr)
library(scCustomize)
library(patchwork)
library(loomR)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 5000 * 1024^2)

# 2.0 Load Gene Lists for Module Scoring ----------------------------------
    # Gene lists available in SI Tables XX-XX
shared_sig <- "Load Microglia Meta Cell Score"

shared_sig_ensembl <- "Convert gene names to ensembl"

homeostatic_mg <- "Homeostatic microglia gene list"

homeostatic_mg_ensembl <- "Convert gene names to ensembl"