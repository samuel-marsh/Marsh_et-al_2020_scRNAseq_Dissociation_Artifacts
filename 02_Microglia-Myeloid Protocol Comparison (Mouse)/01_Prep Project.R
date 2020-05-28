# Start packrat
library(packrat)
init()

# Install versions
install.packages("versions")

# Load versions and install properly dataed packages
library(versions)
install.dates("tidyverse", "2019-02-28")

install.dates("BiocManager", "2019-02-28")

install.dates(c("lars", "fpc", "dtw", "Hmisc", "doSNOW", "foreach", "hdf5r"), "2019-02-28")

# Install Seurat from CRAN archive
    #wget https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.3.4.tar.gz
install.packages("Seurat_2.3.4.tar.gz", repos = NULL, type = "source")


