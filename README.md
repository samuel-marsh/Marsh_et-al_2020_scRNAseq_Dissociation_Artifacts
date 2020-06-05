# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts
Code to reproduce analysis objects for the data contained in (CITATION AND LINK HERE)

Included is the code necessary to replicate the Seurat object used for analysis and plotting.
- Each R file specifies version of Seurat used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv](https://cran.r-project.org/web/packages/renv/index.html) packages followed by date-specific version installation using [versions](https://cran.r-project.org/web/packages/versions/index.html) package is recommended.

## Data
All the raw data (fastqs) and expression matrices are available at the Gene Expression Omnibus (GEO) under [GSEXXXXXXX](GSEXXXXXXX) (Coming soon). The data in this project can be broadly divided into 5 subprojects:
- Single Cell Sequencing of Microglia with 4 different dissociation protocols; 12 samples (n=3 per group) ([GSEXXXXXXX](GSEXXXXXXX))

### Processed Data
All processed data files from Cell Ranger output are available via NCBI GEO.  There are 3 files per library:
  1. Sample_Name_barcodes.tsv.gz: corresponds to the cell barcodes (column names).
  2. Sample_Name_features.tsv.gz: corresponds to the gene identifiers.
  3. Sample_Name_matrix*mtx.gz: expression matrix in sparse format.

### fastq Files
All raw data fastq files can be downloaded from SRA linked from NCBI GEO records.

## ** Scripts to Figures Guides **
- Figures 1, SI Figures
