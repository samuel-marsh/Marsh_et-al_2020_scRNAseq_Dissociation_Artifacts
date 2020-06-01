# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts
Code to reproduce analysis objects for the data contained in (CITATION AND LINK HERE)

Included is the code necessary to replicate the Seurat object used for analysis and plotting.
- Each R file specifies version of Seurat used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv] packages followed by date-specific version installation using [versions](https://cran.r-project.org/web/packages/versions/index.html) package is recommended.
