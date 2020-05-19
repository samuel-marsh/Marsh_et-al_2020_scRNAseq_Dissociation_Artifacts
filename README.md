# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts
Code to reproduce analysis objects for the data contained in (CITATION AND LINK HERE)

Included is the code necessary to replicate the Seurat object used for analysis and plotting.
- Each R file specifies version of Seurat used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.
- Each R File also specifies which Figures and SI Tables analysis corresponds to.

- Where possible date of analysis performed prior to version control is specified.  To replicate analyses performed on specific date use of the [versions](https://cran.r-project.org/web/packages/versions/index.html) package is recommended.
