# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts
Code to reproduce analysis objects for the data contained in [Tissue Processing and Enzymatic Dissociation Induce Artifactual Gene Expression of Brain Cell Populations in Mice and Humans](CITATION)

## Code
Included is the code necessary to replicate the Seurat or LIGER (or both) objects used for analysis and plotting.
- Each R file specifies version of Seurat/LIGER used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv](https://cran.r-project.org/web/packages/renv/index.html) packages followed by date-specific version installation using [versions](https://cran.r-project.org/web/packages/versions/index.html) package is recommended.

## Data
The data in this project can be broadly divided into 2 categories (5 sub-projects).  Please see [SI Table 1](LINK) for breakdown by sample and more information.

***For experiments 1-4***, all the raw data (fastqs) and expression matrices are available at the NCBI Gene Expression Omnibus (GEO) under [GSE152184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152184).   

***For experiment 5***, the raw data (fastqs) and expression matrices are available through the European phenome-Genome Archive (EGA) (Accession ID: [EGADXXXXXXX](EGADXXXXXXX))

**scRNAseq of cells isolated from mouse brain (all cell types and sorted myeloid cells)**

  1. scRNAseq of microglia with 4 different dissociation protocols | 10X 3' V2 | 12 samples (n=3 per group) | ([GSE152183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152183))
  2. scRNAseq of all CNS cells with or without inhibitors | 10X 3' V2 | 4 samples (n=2 per group) | ([GSE152182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152182))
  3. scRNAseq of microglia (tail vein PBS injection) | 10X 3' V2 | 4 samples | ([GSE152210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152210))

**snRNAseq of nuclei isolated from frozen human brain tissue**

  4. snRNAseq of post-mortem brain tissue | 10X 3' V3 | 3 samples | ([GSEXXXXXX](GSEXXXXXX))  
  5. snRNAseq of surgically resected brain tissue with or without freezing time delay | 10X 3' V3 | 4 samples (n= 2; x2 timepoints per sample) | ([GSEXXXXXX](GSEXXXXXX))

### Processed Data
_**Experiments 1-4**_  
All processed data files from Cell Ranger `count` outputs are available via NCBI GEO.  Information on Cell Ranger version and Genome/Annotation version can be found in each GEO record.
There are 3 processed data files per library:
  1. GSM\*\_Sample_Name_barcodes.tsv.gz: corresponds to the cell barcodes (i.e. column names).
  2. GSM\*\_Sample_Name_features.tsv.gz: corresponds to the gene identifiers (i.e. row names).
  3. GSM\*\_Sample_Name_matrix.mtx.gz: expression matrix in sparse format.

### Raw fastq Files
All raw data fastq files can be downloaded from SRA linked from NCBI GEO records or from EGA.

## ** Scripts to Figures Guides **
- Add when final version submitted
