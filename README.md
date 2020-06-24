# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts  
#### Code to reproduce analysis objects for the data contained in:  
[Tissue Processing and Enzymatic Dissociation Induce Artifactual Gene Expression of Brain Cell Populations in Mice and Humans](CITATION)  
Samuel E. Marsh, Tushar Kamath, Timothy R. Hammond, Alec J. Walker, Lasse Dissing-Olesen, Adam M.H. Young, Abdul Abdul, Naeem Nadaf, Connor Dufort, Sarah Murphy, Velina Kozareva, Charles Vanderburg, Soyon Hong, Harry Bulstrode, Peter J. Hutchinson, Robin J.M. Franklin, Evan Macosko, & Beth Stevens

## Code
Included is the code necessary to replicate the Seurat or LIGER (or both) objects used for analysis and plotting.
- Each R file specifies version of Seurat/LIGER used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.
    - Seurat V2.3.4 source package can be downloaded here from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source.

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date the following actions are recommended or described in code:
  - Use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv](https://cran.r-project.org/web/packages/renv/index.html) packages. Followed by date-specific version installation of CRAN packages using [versions](https://cran.r-project.org/web/packages/versions/index.html) package.
  - Archived source versions of specific packages may also be needed depending on version of R and can be downloaded from CRAN archives and installed from local source.

- LIGER analyses were performed using the in development ["online"](https://github.com/MacoskoLab/liger/tree/online) branch, updating throughout analysis to accommodate bug fixes.  As such analysis code may not fully reproduce identical figures as presented in the paper.  Instances such as this are denoted in code and .RDS objects will be linked in code document.
  - LIGER analyses also utilize multiple versions of Seurat as specified in code for some of the following situations:
    - Seurat V3 used used for data import, QC filtering (genes, UMIs, % mito), and final plotting due to more advanced plotting features and patchwork compatibility
    - Seurat V2 was used during LIGER workflow to accommodate use of `clusterLouvainJaccard` function which relied on Seurat V2 object structure
    - Conversion between Seurat and LIGER objects was performed using built in LIGER functions `seuratToLiger` and `ligerToSeurat`

## Data  
### Original Data
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
**Experiments 1-4**  
All processed data files from Cell Ranger `count` outputs are available via NCBI GEO.  Information on Cell Ranger version and Genome/Annotation version can be found in each GEO record.
There are 3 processed data files per library:
  1. GSM\*\_Sample_Name_barcodes.tsv.gz: corresponds to the cell barcodes (i.e. column names).
  2. GSM\*\_Sample_Name_features.tsv.gz: corresponds to the gene identifiers (i.e. row names).
  3. GSM\*\_Sample_Name_matrix.mtx.gz: expression matrix in sparse format.

### Raw fastq Files
All raw data fastq files can be downloaded from SRA linked from NCBI GEO records or from EGA.

### Literature Reanalysis
Reanalyzed data from literature is summarized detailed in table below.
| Dataset | Species | Seq Used | Raw/Processed Data | Publication |
| :-----: | :-----: | :------: | :----------------: | :---------: |
| Mathys | Mouse | scRNAseq (Smart-seq2) | [GEO103334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103334) & Authors | [Mathys et al., 2017 (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31314-1?) |
| Plemel | Mouse | scRNAseq (10X 3' V2) | [GSE115803](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115803) | [Plemel et al., 2020 (Science Advances)](https://advances.sciencemag.org/content/6/3/eaay6324) |
| Zywitza | Mouse | scRNAseq (Drop-Seq) | [GSE111527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527) | [Zywitza et al., 2018 (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31732-7?) |
| Mizrak | Mouse | scRNAseq (Microwell Seq) | [GSE109447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109447) | [Mizrak et al., 2019 (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31974-0?) |
| Zeisel | Mouse | scRNAseq (10X 3' V1 & V2) | [mousebrain.org](mousebrain.org) | [Zeisel et al., 2018 (Cell)](https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X?) |
| Hammond | Mouse | scRNAseq (10X 3' V1 & V2) | [GSE121654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654) | [Hammond et al., 2019 (Immunity)](https://www.cell.com/immunity/fulltext/S1074-7613(18)30485-0?) |
| Zhou | Human | snRNAseq (10X 5' V1) | [syn21670836](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn21670836) | [Zhou et al., 2020 (Nature Medicine)](https://www.nature.com/articles/s41591-019-0695-9?) |
| Mathys (Hu) | Human | snRNAseq (10X 3' V2) | [syn18485175](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn18485175) | [Mathys et al., 2019 (Nature)](https://www.nature.com/articles/s41586-019-1195-2) |
| Morabito | Human | snRNAseq (10X 3' V3) | [syn18915937](https://www.synapse.org/#!Synapse:syn18915937/wiki/592740) | [Morabito et al., 2019 (biorxiv)](https://www.biorxiv.org/content/10.1101/695221v1) |
| Dataset | Species | Seq | [Data](link) | [Publication](link) |
| Dataset | Species | Seq | [Data](link) | [Publication](link) |



## ** Scripts to Figures Guides **
- Add when final version submitted
