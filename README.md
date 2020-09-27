# Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts  
#### Code to reproduce analysis objects for the data contained in:  
[**Tissue Processing and Enzymatic Dissociation Induce Artifactual Gene Expression in Brain Cell Populations of Mice and Humans**](CITATION)  
Samuel E. Marsh<sup>1,\* </sup>, Tushar Kamath<sup>1</sup>, Timothy R. Hammond<sup>2</sup>, Alec J. Walker, Lasse Dissing-Olesen, Adam M.H. Young, Abdul Abdul, Naeem Nadaf, Connor Dufort, Sarah Murphy, Velina Kozareva<sup>2</sup>, Charles Vanderburg, Soyon Hong, Harry Bulstrode, Peter J. Hutchinson, Robin J.M. Franklin, Evan Macosko, & Beth Stevens

<sup><sup>1</sup>Performed analysis</sup>   
<sup><sup>2</sup>Assisted analysis</sup>  
<sup><sup>\*</sup>Analysis lead (contact: samuel.marsh@childrens.harvard.edu)</sup>  

## Code
Included is the code necessary to replicate the Seurat or LIGER (or both) objects used for analysis and plotting.
- Each R file specifies version of Seurat/LIGER used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.
    - Seurat V2.3.4 source package can be downloaded here from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source.
    - Seurat V3.2+ was released near the end of analysis.  To maintain consistency Seurat V3.1.5 was downloaded from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source when switching between V2 and V3 was necessary.  

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date the following actions are recommended or described in code:
  - Use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv](https://cran.r-project.org/web/packages/renv/index.html) packages. Followed by date-specific version installation of CRAN packages using [versions](https://cran.r-project.org/web/packages/versions/index.html) package.
  - Archived source versions of specific packages may also be needed depending on version of R and can be downloaded from CRAN archives and installed from local source.

- LIGER analyses were performed using the in development ["online"](https://github.com/MacoskoLab/liger/tree/online) branch, updating throughout analysis to accommodate bug fixes.  
  - LIGER analyses also utilize multiple versions of Seurat as specified in code for some of the following situations:
    - Seurat V3 used used for data import, QC filtering (genes, UMIs, % mito), and majority of plotting.
    - Seurat V2 was used during LIGER analysis workflow to accommodate use of now deprecated [`clusterLouvainJaccard` function](https://github.com/samuel-marsh/Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts/tree/master/08_Misc) which relied on Seurat V2 object structure.
    - Conversion between Seurat and LIGER objects was performed using built in LIGER functions `seuratToLiger` and `ligerToSeurat`.

## Data  
### Original Data
The data in this project can be broadly divided into 2 categories (5 sub-projects).  Please see [SI Table 1](LINK) (Mouse; Experiments 1-3) and [SI Table 2](LINK) (Human; Experiments 4-5) for breakdown by sample and more information.

***For experiments 1-4***, all the raw data (fastqs) and expression matrices are available at the NCBI Gene Expression Omnibus (GEO) under SuperSeries [GSE152184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152184).   

***For experiment 5***, the raw data (fastqs) and expression matrices are available through the European phenome-Genome Archive (EGA) (Accession ID: [EGADXXXXXXX](EGADXXXXXXX))

**scRNAseq of cells isolated from mouse brain (all cell types and sorted myeloid cells)**

  1. scRNAseq of microglia with 4 different dissociation protocols | 10X 3' V2 | 12 samples (n=3 per group) | ([GSE152183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152183))
  2. scRNAseq of all CNS cells with or without inhibitors | 10X 3' V2 | 4 samples (n=2 per group) | ([GSE152182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152182))
  3. scRNAseq of microglia (tail vein PBS injection) | 10X 3' V2 | 4 samples | ([GSE152210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152210))

**snRNAseq of nuclei isolated from frozen human brain tissue**

  4. snRNAseq of post-mortem brain tissue | 10X 3' V3.0 | 3 samples | ([GSE157760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157760))  
  5. snRNAseq of surgically resected brain tissue with or without freezing time delay | 10X 3' V3.0 | <br> 4 samples (n= 2; x2 timepoints per sample) | ([EGAXXXXXXX](EGAXXXXXXX))

### Processed Data
**Experiments 1-4**  
All processed data files from Cell Ranger `count` outputs are available via NCBI GEO.  Information on Cell Ranger version and Genome/Annotation for each experiment can be found in each GEO record and [SI Table 1](LINK) or [2](LINK).
There are 3 processed data files per library:
  1. GSM\*\_*Sample-Name*_barcodes.tsv.gz: corresponds to the cell barcodes (i.e. column names).
  2. GSM\*\_*Sample-Name*_features.tsv.gz: corresponds to the gene identifiers (i.e. row names).
  3. GSM\*\_*Sample-Name*_matrix.mtx.gz: expression matrix in sparse format.

### Raw fastq Files
All raw data fastq files can be downloaded from SRA linked from NCBI GEO records or from EGA.

### Literature Reanalysis
Reanalyzed data from literature is summarized detailed in table below.
| Dataset | Species | Seq Used | Raw/Count Data | Publication |
| :-----: | :-----: | :------: | :------------: | :---------: |
| Mathys | Mouse | scRNAseq (Smart-seq2) | [GEO103334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103334) <br> & Authors<sup>a</sup> | [Mathys et al., 2017 <br> (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31314-1?) |
| Plemel | Mouse | scRNAseq (10X 3' V2) | [GSE115803](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115803) | [Plemel et al., 2020 <br> (Science Advances)](https://advances.sciencemag.org/content/6/3/eaay6324) |
| Zywitza | Mouse | scRNAseq (Drop-Seq) | [GSE111527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527) | [Zywitza et al., 2018 <br> (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31732-7?) |
| Mizrak | Mouse | scRNAseq (Microwell Seq) | [GSE109447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109447) | [Mizrak et al., 2019 <br> (Cell Reports)](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31974-0?) |
| Zeisel | Mouse | scRNAseq (10X 3' V1 & V2) | [mousebrain.org](mousebrain.org) | [Zeisel et al., 2018 <br> (Cell)](https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X?) |
| Hammond | Mouse | scRNAseq (10X 3' V1 & V2) | [GSE121654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654) | [Hammond et al., 2019 <br> (Immunity)](https://www.cell.com/immunity/fulltext/S1074-7613(18)30485-0?) |
| Zhou | Human | snRNAseq (10X 5' V1) | [syn21670836](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn21670836) | [Zhou et al., 2020 <br> (Nature Medicine)](https://www.nature.com/articles/s41591-019-0695-9?) |
| Morabito<sup>i</sup> | Human | snRNAseq (10X 3' V3.0) | [syn18915937](https://www.synapse.org/#!Synapse:syn18915937/wiki/592740) | [Morabito et al., 2019 <br> (biorxiv)](https://www.biorxiv.org/content/10.1101/695221v1) |
| Leng & Li<sup>b</sup> | Human | snRNAseq (10X 3' V2) | [syn21788402](https://www.synapse.org/#!Synapse:syn21788402/wiki/601825)<sup>c</sup> <br> & [GSE147528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147528) | [Leng & Li et al., 2020 <br> (biorxiv)](https://www.biorxiv.org/content/10.1101/2020.04.04.025825v2) |

<sup><sup>a</sup>FPKM data and raw fastq files are available via GEO.  Raw count matrix was obtained via personal communication with authors.</sup>  
<sup><sup>b</sup>Additional metadata obtained via personal communication with the authors</sup>     
<sup><sup>c</sup>Data on synapse are post-QC and were used for re-analysis.  GEO records contain the all barcodes (unfiltered) HDF5 cellranger output files and fastqs.</sup>  
<sup><sup>i</sup>Reanalysis of Morabito et al., was also used for calculation of cell type proportions in [Liddelow, Marsh, & Stevens et al., 2020 (Trends in Immunology)](https://www.cell.com/trends/immunology/fulltext/S1471-4906(20)30155-1)</sup>

#### Human Data Reanalysis Meta Data
Meta data for human data was assembled from published SI Tables, public data on synapse, or restricted access data on synapse
  - Compiled publicly available meta data variables for each human dataset can be found in [SI Table 2](LINK).
  - "DUC" in the table indicates data available from synapse following submission and approval of Data Use Certificate.

### Acknowledgements:
This study was supported by funding from Cure Alzheimer's Fund (B.S.).  Special thanks to authors Tushar Kamath, Tim Hammond, Alec Walker, Lasse-Dissing-Olesen, Velina Kozareva, Evan Macosko, as well other members of Stevens and Macosko labs for helpful discussions and assistance during the analysis of this project.  

The analysis and results published here from Zhou et al., 2020 in whole or in part are based on data obtained from the [AMP-AD Knowledge Portal](https://adknowledgeportal.synapse.org/). Samples for this study were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute. Raw data used in analysis here are available from AMP-AD/Synapse database through links provided in table above.  Additional ROSMAP data can be requested at [https://www.radc.rush.edu](https://www.radc.rush.edu).

The analysis and results published here for Morabito et al., 2019 are based on reanalysis of study data downloaded from Synapse as provided by Dr. Vivek Swarup, Institute for Memory Impairments and Neurological Disorders, University of California, Irvine.  Data collection was supported through funding UCI Startup funds and American Federation of Aging Research.  Raw data used in analysis here are available from the Synapse database through link provided in table above.
