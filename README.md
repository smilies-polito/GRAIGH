# Gene Regulation Accessibility Integrating GeneHancer (GRAIGH)

This repository contains the code of the Gene Regulation Accessibility Integrating GeneHancer (GRAIGH), a simple apporoach to integrate the GeneHancer elements with the scATAC-seq data. 

This methosd aims to analysize epigenetic single-cell data through uniquely identified features.
## Release notes

GRAIGH v1.0: 

* First release of GRAIGH.

## How to cite

### NSS Primary publication

* Martini, L., Amprimo, G., Di Carlo, S., Olmo, G., Ferraris, C., Savino, A. and Bardini, R., Neuronal Spike Shapes (NSS): a simple approach to study electrophysiological data for heterogeneity investigation, 2023 (submitted to Elsevier Computers in Biology and Medicine).

NSS validation relies on two datasets of murine cortical interneurons.
## Experimental setup

Follow these steps to setup for reproducing the experiments provided in _Martini et al., 2023_.
1) Install `Singularity` from https://docs.sylabs.io/guides/3.0/user-guide/installation.html:
	* Install `Singularity` release 3.10.2, with `Go` version 1.18.4
	* Suggestion: follow instructions provided in _Download and install singularity from a release_ section after installing `Go`
	* Install dependencies from: https://docs.sylabs.io/guides/main/admin-guide/installation.html
2) Clone the GRAIGH repository in your home folder
```
git clone https://github.com/smilies-polito/GRAIGH.git
```
3) Move to the GRAIGH source subfolder, and build the singularity container with 
```
mv GRAIGH/source
sudo singularity build GRAIGH.sif GRAIGH.def
```
or
```
mv GRAIGH/source
singularity build --fakeroot GRAIGH.sif GRAIGH.def
```

# Reproducing GRAIGH analysis

## Data required

In order to reproduce the analysis, it is necessary to gather required data files, organizing them in the repository folders as provided below.

### GeneHancer database
From Genecard, one can directly download only the older 2017 version https://www.genecards.org/GeneHancer_Version_4-4, but one can quickly request the authors' access to the latest versions from their online platform https://www.genecards.org/Guide/DatasetRequest. 
Specifically this work employs the files from the last version:
 *   `GeneHancer_v5.15.gff`
 *   `GeneHancer_AnnotSV_elements_v5.15.txt`
 *   `GeneHancer_AnnotSV_gene_association_scores_v5.15.txt`

They must be put in the `DATA/Genehancer` folder.

### scATAC-seq data
The 10X genomics dataset is freely available at https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard.
Specifically you need to download the files:
*  `Peak by cell matrix HDF5 (filtered)`
*  `Fragments (TSV)`
*  `Fragments index (TBI)`

All this files must go in `DATA/Human_PBMC` folder. The `filtered_peak_bc_matrix` must also be unzipped.

### Reference data

The reference data for the Seurat cell-type integration is available at https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds.
The file must be put in `DATA` folder.

## Reproducing the analysis running the Singularity container

To reproduce the analyses from _Martini et al., 2023_, run the `NSS.sif` container.

```
singularity run --no-home --bind  /local/path/to/GRAIGH:/local/path/to/home/ GRAIGH.sif
```

## Repository structure

```
|
├── DATA                                        // Data files
|    ├── GeneHancer                             // Data files for the PatchSeqDataset
|    └── Human_PBMC                             // Data files for the PatchClampDataset   
|   
├── Results                                     // Results of the analysis
|    ├── IMAGES                                 // IMAGES produced by the code
|    |    ├── ATAC_embedding.pdf                // ATAC 2D cells embedding 
|    |    ├── CD4_accessiblity.pdf              // CD4 accessibility og the GH element
|    |    ├── CD4_activity.pdf                  // CD4 activity
|    |    ├── CD14_accessibility.pdf            // CD14 accessibility og the GH element
|    |    ├── CD14_activity.pdf                 // CD14 activity
|    |    ├── Cell_types.pdf                    // Cell-type labels obtaine with Seurat integration
|    |    ├── GH_embedding.pdf                  // GH 2D cells embedding 
|    |    ├── GRAIGH WORKFLOW.png               // GRAIGH WORKFLOW and grafical abstract
|    |    └── peaks_per_GH.pdf                  // bar plot of peaks per GH element 
|    | 
|    └── Tables                                 // Folder for all the csv Tables
|         ├── cell_type_table.csv               // Cell-types labels for all the cells
|         └── gh_markers_ct.csv                 // Results of DA analysis of the GH matrix
| 
├── Source_Code                                 // Folder for the R code
|    └── GRAIGH.R                               // R script for the analysis
| 
├── renv.lock                                   // R virtual environment lock file
|
└── README.md                                   // This README file          
```
  
