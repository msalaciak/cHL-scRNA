# Longitudinal single-cell RNA-seq and TCR analysis in classical Hodgkin lymphoma

This repository contains the R analysis code used for the manuscript:

**Immunophenotyping of peripheral CD8+ T cells from classical Hodgkin Lymphoma patients treated with PD-1 blockade**

The analysis follows longitudinal peripheral blood mononuclear cell (PBMC) samples collected across six time points and includes PBMC preprocessing, differential expression analysis, T-cell reclustering and annotation, TCR clonotype analysis, exhaustion-signature analysis, CD8+ T-cell analysis, WGCNA, UCell scoring, and generation of manuscript figures.

All raw and processed scRNA-seq, scTCR-seq and TCR-seq data is publicly available in Gene Expression Omnibus (GEO) at (GSE337135, GSE337136). 
## Repository contents

```text
.
├── README.md
├── 01_PBMC_preprocessing.R
├── 02_PBMC_DEG_figures.R
└── 03_Tcell_TCR_CD8_analysis_figures.R
```

### `01_PBMC_preprocessing.R`

Performs the initial PBMC preprocessing workflow, including:

- loading six longitudinal 10x Genomics gene-expression datasets
- optional SoupX ambient-RNA correction
- creation of Seurat objects
- addition of clinical/time-point metadata
- removal of Scrublet-predicted doublets
- mitochondrial, ribosomal, and hemoglobin QC metrics
- miQC filtering
- merging of the six time points
- saving the QC-filtered PBMC object and QC summaries

The historical analysis used the original count matrices for downstream Seurat analysis; therefore, SoupX correction is retained as an optional documented step and is disabled by default.

### `02_PBMC_DEG_figures.R`

Uses the historical processed PBMC Seurat object to reproduce PBMC-level analyses and manuscript figures, including:

- PBMC UMAP
- canonical cell-type marker DotPlot
- cell-type proportions across time points
- differential expression across longitudinal time points
- Azimuth-based annotation validation

### `03_Tcell_TCR_CD8_analysis_figures.R`

Contains the downstream T-cell and CD8+ T-cell analysis, including:

- historical T-cell normalization, variable-feature selection, scaling, PCA, Harmony integration, UMAP, and clustering
- manual T-cell subset annotation
- T-cell manuscript figures
- Scirpy-derived single-cell TCR integration
- MiXCR bulk/tissue TCR matching
- blood/tissue clonotype tracking
- longitudinal expanded-clonotype analysis
- exhaustion-marker analysis
- Beltra intermediate- and terminal-exhaustion signature scoring
- CD8+ T-cell differential-expression analysis
- WGCNA metacell construction and network analysis
- UCell scoring of WGCNA modules
- manuscript and supplemental figure generation

## Data availability

**All data files required to run these scripts will be deposited on Zenodo.**

Zenodo record:

**[Zenodo DOI/link will be added here after publication]**

The Zenodo deposit will contain the processed Seurat objects and auxiliary input files required by the analysis, including the historical checkpoint objects used to reproduce the manuscript figures.

Expected data files include:

```text
data/
├── original_PBMC_final.rds
├── Tcells_original.rds
├── Tcells_processed.rds
├── CD8_processed.rds
├── tcr_info.csv
├── sig_beltra.csv
├── f15.csv
├── azmiuth_jan232022.csv
└── trb_mixcr/
```

The raw/initial PBMC preprocessing inputs used by `01_PBMC_preprocessing.R`, including the six 10x Genomics datasets and Scrublet doublet calls, will also be provided through the associated Zenodo data record where applicable.

After downloading the Zenodo data, arrange the project so that the required files are available under the paths specified at the beginning of each script.

## Historical processed-object checkpoints

Several downstream analyses use saved historical Seurat objects rather than regenerating the manuscript objects from scratch.

The main checkpoints are:

- `original_PBMC_final.rds` — processed PBMC object used for PBMC-level manuscript analyses
- `Tcells_original.rds` — T-cell subset before T-cell-specific reclustering and annotation
- `Tcells_processed.rds` — historical processed T-cell object used for downstream manuscript analyses and figures
- `CD8_processed.rds` — historical processed CD8+ T-cell subset used for downstream CD8 analyses

The processing steps that generated these objects are documented in the scripts, but the saved historical checkpoints are treated as the authoritative object state for figure reproduction. This avoids small differences that can arise when stochastic dimensional-reduction and clustering procedures are rerun under different software environments.

## Analysis workflow

The scripts are intended to be read/run in numerical order:

```text
01_PBMC_preprocessing.R
        ↓
02_PBMC_DEG_figures.R
        ↓
03_Tcell_TCR_CD8_analysis_figures.R
```

For exact reproduction of the manuscript figures, the processed historical objects supplied through Zenodo can be used directly by scripts 02 and 03.

## Directory structure

The scripts currently assume a project directory similar to:

```text
/home/rstudio/CIR_manuscript/
├── data/
├── figures/
├── objects/
└── results/
```

The main project path is defined near the beginning of each script and can be changed to match the user's local environment.

`01_PBMC_preprocessing.R` additionally expects locations for raw 10x data and Scrublet doublet calls:

```r
raw_dir <- "/data/10x"
doublet_dir <- "/data/doublets"
```

These paths should be updated as needed.

## Software

The analysis was originally performed using an R/Seurat 4-era workflow. The code was cleaned for public release while preserving the historical analytical parameters wherever possible.

Major R packages used across the scripts include:

- Seurat
- dplyr
- tidyr
- tibble
- readr
- ggplot2
- patchwork
- harmony
- miQC
- SoupX
- SeuratWrappers
- immunarch
- scWGCNA
- WGCNA
- UCell
- ComplexHeatmap
- circlize
- viridis
- pheatmap
- RColorBrewer

The downstream T-cell analysis was developed around **R 4.1.x / Seurat 4.0.x**. For the closest reproduction of the historical analysis, using an equivalent software environment is recommended.

Each script writes `sessionInfo()` output so that the package environment used for a given run can be recorded.

## TCR analysis

Single-cell TCR annotations were generated externally and imported through `tcr_info.csv`.

The analysis matches single-cell clonotypes to MiXCR repertoire data by TRB CDR3 amino-acid sequence and uses the resulting data to examine:

- clonotype expansion
- longitudinal clonotype persistence
- overlap between single-cell PBMC, bulk PBMC TCR, and tumor-tissue TCR repertoires
- blood/tissue-associated clonotypes

The MiXCR repertoire files used by these analyses are distributed in the Zenodo data deposit under `trb_mixcr/`.

## WGCNA and UCell

CD8+ T-cell co-expression analysis was performed using a metacell-based scWGCNA workflow.

Historical parameters include:

- metacell `k = 25`
- signed network
- bicor soft-threshold analysis
- selected soft-thresholding power = 18
- minimum module size = 50
- deep split = 4

WGCNA modules were subsequently scored at the single-cell level using UCell. The pink and yellow modules correspond to the memory-stem-like and terminal-effector-like signatures used in the manuscript.

## Outputs

The scripts write generated files to the following project directories:

```text
figures/
results/
objects/
```

Outputs include:

- manuscript figure panels
- supplemental figure panels
- differential-expression tables
- TCR metadata tables
- WGCNA module assignments
- QC summaries
- saved intermediate/result objects
- `sessionInfo()` files

## Notes on reproducibility

This repository is intended to document the analysis used for the manuscript as faithfully as possible.

Because the study analysis was developed over time, some final manuscript analyses rely on historical saved Seurat objects. The public scripts therefore distinguish between:

1. **documenting the historical processing workflow**, and
2. **using the saved historical processed object as the authoritative state for downstream manuscript figure reproduction**.

This is particularly relevant for stochastic procedures such as UMAP, Harmony, nearest-neighbor graph construction, and clustering, for which exact reruns can vary across package versions or computing environments.

## Citation

If you use this code or the associated dataset, please cite the corresponding manuscript:

> **Immunophenotyping of peripheral CD8+ T cells from classical Hodgkin Lymphoma patients treated with PD-1 blockade.**

Full publication details and DOI will be added once available.

## Data citation

The Zenodo dataset DOI will be added here once the data record is publicly released:

> **Zenodo:** [DOI/link pending]

## Contact

For questions regarding the analysis or repository, please open an issue in this GitHub repository.
