# PBMC preprocessing
# Longitudinal cHL scRNA-seq analysis


# Load packages
library(Seurat)
library(dplyr)
library(tibble)
library(miQC)
library(SoupX)
library(SeuratWrappers)


# Set directories
project_dir <- "/home/rstudio/CIR_manuscript"
raw_dir <- "/data/10x"
doublet_dir <- "/data/doublets"

dir.create(file.path(project_dir, "objects"),
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(file.path(project_dir, "results"),
           recursive = TRUE,
           showWarnings = FALSE)


# Sample IDs
sample_ids <- c(
  "151231",
  "160411",
  "161962",
  "171094",
  "171642",
  "180251"
)

names(sample_ids) <- paste0("TP", 1:6)


# Load filtered 10x matrices
timepoint_data <- lapply(sample_ids, function(sample_id) {
  
  Read10X(
    data.dir = file.path(
      raw_dir,
      sample_id,
      "outs",
      "filtered_gene_bc_matrices",
      "GRCh38"
    )
  )
  
})

names(timepoint_data) <- paste0("TP", 1:6)


# Check matrix dimensions
lapply(timepoint_data, dim)


# Optional SoupX correction
# The original script calculated SoupX-corrected matrices,
# but later created Seurat objects from the original matrices.
# Keep FALSE until we confirm which counts match the original analysis.

USE_SOUPX <- FALSE


run_soupx <- function(sample_id) {
  
  sc <- load10X(
    file.path(raw_dir, sample_id, "outs")
  )
  
  counts <- Read10X(
    data.dir = file.path(
      raw_dir,
      sample_id,
      "outs",
      "filtered_gene_bc_matrices",
      "GRCh38"
    )
  )
  
  seu <- CreateSeuratObject(counts = counts)
  
  seu <- SCTransform(seu, verbose = TRUE)
  seu <- RunPCA(seu, verbose = TRUE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = TRUE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = TRUE)
  seu <- FindClusters(seu, resolution = 0.8, verbose = TRUE)
  
  sc <- setClusters(
    sc,
    seu$seurat_clusters
  )
  
  platelet_genes <- c("PPBP", "PF4")
  
  non_expressing <- estimateNonExpressingCells(
    sc,
    nonExpressedGeneList = list(
      platelet = platelet_genes
    ),
    clusters = FALSE
  )
  
  sc <- calculateContaminationFraction(
    sc,
    list(
      platelet = platelet_genes
    ),
    useToEst = non_expressing
  )
  
  adjustCounts(sc)
}


# Run SoupX if selected
if (USE_SOUPX) {
  
  counts_for_analysis <- lapply(
    sample_ids,
    run_soupx
  )
  
  names(counts_for_analysis) <- paste0("TP", 1:6)
  
} else {
  
  counts_for_analysis <- timepoint_data
  
}


# Create Seurat objects
seurat_list <- lapply(seq_along(counts_for_analysis), function(i) {
  
  CreateSeuratObject(
    counts = counts_for_analysis[[i]],
    project = paste0("PBMC timepoint ", i),
    min.cells = 3,
    min.features = 200
  )
  
})

names(seurat_list) <- paste0("TP", 1:6)


# Add clinical metadata
clinical_metadata <- data.frame(
  
  Timepoint = as.character(1:6),
  
  Disease_status = c(
    "HL",
    "Remission",
    "Remission",
    "Relapse",
    "Relapse",
    "Relapse"
  ),
  
  Pembrolizumab = c(
    "No",
    "Yes",
    "Yes",
    "No",
    "Yes",
    "No"
  ),
  
  iRAE = c(
    "No",
    "Yes",
    "No",
    "No",
    "Yes",
    "No"
  )
  
)


for (i in seq_along(seurat_list)) {
  
  seurat_list[[i]]$Timepoint <-
    clinical_metadata$Timepoint[i]
  
  seurat_list[[i]]$Disease_status <-
    clinical_metadata$Disease_status[i]
  
  seurat_list[[i]]$Pembrolizumab <-
    clinical_metadata$Pembrolizumab[i]
  
  seurat_list[[i]]$iRAE <-
    clinical_metadata$iRAE[i]
  
}


# Load Scrublet doublet calls
doublet_files <- file.path(
  doublet_dir,
  paste0("tp", 1:6, "-doublet.txt")
)


# Remove Scrublet-predicted doublets
doublet_summary <- data.frame()

for (i in seq_along(seurat_list)) {
  
  doublets <- read.table(
    doublet_files[i],
    sep = ",",
    header = TRUE
  )
  
  doublets <- doublets[
    doublets$predicted_doublets == "True",
  ]
  
  before <- ncol(seurat_list[[i]])
  
  cells_to_keep <- colnames(seurat_list[[i]])[
    !colnames(seurat_list[[i]]) %in% doublets$barcode
  ]
  
  seurat_list[[i]] <- subset(
    seurat_list[[i]],
    cells = cells_to_keep
  )
  
  after <- ncol(seurat_list[[i]])
  
  doublet_summary <- rbind(
    doublet_summary,
    data.frame(
      Timepoint = paste0("TP", i),
      Cells_before = before,
      Predicted_doublets = nrow(doublets),
      Cells_after = after
    )
  )
  
}


# View doublet-removal summary
print(doublet_summary)


# Calculate QC metrics
for (i in seq_along(seurat_list)) {
  
  seurat_list[[i]] <- PercentageFeatureSet(
    seurat_list[[i]],
    pattern = "^MT-",
    col.name = "percent.mt"
  )
  
  seurat_list[[i]] <- PercentageFeatureSet(
    seurat_list[[i]],
    pattern = "^RP[SL]",
    col.name = "percent_ribo"
  )
  
  seurat_list[[i]] <- PercentageFeatureSet(
    seurat_list[[i]],
    pattern = "^HB[^(P)]",
    col.name = "percent_hb"
  )
  
}


# Run miQC on TP1-TP5
for (i in 1:5) {
  
  seurat_list[[i]] <- RunMiQC(
    seurat_list[[i]],
    percent.mt = "percent.mt",
    nFeature_RNA = "nFeature_RNA",
    posterior.cutoff = 0.75,
    model.slot = "flexmix_model"
  )
  
}


# Run miQC on TP6
seurat_list[[6]] <- RunMiQC(
  seurat_list[[6]],
  percent.mt = "percent.mt",
  nFeature_RNA = "nFeature_RNA",
  posterior.cutoff = 0.90,
  model.slot = "flexmix_model",
  backup.option = "percentile",
  backup.percentile = 0.95
)


# Keep cells passing miQC
miqc_summary <- data.frame()

for (i in seq_along(seurat_list)) {
  
  before <- ncol(seurat_list[[i]])
  
  seurat_list[[i]] <- subset(
    seurat_list[[i]],
    miQC.keep == "keep"
  )
  
  after <- ncol(seurat_list[[i]])
  
  miqc_summary <- rbind(
    miqc_summary,
    data.frame(
      Timepoint = paste0("TP", i),
      Cells_before = before,
      Cells_after = after
    )
  )
  
}


# View miQC summary
print(miqc_summary)


# Save QC summaries
write.csv(
  doublet_summary,
  file.path(
    project_dir,
    "results",
    "doublet_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  miqc_summary,
  file.path(
    project_dir,
    "results",
    "miQC_summary.csv"
  ),
  row.names = FALSE
)


# Merge the six filtered timepoints
pbmc.raw <- merge(
  seurat_list[[1]],
  y = seurat_list[2:6],
  add.cell.ids = c(
    "Timepoint_1",
    "Timepoint_2",
    "Timepoint_3",
    "Timepoint_4",
    "Timepoint_5",
    "Timepoint_6"
  )
)


# Check total cell count
print(pbmc.raw)

cat(
  "\nTotal cells after doublet removal and miQC:",
  ncol(pbmc.raw),
  "\n"
)


# Save the merged QC-filtered PBMC object
saveRDS(
  pbmc.raw,
  file.path(
    project_dir,
    "objects",
    "pbmc_raw_QC_filtered.rds"
  )
)


# Save package information
writeLines(
  capture.output(sessionInfo()),
  file.path(
    project_dir,
    "sessionInfo.txt"
  )
)