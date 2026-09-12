# T-cell, TCR, CD8, WGCNA, UCell and manuscript figures
#
# Analysis outline
#
# 1. Load data and historical checkpoints
#
# 2. Reproduce T-cell processing
#    2.1 Normalization and variable features
#    2.2 PCA and Harmony
#    2.3 UMAP and clustering
#    2.4 T-cell annotation
#
# 3. T-cell figures
#    Figure 2a - T-cell UMAP
#    Figure 2b - Canonical T-cell marker expression
#
# 4. TCR analysis
#    4.1 Scirpy clonotypes
#    4.2 MiXCR bulk/tissue matching
#    Figure 2c - Blood/tissue-associated clonotypes
#    Figure 2d - TCR clonotype tracking
#    Figure 2e - Expanded scRNA clonotypes
#    Figure 2f - Longitudinal clonotypes
#
# 5. Exhaustion analysis
#    Figure 3a - Exhaustion marker heatmap
#    Figure 3b - Exhaustion FeaturePlots
#    Supplemental Figure 7 - Beltra Tex-int/Tex-term signatures
#
# 6. CD8 analysis
#    Supplemental Figure 4 - CD8 TEM-3 DEG heatmap
#    Supplemental Figure 5 - CD8 subset DEG plot
#
# 7. WGCNA
#    7.1 CD8 metacell construction
#    7.2 Network construction
#    7.3 Module identification
#
# 8. UCell and module-associated analysis
#    Figure 4a - Memory-stem-like vs terminal-effector-like DEG heatmap
#    Supplemental Figure 6 - WGCNA/UCell module scores
#
# Saved checkpoints:
# tcells.original = original T-cell subset from PBMC analysis
# tcells          = historical processed T-cell object
# cd8.sub         = historical processed CD8 subset
#
# Historical checkpoints are never overwritten.

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(harmony)
library(immunarch)
library(scWGCNA)
library(WGCNA)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(pheatmap)

options(stringsAsFactors = FALSE)

project_dir <- "/home/rstudio/CIR_manuscript"
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")
result_dir <- file.path(project_dir, "results")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

tcell_original_file <- file.path(data_dir, "Tcells_original.rds")
tcell_processed_file <- file.path(data_dir, "Tcells_processed.rds")
cd8_processed_file <- file.path(data_dir, "CD8_processed.rds")
tcr_file <- file.path(data_dir, "tcr_info.csv")
mixcr_dir <- file.path(data_dir, "trb_mixcr")
beltra_file <- file.path(data_dir, "sig_beltra.csv")
f15_file <- file.path(data_dir, "f15.csv")

stopifnot(
  file.exists(tcell_original_file),
  file.exists(tcell_processed_file),
  file.exists(cd8_processed_file),
  file.exists(tcr_file),
  dir.exists(mixcr_dir),
  file.exists(beltra_file)
)

tcells.original <- readRDS(tcell_original_file)
tcells <- readRDS(tcell_processed_file)
cd8.sub <- readRDS(cd8_processed_file)
tcr_info <- as.data.frame(read_csv(tcr_file))
immdata_mixcr <- repLoad(mixcr_dir)
sig_beltra <- read_csv(beltra_file)

print(tcells.original)
print(tcells)
print(cd8.sub)

# 2. Reproduce T-cell processing
tcells.reprocessed <- tcells.original
tcells.reprocessed <- NormalizeData(tcells.reprocessed)
tcells.reprocessed <- FindVariableFeatures(tcells.reprocessed)
tcells.reprocessed <- ScaleData(tcells.reprocessed, features = rownames(tcells.reprocessed))
tcells.reprocessed <- RunPCA(tcells.reprocessed, features = VariableFeatures(tcells.reprocessed))
tcells.reprocessed <- RunHarmony(tcells.reprocessed, "Timepoint", plot_convergence = FALSE)
tcells.reprocessed <- RunUMAP(tcells.reprocessed, reduction = "harmony", dims = 1:20)
tcells.reprocessed <- FindNeighbors(tcells.reprocessed, reduction = "harmony", dims = 1:20)
tcells.reprocessed <- FindClusters(tcells.reprocessed, resolution = 0.7)

tcell.markers <- FindAllMarkers(
  tcells.reprocessed,
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.09
) %>% mutate(Difference = pct.1 - pct.2)

write.csv(tcell.markers, file.path(result_dir, "Tcell_cluster_markers.csv"), row.names = FALSE)

tcell.ids <- c(
  `0`="CD4 TH2", `1`="CD4 TCM-1", `2`="CD8 TEM-1", `3`="CD8 TEM-2",
  `4`="CD8 TEM-3", `5`="CD8 Exhausted-1", `6`="CD8 TEM-4", `7`="CD4 Tregs",
  `8`="CD4 TH17", `9`="CD8 TEM-5", `10`="CD8 Exhausted-2", `11`="Gd T cells"
)
tcells.reprocessed <- RenameIdents(tcells.reprocessed, tcell.ids)
tcells.reprocessed$ident <- Idents(tcells.reprocessed)

tcell.bulk <- c(
  "CD4 TH2"="CD4","CD4 TCM-1"="CD4","CD4 Tregs"="CD4","CD4 TH17"="CD4",
  "CD8 TEM-1"="CD8","CD8 TEM-2"="CD8","CD8 TEM-3"="CD8","CD8 TEM-4"="CD8",
  "CD8 TEM-5"="CD8","CD8 Exhausted-1"="CD8","CD8 Exhausted-2"="CD8",
  "Gd T cells"="Gd T cells"
)
tcells.reprocessed$bulk <- unname(tcell.bulk[as.character(tcells.reprocessed$ident)])

tcells.reprocessed$Groups <- case_when(
  tcells.reprocessed$Timepoint == 1 ~ "1",
  tcells.reprocessed$Timepoint == 2 ~ "9",
  tcells.reprocessed$Timepoint == 3 ~ "2",
  tcells.reprocessed$Timepoint == 4 ~ "3",
  tcells.reprocessed$Timepoint == 5 ~ "9",
  tcells.reprocessed$Timepoint == 6 ~ "4"
)
rm(tcells.reprocessed)

# 3. T-cell figures

# Figure 2a
Idents(tcells) <- "ident"
p_fig2a <- DimPlot(tcells, label = FALSE, group.by = "ident") + theme_void()
ggsave(file.path(figure_dir, "Figure2a_Tcell_UMAP.png"), p_fig2a, width=6, height=5, dpi=600)

# Figure 2b
tcell.genes <- c(
  "TCF7","IL7R","SELL","GATA3","PRF1","KLRG1","GZMB","GZMK","GZMH",
  "FOXP3","CTLA4","KLRB1","PDCD1","LAG3","TIGIT","HAVCR2","STMN1","TRDV3"
)
tcell.genes <- tcell.genes[tcell.genes %in% rownames(tcells)]
p_fig2b <- DotPlot(tcells, features=tcell.genes, cols=c("lightgrey","red")) +
  coord_flip() + RotatedAxis()
ggsave(file.path(figure_dir, "Figure2b_Tcell_markers.png"), p_fig2b, width=8, height=4.75, dpi=600)

# 4. TCR analysis

# 4.1 Scirpy clonotypes and MiXCR matching
tcr <- tcr_info[,c(1,16,17,22:25,32,33,38:53,58:62,67,68)]
tcr$clone_id <- paste0("TCR_", as.character(tcr$clone_id + 1))
tcr$CDR3.aa <- tcr$IR_VDJ_1_junction_aa

join_mixcr <- function(x, sample_name, remove_cols, new_name) {
  x <- left_join(x, immdata_mixcr$data[[sample_name]], by="CDR3.aa")
  x <- x[, -remove_cols]
  x <- distinct(x, ...1, .keep_all=TRUE)
  rename(x, !!new_name := Clones)
}

tcr <- join_mixcr(tcr, "151389", 35:51, "Clones_151389")
tcr <- join_mixcr(tcr, "151403-1a", 36:52, "Clones_151403_1a")
tcr <- join_mixcr(tcr, "151403-2a", 37:53, "Clones_151403_2a")
tcr <- join_mixcr(tcr, "171094", 38:54, "Clones_171094")
tcr <- join_mixcr(tcr, "180387", 39:55, "Clones_180387")
tcr <- join_mixcr(tcr, "f1841", 40:56, "Clones_f1841")

clone_columns <- c(
  "Clones_151389","Clones_151403_1a","Clones_151403_2a",
  "Clones_171094","Clones_180387","Clones_f1841"
)
tcr[clone_columns][is.na(tcr[clone_columns])] <- 0
tcr <- rename(tcr, barcode = ...1)

tcr <- tcr %>%
  rowwise() %>%
  mutate(Bulk_Clones_Total = sum(c_across(all_of(clone_columns)))) %>%
  ungroup() %>%
  mutate(
    in_bulk = ifelse(Bulk_Clones_Total >= 1, "Yes", "No"),
    in_tissue_relapse = ifelse(Clones_f1841 >= 1, "Yes", "No"),
    in_tissue = ifelse(
      Clones_f1841 >= 1 | Clones_151403_1a >= 1 | Clones_151403_2a >= 1,
      "Yes","No"
    )
  )

tcells.metadata <- tcells@meta.data
tcells.metadata$barcode <- rownames(tcells.metadata)
tcells.metadata <- left_join(tcells.metadata, tcr, by="barcode")
tcells.metadata <- tcells.metadata %>%
  mutate(Clonotype_Size = case_when(
    clone_id_size >= 20 ~ "Large",
    clone_id_size > 1 & clone_id_size < 20 ~ "Small",
    clone_id_size == 1 ~ "Single",
    TRUE ~ "0"
  ))

# Historical pre-pembrolizumab tissue file
if (file.exists(f15_file)) {
  f15 <- as.data.frame(read_csv(f15_file))
  tcells.metadata <- left_join(tcells.metadata, f15, by="CDR3.aa", suffix=c("",".f15"))
  if ("Clones" %in% colnames(tcells.metadata)) {
    tcells.metadata <- rename(tcells.metadata, Clones_f15244 = Clones)
  } else if ("Clones.f15" %in% colnames(tcells.metadata)) {
    tcells.metadata <- rename(tcells.metadata, Clones_f15244 = Clones.f15)
  }
  if ("Clones_f15244" %in% colnames(tcells.metadata)) {
    tcells.metadata$Clones_f15244[is.na(tcells.metadata$Clones_f15244)] <- 0
  }
}

write.csv(tcells.metadata, file.path(result_dir, "TCR_metadata_processed.csv"), row.names=FALSE)

get_tcr_cells <- function(id) filter(tcells.metadata, clone_id == id)$barcode
TCR_1 <- get_tcr_cells("TCR_1"); TCR_3 <- get_tcr_cells("TCR_3")
TCR_5 <- get_tcr_cells("TCR_5"); TCR_14 <- get_tcr_cells("TCR_14")
TCR_27 <- get_tcr_cells("TCR_27"); TCR_28 <- get_tcr_cells("TCR_28")
TCR_51 <- get_tcr_cells("TCR_51"); TCR_74 <- get_tcr_cells("TCR_74")
TCR_84 <- get_tcr_cells("TCR_84"); TCR_102 <- get_tcr_cells("TCR_102")
TCR_135 <- get_tcr_cells("TCR_135"); TCR_259 <- get_tcr_cells("TCR_259")
TCR_312 <- get_tcr_cells("TCR_312"); TCR_532 <- get_tcr_cells("TCR_532")
TCR_555 <- get_tcr_cells("TCR_555"); TCR_617 <- get_tcr_cells("TCR_617")
TCR_1669 <- get_tcr_cells("TCR_1669")


# Figure 2c - Blood/tissue-associated clonotypes
pre.blood <- tcells.metadata %>%
  filter(Clones_151389 != 0, has_ir == "True") %>%
  distinct(clone_id, .keep_all=TRUE) %>%
  slice_max(n=19, order_by=Clones_151389)
pre.blood <- filter(tcells.metadata, clone_id %in% pre.blood$clone_id)$barcode

relapse.blood <- tcells.metadata %>%
  filter(Clones_171094 != 0, has_ir == "True") %>%
  distinct(clone_id, .keep_all=TRUE) %>%
  slice_max(n=19, order_by=Clones_171094)
relapse.blood <- filter(tcells.metadata, clone_id %in% relapse.blood$clone_id)$barcode

post.blood <- tcells.metadata %>%
  filter(Clones_180387 != 0, has_ir == "True") %>%
  distinct(clone_id, .keep_all=TRUE) %>%
  slice_max(n=19, order_by=Clones_180387)
post.blood <- filter(tcells.metadata, clone_id %in% post.blood$clone_id)$barcode

bulk.blood <- unique(c(post.blood, pre.blood, relapse.blood))

post.tissue <- tcells.metadata %>%
  filter(Clones_f1841 != 0, has_ir == "True") %>%
  distinct(clone_id, .keep_all=TRUE) %>%
  slice_max(n=15, order_by=Clones_f1841)
post.tissue <- filter(tcells.metadata, clone_id %in% post.tissue$clone_id)$barcode

if ("Clones_f15244" %in% colnames(tcells.metadata)) {
  pre.tissue <- tcells.metadata %>%
    filter(Clones_f15244 != 0, has_ir == "True") %>%
    distinct(clone_id, .keep_all=TRUE) %>%
    slice_max(n=15, order_by=Clones_f15244)
  pre.tissue <- filter(tcells.metadata, clone_id %in% pre.tissue$clone_id)$barcode

  p_fig2c <- DimPlot(
    tcells,
    label=FALSE,
    repel=TRUE,
    cells.highlight=list(unique(bulk.blood), post.tissue, pre.tissue)
  ) +
    scale_color_manual(
      labels=c(
        "unselected",
        "Top Pre-Pembro Tissue Clones",
        "Top Pembro Progression Tissue Clones",
        "Top Blood Clones"
      ),
      values=c("lightgrey","coral","orchid4","skyblue3")
    ) +
    theme_void()

  ggsave(
    file.path(figure_dir, "Figure2c_TCR_blood_tissue_UMAP.png"),
    p_fig2c, width=6, height=3.5, dpi=600
  )
} else {
  warning("Figure 2c is missing the pre-pembrolizumab tissue component because data/f15.csv is not present.")
}


# Figure 2d - TCR clonotype tracking

clones <- c(
  trackClonotypes(
    immdata_mixcr$data,
    list("151389", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("f15244", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("f1841", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("180387", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("151231_sc", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("160411_sc", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("161962_sc", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("171094_sc", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("171642_sc", 2),
    .col = "aa"
  )$CDR3.aa,

  trackClonotypes(
    immdata_mixcr$data,
    list("180251_sc", 2),
    .col = "aa"
  )$CDR3.aa
)

clones <- unique(
  clones
)

p_fig2d <- trackClonotypes(
  immdata_mixcr$data,
  clones,
  .col = "aa"
) %>%
  vis(
    .order = c(
      1,
      10,
      2,
      3,
      4,
      5,
      6,
      7,
      8,
      9,
      11
    )
  ) +
  theme_classic() +
  scale_x_discrete(
    labels = c(
      "pbmc scRNA-tcr",
      "tissue bulk-tcr",
      "pbmc bulk-tcr",
      "pbmc scRNA-tcr",
      "pbmc scRNA-tcr",
      "pbmc scRNA-tcr",
      "pbmc bulk-tcr",
      "pbmc scRNA-tcr",
      "pbmc scRNA-tcr",
      "pbmc bulk-tcr",
      "tissue bulk-tcr"
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2d_TCR_clonotype_tracking.png"
  ),
  p_fig2d,
  width = 8,
  height = 4.5,
  dpi = 600
)


# Figure 2e - most expanded scRNA clonotypes
p_fig2e <- DimPlot(
  tcells,
  label=FALSE,
  repel=TRUE,
  cells.highlight=list(TCR_3,TCR_51,TCR_14,TCR_28,TCR_27)
) +
  scale_color_manual(
    labels=c("unselected","TCR_27","TCR_28","TCR_14","TCR_51","TCR_3"),
    values=c("lightgrey","black","deeppink2","blue","purple","coral3")
  ) +
  theme_void()

ggsave(
  file.path(figure_dir, "Figure2e_top_scRNA_clonotypes.png"),
  p_fig2e, width=6, height=3.5, dpi=600
)


# Figure 2f - clonotypes across timepoints
# Historical split-by-timepoint clone plot found.
p_fig2f_candidate <- DimPlot(
  tcells,
  label=FALSE,
  repel=TRUE,
  cells.highlight=list(
    TCR_259,TCR_135,TCR_532,TCR_312,TCR_3,TCR_51,TCR_14,TCR_28,TCR_27
  ),
  split.by="Timepoint"
)

ggsave(
  file.path(figure_dir, "Figure2f_candidate_clonotypes_by_timepoint.png"),
  p_fig2f_candidate, width=12, height=7, dpi=600
)

# Another historical Figure 2f-like plot used metadata named TCR_Groups.
# Its definition was not found in the supplied scripts.
if ("TCR_Groups" %in% colnames(tcells@meta.data)) {
  p_fig2f_tcrgroups <- DimPlot(
    subset(tcells, subset=has_ir=="True"),
    group.by="TCR_Groups"
  ) + theme_classic()

  ggsave(
    file.path(figure_dir, "Figure2f_TCR_Groups_historical.png"),
    p_fig2f_tcrgroups, width=6, height=4, dpi=600
  )
}


# 5. Exhaustion analysis

# Figure 3a - Exhaustion average-expression heatmap
Idents(tcells) <- "bulk"

cd8.heatmap.timepoint <- AverageExpression(
  subset(tcells, idents="CD8"),
  verbose=TRUE,
  group.by="Timepoint"
)

cd8.heat <- log1p(cd8.heatmap.timepoint$RNA)

pre.exhaust <- c(
  "HAVCR2","CD244","PDCD1","ENTPD1","LAG3","CTLA4","CD160","SLAMF6",
  "CXCR3","IL7R","SELL","CD28","CD69","ID2","PRDM1","TOX","TBX21",
  "EOMES","TIGIT","NR4A2","CX3CR1","TCF7"
)
pre.exhaust <- pre.exhaust[pre.exhaust %in% rownames(cd8.heat)]
cd8.heat <- cd8.heat[pre.exhaust,,drop=FALSE]
cd8.heat <- t(scale(t(cd8.heat)))

png(
  file.path(figure_dir, "Figure3a_exhaustion_heatmap.png"),
  height=10, width=6, units="in", res=1200
)
Heatmap(
  cd8.heat,
  col=colorRamp2(seq(-1,2,0.1), viridis(length(seq(-1,2,0.1)))),
  cluster_rows=TRUE,
  cluster_columns=FALSE,
  row_names_side="left",
  column_names_rot=0,
  column_names_side="top",
  row_names_gp=gpar(fontsize=10),
  row_dend_reorder=TRUE,
  show_row_dend=FALSE,
  column_dend_side="bottom"
)
dev.off()


# Figure 3b - TP3 and TP4 exhaustion FeaturePlots
p_fig3b_tp3 <- FeaturePlot(
  subset(tcells, subset=Timepoint==3),
  features=c("TBX21","CX3CR1","PRDM1","BATF","ZEB2","TIGIT"),
  order=TRUE, label=FALSE, cols=c("lightgrey","red"),
  combine=TRUE, max.cutoff="q99", ncol=2
) & theme_void()

ggsave(
  file.path(figure_dir, "Figure3b_Timepoint3.png"),
  p_fig3b_tp3, width=8, height=6, dpi=600
)

p_fig3b_tp4 <- FeaturePlot(
  subset(tcells, subset=Timepoint==4),
  features=c("TBX21","CX3CR1","PRDM1","BATF","ZEB2","TIGIT"),
  order=TRUE, label=FALSE, cols=c("lightgrey","red"),
  combine=TRUE, max.cutoff="q99", ncol=2
) & theme_void()

ggsave(
  file.path(figure_dir, "Figure3b_Timepoint4.png"),
  p_fig3b_tp4, width=8, height=6, dpi=600
)


# Beltra exhaustion signatures
sig_beltra <- sig_beltra %>%
  mutate(across(where(is.character), toupper))
gene_sets <- lapply(sig_beltra, function(x) x[!is.na(x)])

tcells.beltra <- AddModuleScore(
  tcells,
  features=list(
    gene_sets$TextInt,
    gene_sets$TextTerm1,
    gene_sets$TextTerm2,
    gene_sets$Textprog1,
    gene_sets$Textprog2,
    gene_sets$Textprog3,
    gene_sets$Textprog4
  ),
  name=c("TexInt","TextTerm1","TextTerm2","TextProg1","TextProg2","TextProg3","TextProg4")
)


# Supplemental Figure 7 - Beltra signatures
tcells_scores <- FetchData(
  tcells,
  vars=c("ident","TexInt1","TextTerm12","TextTerm23","TextProg14","TextProg25","TextProg36","TextProg47")
)

tcells_scores_long <- tcells_scores %>%
  pivot_longer(
    cols=c(TexInt1,TextTerm12,TextTerm23,TextProg14,TextProg25,TextProg36,TextProg47),
    names_to="Signature",
    values_to="Score"
  )

cd8_idents <- c(
  "CD8 TEM-1","CD8 TEM-2","CD8 TEM-3","CD8 TEM-4","CD8 TEM-5",
  "CD8 Exhausted-1","CD8 Exhausted-2"
)

tcells_scores_filtered <- tcells_scores_long %>%
  filter(
    ident %in% cd8_idents,
    Signature %in% c("TexInt1","TextTerm23")
  ) %>%
  mutate(Signature=recode(Signature, TexInt1="Tex-int", TextTerm23="Tex-term"))

beltra_mat <- tcells_scores_filtered %>%
  group_by(ident, Signature) %>%
  summarise(mean_score=mean(Score, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=Signature, values_from=mean_score) %>%
  column_to_rownames("ident") %>%
  as.matrix()

png(
  file.path(figure_dir, "Supplemental_Figure7_Beltra.png"),
  width=2400, height=1800, res=300
)
pheatmap(beltra_mat, scale="row")
dev.off()


# 6. CD8 analysis

# CD8 setup
cd8.sub$ident <- droplevels(cd8.sub$ident)
Idents(cd8.sub) <- "ident"


# CD8 DEGs across timepoints
Idents(cd8.sub) <- "Timepoint"
cd8.groups.markers <- FindAllMarkers(
  cd8.sub,
  only.pos=FALSE,
  min.diff.pct=-Inf,
  logfc.threshold=0.10,
  min.pct=0.05
) %>% mutate(Difference=pct.1-pct.2)

cd8.groups.markers.pval <- filter(cd8.groups.markers, p_val_adj<=0.05)

write.csv(
  cd8.groups.markers,
  file.path(result_dir, "CD8_timepoint_DEGs_all.csv"),
  row.names=FALSE
)
write.csv(
  cd8.groups.markers.pval,
  file.path(result_dir, "CD8_timepoint_DEGs_significant.csv"),
  row.names=FALSE
)


# Per-CD8-subset timepoint DEGs
Idents(tcells) <- "ident"
cd8_cluster_names <- c(
  "CD8 TEM-1","CD8 TEM-2","CD8 TEM-3","CD8 TEM-4","CD8 TEM-5",
  "CD8 Exhausted-1","CD8 Exhausted-2"
)

cd8_cluster_timepoint_markers <- list()

for (cluster_name in cd8_cluster_names) {
  tmp <- subset(tcells, idents=cluster_name)
  tmp <- NormalizeData(tmp, assay="RNA")
  Idents(tmp) <- "Timepoint"

  tmp_markers <- FindAllMarkers(
    tmp,
    only.pos=FALSE,
    min.diff.pct=-Inf,
    logfc.threshold=0.10,
    min.pct=0.05
  ) %>% mutate(Difference=pct.1-pct.2)

  cd8_cluster_timepoint_markers[[cluster_name]] <- tmp_markers
}
rm(tmp, tmp_markers)


# Supplemental Figure 4 - CD8 TEM-3 heatmap
Idents(tcells) <- "ident"

cd8tem3.relapse.markers <- FindMarkers(
  subset(tcells, subset=ident=="CD8 TEM-3"),
  ident.1=4,
  group.by="Timepoint"
)
cd8tem3.relapse.markers$gene <- rownames(cd8tem3.relapse.markers)
cd8tem3.relapse.markers <- cd8tem3.relapse.markers %>%
  mutate(Difference=pct.1-pct.2)

cd8tem3.relapse.markers.pval <- cd8tem3.relapse.markers %>%
  filter(p_val_adj<=0.05, abs(Difference)>=0.1)

cluster.averages <- subset(tcells, idents="CD8 TEM-3")
Idents(cluster.averages) <- "Timepoint"
cluster.averages <- AverageExpression(
  cluster.averages,
  verbose=TRUE,
  return.seurat=TRUE,
  group.by="Timepoint"
)
cluster.averages$orig.ident <- rownames(cluster.averages@meta.data)

p_supp4 <- DoHeatmap(
  cluster.averages,
  features=unique(cd8tem3.relapse.markers.pval$gene),
  size=3,
  draw.lines=FALSE,
  group.by="orig.ident",
  disp.min=-1
) +
  scale_fill_viridis_c() +
  theme(axis.text.y=element_text(size=8))

ggsave(
  file.path(figure_dir, "Supplemental_Figure4_CD8_TEM3_DEG_heatmap.png"),
  p_supp4, width=6, height=10, dpi=600
)


# Supplemental Figure 5 - top DEGs from each CD8 subset across timepoints
sig_cluster_markers <- lapply(
  cd8_cluster_timepoint_markers,
  function(x) filter(x, p_val_adj<=0.05)
)

top_cluster_genes <- lapply(sig_cluster_markers, function(x) {
  top <- rbind(
    x %>% group_by(cluster) %>% slice_max(n=10, order_by=avg_log2FC),
    x %>% group_by(cluster) %>% slice_min(n=10, order_by=avg_log2FC)
  )

  unwanted <- paste0(
    "^IGHV|^IGHJ|^IGHD|^IGKV|^IGLV|^TRBV|^TRBD|^TRBJ|",
    "^TRDV|^TRDD|^TRDJ|^TRAV|^TRAJ|^TRGV|^TRGJ|^TRD"
  )
  remove <- top$gene[grep(unwanted, top$gene)]
  top[!top$gene %in% remove,]
})

supp5_genes <- unique(unlist(lapply(top_cluster_genes, function(x) x$gene)))

p_supp5 <- DotPlot(
  subset(tcells, subset=bulk=="CD8"),
  features=supp5_genes,
  group.by="ident",
  col.min=-1.5,
  col.max=1,
  cluster.idents=TRUE,
  cols=c("lightgrey","red")
) +
  coord_flip() +
  RotatedAxis()

ggsave(
  file.path(figure_dir, "Supplemental_Figure5_CD8_top_DEGs.png"),
  p_supp5, width=8, height=12, dpi=600
)


# 7. WGCNA
cd8.wgcna <- NormalizeData(cd8.sub)
cd8.wgcna <- FindVariableFeatures(cd8.wgcna, nfeatures=2000)
Idents(cd8.wgcna) <- "ident"

seurat_list <- list()
for (group in unique(as.character(cd8.wgcna$ident))) {
  print(group)
  cur_seurat <- subset(cd8.wgcna, subset=ident==group)

  cur_metacell_seurat <- scWGCNA::construct_metacells(
    cur_seurat,
    name=group,
    k=25,
    reduction="umap",
    assay="RNA",
    slot="data"
  )

  seurat_list[[group]] <- cur_metacell_seurat
}

metacell_seurat <- merge(
  seurat_list[[1]],
  seurat_list[2:length(seurat_list)]
)

rm(cur_seurat, cur_metacell_seurat)

genes.use <- rownames(metacell_seurat)

datExpr <- as.data.frame(
  GetAssayData(
    metacell_seurat,
    assay="RNA",
    slot="data"
  )[genes.use,]
)

datExpr <- as.data.frame(t(datExpr))
datExpr <- datExpr[,goodGenes(datExpr)]

powers <- c(seq(1,10,by=1), seq(12,30,by=2))

powerTable <- list(
  data=pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose=100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
)

softPower <- 18

multiExpr <- list()
multiExpr[["ODC"]] <- list(data=datExpr)
checkSets(multiExpr)

net <- blockwiseConsensusModules(
  multiExpr,
  blocks=NULL,
  maxBlockSize=30000,
  randomSeed=12345,
  corType="pearson",
  power=softPower,
  consensusQuantile=0.3,
  networkType="signed",
  TOMType="unsigned",
  TOMDenom="min",
  scaleTOMs=TRUE,
  scaleQuantile=0.8,
  sampleForScaling=TRUE,
  sampleForScalingFactor=1000,
  useDiskCache=TRUE,
  chunkSize=NULL,
  deepSplit=4,
  pamStage=FALSE,
  detectCutHeight=0.995,
  minModuleSize=50,
  mergeCutHeight=0.2,
  saveConsensusTOMs=TRUE,
  consensusTOMFilePattern=file.path(
    result_dir,
    "ConsensusTOM-block.%b.rda"
  )
)

moduleLabels <- net$colors
moduleColors <- as.character(moduleLabels)

MEs <- moduleEigengenes(
  multiExpr[[1]]$data,
  colors=moduleColors,
  nPC=1
)$eigengenes
MEs <- orderMEs(MEs)

KMEs <- signedKME(
  datExpr,
  MEs,
  outputColumnName="kME",
  corFnc="bicor"
)

geneInfo <- as.data.frame(
  cbind(
    colnames(datExpr),
    moduleColors,
    KMEs
  )
)
colnames(geneInfo)[1:2] <- c(
  "GeneSymbol",
  "Initially.Assigned.Module.Color"
)

write.csv(
  geneInfo,
  file.path(result_dir, "WGCNA_gene_modules.csv"),
  row.names=FALSE
)


# 8. UCell and module-associated analysis

# UCell signatures from WGCNA modules
markers <- list()
markers$black <- colnames(datExpr)[net$colors=="black"]
markers$blue <- colnames(datExpr)[net$colors=="blue"]
markers$brown <- colnames(datExpr)[net$colors=="brown"]
markers$green <- colnames(datExpr)[net$colors=="green"]
markers$pink <- colnames(datExpr)[net$colors=="pink"]
markers$red <- colnames(datExpr)[net$colors=="red"]
markers$turquoise <- colnames(datExpr)[net$colors=="turquoise"]
markers$yellow <- colnames(datExpr)[net$colors=="yellow"]
markers$ex <- c(
  "PDCD1","CTLA4","HAVCR2","ENTPD1","CD160","CD244","TIGIT","TBX21","LAG3"
)

test <- AddModuleScore_UCell(
  tcells,
  features=markers
)

signature.names <- paste0(names(markers), "_UCell")


# Supplemental Figure 6
p_supp6 <- VlnPlot(
  subset(test, subset=bulk=="CD8"),
  features=signature.names[1:8],
  pt.size=0,
  same.y.lims=FALSE
) &
  geom_boxplot()

ggsave(
  file.path(figure_dir, "Supplemental_Figure6_WGCNA_UCell.png"),
  p_supp6, width=10, height=9, dpi=600
)


# Memory-stem-like and terminal-effector-like module maps
p_memory <- FeaturePlot(
  subset(test, subset=bulk=="CD8"),
  features="pink_UCell",
  order=TRUE,
  label=FALSE,
  cols=c("lightgrey","red"),
  min.cutoff="q45",
  max.cutoff="q99"
) &
  theme_void() &
  NoLegend()

p_terminal <- FeaturePlot(
  subset(test, subset=bulk=="CD8"),
  features="yellow_UCell",
  order=TRUE,
  label=FALSE,
  cols=c("lightgrey","red"),
  min.cutoff="q35",
  max.cutoff="q99"
) &
  theme_void()

ggsave(
  file.path(figure_dir, "WGCNA_memory_terminal_modules.png"),
  p_memory | p_terminal,
  width=10, height=4, dpi=600
)


# Mean WGCNA/UCell signature scores by original CD8 subset
test.meta <- test@meta.data

mean_score <- function(cluster, score) {
  mean(
    filter(test.meta, ident==cluster)[[score]],
    na.rm=TRUE
  )
}

signat.heatmap <- data.frame(
  `CD8 TEM-1`=c(mean_score("CD8 TEM-1","pink_UCell"),mean_score("CD8 TEM-1","yellow_UCell")),
  `CD8 TEM-2`=c(mean_score("CD8 TEM-2","pink_UCell"),mean_score("CD8 TEM-2","yellow_UCell")),
  `CD8 TEM-3`=c(mean_score("CD8 TEM-3","pink_UCell"),mean_score("CD8 TEM-3","yellow_UCell")),
  `CD8 TEM-4`=c(mean_score("CD8 TEM-4","pink_UCell"),mean_score("CD8 TEM-4","yellow_UCell")),
  `CD8 TEM-5`=c(mean_score("CD8 TEM-5","pink_UCell"),mean_score("CD8 TEM-5","yellow_UCell")),
  `CD8 Exhausted-1`=c(mean_score("CD8 Exhausted-1","pink_UCell"),mean_score("CD8 Exhausted-1","yellow_UCell")),
  `CD8 Exhausted-2`=c(mean_score("CD8 Exhausted-2","pink_UCell"),mean_score("CD8 Exhausted-2","yellow_UCell"))
)
rownames(signat.heatmap) <- c("Memory Stem Like","Terminal Effector Like")


# Figure 4a
# Historical CD8 module-based reclustering and DEG heatmap
Idents(tcells) <- "bulk"
cd16.tcells <- subset(tcells, idents="CD8")
Idents(cd16.tcells) <- "ident"

cd16.tcells <- NormalizeData(cd16.tcells)
cd16.tcells <- FindVariableFeatures(cd16.tcells)
cd16.tcells <- ScaleData(cd16.tcells, rownames(cd16.tcells))

wgcna_features <- unique(c(markers$yellow, markers$pink))
wgcna_features <- wgcna_features[wgcna_features %in% rownames(cd16.tcells)]

cd16.tcells <- RunPCA(
  cd16.tcells,
  features=wgcna_features
)
cd16.tcells <- RunHarmony(
  cd16.tcells,
  "Timepoint",
  plot_convergence=FALSE
)
cd16.tcells <- RunUMAP(
  cd16.tcells,
  reduction="harmony",
  dims=1:10
)
cd16.tcells <- FindNeighbors(
  cd16.tcells,
  reduction="harmony",
  dims=1:10
)
cd16.tcells <- FindClusters(
  cd16.tcells,
  resolution=0.5
)

Idents(cd16.tcells) <- "RNA_snn_res.0.5"

term.v.stem <- FindMarkers(
  cd16.tcells,
  ident.1=1,
  ident.2=c(2,3)
)
term.v.stem$gene <- rownames(term.v.stem)
term.v.stem.pval <- filter(term.v.stem, p_val_adj<=0.05)

unwanted_fig4 <- paste0(
  "^IGHV|^IGHJ|^IGHD|^IGKV|^IGLV|^TRBV|^TRBD|^TRBJ|",
  "^TRDV|^TRDD|^TRDJ|^TRAV|^TRAJ|^TRGV|^TRGJ|^TRD|",
  "^MT-|^RP|^MTRN"
)
remove_fig4 <- term.v.stem.pval$gene[
  grep(unwanted_fig4, term.v.stem.pval$gene)
]
term.v.stem.pval <- term.v.stem.pval[
  !term.v.stem.pval$gene %in% remove_fig4,
]

top.term.v.stem <- rbind(
  term.v.stem.pval %>% slice_max(n=55, order_by=avg_log2FC),
  term.v.stem.pval %>% slice_min(n=55, order_by=avg_log2FC)
)

Idents(cd16.tcells) <- "ident"

cd16.tcells.avg <- AverageExpression(
  cd16.tcells,
  verbose=TRUE,
  group.by="ident"
)
cd16.tcells.avg <- log1p(cd16.tcells.avg$RNA)

fig4_genes <- unique(top.term.v.stem$gene)
fig4_genes <- fig4_genes[
  fig4_genes %in% rownames(cd16.tcells.avg)
]
cd16.tcells.avg <- cd16.tcells.avg[
  fig4_genes,
  ,
  drop=FALSE
]
cd16.tcells.avg <- t(scale(t(cd16.tcells.avg)))

term.color <- colorRamp2(c(0.3,0.5), c("white","#fc4dcf"))
stem.color <- colorRamp2(c(0.2,0.4), c("white","#0070fc"))

ha_fig4 <- HeatmapAnnotation(
  stem=t(signat.heatmap[1,]),
  term=t(signat.heatmap[2,]),
  show_legend=TRUE,
  col=list(
    stem=stem.color,
    term=term.color
  )
)

fig4_column_order <- c(
  "CD8 Exhausted-1","CD8 TEM-2","CD8 TEM-4",
  "CD8 TEM-3","CD8 TEM-1","CD8 TEM-5","CD8 Exhausted-2"
)
fig4_column_order <- fig4_column_order[
  fig4_column_order %in% colnames(cd16.tcells.avg)
]

png(
  file.path(figure_dir, "Figure4a_memory_terminal_DEG_heatmap.png"),
  height=10, width=4.5, units="in", res=600
)
Heatmap(
  cd16.tcells.avg,
  col=colorRamp2(seq(-1,1.5,0.1), viridis(length(seq(-1,1.5,0.1)))),
  cluster_rows=TRUE,
  cluster_columns=FALSE,
  row_names_side="left",
  column_names_rot=90,
  column_names_side="top",
  column_names_gp=gpar(fontsize=8),
  row_names_gp=gpar(fontsize=6),
  row_labels=rownames(cd16.tcells.avg),
  row_dend_reorder=TRUE,
  show_row_dend=FALSE,
  column_dend_reorder=FALSE,
  column_dend_side="bottom",
  show_column_names=TRUE,
  cluster_row_slices=FALSE,
  top_annotation=ha_fig4,
  column_order=fig4_column_order
)
dev.off()


# Save WGCNA objects and session information
saveRDS(
  list(
    net=net,
    MEs=MEs,
    KMEs=KMEs,
    geneInfo=geneInfo,
    markers=markers,
    powerTable=powerTable
  ),
  file.path(result_dir, "WGCNA_results.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "sessionInfo_03_Tcell_analysis.txt")
)
