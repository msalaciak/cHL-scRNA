# PBMC DEGs and figures

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# Set directories
project_dir <- "/home/rstudio/CIR_manuscript"
figure_dir <- file.path(project_dir, "figures")
result_dir <- file.path(project_dir, "results")

dir.create(figure_dir, showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)


# Load original PBMC object
pbmc <- readRDS(
  file.path(project_dir, "data", "original_PBMC_final.rds")
)


# Set cell type order
cell_types <- c(
  "CD4 T cells",
  "CD8 T cells",
  "NK cells",
  "B cells",
  "CD14 Monocytes",
  "CD14 Monocytes-2",
  "CD16 Monocytes",
  "cDCs",
  "pDCs"
)

pbmc$bulk <- factor(
  pbmc$bulk,
  levels = cell_types
)


# Figure colors
pbmc_colors <- c(
  "#66C2A5",
  "#E9936A",
  "#83798A",
  "#86A37A",
  "#C6B18B",
  "#C7D846",
  "#F8D348",
  "#DD97DE",
  "#805A24"
)


# Figure 1b - PBMC UMAP
p_umap <- DimPlot(
  pbmc,
  group.by = "bulk",
  cols = pbmc_colors,
  label = FALSE
) +
  theme_void() +
  NoLegend()

ggsave(
  file.path(figure_dir, "Figure1b_PBMC_UMAP.png"),
  p_umap,
  width = 3.5,
  height = 3.5,
  dpi = 600
)


# Figure 1c - canonical markers
pbmc_genes <- c(
  "CD3D",
  "CD8A",
  "IL7R",
  "CD40LG",
  "KLRD1",
  "NCAM1",
  "CD79A",
  "MS4A1",
  "CD14",
  "PF4",
  "FCGR3A",
  "CD1C",
  "CLEC10A",
  "LILRA4"
)

p_dot <- DotPlot(
  pbmc,
  features = pbmc_genes,
  group.by = "bulk",
  cols = c("lightgrey", "red")
) +
  coord_flip() +
  RotatedAxis()

ggsave(
  file.path(figure_dir, "Figure1c_PBMC_markers.png"),
  p_dot,
  width = 8,
  height = 4.5,
  dpi = 600
)


# Figure 1d - cell type proportions
Idents(pbmc) <- "bulk"

pbmc_prop <- table(
  Idents(pbmc),
  pbmc$Timepoint
)

pbmc_prop <- as.data.frame(pbmc_prop)

colnames(pbmc_prop) <- c(
  "Cell_type",
  "Timepoint",
  "Count"
)

pbmc_prop$Cell_type <- factor(
  pbmc_prop$Cell_type,
  levels = cell_types
)

p_prop <- ggplot(
  pbmc_prop,
  aes(
    x = Count,
    y = Timepoint,
    fill = Cell_type
  )
) +
  geom_col(
    position = position_fill(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = pbmc_colors
  ) +
  labs(
    x = "Proportion",
    y = "Time point",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "Figure1d_PBMC_proportions.png"),
  p_prop,
  width = 8,
  height = 4.5,
  dpi = 600
)


# Save cell proportions
pbmc_prop <- pbmc_prop %>%
  group_by(Timepoint) %>%
  mutate(
    Proportion = Count / sum(Count)
  )

write.csv(
  pbmc_prop,
  file.path(result_dir, "PBMC_celltype_proportions.csv"),
  row.names = FALSE
)


# Figure 1e - DEGs across timepoints
Idents(pbmc) <- "Timepoint"


# CD4 subclusters
cd4_1 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD4 T cells-1"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cd4_2 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD4 T cells-2"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)


# CD8 subclusters
cd8_1 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD8 T cells-1"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cd8_2 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD8 T cells-2"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cd8_3 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD8 T cells-3"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cd8_4 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD8 T cells-4"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)


# CD14 monocyte subclusters
cd14_1 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD14 Monocytes-1"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cd14_2 <- FindAllMarkers(
  subset(pbmc, subset = non.bulk.id == "CD14 Monocytes-2"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)


# Other PBMC populations
cd16 <- FindAllMarkers(
  subset(pbmc, subset = bulk == "CD16 Monocytes"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

bcells <- FindAllMarkers(
  subset(pbmc, subset = bulk == "B cells"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

nkcells <- FindAllMarkers(
  subset(pbmc, subset = bulk == "NK cells"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

cdcs <- FindAllMarkers(
  subset(pbmc, subset = bulk == "cDCs"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)

pdcs <- FindAllMarkers(
  subset(pbmc, subset = bulk == "pDCs"),
  only.pos = FALSE,
  min.diff.pct = -Inf,
  logfc.threshold = 0.10,
  min.pct = 0.05
)


# Combine CD4 results
cd4_degs <- bind_rows(
  cd4_1,
  cd4_2
) %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "CD4 T cells")


# Combine CD8 results
cd8_degs <- bind_rows(
  cd8_1,
  cd8_2,
  cd8_3,
  cd8_4
) %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "CD8 T cells")


# Combine CD14 results
cd14_degs <- bind_rows(
  cd14_1,
  cd14_2
) %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "CD14 Monocytes")


# Format remaining populations
cd16_degs <- cd16 %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "CD16 Monocytes")

bcell_degs <- bcells %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "B cells")

nk_degs <- nkcells %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "NK cells")

cdc_degs <- cdcs %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "cDCs")

pdc_degs <- pdcs %>%
  filter(p_val_adj <= 0.05) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cell_type = "pDCs")


# Combine all DEGs
pbmc_degs <- bind_rows(
  cd4_degs,
  cd8_degs,
  nk_degs,
  bcell_degs,
  cd14_degs,
  cd16_degs,
  cdc_degs,
  pdc_degs
)

pbmc_degs$Cell_type <- factor(
  pbmc_degs$Cell_type,
  levels = cell_types
)


# Save DEG table
write.csv(
  pbmc_degs,
  file.path(result_dir, "PBMC_DEGs_across_timepoints.csv"),
  row.names = FALSE
)


# Plot DEG fold changes
p_deg <- ggplot(
  pbmc_degs,
  aes(
    x = Cell_type,
    y = avg_log2FC,
    color = Cell_type
  )
) +
  geom_jitter(
    width = 0.45
  ) +
  geom_hline(
    yintercept = c(-1.5, 1.5),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = pbmc_colors
  ) +
  labs(
    x = NULL,
    y = "Average Log Fold Change"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "Figure1e_PBMC_DEGs.png"),
  p_deg,
  width = 6,
  height = 4.5,
  dpi = 600
)

# Supplemental Figure 2 - Azimuth validation

azimuth_file <- file.path(
  project_dir,
  "data",
  "azmiuth_jan232022.csv"
)

if (file.exists(azimuth_file)) {
  
  # Load Azimuth prediction scores
  azimuth <- read.csv(
    azimuth_file,
    row.names = 1,
    check.names = FALSE
  )
  
  # Scale prediction scores
  azimuth_matrix <- scale(
    azimuth,
    center = FALSE
  )
  
  azimuth_matrix[is.na(azimuth_matrix)] <- 0
  
  # Build heatmap
  ht <- Heatmap(
    azimuth_matrix,
    col = colorRamp2(
      seq(0, 3, 0.1),
      viridis(31)
    ),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_dend_side = "bottom"
  )
  
  # Display heatmap in RStudio
  draw(ht)
  
  # Save heatmap
  png(
    file.path(
      figure_dir,
      "Supplemental_Figure2_Azimuth.png"
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 600
  )
  
  draw(ht)
  
  dev.off()
  
} else {
  
  message(
    "Azimuth file not found: ",
    azimuth_file
  )
  
}

# Save session information
writeLines(
  capture.output(sessionInfo()),
  file.path(
    project_dir,
    "sessionInfo_PBMC_figures.txt"
  )
)