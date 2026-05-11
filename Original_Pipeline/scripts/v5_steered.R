# Script for Seurat-based scRNA-seq analysis and figure generation
# Last edited: 27/04/2026. Rohit
# Uses FPKM expression matrix and SRA metadata to perform PCA, clustering,
# marker gene identification, and generate Figures 3D, 3E and 3F

library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(ggdendro)

# load expression matrix
fpkm_data <- fread("~/Steered_Project/expression_matrix_full.csv", sep = ",")
fpkm_df <- as.data.frame(fpkm_data)

rownames(fpkm_df) <- make.unique(fpkm_df$gene_name)
fpkm_df <- fpkm_df[, !colnames(fpkm_df) %in% c("gene_id", "gene_name")]
fpkm_matrix <- as.matrix(fpkm_df)
mode(fpkm_matrix) <- "numeric"

cat("Genes:", nrow(fpkm_matrix), "\n")
cat("Cells:", ncol(fpkm_matrix), "\n")

# load metadata
sra_metadata <- read.csv("~/Steered_Project/SraRunTable.csv",
                         stringsAsFactors = FALSE)

cell_metadata <- sra_metadata[match(colnames(fpkm_matrix), sra_metadata$Run), ]

cat("Cells matched:", sum(!is.na(cell_metadata$Run)), "\n")

# log2 transform FPKM
log2_data <- log2(fpkm_matrix + 1)

# housekeeping gene check
if ("ACTB" %in% rownames(fpkm_matrix)) {
  cat("ACTB expressing cells:", sum(fpkm_matrix["ACTB", ] > 0), "\n")
}

if ("GAPDH" %in% rownames(fpkm_matrix)) {
  cat("GAPDH expressing cells:", sum(fpkm_matrix["GAPDH", ] > 0), "\n")
}

# filter genes
gene_var <- apply(log2_data, 1, var)
cells_expr <- apply(fpkm_matrix, 1, function(x) sum(x > 1))

variable_genes <- rownames(fpkm_matrix)[gene_var > 0.5 & cells_expr > 2]

cat("Variable genes:", length(variable_genes), "\n")

# create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = fpkm_matrix,
  project = "Camp2015",
  min.cells = 0,
  min.features = 0
)

# put log2 data into Seurat data slot
seurat_obj <- SetAssayData(
  object = seurat_obj,
  slot = "data",
  new.data = log2_data
)

# add metadata
seurat_obj$tissue <- cell_metadata$tissue
seurat_obj$stage <- cell_metadata$Stage
seurat_obj$gsm <- cell_metadata$Sample.Name
seurat_obj$source <- cell_metadata$source_name

# fetal vs organoid
seurat_obj$sample_type <- ifelse(
  grepl("Fetal", seurat_obj$tissue, ignore.case = TRUE),
  "fetal",
  "organoid"
)

print(table(seurat_obj$sample_type))

# region grouping for figure 3F
seurat_obj$region_group <- case_when(
  seurat_obj$gsm >= "GSM1957381" & seurat_obj$gsm <= "GSM1957428" ~ "r1",
  seurat_obj$gsm >= "GSM1957429" & seurat_obj$gsm <= "GSM1957476" ~ "r2",
  seurat_obj$gsm >= "GSM1957477" & seurat_obj$gsm <= "GSM1957536" ~ "r3",
  seurat_obj$gsm >= "GSM1957537" & seurat_obj$gsm <= "GSM1957572" ~ "r4",
  seurat_obj$sample_type == "fetal" ~ "fetal",
  TRUE ~ NA_character_)

print(table(seurat_obj$region_group, useNA = "ifany"))

# use organoid cells for figure 3D and 3E
seurat_organoid <- subset(seurat_obj, subset = sample_type == "organoid")

cat("Organoid cells:", ncol(seurat_organoid), "\n")

# recalculate variable genes for organoid cells
log2_organoid <- log2_data[, colnames(seurat_organoid)]

gene_var_org <- apply(log2_organoid, 1, var)
cells_expr_org <- apply(fpkm_matrix[, colnames(seurat_organoid)], 1,
                        function(x) sum(x > 1))

variable_genes_org <- rownames(fpkm_matrix)[gene_var_org > 0.5 & cells_expr_org > 2]

cat("Organoid variable genes:", length(variable_genes_org), "\n")

VariableFeatures(seurat_organoid) <- variable_genes_org

# PCA
seurat_organoid <- ScaleData(
  seurat_organoid,
  features = variable_genes_org,
  verbose = FALSE
)

seurat_organoid <- RunPCA(
  seurat_organoid,
  features = variable_genes_org,
  npcs = 50,
  verbose = FALSE
)

# JackStraw permutation test
seurat_organoid <- JackStraw(
  seurat_organoid,
  num.replicate = 200,
  prop.freq = 0.01,
  dims = 50,
  verbose = TRUE
)

seurat_organoid <- ScoreJackStraw(seurat_organoid, dims = 1:50)

js_data <- JS(seurat_organoid[["pca"]])
pvals <- js_data@overall.p.values

sig_pcs <- which(pvals[, 2] < 1e-20)

cat("Significant PCs at p < 1e-20:", length(sig_pcs), "\n")
print(sig_pcs)

# if the strict threshold gives too few PCs, use first 15 to keep biological structure
if (length(sig_pcs) < 10) {
  sig_pcs <- 1:15
  cat("Using first 15 PCs because strict threshold gave too few PCs\n")
}

# tSNE and clustering
set.seed(42)

seurat_organoid <- RunTSNE(
  seurat_organoid,
  dims = sig_pcs,
  perplexity = 5,
  seed.use = 42
)

seurat_organoid <- FindNeighbors(
  seurat_organoid,
  dims = sig_pcs,
  verbose = FALSE
)

# test resolutions to find one close to 11 clusters
cluster_test <- data.frame(resolution = numeric(), clusters = numeric())

for (res in seq(0.1, 3.0, by = 0.1)) {
  seurat_organoid <- FindClusters(
    seurat_organoid,
    resolution = res,
    verbose = FALSE
  )
  
  n_clusters <- length(unique(seurat_organoid$seurat_clusters))
  cluster_test <- rbind(cluster_test,
                        data.frame(resolution = res, clusters = n_clusters))
  
  cat("resolution:", res, "| clusters:", n_clusters, "\n")
}

write.csv(cluster_test, "resolution_cluster_test.csv", row.names = FALSE)

# choose final resolution
final_resolution <- cluster_test$resolution[which.min(abs(cluster_test$clusters - 11))]

cat("Selected resolution:", final_resolution, "\n")

seurat_organoid <- FindClusters(
  seurat_organoid,
  resolution = final_resolution,
  verbose = FALSE
)

cat("Final clusters:", length(unique(seurat_organoid$seurat_clusters)), "\n")

# marker detection using ROC
cluster_markers <- FindAllMarkers(
  seurat_organoid,
  test.use = "roc",
  only.pos = TRUE,
  min.pct = 0.25,
  verbose = FALSE
)

write.csv(cluster_markers, "roc_markers.csv", row.names = FALSE)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = myAUC, n = 10)

write.csv(top_markers, "top_roc_markers.csv", row.names = FALSE)

# marker validation for cluster identity
markers_check <- c("FOXG1", "OTX2", "NEUROD6", "RSPO2",
                   "DCN", "PAX6", "SOX2", "LIN28A",
                   "ASPM", "MYT1L", "WNT2B", "NFIA", "NFIB",
                   "GAD1", "GAD2", "VIM")

markers_check <- markers_check[markers_check %in% rownames(seurat_organoid)]

avg_expr <- AverageExpression(
  seurat_organoid,
  features = markers_check,
  slot = "data"
)$RNA

print(round(avg_expr, 2))
write.csv(round(avg_expr, 2), "cluster_marker_expression.csv")

# manually assign cluster labels after checking marker expression
# edit these labels if your marker table shows a different mapping
Idents(seurat_organoid) <- seurat_organoid$seurat_clusters

seurat_organoid <- RenameIdents(seurat_organoid,
                                "0"  = "c4 Dorsal forebrain Ns",
                                "1"  = "c3 Dorsal forebrain NPCs",
                                "2"  = "c2 Dorsal forebrain NPCs",
                                "3"  = "c10 Mesenchymal",
                                "4"  = "c1 Dorsal forebrain NPCs",
                                "5"  = "c5 Cycling NPCs",
                                "6"  = "c8 RSPO+",
                                "7"  = "c7 Ventral forebrain",
                                "8"  = "c9 Dorsal forebrain Ns",
                                "9"  = "c6 Dorsal forebrain Ns",
                                "10" = "c11 Dorsal forebrain Ns")



# dendrogram
seurat_organoid <- BuildClusterTree(seurat_organoid, verbose = FALSE)

seurat_organoid <- BuildClusterTree(seurat_organoid, verbose = FALSE)

cluster_tree <- Tool(seurat_organoid, slot = "BuildClusterTree")
hc <- as.hclust(cluster_tree)

p_tree <- ggdendrogram(hc, rotate = FALSE, theme_dendro = FALSE) +
  labs(title = "Phylogenetic tree of organoid cell clusters",
       x = "Cluster", y = "Height") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 8))

print(p_tree)
ggsave("Figure1_dendrogram.png", p_tree, width = 10, height = 6, dpi = 300)

# figure 3D

# the paper splits points by both time point and microdissected region
# so stage alone wasn't enough - made a new column to handle this
seurat_organoid$shape_group <- case_when(
  seurat_organoid$gsm >= "GSM1957381" &
    seurat_organoid$gsm <= "GSM1957428" ~ "r1 53d",
  seurat_organoid$gsm >= "GSM1957429" &
    seurat_organoid$gsm <= "GSM1957476" ~ "r2 53d",
  seurat_organoid$gsm >= "GSM1957477" &
    seurat_organoid$gsm <= "GSM1957536" ~ "r3 58d",
  seurat_organoid$gsm >= "GSM1957537" &
    seurat_organoid$gsm <= "GSM1957572" ~ "r4 58d",
  
  seurat_organoid$stage == "33 days" ~ "33d",
  seurat_organoid$stage == "35 days" ~ "35d",
  seurat_organoid$stage == "37 days" ~ "37d",
  seurat_organoid$stage == "41 days" ~ "41d",
  seurat_organoid$stage == "65 days" ~ "65d",
  
  TRUE ~ "other"
)

print(table(seurat_organoid$shape_group))

# save final object
saveRDS(seurat_organoid, "Camp2015_organoid_final.rds")

fig3d <- DimPlot(seurat_organoid,
                 reduction = "tsne",
                 label = TRUE,
                 label.size = 3,
                 pt.size = 1.5,
                 repel = TRUE,
                 shape.by = "shape_group") +
  scale_shape_manual(values = c(
    "33d" = 16, "35d" = 17, "37d" = 25,
    "65d" = 15, "41d" = 18,
    "r1 53d" = 0, "r2 53d" = 5,
    "r3 58d" = 1, "r4 58d" = 2)) +
  theme_classic() +
  labs(title = "Figure 3D. Organoid cell clusters",
       x = "tSNE 1", y = "tSNE 2",
       shape = "Stage / Region")

print(fig3d)
ggsave("Figure3D.png", fig3d, width = 12, height = 8, dpi = 300)

# figure 3E
genes_3e <- c("FOXG1", "OTX2", "RSPO2", "DCN",
              "ASPM", "LIN28A", "MYT1L", "NEUROD6")

genes_3e <- genes_3e[genes_3e %in% rownames(seurat_organoid)]

fig3e <- FeaturePlot(
  seurat_organoid,
  features = genes_3e,
  reduction = "tsne",
  cols = c("lightgrey", "red"),
  pt.size = 1,
  order = TRUE,
  ncol = 4,
  slot = "data"
) &
  theme_classic() &
  theme(plot.title = element_text(face = "italic", hjust = 0.5))

print(fig3e)
ggsave("Figure3E.png", fig3e, width = 16, height = 8, dpi = 300)

# figure 3F uses organoid regions plus fetal cells
seurat_3f <- subset(seurat_obj,
                    subset = region_group %in% c("r1","r2","r3","r4","fetal"))

print(paste("cells for fig 3f:", ncol(seurat_3f)))
print(table(seurat_3f$region_group))

genes_3f <- c("FOXG1", "NEUROD6", "OTX2")
genes_3f <- genes_3f[genes_3f %in% rownames(seurat_3f)]

Idents(seurat_3f) <- factor(seurat_3f$region_group,
                            levels = c("r1","r2","r3","r4","fetal"))

fig3f <- VlnPlot(seurat_3f,
                 features = genes_3f,
                 ncol = 3,
                 pt.size = 0.1,
                 layer = "data") &
  theme_classic() &
  labs(y = "log2(FPKM + 1)") &
  theme(legend.position = "none")

print(fig3f)
ggsave("Figure3F.png", fig3f, width = 10, height = 8, dpi = 300)

# save objects for website use
saveRDS(seurat_3f, "Camp2015_fig3f.rds")

