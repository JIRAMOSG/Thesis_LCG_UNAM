library(Seurat)
library(patchwork)
library(ggplot2)
setwd("/home/johana/Downloads")
load("0- Lumaquin scRNAseq fish melanoma tumours - SeuratObject.rda")
getwd()


############################################################3#UMAP 
seurat.object<-scLumaquin
DimPlot(seurat.object, group.by = "tumor_type")
seurat.object <- RunUMAP(seurat.object, dims = 1:40, reduction = "pca")

library(Seurat)
library(ggplot2)


get_palette <- function(n){
  base <- c("#1b9e77","#d95f02","#7570b3","#e7298a",
            "#66a61e","#e6ab02","#a6761d","#666666",
            "#1f78b4","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
  if(n <= length(base)) base[1:n] else grDevices::colorRampPalette(base)(n)
}

ngroups <- length(unique(seurat.object@meta.data$tumor_type))


p_umap <- DimPlot(
  seurat.object,
  reduction = "umap",
  group.by = "tumor_type",
  label = TRUE,          # etiqueta el centroide de cada grupo
  repel = TRUE,          # que no se encimen las etiquetas
  raster = TRUE,         # puntos en raster: más rápido y sin “aliasing”
  pt.size = 1.5,         # tamaño de punto
  cols = get_palette(ngroups)
) +
  coord_fixed() +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ggtitle("UMAP by Tumor Type")

# Mostrar
p_umap

# Guardar 
ggsave("UMAP_tumor_type.png", p_umap, width = 7, height = 6, dpi = 600, bg = "white")
ggsave("UMAP_tumor_type.pdf", p_umap, width = 7, height = 6, useDingbats = FALSE)  # vectorial



############################################3 HEATMAP
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

obj <- scLumaquin

# Estados y colores (sin TME)
estado_orden <- c("melanocytic","stressed","proliferative","invasive","inflammatory")
Idents(obj) <- obj@active.ident
obj_sub <- subset(obj, idents = estado_orden)
obj_sub$state <- factor(Idents(obj_sub), levels = estado_orden)
Idents(obj_sub) <- obj_sub$state
cols_estado <- c(
  melanocytic   = "#7570b3",
  stressed      = "#66a61e",
  proliferative = "#e7298a",
  invasive      = "#d95f02",
  inflammatory  = "#1b9e77"
)

# Tus genes (y a qué estado esperas que pertenezcan)
genes_label <- c("tyrp1b","pmela","tyrp1a","ddit3","junba","hmgn2","cdk1","krt18a.1","pdgfbb","defbl1","apoda.1")
map_esperado <- c(
  tyrp1b="melanocytic", pmela="melanocytic", tyrp1a="melanocytic",
  ddit3="stressed", junba="stressed", hmgn2="proliferative",
  cdk1="proliferative",
  krt18a.1="invasive",
  pdgfbb="inflammatory", defbl1="inflammatory", apoda.1="inflammatory"
)

genes_label_present <- intersect(genes_label, rownames(obj_sub[["SCT"]]))
if (length(setdiff(genes_label, genes_label_present)))
  message("No están en SCT: ", paste(setdiff(genes_label, genes_label_present), collapse=", "))

# 1) MARCADORES con los parámetros "estándar"
N <- 15
markers <- FindAllMarkers(
  obj_sub, only.pos = FALSE, assay = "SCT",
  logfc.threshold = 0.1,   # ~1.19x en escala lineal
  min.pct = 0.01             # si esto te está filtrando demasiado, bájalo a 0.05 o 0
)
top15 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = N) %>%
  ungroup()

# --- Diagnóstico: ¿tus genes están en el top15 de SU estado esperado?
diag_tab <- tibble(
  gene = genes_label_present,
  estado_esperado = unname(map_esperado[genes_label_present])
) %>%
  left_join(top15 %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj),
            by = "gene") %>%
  mutate(
    en_top15_del_esperado = ifelse(!is.na(cluster) & cluster == estado_esperado, TRUE, FALSE)
  ) %>%
  group_by(estado_esperado) %>%
  mutate(rank_en_cluster = ifelse(en_top15_del_esperado,
                                  rank(-avg_log2FC, ties.method="first"),
                                  NA)) %>%
  ungroup()

print(diag_tab)  # Mira qué genes no cayeron en el top15 de su estado

# --- 2) Asegurar top-15 por estado, incluyendo tus genes si faltan
top15_por_estado <- lapply(estado_orden, function(st){
  g_top <- top15 %>% filter(cluster == st) %>% pull(gene)
  g_need <- genes_label_present[map_esperado[genes_label_present] == st]
  # si alguno de tus genes no está en el top15, forzamos su inclusión
  g_final <- unique(c(g_need, g_top))     # anteponer tus genes
  g_final[seq_len(min(length(g_final), N))]   # recorta a 15 manteniendo los tuyos
})
names(top15_por_estado) <- estado_orden

genes_all <- unique(unlist(top15_por_estado))

# --- Construir matriz y splits por estado del gen
obj_sub <- ScaleData(obj_sub, features = genes_all, assay = "SCT", verbose = FALSE)
M <- GetAssayData(obj_sub, assay = "SCT", layer = "scale.data")[genes_all, , drop = FALSE]

# Orden columnas por estado
ordc <- order(obj_sub$state)
M <- M[, ordc, drop = FALSE]
cell_state <- obj_sub$state[ordc]

# Row split: asigna cada gen a su estado (el bloque del que salió)
gene_state <- unlist(lapply(names(top15_por_estado), function(st){
  setNames(rep(st, length(top15_por_estado[[st]])), top15_por_estado[[st]])
}))
gene_state <- factor(gene_state[rownames(M)], levels = estado_orden)

# Orden de filas: dentro de cada estado, deja primero tus genes (en tu orden), luego los demás por logFC
ordered_genes <- unlist(lapply(estado_orden, function(st){
  g_block <- names(gene_state)[gene_state == st]
  g_lab   <- genes_label_present[map_esperado[genes_label_present] == st]
  g_lab   <- g_lab[g_lab %in% g_block]             # sólo los que están en este bloque
  g_rest  <- setdiff(g_block, g_lab)
  if (length(g_rest)) {
    ord_rest <- top15 %>% filter(cluster == st, gene %in% g_rest) %>%
      arrange(desc(avg_log2FC)) %>% pull(gene)
    c(g_lab, ord_rest)
  } else g_lab
}), use.names = FALSE)

M <- M[ordered_genes, , drop = FALSE]
gene_state <- gene_state[ordered_genes]

# Etiquetas SOLO de tus genes (derecha, itálicas, en tu orden original)
sel <- match(genes_label_present, rownames(M))
sel <- sel[!is.na(sel)]
labs <- genes_label_present[!is.na(match(genes_label_present, rownames(M)))]

mm_per_gene <- 1.5
H_px <- ceiling((nrow(M) * mm_per_gene / 25.4) * 300)

ht <- Heatmap(
  M,
  name = "scaled\nexpression",
  col = colorRamp2(c(-2,0,2), c("#2C5AA0","#F7F7F7","#C0392B")),
  heatmap_legend_param = list(at = c(-2,0,2)),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster = TRUE, rect_gp = gpar(col = NA),
  column_split = cell_state,
  row_split = gene_state,
  column_title_rot = 30,                     # títulos de los splits en diagonal
  column_title_gp = gpar(fontface = "bold"),
  top_annotation = HeatmapAnnotation(
    state = cell_state, col = list(state = cols_estado),
    annotation_name_side = "left"
  ),
  height = unit(nrow(M) * mm_per_gene, "mm")
)

row_anno <- rowAnnotation(
  labels = anno_mark(
    at = sel, labels = labs,
    side = "right",
    labels_gp = gpar(fontsize = 9, fontface = "italic"),
    link_width = unit(4, "mm"), link_height = unit(2, "mm"),
    extend = unit(10, "mm")
  ),
  width = unit(40, "mm")
)

png("heatmap_TOP15_con_tus_genes_y_titulos_diagonales.png", width = 3200, height = max(2400, H_px), res = 300)
draw(ht + row_anno, heatmap_legend_side = "left", annotation_legend_side = "right")
dev.off()
##########################FINAL

obj <- scLumaquin

# --- Estados y colores (sin TME) ---
estado_orden <- c("melanocytic","stressed","proliferative","invasive","inflammatory")
Idents(obj) <- obj@active.ident
obj_sub <- subset(obj, idents = estado_orden)
obj_sub$state <- factor(Idents(obj_sub), levels = estado_orden)
Idents(obj_sub) <- obj_sub$state

cols_estado <- c(
  melanocytic   = "#7570b3",
  stressed      = "#66a61e",
  proliferative = "#e7298a",
  invasive      = "#d95f02",
  inflammatory  = "#1b9e77"
)

# --- Top 15 por estado (log2FC>=0.50) ---
N <- 15
markers <- FindAllMarkers(obj_sub, only.pos = FALSE, assay = "SCT",
                          logfc.threshold = 1, min.pct = 0.10)

top15 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = N) %>%
  ungroup()

# --- Tus 11 genes y estado esperado ---
genes_label <- c("tyrp1b","pmela","tyrp1a","ddit3","junba","hmgn2",
                 "cdk1","krt18a.1","pdgfbb","apoda.1", "defbl1")
map_esperado <- c(
  tyrp1b="melanocytic", pmela="melanocytic", tyrp1a="melanocytic",
  ddit3="stressed", junba="stressed", hmgn2="proliferative",
  cdk1="proliferative",
  krt18a.1="invasive",
  pdgfbb="invasive", apoda.1="inflammatory", defbl1="inflammatory"
)
genes_label_present <- intersect(genes_label, rownames(obj_sub[["SCT"]]))

# Asegurar que estén dentro de los 15 por estado (si falta alguno, sustituye el #15)
top15_por_estado <- lapply(estado_orden, function(st){
  g_top  <- top15 %>% filter(cluster == st) %>% pull(gene)
  g_need <- genes_label_present[map_esperado[genes_label_present] == st]
  g_final <- unique(c(g_need, g_top))
  g_final[seq_len(min(length(g_final), N))]
})
names(top15_por_estado) <- estado_orden
genes_all <- unique(unlist(top15_por_estado))

# --- Matriz escalada y splits ---
obj_sub <- ScaleData(obj_sub, features = genes_all, assay = "SCT", verbose = FALSE)
M <- GetAssayData(obj_sub, assay = "SCT", layer = "scale.data")[genes_all, , drop = FALSE]

# Orden columnas por estado
ordc <- order(obj_sub$state)
M <- M[, ordc, drop = FALSE]
cell_state <- obj_sub$state[ordc]

# Estado de cada gen (para row_split) + orden de filas
gene_state <- unlist(lapply(names(top15_por_estado), function(st){
  setNames(rep(st, length(top15_por_estado[[st]])), top15_por_estado[[st]])
}))
gene_state <- factor(gene_state[rownames(M)], levels = estado_orden)

ordered_genes <- unlist(lapply(estado_orden, function(st){
  g_block <- names(gene_state)[gene_state == st]
  g_lab   <- genes_label_present[map_esperado[genes_label_present] == st]
  g_lab   <- g_lab[g_lab %in% g_block]
  g_rest  <- setdiff(g_block, g_lab)
  if (length(g_rest)) {
    ord_rest <- top15 %>% filter(cluster == st, gene %in% g_rest) %>%
      arrange(desc(avg_log2FC)) %>% pull(gene)
    c(g_lab, ord_rest)
  } else g_lab
}), use.names = FALSE)

M <- M[ordered_genes, , drop = FALSE]
gene_state <- gene_state[ordered_genes]

# --- Anotación superior SIN el texto "state" ---
top_ann <- HeatmapAnnotation(
  state = cell_state,
  col = list(state = cols_estado),
  annotation_name_gp = gpar(fontsize = 0)  # <-- oculta la palabra "state"
)

# --- Heatmap (títulos de columnas en diagonal; SIN títulos de filas/splits) ---
mm_per_gene <- 1.5
H_px <- ceiling((nrow(M) * mm_per_gene / 25.4) * 300)

ht <- Heatmap(
  M,
  name = "scaled\nexpression",
  col = colorRamp2(c(-2,0,2), c("#2C5AA0","#F7F7F7","#C0392B")),
  heatmap_legend_param = list(at = c(-2,0,2)),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster = TRUE, rect_gp = gpar(col = NA),
  column_split = cell_state,
  row_split = gene_state,
  row_title = NULL,                 # <-- quita nombres de estados a la izquierda
  column_title_rot = 30,            # títulos de columnas en diagonal
  column_title_gp = gpar(fontface = "bold"),
  top_annotation = top_ann,
  height = unit(nrow(M) * mm_per_gene, "mm")
)

# --- Etiquetas SOLO de tus 11 genes (derecha) ---
sel  <- match(genes_label_present, rownames(M))
labs <- genes_label_present
row_anno <- rowAnnotation(
  labels = anno_mark(
    at = sel, labels = labs,
    side = "right",
    labels_gp = gpar(fontsize = 9, fontface = "italic"),
    link_width = unit(4, "mm"), link_height = unit(2, "mm"),
    extend = unit(10, "mm")
  ),
  width = unit(40, "mm")
)

png("heatmap_top15_F1025.png",
    width = 3200, height = max(2400, H_px), res = 300)
draw(ht + row_anno, heatmap_legend_side = "left", annotation_legend_side = "right")
dev.off()

# --- Diagnóstico: ¿tus genes están en el top15 de SU estado esperado?
diag_tab <- tibble(
  gene = genes_label_present,
  estado_esperado = unname(map_esperado[genes_label_present])
) %>%
  left_join(top15 %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj),
            by = "gene") %>%
  mutate(
    en_top15_del_esperado = ifelse(!is.na(cluster) & cluster == estado_esperado, TRUE, FALSE)
  ) %>%
  group_by(estado_esperado) %>%
  mutate(rank_en_cluster = ifelse(en_top15_del_esperado,
                                  rank(-avg_log2FC, ties.method="first"),
                                  NA)) %>%
  ungroup()

print(diag_tab) 