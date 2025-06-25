# Instalar y cargar Monocle3
#if (!requireNamespace("monocle3", quietly = TRUE)) {
#  install.packages("devtools")
#  devtools::install_github("cole-trapnell-lab/monocle3")
#}

library(monocle3)
library(Seurat)
library(Matrix)

# Leer los datos de conteo
count_data <- read.csv("C:/Users/jramosgalguera/Downloads/GSE201378_scRNAseq_counts.csv", row.names = 1)
cell_metadata <- read.csv("C:/Users/jramosgalguera/Downloads/GSE201378_scRNAseq_metadata.csv", row.names = 1)
rownames(cell_metadata)<- gsub('-', '.', rownames(cell_metadata))
# Crear un objeto Seurat
# Verificar que los datos están en formato de matriz
count_matrix <- as.matrix(count_data)

# Crear un objeto CellDataSet
# Se necesitan datos de células y características (genes)
#cell_metadata <- data.frame(row.names = colnames(count_matrix))  # Metadata de células
gene_metadata <- data.frame(gene_short_name = rownames(count_matrix))  # Metadata de genes
rownames(gene_metadata) <- rownames(count_matrix)

cds <- new_cell_data_set(
  expression_data = count_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
# Preprocesar los datos
cds <- preprocess_cds(cds, num_dim = 50)
# 
print(cds)
# Verificar la estructura del objeto cds después de reducir la dimensionalidad
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- reduce_dimension(cds, preprocess_method = "UMAP")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
plot_cells(cds)
#Verificar la estructura del objeto cds
print(cds)
# Aprender el grafo de células para pseudotiempo
cds <- learn_graph(cds)
# # Visualizar las trayectorias
plot_cells(cds, color_cells_by = "pseudotime")
# plot_cells(cds, color_cells_by = "cluster")
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
          label_branch_points=FALSE)
cds <- order_cells(cds)



plot_cells(cds,
           color_cells_by = "tumor_type",
           label_groups_by_cluster = F, 
           group_label_size = 4,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=7)
