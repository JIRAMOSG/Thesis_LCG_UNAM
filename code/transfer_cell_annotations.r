library(Seurat)

# 1. SCTransform normalization and dimensionality reduction for both objects
scLumaquin <- SCTransform(seurat1)
scLumaquin <- RunPCA(seurat1)
scLumaquin <- RunUMAP(seurat1, dims = 1:30)

Weiss <- SCTransform(seurat2)
Weiss <- RunPCA(seurat2)
Weiss <- RunUMAP(seurat2, dims = 1:30)

# 2. Find transfer anchors between reference and query using SCT normalization
anchors <- FindTransferAnchors(
  reference = scLumaquin,
  query = Weiss,
  normalization.method = "SCT",
  dims = 1:25
)

# 3. Transfer cell state/cluster labels from seurat1 to seurat2
# Using the name of the column with your cluster or cell type annotation in seurat1@meta.data
predictions <- TransferData(
  anchorset = anchors,
  refdata = scLumaquin$tumor_type,
  dims = 1:25
)

# 4. Add the predicted labels as metadata to seurat2
seurat2 <- AddMetaData(Weiss, metadata = predictions)

# 5. The column seurat2$predicted.id now contains the transferred cell state labels
table(Weiss$predicted.id)
