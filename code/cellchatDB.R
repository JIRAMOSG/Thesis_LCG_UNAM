# Load necessary libraries
library(CellChat)
library(patchwork)
library(dplyr)

# Load the Seurat object
load("~/Lumaquin/scRNA-seq/0- Lumaquin scRNAseq fish melanoma tumours - SeuratObject.rda")

# Extract labels and metadata
labels <- Idents(scLumaquin)
meta <- data.frame(labels = labels, row.names = names(labels))

# Access normalized data matrix from the 'SCT' assay
data.input <- scLumaquin[["SCT"]]@data

# Create a CellChat object
cellchat <- createCellChat(object = scLumaquin, group.by = "ident", assay = "SCT")

# Load the CellChat database for zebrafish
CellChatDB <- CellChatDB.zebrafish

# Show the categories in the database
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# Subset the database to use only "Secreted Signaling" interactions
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

# Assign the subsetted database to the CellChat object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)

# Subset the expression data of signaling genes to save computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

# Print a summary of the CellChat object to verify
print(cellchat)

