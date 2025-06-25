# for some details you can refer to the Seurat vignettes https://satijalab.org/seurat/

#### ONLY FIRST USED AND IF REQUIRED ####
install.packages("readr")
install.packages('Seurat')
install.packages('SeuratObject')


#### END OF ONLY FIRST USED AND IF REQUIRED ####

#### Load the data - to adapt according to what are your inputs
library(Seurat)
library(SeuratObject)


#### 1) if input = matrix (a table with all cells in column and genes in row)

library(readr)
Matrix <- as.data.frame(read_csv("~path/matrix.csv")) #load data, might take a bit of time
#set gene names as row names. To adapt according to your matrix
#if error occurs because of duplicate, wrap the argument in make.unique()
rownames(Matrix) <- Matrix[,which(colnames(Matrix)== "name of the column containing the gene names"] 
# delete the column containing the gene names as the information is now stored as row names
Matrix <- Matrix[,-which(colnames(Matrix)== "name of the column containing the gene names"]
# Potentially, you may have to remove some other annotation columns, follow the procedure done in the previous line.
# If you want to keep these information, create a data frame. # Gene.annotations <- Matrix[, c("column to keep 1", "column to keep2", etc.)]

#if you have some metadata available
Metadata <- as.data.frame(read_csv("~/path/_Metadata.csv")) #load the metadata, and transform them as data.frame
rownames(Metadata) <- Metadata$Cell_names #put as rownames the exact same names as used in the Matrix column names

# Create Seurat Object, 
seurat <- CreateSeuratObject(counts = raw_matrix, meta.data = Metadata, 
                             project = "project_name" # adding the project name is useful when merging seurat objects
                             min.features = 500, # keep only cells expressing  at least y different genes. Important,
                                        # to only retain cells that are informative, expressing sufficient genes.
                             min.cells = 3) # keep only genes expressed in at least x cells. You can adapt this value at your convenience

  ##go to the next section, skip 2 and 3)

#### 2) if input = 10X files
# you need to have in a folder 3 files: 1) barcodes.tsv.gz (list of cell barcodes), 
      # 2) features.tsv.gz (list of genes), and 3) matrix.mtx.gz (count data). Do not unzip your files!
seurat <- Read10X(data.dir = "path/folder containing the data/") # load the data

#### 3) if input = h5 files
seurat <- Read10X_h5("path/folder containing the data/file.h5", use.names = TRUE, unique.features = TRUE)


# for 2 and 3)
seurat <- CreateSeuratObject(counts = seurat, project = "project_name", # adding the project name is useful when merging seurat objects
                             min.cells = 3,  # keep only genes expressed in at least x cells. You can adapt this value at your convenience
                             min.features = 500) # keep only cells expressing  at least y different genes. Important,
      # to only retain cells that are informative, expressing sufficient genes. 500 might be low depending on the situation



#### Quality controls, have to be adapted to your situation
### 1) calculate the percentage of mitochondrial reads 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

### 2) excluding the cell doublets
#A good idea would be to look at your data 
hist(log10(seurat$nFeature_RNA), n =100) # look at the distribution of the number of different genes expressed per cell
summary(seurat$nFeature_RNA) #give some statistical characteristics of the distribution
hist(log10(seurat$nCount_RNA), n =100) # look at the distribution of the number of reads associated to each cell
summary(seurat$nCount_RNA) #give some statistical characteristics of the distribution

#If you see a tail in the distribution, it might be good to remove it. Also we tend to remove outliers (X% highest value)
seurat <- subset(x = seurat, subset = nCount_RNA > 50000) # number here is random (other option = X*median(seurat$nCount_RNA) )
seurat <- subset(x = seurat, subset = nFeature_RNA < 110000 )  #number here is random (other option = X*median(seurat$nFeature_RNA) )

### 3) excluding the cell with a too high percentage of mitochondrial genes
# the idea is that low-quality / dying cells often exhibit extensive mitochondrial contamination
# but this concept is tricky especially in the melanocyte/melanoma field as mitochondrial gene expression varies quite a lot according to the cell state
# according to the TCGA, a melanoma tumour has a median of 18.8% of its expression coming from the mitochondria (25-75% range = 13.7 to 26.9)
# and according to the CCLE, a melanoma cell have a median of 5.1% of its expression coming from the mitochondria (25-75% range = 4 to 7.2)
# This subseting might be optional or with a quite high treshold
seurat <- subset(x = seurat, subset = percent.mt < 20 )  #number here is random (other option = X*median(seurat$nFeature_RNA) )



#### Normalization and Clustering

### 1) Normalization (log transformation)
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

### 2) Find the variable features, 2000 is a quite accepted number but you can use whatever number you want
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(seurat) # you can visualize how the features are variating. It may help to set the treshold

### 3) scale the data (important to give each gene the same weight)
seurat <- ScaleData(seurat, features = rownames(seurat))

### 4) Run the PCA (linear dimensional reduction)
seurat <- RunPCA(seurat, 
                 features = VariableFeatures(object = seurat),
                 verbose = F)
DimPlot(seurat, reduction = "pca") + NoLegend() # if you want to have a look at the PCA

### 5) Identify the meaningfull dimensions
ElbowPlot(seurat, n =50) # multiple methods exist but the most used is the elbowplot. 
  # you keep dimension until you reach the elbow (i.e. the last PC before the curve start to look linear)
  # feel free to play a bit with the data

### 6) Cluster the cells
seurat <- FindNeighbors(seurat, dims = 1:30) # use the knn methods to establish the distance in between cells
seurat <- FindClusters(seurat, resolution = 0.8 ) # identify the cluster. Usually 0.4-1.2 are good resolution treshold but
    # play with it! 

### 7) Run the non-linear reduction to produce the UMAP
seurat <- RunUMAP(seurat, dims = 1:30)
DimPlot(seurat, reduction = "umap", # to visualize the umap
        label = T) + NoLegend() #you can keep the legend and/or remove the labels

# a good practice would be to look at the QC parameter per cluster 
# to check if some cluster might be linked to bad quality/duplicate
VlnPlot( object = seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)



#### Represent gene expression 
### 1) Gene expression per cluster as violin plot
VlnPlot( 
  object = seurat,
  features = c("Gene1", "Gene2"),
  ncol = 2 # to set the number of graph per row 
)

### 2) Gene expression per cluster as umap
FeaturePlot( 
  object = seurat,
  features = c("Gene1", "Gene2"),
  ncol = 2 # to set the number of graph per row 
)

#to see side by side the clustering and the gene expression
FeaturePlot( object = seurat,features = c("Gene1")) + (DimPlot(seurat, reduction = "umap",label = T) + NoLegend())



#### Find cluster markers
markers <- FindAllMarkers(seurat, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)






