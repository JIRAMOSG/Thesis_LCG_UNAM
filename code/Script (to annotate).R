fish_mel_sig <- list(
  Melanocytic = DR2HS_genes[DR2HS_genes$`Human_Symbol_DIOPT≥2` %in% subset(Tsoi_mel_sig,Signature == "Melanocytic")$Gene,]$Zebrafish_Symbol,
  Transitory = DR2HS_genes[DR2HS_genes$`Human_Symbol_DIOPT≥2` %in% subset(Tsoi_mel_sig,Signature == "Transitory")$Gene,]$Zebrafish_Symbol,
  NCSC = DR2HS_genes[DR2HS_genes$`Human_Symbol_DIOPT≥2` %in% subset(Tsoi_mel_sig,Signature == "Neural crest-like")$Gene,]$Zebrafish_Symbol,
  Undifferentiated = DR2HS_genes[DR2HS_genes$`Human_Symbol_Unique≥2` %in% subset(Tsoi_mel_sig,Signature == "Undifferentiated")$Gene,]$Zebrafish_Symbol
  )
DR2HS_genes$
library(escape)

scLumaquin <- runEscape(scLumaquin,
                        gene.sets = fish_mel_sig,
                        method = "ssGSEA",
                        groups = 1000,
                        min.size = 5, 
                        normalize = T)

FeaturePlot(scLumaquin, features = c("Melanocytic", "Transitory", "NCSC","Undifferentiated"))


MelLum <- subset(scLumaquin, subset = tumor_type %in% c("melanocytic","invasive","stressed","proliferative" ))


MelLum <- SCTransform(MelLum)
MelLum <- RunPCA(MelLum, npcs = 100)
ElbowPlot(MelLum, n = 50)



#non linear reduction using the UMAP 
MelLum <- RunUMAP(MelLum, reduction = "pca" ,dims = 1:30, verbose = T,reduction.name = "umap")

#build the map by defining the cell proximity. Using the nearest-neighbour method
MelLum <- FindNeighbors(MelLum, dims = 1:30, verbose = T,
                                       reduction = "pca")

#define the cluster. Start with a low resolution. Use the Louvain algorithm
MelLum <- FindClusters(MelLum, verbose = T, resolution = .3)

#plot the UMAP
DimPlot(MelLum, label = T, label.size = 8,
        reduction = "umap", group.by = "tumor_type",
        pt.size = 1.5)+NoLegend()+ggtitle("Cluster Annotation")+
#plot the UMAP
DimPlot(MelLum, label = T, label.size = 8,
        reduction = "umap", 
        pt.size = 1.5)+NoLegend()+ggtitle("Cluster Annotation")


Markers <- FindAllMarkers(MelLum, only.pos = T)
MarkersROC <- FindAllMarkers(MelLum, only.pos = T, test.use = "roc")

FeaturePlot(MelLum, features = "hsBRAF")
VlnPlot(MelLum, features = "hsBRAF")



NewIDs <- c(
  "Patient specific","Melanocytic","Mitotic","Mitochondrial","Neural Crest like","Antigen Presenting / Stressed", "Mesenchymal")
names(NewIDs) <- levels(MelLum)
MelLum <- RenameIdents(MelLum, NewIDs)
DimPlot(MelLum, reduction = "umap", label = TRUE, pt.size = 2, label.size = 6) + NoLegend()

MelLum$Cell_State <- Idents(MelLum)




Idents(scLumaquin) <- scLumaquin$seurat_clusters
DimPlot(scLumaquin, label = T, label.size = 10)
Marks <- FindAllMarkers(scLumaquin)

2 


NewIDs <- c(
  "Melanoma","Irido(leuko?)phores","Macrophages","Melanoma","Fibroblasts","Melanoma",
  "Melanoma","Suprabasal Cells","Basal Cells","Melanoma")
names(NewIDs) <- levels(scLumaquin)
scLumaquin <- RenameIdents(scLumaquin, NewIDs)
DimPlot(scLumaquin, reduction = "umap", label = TRUE, pt.size = 2, label.size = 4) + NoLegend()

scLumaquin$Cell_Type <- Idents(scLumaquin)

CellStates <- MelLum[["Cell_State"]]
SeuratID <- scLumaquin[["Cell_Type"]]
SeuratID$Cell_State <- as.character(SeuratID$Cell_Type)
SeuratID[rownames(subset(CellStates,Cell_State == "Melanocytic")),2] <-  "Melanocytic"
SeuratID[rownames(subset(CellStates,Cell_State == "Antigen Presenting / Stressed")),2] <-  "Antigen Presenting / Stressed"
SeuratID[rownames(subset(CellStates,Cell_State == "Neural Crest like")),2] <-  "Neural Crest like"
SeuratID[rownames(subset(CellStates,Cell_State == "Mesenchymal")),2] <-  "Mesenchymal"
SeuratID[rownames(subset(CellStates,Cell_State == "Mitotic")),2] <-  "Mitotic"
SeuratID[rownames(subset(CellStates,Cell_State == "Patient specific")),2] <-  "Patient specific"
SeuratID[rownames(subset(CellStates,Cell_State == "Mitochondrial")),2] <-  "Mitochondrial"

scLumaquin <- AddMetaData(scLumaquin, metadata = SeuratID)

Idents(scLumaquin) <- scLumaquin$Cell_State

DimPlot(scLumaquin, label = T)+NoLegend()


VlnPlot(scLumaquin, features = "mitfa")


VlnPlot(Srt_Skin, features = "ENSDARG00000012405", group.by = "cell_type_broad")+NoLegend()

Idents(Srt_Skin) <- Srt_Skin$cell_type_broad
ma <- FindAllMarkers(Srt_Skin)

save(scLumaquin, file = "0- Lumaquin scRNAseq fish melanoma tumours - SeuratObject.rda")
setwd("Lumaquin/scRNA-seq/")

library(Seurat)
library(ShinyCell)


scConf = createConfig(scLumaquin)
makeShinyApp(scLumaquin, scConf, gene.mapping = F,
             shiny.title = "Lumaquin Melanoma dataset",  
             shiny.footnotes = "Lumaquin et al. Nature Communications 2023",
             gex.assay = "SCT",
             shiny.dir = "shinyAppLumaq",
             default.gene1 = "mitfa",
             default.gene2 = "ubb") 

