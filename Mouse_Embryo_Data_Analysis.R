### ========================
### 1st part: Global setting ----
### ========================
setwd("/home/yhw/bioinfo/project-mine/MultiOmics")
load("R/CodeData/Multiome_Data_Analysis_Mouse.RData")
# - species
species <- "mouse"
# - library path
pkg.lib <- "/home/yhw/software/anaconda3/envs/rstudio-4.2.1/lib/R/library"
.libPaths(pkg.lib)
# - human and mouse homologous genes
hs.mm.hm <- read.table("/home/yhw/document/ensembl/Homologous_Genes/release108_GRCh38_GRCm38_mart_export.txt", header = T, stringsAsFactors = F, sep = "\t")
hs.mm.hm <- subset(hs.mm.hm, Mouse.homology.type == "ortholog_one2one")
hs.mm.hm <- hs.mm.hm[hs.mm.hm$Gene.name != "",]
hs.mm.hm <- hs.mm.hm[!duplicated(hs.mm.hm$Gene.name),]
rownames(hs.mm.hm) <- hs.mm.hm$Mouse.gene.name
# - human marker genes
marker <- list()
marker$super.cluster <- c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3","TP63","PERP","KRT18","FOXD3","MPZ","PLP1","S100B","ALX1","ALX3",
                          "DLX1","DLX2","SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9","PAX8","WT1","TBX5","PITX1","FOXF1","GATA6",
                          "GATA5","NKX2-5","PLVAP","CD53","PLAC8","CD68","CORO1A","FOXA2","CLDN6","COL1A1","COL1A2","COL3A1","POSTN","DCN","PAX2","OSR1")
# - dotplot theme setting
theme_dp <- theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
                  axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
                  axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
                  axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
pd.col <- c("#A6CEE3","#1F78B4","#2CA9FD",
            "#B2DF8A","#33A02C","#1DF310",
            "#FB9A99","#A41617","#E31A1C",
            "#FDBF6F","#B15928","#FF7F00",
            "#CAB2D6","#6A3D9A","#8717FF",
            "#FFFF99","#C4C42D","#F3F310",
            "#7DE4E6","#00FBFF","#00BCBF",
            "#D277FD","#7D23A7","#AE00FF",
            "#7789FD","#2234A8","#203EFF",
            "#FC8ACD","#B01E74","#FF20A4",
            "#8AFCE9","#21AD96","#1CFFD9",
            "#97FF99","#F912FD",
            "#BFBFBF","#7B7B7B")
pd.col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
            "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
            "#00FBFF","#00BCBF","#D277FD","#AE00FF","#7789FD","#203EFF",
            "#FC8ACD","#FF20A4","#8AFCE9","#1CFFD9","#97FF99","#1CFF20",
            "#7B7B7B","#BFBFBF")


### =================
### 2nd part: Library ----
### =================
library(devtools)
library(BiocManager)

### >>> 1. CRAN packages
cran.pks <- c("Seurat", "harmony", "rliger", "stringr", "circlize", "ggplot2", "ggalluvial", "RCurl",
              "ggrepel", "forcats", "dplyr", "tidyr", "tidyverse", "stringr", "Signac")
for (pks in cran.pks) { library(pks, character.only = T) }

### >>> 2. Bioconductor packages
bioc.pks <- c("BSgenome.Hsapiens.UCSC.hg38", "SingleR", "scater")
for (pks in bioc.pks) { library(pks, character.only = T) }

### >>> 3. GitHub packages
if (!require("SeuratDisk", character.only = TRUE)) {
  devtools::install_github('mojaveazure/seurat-disk')
  library("SeuratDisk", character.only = T)
} else {
  library("SeuratDisk", character.only = T)
}



### ============================================
### 3nd part: Process mouse embryo scRNAseq data ----
### ============================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Mm_embryo")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load raw data
# - human seurat object
hs.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.new.final.rds")
hs.embryo$SuperCluster <- str_split_fixed(hs.embryo$CellType, ":", 2)[,1]
hs.embryo.sub <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.merge2.sub1.rds")
# - mouse data
mm.embryo <- readRDS("/home/yhw/bioinfo/project-temp/jobs/proj6th/GSE119945_RAW_Mm/gene_count_cleaned.RDS")
cell.anno <- read.csv("/home/yhw/bioinfo/project-temp/jobs/proj6th/GSE119945_RAW_Mm/cell_annotate.csv", row.names = 1)
if (all(colnames(mm.embryo) %in% rownames(cell.anno))) {
  cell.anno <- cell.anno[colnames(mm.embryo),]
}
gene.anno <- read.csv("/home/yhw/bioinfo/project-temp/jobs/proj6th/GSE119945_RAW_Mm/gene_annotate.csv", row.names = 1)
if (all(rownames(mm.embryo) == rownames(gene.anno))) {
  rownames(mm.embryo) <- gene.anno$gene_short_name
}
mm.embryo <- mm.embryo[intersect(rownames(mm.embryo),hs.mm.hm$Mouse.gene.name),]
hs.mm.hm <- hs.mm.hm[rownames(mm.embryo),]
rownames(mm.embryo) <- hs.mm.hm$Gene.name
tail(rownames(mm.embryo))
tail(hs.mm.hm$Gene.name)


### >>> 3. Add metadata information and filter data (E9.5-E12.5)
# create seurat objective
if (all(colnames(mm.embryo) == rownames(cell.anno))) {
  mm.embryo <- CreateSeuratObject(counts = mm.embryo, meta.data = cell.anno, project = "Mm.embryo")
  mm.embryo@meta.data$development_stage <- paste0("E", mm.embryo@meta.data$development_stage)
}
mm.embryo[["percent.mt"]] <- PercentageFeatureSet(mm.embryo, pattern = "^mt-")
pdf(file.path(res.out, "Quality_control_before_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = mm.embryo, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3, split.by = "development_stage", pt.size = 0, group.by = "development_stage")
dev.off()
# filter cells with low quality
mm.embryo <- subset(
  x = mm.embryo,
  subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000 & percent.mt <= 30
)
# filter development stage
mm.embryo <- subset(mm.embryo, development_stage %in% c("E9.5", "E10.5", "E11.5", "E12.5"))
#keep.cell <- sample(1:ncol(mm.embryo), 2*ncol(hs.embryo))
#mm.embryo <- mm.embryo[,keep.cell]
pdf(file.path(res.out, "Quality_control_after_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = mm.embryo, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3, split.by = "development_stage", pt.size = 0, group.by = "development_stage")
dev.off()


### >>> 4. Dimension reduction and clustering
# - post-processing
DefaultAssay(mm.embryo) <- "RNA"
mm.embryo <- NormalizeData(mm.embryo) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mm.embryo <- ScaleData(mm.embryo, verbose = FALSE, features = rownames(mm.embryo))
mm.embryo <- RunPCA(mm.embryo, npcs = 50, verbose = FALSE, features = VariableFeatures(mm.embryo))
ElbowPlot(mm.embryo, ndims = 30)
dim.n <- 20
# - dimension reduction and clustering
mm.embryo <- FindNeighbors(mm.embryo, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo <- RunUMAP(mm.embryo, reduction = "pca", dims = 1:dim.n, return.model = TRUE)
# - visualize clustering
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_uncorrected_clustering.pdf"), height = 5, width = 12.5)
DimPlot(mm.embryo, reduction = "umap", group.by = c("development_stage", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# - remove batch effect
library(harmony)
mm.embryo <- RunHarmony(mm.embryo, group.by.vars = "development_stage")
mm.embryo <- RunUMAP(mm.embryo, reduction = "harmony", dims = 1:dim.n, return.model = TRUE)
mm.embryo <- FindNeighbors(mm.embryo, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_corrected_clustering.pdf"), height = 5, width = 12.5)
DimPlot(mm.embryo, reduction = "umap", group.by = c("development_stage", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.25)
dev.off()
# - find marker genes
#mm.embryo.markers <- FindAllMarkers(mm.embryo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


### >>> 5. Annotating cell types by Seurat (R)
# - filtering non-homologous genes
hs.mm.hm.tar <- unique(hs.mm.hm[,c(1,8)]) %>%
  dplyr::filter(Gene.name %in% rownames(hs.embryo)) %>%
  dplyr::filter(Gene.name %in% rownames(mm.embryo))
sr.predict <- list()
sr.predict$hs <- hs.embryo[rownames(hs.embryo) %in% hs.mm.hm.tar$Gene.name,]
sr.predict$mm <- mm.embryo[rownames(mm.embryo) %in% hs.mm.hm.tar$Gene.name,]
# - annotate cells using human dataset as reference (Seurat 4.0.6)
anchors <- FindTransferAnchors(reference = sr.predict$hs, query = sr.predict$mm, dims = 1:30, reference.reduction = "pca")
sr.predict$map <- MapQuery(anchorset = anchors, reference = sr.predict$hs, query = sr.predict$mm,
                           refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
sr.predict$hs$SuperCluster <- str_split_fixed(sr.predict$hs$CellType, ":", 2)[,1]
sr.predict$map$CellType.Seurat <- sr.predict$map$predicted.celltype
sr.predict$map$SuperCluster.Seurat <- str_split_fixed(sr.predict$map$predicted.celltype, ":", 2)[,1]
pd.col.fix <- pd.col[1:length(unique(sr.predict$hs$SuperCluster))]
names(pd.col.fix) <- sort(unique(sr.predict$hs$SuperCluster))
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map, reduction = "ref.umap", group.by = "SuperCluster.Seurat", label = F, pt.size = 0.15) +
  ggtitle("Mouse embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_seurat.pdf"), height = 5, width = 15)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, anchors)


### >>> 6. Annotating cell types by SingleR (R)
# - transform Seurat object into SingleCellExperiment object
hs.sce <- as.SingleCellExperiment(sr.predict$hs)
mm.sce <- as.SingleCellExperiment(sr.predict$mm)
# - cell type annotation
library(SingleR)
predictions.sgr <- SingleR(test=mm.sce, assay.type.test=1, ref=hs.sce, labels=hs.sce$CellType, num.threads=16)
# - visualization
if (all(colnames(sr.predict$mm) == rownames(predictions.sgr)) & all(colnames(sr.predict$map) == rownames(predictions.sgr))) {
  # pig object
  sr.predict$mm$CellType.SingleR <- predictions.sgr$labels
  sr.predict$mm$SuperCluster.SingleR <- str_split_fixed(sr.predict$mm$CellType.SingleR, ":", 2)[,1]
  # mapped object
  sr.predict$map$CellType.SingleR <- predictions.sgr$labels
  sr.predict$map$SuperCluster.SingleR <- str_split_fixed(sr.predict$map$CellType.SingleR, ":", 2)[,1]
}
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map, reduction = "ref.umap", group.by = "SuperCluster.SingleR", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_singleR.pdf"), height = 5, width = 15)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, hs.sce, mm.sce)


### >>> 7. Annotating cell types by SciBet (R)
library("Rcpp")
library("RcppEigen")
library("ggsci")
library("viridis")
library("tidyverse")
library("scibet")
# - prepare data
ref.set <- GetAssayData(hs.embryo, slot = "data") %>% t() %>% as.data.frame()
ref.set$label <- hs.embryo$CellType
query.set <- GetAssayData(mm.embryo, slot = "data") %>% t() %>% as.data.frame()
# - predict cell type
predictions.scibet.5k <- SciBet(ref.set, query.set, k = 5000)
predictions.scibet.4k <- SciBet(ref.set, query.set, k = 4000)
predictions.scibet.3k <- SciBet(ref.set, query.set, k = 3000)
predictions.scibet.2k <- SciBet(ref.set, query.set, k = 2000)
predictions.scibet.1k <- SciBet(ref.set, query.set, k = 1000)
# - plot
sr.predict$map$CellType.SciBet <- predictions.scibet.1k
sr.predict$map$SuperCluster.SciBet <- str_split_fixed(sr.predict$map$CellType.SciBet, ":", 2)[,1]
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map, reduction = "ref.umap", group.by = "SuperCluster.SciBet", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_SciBet.pdf"), height = 5, width = 15)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, ref.set, query.set)


### >>> 9. Compare the annotation between SingleR, Seurat and sciBet
# majority-rule
cell.type <- sr.predict$map@meta.data[,(ncol(sr.predict$map@meta.data)-5):ncol(sr.predict$map@meta.data)]
cell.type$CellType.Stat <- apply(cell.type[,seq(1,5,2)], 1, function(x){length(sort(table(as.character(x)),decreasing=TRUE))})
cell.type$CellType.Final <- apply(cell.type[,seq(1,5,2)], 1, function(x){names(sort(table(as.character(x)),decreasing=TRUE)[1])})
cell.type$CellType.Final[cell.type$CellType.Stat == 3] <- cell.type$CellType.Seurat[cell.type$CellType.Stat == 3]
if (identical(colnames(sr.predict$map), rownames(cell.type)) & identical(colnames(mm.embryo), rownames(cell.type))) {
  sr.predict$map$CellType.Final <- cell.type$CellType.Final
  sr.predict$map$SuperCluster.Final <- str_split_fixed(sr.predict$map$CellType.Final, ":", 2)[,1]
  mm.embryo$CellType <- cell.type$CellType.Final
  mm.embryo$SuperCluster <- str_split_fixed(mm.embryo$CellType, ":", 2)[,1]
}
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map[,sample(1:ncol(sr.predict$map),100000)], reduction = "ref.umap",
              group.by = "SuperCluster.Final", label = F, pt.size = 0.15) +
  ggtitle("Mouse embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_aggregation.pdf"), height = 5, width = 15)
p1+p2
dev.off()
png(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_aggregation.png"),
    height = 5, width = 15, res = 1000, units = "in")
p1+p2
dev.off()
# - proportion of cell types
# human
hs.embryo@meta.data <- hs.embryo@meta.data %>%
  mutate(SuperCluster=str_split_fixed(CellType, ":", 2)[,1],
         CellType.sub=str_split_fixed(CellType, ":", 2)[,2])
pd <- hs.embryo@meta.data[,c("SuperCluster","CellType.sub")]
library(webr)
pdf(file.path(res.out, "PieDonut_plot_to_show_human_cell_type_proportion.pdf"), height = 6, width = 6)
PieDonut(pd, aes(pies = SuperCluster), addDonutLabel = TRUE, showRatioDonut = T,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2, r0 = 0.5) +
  scale_fill_manual(values = pd.col.fix)
dev.off()
pd <- table(pd$SuperCluster) %>% as.data.frame() %>% mutate(Species = "human")
pd.col.fix <- pd.col[1:length(unique(sr.predict$hs$SuperCluster))]
names(pd.col.fix) <- sort(unique(sr.predict$hs$SuperCluster))
pdf(file.path(res.out, "Barplot_plot_to_show_human_cell_type_proportion.pdf"), height = 10, width = 5)
ggplot(pd, aes(x = Species, y = Freq, fill = Var1)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = Var1), position = position_stack(vjust= 0.5), colour = "white", size = 3) +
  scale_fill_manual(values = pd.col.fix) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
dev.off()
# mouse
pd <- mm.embryo@meta.data[,c("SuperCluster","CellType")]
library(webr)
pdf(file.path(res.out, "PieDonut_plot_to_show_mouse_cell_type_proportion.pdf"), height = 6, width = 6)
PieDonut(pd, aes(pies = SuperCluster), addDonutLabel = TRUE, showRatioDonut = T,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2, r0 = 0.5) +
  scale_fill_manual(values = pd.col.fix)
dev.off()
pd <- table(pd$SuperCluster) %>% as.data.frame() %>% mutate(Species = "human")
pd.col.fix <- pd.col[1:length(unique(sr.predict$hs$SuperCluster))]
names(pd.col.fix) <- sort(unique(sr.predict$hs$SuperCluster))
pdf(file.path(res.out, "Barplot_plot_to_show_mouse_cell_type_proportion.pdf"), height = 10, width = 5)
ggplot(pd, aes(x = Species, y = Freq, fill = Var1)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = Var1), position = position_stack(vjust= 0.5), colour = "white", size = 3) +
  scale_fill_manual(values = pd.col.fix) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
dev.off()
# - heatmap
library(viridis)
predictions <- table(sr.predict$map$CellType.Final, sr.predict$map$CellType.SciBet)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (SciBet)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
predictions <- table(sr.predict$map$CellType.Final, sr.predict$map$CellType.Seurat)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p2 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (Seurat)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
predictions <- table(sr.predict$map$CellType.Final, sr.predict$map$CellType.SingleR)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p3 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (SingleR)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
pdf(file.path(res.out, "comparing_of_cell_type_predicted_by_Seurat_SingleR_SciBet_heatmap.pdf"), width = 60, height = 19)
p1 + p2 + p3
dev.off(); rm(p1, p2, p3)
# output mouse embryo data
if (identical(colnames(mm.embryo), colnames(sr.predict$map))) {
  mm.embryo@meta.data <- cbind(mm.embryo@meta.data, sr.predict$map@meta.data[,(ncol(sr.predict$map@meta.data)-7):ncol(sr.predict$map@meta.data)])
}
# - create new mouse data with mouse gene names
mm.embryo.raw <- readRDS("/home/yhw/bioinfo/project-temp/jobs/proj6th/GSE119945_RAW_Mm/gene_count_cleaned.RDS")
if (all(colnames(mm.embryo.raw) %in% rownames(cell.anno))) {
  cell.anno <- cell.anno[colnames(mm.embryo.raw),]
}
if (all(rownames(mm.embryo.raw) == rownames(gene.anno))) {
  rownames(mm.embryo.raw) <- gene.anno$gene_short_name
}
if (all(colnames(mm.embryo.raw) == rownames(cell.anno))) {
  mm.embryo.raw <- CreateSeuratObject(counts = mm.embryo.raw, meta.data = cell.anno, project = "mm.embryo.raw")
  mm.embryo.raw@meta.data$development_stage <- paste0("E", mm.embryo.raw@meta.data$development_stage)
}
mm.embryo.raw <- mm.embryo.raw[,colnames(mm.embryo)]
if (all(colnames(mm.embryo.raw) == colnames(mm.embryo))) {
  mm.embryo.raw@meta.data <- mm.embryo@meta.data
}
DefaultAssay(mm.embryo.raw) <- "RNA"
mm.embryo.raw <- NormalizeData(mm.embryo.raw) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mm.embryo.raw$CellType.sub <- str_split_fixed(mm.embryo$CellType, ":", "2")[,2]
#
Idents(mm.embryo) <- mm.embryo$SuperCluster
pd <- subset(mm.embryo, SuperCluster %in% unique(setdiff(mm.embryo$SuperCluster, c("Undefined", "Fibroblast"))))
DefaultAssay(pd) <- "RNA"
pd <- NormalizeData(pd) %>% ScaleData()
pd@meta.data %>%
  mutate(SuperCluster=factor(SuperCluster,
                             levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", "Craniofacial",
                                        "Head mesoderm", "Somite", "Intermediate mesoderm", "Limb","Lateral plate mesoderm",
                                        "Heart", "Endothelium", "Blood", "Endoderm"))) -> pd@meta.data
pdf(file.path(res.out, "Dotplot_to_show_mouse_super_cluster_markers_grouped_by_cell_type_rna.pdf"),
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39], scale.max = 40,
        assay = "RNA", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()
pd <- SeuratWrappers::RunALRA(pd)
pd <- ScaleData(pd, features = rownames(pd))
pdf(file.path(res.out, "Dotplot_to_show_mouse_super_cluster_markers_grouped_by_cell_type_alra.pdf"),
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39], scale.max = 100,
        assay = "alra", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()



### =====================================================
### 5th step: Re-clustering to check cell type annotation ----
### =====================================================

### >>> 1. re-clustering
res.out <- file.path(getwd(), "R/Graphs/Mm_embryo/cell_annotation_final")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
mm.embryo.sub <- list()
pd.col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
            "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
            "#00FBFF","#00BCBF","#D277FD","#AE00FF","#7789FD","#203EFF",
            "#FC8ACD","#FF20A4","#8AFCE9","#1CFFD9","#97FF99","#1CFF20",
            "#7B7B7B","#BFBFBF")
mm.embryo$Stage <- mm.embryo$development_stage
# - Craniofacial
mm.embryo.sub$Craniofacial <- mm.embryo[,grep("Craniofacial", mm.embryo$CellType)]
mm.embryo.sub$Craniofacial <- mm.embryo.sub$Craniofacial %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Craniofacial)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Craniofacial <- RunHarmony(mm.embryo.sub$Craniofacial, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Craniofacial, reduction = "harmony", ndims = 30)
dim.n <- 20
mm.embryo.sub$Craniofacial <- RunUMAP(mm.embryo.sub$Craniofacial, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Craniofacial$CellType.sub <- str_split_fixed(mm.embryo.sub$Craniofacial$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Craniofacial) <- mm.embryo.sub$Craniofacial$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Craniofacial_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Craniofacial, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Head mesoderm
mm.embryo.sub$HeadMesoderm <- mm.embryo[,grep("Head mesoderm", mm.embryo$CellType)]
mm.embryo.sub$HeadMesoderm <- mm.embryo.sub$HeadMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$HeadMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$HeadMesoderm <- RunHarmony(mm.embryo.sub$HeadMesoderm, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$HeadMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 15
mm.embryo.sub$HeadMesoderm <- RunUMAP(mm.embryo.sub$HeadMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$HeadMesoderm$CellType.sub <- str_split_fixed(mm.embryo.sub$HeadMesoderm$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$HeadMesoderm) <- mm.embryo.sub$HeadMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_HeadMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$HeadMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Lateral plate mesoderm
mm.embryo.sub$LateralPlateMesoderm <- mm.embryo[,grep("Lateral plate mesoderm", mm.embryo$CellType)]
table(mm.embryo.sub$LateralPlateMesoderm$CellType)
mm.embryo.sub$LateralPlateMesoderm <- mm.embryo.sub$LateralPlateMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$LateralPlateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$LateralPlateMesoderm <- RunHarmony(mm.embryo.sub$LateralPlateMesoderm, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$LateralPlateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$LateralPlateMesoderm <- RunUMAP(mm.embryo.sub$LateralPlateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$LateralPlateMesoderm$CellType.sub <- str_split_fixed(mm.embryo.sub$LateralPlateMesoderm$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$LateralPlateMesoderm) <- mm.embryo.sub$LateralPlateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_LateralPlateMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$LateralPlateMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Limb
mm.embryo.sub$Limb <- mm.embryo[,grep("Limb", mm.embryo$CellType)]
mm.embryo.sub$Limb <- mm.embryo.sub$Limb %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Limb)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Limb <- RunHarmony(mm.embryo.sub$Limb, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Limb, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$Limb <- RunUMAP(mm.embryo.sub$Limb, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Limb$CellType.sub <- str_split_fixed(mm.embryo.sub$Limb$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Limb) <- mm.embryo.sub$Limb$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Limb_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Limb, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Somite
mm.embryo.sub$Somite <- mm.embryo[,grep("Somite", mm.embryo$CellType)]
mm.embryo.sub$Somite <- mm.embryo.sub$Somite %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Somite)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Somite <- RunHarmony(mm.embryo.sub$Somite, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Somite, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$Somite <- RunUMAP(mm.embryo.sub$Somite, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Somite$CellType.sub <- str_split_fixed(mm.embryo.sub$Somite$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Somite) <- mm.embryo.sub$Somite$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Somite_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Somite, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Endothelium
mm.embryo.sub$Endothelium <- mm.embryo[,grep("Endothelium", mm.embryo$CellType)]
mm.embryo.sub$Endothelium <- mm.embryo.sub$Endothelium %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Endothelium)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Endothelium <- RunHarmony(mm.embryo.sub$Endothelium, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Endothelium, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$Endothelium <- RunUMAP(mm.embryo.sub$Endothelium, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Endothelium$CellType.sub <- str_split_fixed(mm.embryo.sub$Endothelium$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Endothelium) <- mm.embryo.sub$Endothelium$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endothelium_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Endothelium, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Blood
mm.embryo.sub$Blood <- mm.embryo[,grep("Blood", mm.embryo$CellType)]
mm.embryo.sub$Blood <- mm.embryo.sub$Blood %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Blood)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Blood <- RunHarmony(mm.embryo.sub$Blood, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Blood, reduction = "harmony", ndims = 30)
dim.n <- 8
mm.embryo.sub$Blood <- RunUMAP(mm.embryo.sub$Blood, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Blood$CellType.sub <- str_split_fixed(mm.embryo.sub$Blood$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Blood) <- mm.embryo.sub$Blood$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Blood_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Blood, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Intermediate mesoderm
mm.embryo.sub$IntermediateMesoderm <- mm.embryo[,grep("Intermediate mesoderm", mm.embryo$CellType)]
mm.embryo.sub$IntermediateMesoderm <- mm.embryo.sub$IntermediateMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$IntermediateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$IntermediateMesoderm <- RunHarmony(mm.embryo.sub$IntermediateMesoderm, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$IntermediateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 8
mm.embryo.sub$IntermediateMesoderm <- RunUMAP(mm.embryo.sub$IntermediateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$IntermediateMesoderm$CellType.sub <- str_split_fixed(mm.embryo.sub$IntermediateMesoderm$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$IntermediateMesoderm) <- mm.embryo.sub$IntermediateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_IntermediateMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$IntermediateMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Endoderm
mm.embryo.sub$Endoderm <- mm.embryo[,grep("Endoderm", mm.embryo$CellType)]
mm.embryo.sub$Endoderm <- mm.embryo.sub$Endoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Endoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Endoderm <- RunHarmony(mm.embryo.sub$Endoderm, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Endoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$Endoderm <- RunUMAP(mm.embryo.sub$Endoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Endoderm$CellType.sub <- str_split_fixed(mm.embryo.sub$Endoderm$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Endoderm) <- mm.embryo.sub$Endoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Endoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Neuron
mm.embryo.sub$Neuron <- mm.embryo[,grep("Neuron", mm.embryo$CellType)]
mm.embryo.sub$Neuron <- mm.embryo.sub$Neuron %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Neuron)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Neuron <- RunHarmony(mm.embryo.sub$Neuron, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Neuron, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$Neuron <- RunUMAP(mm.embryo.sub$Neuron, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Neuron$CellType.sub <- str_split_fixed(mm.embryo.sub$Neuron$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Neuron) <- mm.embryo.sub$Neuron$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Neuron_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Neuron, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Neural progenitor
mm.embryo.sub$NeuralProgenitor <- mm.embryo[,grep("Neural progenitor", mm.embryo$CellType)]
mm.embryo.sub$NeuralProgenitor <- mm.embryo.sub$NeuralProgenitor %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$NeuralProgenitor)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$NeuralProgenitor <- RunHarmony(mm.embryo.sub$NeuralProgenitor, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$NeuralProgenitor, reduction = "harmony", ndims = 30)
dim.n <- 10
mm.embryo.sub$NeuralProgenitor <- RunUMAP(mm.embryo.sub$NeuralProgenitor, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$NeuralProgenitor$CellType.sub <- str_split_fixed(mm.embryo.sub$NeuralProgenitor$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$NeuralProgenitor) <- mm.embryo.sub$NeuralProgenitor$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_NeuralProgenitor_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$NeuralProgenitor, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Schwann cells
mm.embryo.sub$Schwann <- mm.embryo[,grep("Schwann cells", mm.embryo$CellType)]
mm.embryo.sub$Schwann <- mm.embryo.sub$Schwann %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Schwann)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Schwann <- RunHarmony(mm.embryo.sub$Schwann, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Schwann, reduction = "harmony", ndims = 30)
dim.n <- 8
mm.embryo.sub$Schwann <- RunUMAP(mm.embryo.sub$Schwann, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Schwann$CellType.sub <- str_split_fixed(mm.embryo.sub$Schwann$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Schwann) <- mm.embryo.sub$Schwann$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Schwann_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Schwann, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Undefined
mm.embryo.sub$Undefined <- mm.embryo[,grep("undefined", mm.embryo$CellType)]
mm.embryo.sub$Undefined <- mm.embryo.sub$Undefined %>% FindVariableFeatures() %>%
  ScaleData(rownames(mm.embryo.sub$Undefined)) %>% RunPCA(verbose = FALSE)
set.seed(100)
mm.embryo.sub$Undefined <- RunHarmony(mm.embryo.sub$Undefined, group.by.vars = "Stage")
ElbowPlot(mm.embryo.sub$Undefined, reduction = "harmony", ndims = 30)
dim.n <- 15
mm.embryo.sub$Undefined <- RunUMAP(mm.embryo.sub$Undefined, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
mm.embryo.sub$Undefined$CellType.sub <- str_split_fixed(mm.embryo.sub$Undefined$CellType, ":", "2")[,2]
Idents(mm.embryo.sub$Undefined) <- mm.embryo.sub$Undefined$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Undefined_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(mm.embryo.sub$Undefined, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Correlation analysis based on specificity score matrix (Genes + TFs)
spe.score.seurat <- list()
for (cluster in names(hs.embryo.sub)) {
  expr.hs = AverageExpression(hs.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
  expr.mm = AverageExpression(mm.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
  deg.hs = VariableFeatures(hs.embryo.sub[[cluster]])
  deg.mm = VariableFeatures(mm.embryo.sub[[cluster]])
  spe.score.seurat[[cluster]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                 ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                                 Species1 = "human", Species2 = "mouse",
                                                 filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", cluster)),
                                                 Height = 12, Width = 12)
}
# all cell types
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.mm = AverageExpression(mm.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.mm = VariableFeatures(mm.embryo)
spe.score.seurat[["all"]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                             ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                             Species1 = "human", Species2 = "mouse",
                                             filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")))
CorrComparePlot.again(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                      ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                      objects.to.return = spe.score.seurat[["all"]],
                      Species1 = "human", Species2 = "mouse",
                      filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")),
                      Height = 70, Width = 70)
# all clusters
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.mm = AverageExpression(mm.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.mm = VariableFeatures(mm.embryo)
spe.score.seurat[["cluster"]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                 ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                                 Species1 = "human", Species2 = "mouse",
                                                 filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", "All_Super_Cluster")))
rm(cluster, expr.hs, expr.mm, deg.hs, deg.mm)
# - Correlation analysis based on specificity score matrix (Genes + LRs)
res.out <- "/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/Mm_embryo/cell_annotation_final/LRs"
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
spe.score.seurat.LR <- list()
for (cluster in names(hs.embryo.sub)) {
  expr.hs = AverageExpression(hs.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
  expr.mm = AverageExpression(mm.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
  deg.hs = VariableFeatures(hs.embryo.sub[[cluster]])
  deg.mm = VariableFeatures(mm.embryo.sub[[cluster]])
  spe.score.seurat.LR[[cluster]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                      ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                                      Species1 = "human", Species2 = "mouse",
                                                      filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", cluster)),
                                                      Height = 12, Width = 12)
}
# all cell types
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.mm = AverageExpression(mm.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.mm = VariableFeatures(mm.embryo)
spe.score.seurat.LR[["all"]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                  ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                                  Species1 = "human", Species2 = "mouse",
                                                  filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")))
CorrComparePlotLR.again(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                        ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                        objects.to.return = spe.score.seurat[["all"]],
                        Species1 = "human", Species2 = "mouse",
                        filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")),
                        Height = 70, Width = 70)
# all clusters
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.mm = AverageExpression(mm.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.mm) = paste("Mm", gsub("RNA.", "", colnames(expr.mm)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.mm = VariableFeatures(mm.embryo)
spe.score.seurat.LR[["cluster"]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                      ExpressionTableSpecies2 = expr.mm, DEgenesSpecies2 = deg.mm,
                                                      Species1 = "human", Species2 = "mouse",
                                                      filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", "All_Super_Cluster")))
rm(cluster, expr.hs, expr.mm, deg.hs, deg.mm)



### ===============================================
### 6th step: Annotate human embryo by mouse embryo ----
### ===============================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Mm_embryo")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Annotating cell types by Seurat (R)
rev.predict <- list()
rev.predict$hs <- hs.embryo[rownames(hs.embryo) %in% hs.mm.hm.tar$Gene.name,]
rev.predict$mm <- mm.embryo[rownames(mm.embryo) %in% hs.mm.hm.tar$Gene.name,]
anchors <- FindTransferAnchors(reference = rev.predict$mm, query = rev.predict$hs, dims = 1:30, reference.reduction = "pca")
rev.predict$map <- MapQuery(anchorset = anchors, reference = rev.predict$mm, query = rev.predict$hs,
                            refdata = list(celltype = "Main_cell_type"), reference.reduction = "pca", reduction.model = "umap")
if (all(colnames(rev.predict$hs) == colnames(rev.predict$map))) {
  # pig object
  rev.predict$hs$Rev.CellType.Seurat <- rev.predict$map$predicted.celltype
  # mapped object
  rev.predict$map$Rev.CellType.Seurat <- rev.predict$map$predicted.celltype
}
pd.col.fix <- pd.col[1:length(unique(rev.predict$mm$Main_cell_type))]
names(pd.col.fix) <- sort(unique(rev.predict$mm$Main_cell_type))
p1 <- DimPlot(rev.predict$mm, reduction = "umap", group.by = "Main_cell_type", label = T, pt.size = 0.15) +
  ggtitle("Mouse embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(rev.predict$map, reduction = "ref.umap", group.by = "Rev.CellType.Seurat", label = T, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "New_integrated_clustering_with_predicted_cell_types_by_seurat.pdf"), height = 5, width = 20.5)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, anchors)


### >>> 6. Annotating cell types by SingleR (R)
# - transform Seurat object into SingleCellExperiment object
hs.sce <- as.SingleCellExperiment(rev.predict$hs)
mm.sce <- as.SingleCellExperiment(rev.predict$mm)
# - cell type annotation
library(SingleR)
predictions.sgr <- SingleR(test=hs.sce, assay.type.test=1, ref=mm.sce, labels=mm.sce$Main_cell_type, num.threads=20)
# - visualization
if (all(colnames(rev.predict$hs) == rownames(predictions.sgr)) & all(colnames(rev.predict$map) == rownames(predictions.sgr))) {
  # pig object
  rev.predict$hs$Rev.CellType.SingleR <- predictions.sgr$labels
  # mapped object
  rev.predict$map$Rev.CellType.SingleR <- predictions.sgr$labels
}
p1 <- DimPlot(rev.predict$mm, reduction = "umap", group.by = "Main_cell_type", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(rev.predict$map, reduction = "ref.umap", group.by = "Rev.CellType.SingleR", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "New_integrated_clustering_with_predicted_cell_types_by_singleR.pdf"), height = 5, width = 20.5)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, hs.sce, mm.sce)


### >>> 7. Annotating cell types by SciBet (R)
library("Rcpp")
library("RcppEigen")
library("ggsci")
library("viridis")
library("tidyverse")
library("scibet")
# - prepare data
ref.set <- GetAssayData(mm.embryo, slot = "data") %>% t() %>% as.data.frame()
ref.set$label <- mm.embryo$Main_cell_type
query.set <- GetAssayData(hs.embryo, slot = "data") %>% t() %>% as.data.frame()
# - predict cell type
predictions.scibet.5k <- SciBet(ref.set, query.set, k = 5000)
predictions.scibet.4k <- SciBet(ref.set, query.set, k = 4000)
predictions.scibet.3k <- SciBet(ref.set, query.set, k = 3000)
predictions.scibet.2k <- SciBet(ref.set, query.set, k = 2000)
predictions.scibet.1k <- SciBet(ref.set, query.set, k = 1000)
# - plot
rev.predict$map$Rev.CellType.SciBet <- predictions.scibet.1k
p1 <- DimPlot(rev.predict$mm, reduction = "umap", group.by = "Main_cell_type", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(rev.predict$map, reduction = "ref.umap", group.by = "Rev.CellType.SciBet", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "New_integrated_clustering_with_predicted_cell_types_by_SciBet.pdf"), height = 5, width = 20.5)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, ref.set, query.set)


### >>> 9. Compare the annotation between SingleR, Seurat and sciBet
# majority-rule
cell.type <- rev.predict$map@meta.data[,(ncol(rev.predict$map@meta.data)-2):ncol(rev.predict$map@meta.data)]
cell.type$CellType.Stat <- apply(cell.type, 1, function(x){length(sort(table(as.character(x)),decreasing=TRUE))})
cell.type$CellType.Final <- apply(cell.type, 1, function(x){names(sort(table(as.character(x)),decreasing=TRUE)[1])})
cell.type$CellType.Final[cell.type$CellType.Stat == 3] <- cell.type$Rev.CellType.Seurat[cell.type$CellType.Stat == 3]
if (identical(colnames(rev.predict$map), rownames(cell.type)) & identical(colnames(hs.embryo), rownames(cell.type))) {
  rev.predict$map$Rev.CellType.Final <- cell.type$CellType.Final
  hs.embryo$Rev.CellType <- cell.type$CellType.Final
}
p1 <- DimPlot(rev.predict$mm, reduction = "umap", group.by = "Main_cell_type", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(rev.predict$map, reduction = "ref.umap", group.by = "Rev.CellType.Final", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "New_integrated_clustering_with_predicted_cell_types_by_aggregation.pdf"), height = 5, width = 20.5)
p1+p2
dev.off()
# - proportion of cell types
# human
pd <- hs.embryo@meta.data[,c("SuperCluster","CellType.sub","Rev.CellType")]
library(webr)
pdf(file.path(res.out, "New_PieDonut_plot_to_show_human_cell_type_proportion.pdf"), height = 6, width = 6)
PieDonut(pd, aes(pies = Rev.CellType), addDonutLabel = TRUE, showRatioDonut = T,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2, r0 = 0.5) +
  scale_fill_manual(values = pd.col.fix)
dev.off()
# mouse
pd <- mm.embryo@meta.data[,c("SuperCluster","CellType","Main_cell_type")]
library(webr)
pdf(file.path(res.out, "New_PieDonut_plot_to_show_mouse_cell_type_proportion.pdf"), height = 6, width = 6)
PieDonut(pd, aes(pies = Main_cell_type), addDonutLabel = TRUE, showRatioDonut = T,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2, r0 = 0.5) +
  scale_fill_manual(values = pd.col.fix)
dev.off()

# - heatmap
library(viridis)
predictions <- table(rev.predict$map$Rev.CellType.Final, rev.predict$map$Rev.CellType.SciBet)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (SciBet)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
predictions <- table(rev.predict$map$Rev.CellType.Final, rev.predict$map$Rev.CellType.Seurat)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p2 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (Seurat)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
predictions <- table(rev.predict$map$Rev.CellType.Final, rev.predict$map$Rev.CellType.SingleR)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p3 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_viridis() +
  xlab("Cell type annotation (Final)") + ylab("Cell type annotation (SingleR)") +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0))
pdf(file.path(res.out, "New_comparing_of_cell_type_predicted_by_Seurat_SingleR_SciBet_heatmap.pdf"), width = 30, height = 9.5)
p1 + p2 + p3
dev.off(); rm(p1, p2, p3)



### ====================
### Last part: Save Data ----
### ====================
saveRDS(mm.embryo.raw, "mm.embryo.final.version.rds")
save.image("Mm_Embryo_Multiome.RData")
rm(list = ls())
