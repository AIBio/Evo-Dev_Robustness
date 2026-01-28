# ========================
# 1st part: Global setting ----
# ========================

### >>> 1. Setting workding directory
setwd("/home/yhw/bioinfo/project-mine/MultiOmics")
dir.create("R/CodeData", recursive = T)
dir.create("R/Graphs", recursive = T)
dir.create("R/Table", recursive = T)


### >>> 2. Setting library
pkg.lib <- "/home/laborer/software/anaconda3/envs/rs-4.2.3/lib/R/library"
.libPaths(pkg.lib)



# =================
# 2nd part: Library ----
# =================

### >>> 1. Packages
cran.pks <- c("SCP", "forcats", "dplyr", "tidyr", "tidyverse", "stringr", "circlize", 
              "ggplot2", "ggalluvial", "ggrepel", "RColorBrewer", "cowplot",
              "Seurat", "harmony", "COSG")
for (pks in cran.pks) {
  library(pks, character.only = T)
}


### >>> 2. Functions
theme_dp <- theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
                  axis.title.y = element_text(face = "plain", colour = "#000000", size = 13, angle = 90),
                  axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
                  axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
cluster.col <- c("#607d8b","#795548","#ff5722","#ffc107","#cddc39","#4caf50","#009688",
                 "#00bcd4","#2196f3","#3f51b5","#673ab7","#9c27b0","#e91e63","#f44336",
                 "#b0bec5","#bcaaa4","#ffab91","#ffe082","#e6ee9c","#a5d6a7","#80cbc4",
                 "#80deea","#90caf9","#9fa8da","#b39ddb","#ce93d8","#f48fb1","#ef9a9a",
                 "#37474f","#4e342e","#d84315","#ff8f00","#9e9d24","#2e7d32","#00695c",
                 "#00838f","#1565c0","#283593","#4527a0","#6a1b9a","#ad1457","#c62828")



# =========================
# 3rd part: Seurat pipeline ----
# =========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Rabbit_Embryos")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
pd.col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
            "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
            "#00FBFF","#00BCBF","#D277FD","#AE00FF","#7789FD","#203EFF",
            "#FC8ACD","#FF20A4","#8AFCE9","#1CFFD9","#97FF99","#1CFF20",
            "#7B7B7B","#BFBFBF")


### >>> 2. Load data
raw.file <- "/home/data/YHW_rabbit_scRNAseq/analysis/results/crrna/count/OryCun2_rs112"
sr.dir <- grep("filtered_feature_bc_matrix$", list.dirs(raw.file, full.names = T, recursive = T), value = T)
sr.list <- list()
for (dir in sr.dir) {
  sr.list[[dir]] <- LoadSCdata(
    mode = "10x", dir = dir,
    proj.name = str_split_fixed(dir, "/", 11)[, 10]
  )
}
names(sr.list) <- str_split_fixed(sr.dir, "/", 11)[, 10]
sr.rabbit <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)), sr.list)
meta <- sr.rabbit@meta.data
rownames(meta) <- gsub("SeuratProject_", "", rownames(meta))
count <- GetAssayData(sr.rabbit, slot = "count")
count <- count[-grep("ENSOCUG", rownames(count)), ]
colnames(count) <- gsub("SeuratProject_", "", colnames(count))
if (all(colnames(count) == rownames(meta))) {
  sr.rabbit <- CreateSeuratObject(counts = count, meta.data = meta)
  rm(count, meta)
}
table(sr.rabbit$orig.ident)
sr.rabbit@meta.data <- sr.rabbit@meta.data %>% 
  mutate(Group = case_when(orig.ident == "Ra-E10-1-10SN" ~ "E10-1", orig.ident == "Ra-E10-2-10SN" ~ "E10-2",
                           orig.ident == "Ra-E12-1-10SN" ~ "E12-1", orig.ident == "Ra-E12-2-10SN" ~ "E12-2",
                           orig.ident == "E14-1-Fetal-1-10SN" ~ "E14-1", orig.ident == "E14-1-Fetal-2-10SN" ~ "E14-2"),
         Stage = case_when(orig.ident %in% c("Ra-E10-1-10SN", "Ra-E10-2-10SN") ~ "E10", 
                           orig.ident %in% c("Ra-E12-1-10SN", "Ra-E12-2-10SN") ~ "E12", 
                           orig.ident %in% c("E14-1-Fetal-1-10SN", "E14-1-Fetal-2-10SN") ~ "E14"))


### >>> 3. Seurat pipeline
# MT ratio
mt.genes <- read.table("/home/data/YHW_rabbit_scRNAseq/analysis/metadata/MT_gene.txt")
sr.rabbit[["percent.mt"]] <- round(PercentageFeatureSet(sr.rabbit, features = mt.genes$V1), 2)
# plot QC
pdf(file.path(res.out, "Quality_control_before_filtering.pdf"), height = 5, width = 15)
FeatureStatPlot(sr.rabbit, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                group.by = "orig.ident", add_box = T)
dev.off()
dim(sr.rabbit)
# filter cells
sr.rabbit <- subset(sr.rabbit, subset = nFeature_RNA >= 500 & nFeature_RNA <= 5000 & nCount_RNA >= 1000 & percent.mt <= 30)
dim(sr.rabbit)
table(sr.rabbit$Group)
# plot QC
pdf(file.path(res.out, "Quality_control_after_filtering.pdf"), height = 5, width = 15)
FeatureStatPlot(sr.rabbit, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                group.by = "orig.ident", add_box = T)
dev.off()
# normalization and PCA
sr.rabbit <- NormalizeData(sr.rabbit, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = rownames(sr.rabbit)) %>% 
  RunPCA(verbose = FALSE)
pdf(file.path(res.out, "Elbowplot_after_filtering.pdf"), height = 5, width = 8)
ElbowPlot(sr.rabbit, ndims = 30)
dev.off()
# clustering
dims <- 15
sr.rabbit <- RunUMAP(sr.rabbit, dims = 1:dims, reduction = "pca")
CellDimPlot(srt = sr.rabbit, group.by = "Group", pt.size = 0.25, 
            reduction = "UMAP", theme_use = "theme_blank")
sr.rabbit <- RunHarmony(sr.rabbit, group.by.vars = "Group", dims = 1:dims)
sr.rabbit <- RunUMAP(sr.rabbit, reduction = "harmony", dims = 1:dims)
sr.rabbit <- FindNeighbors(sr.rabbit, reduction = "harmony", dims = 1:dims) %>%
  FindClusters(resolution = seq(0, 4, 0.25))
CellDimPlot(srt = sr.rabbit, group.by = c("Group", "RNA_snn_res.3"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank")
table(sr.rabbit$RNA_snn_res.3, sr.rabbit$Group)
table(sr.rabbit$orig.ident)
# find markers
Idents(sr.rabbit) <- sr.rabbit$RNA_snn_res.3
marker.genes <- FindAllMarkers(sr.rabbit, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(marker.genes, file.path(res.out, "All_marker_genes.csv"))
marker.genes <- cosg(sr.rabbit, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100)
write.csv(marker.genes, file.path(res.out, "All_marker_genes_COSG.csv"))
grep("FOX", rownames(sr.rabbit), value = T)
dev.off()
FeatureDimPlot(srt = sr.rabbit, features = c("FOXA2", "SPON1", "ARX", "SHH", "TBXT"), reduction = "UMAP", theme_use = "theme_blank")
FeatureDimPlot(srt = sr.rabbit, features = c("LMX1A", "MSX2", "RSPO3", "BMP6"), reduction = "UMAP", theme_use = "theme_blank")
FeatureDimPlot(srt = sr.rabbit, features = c("PAX6", "OTX2", "FOXC1", "FOXC2"), reduction = "UMAP", theme_use = "theme_blank")
FeatureStatPlot(sr.rabbit, stat.by = c("PAX6", "NKX6-1", "TBXT", "FOXA2", "SPON1", "ARX", "SHH", "DHH", "IHH"), 
                group.by = "RNA_snn_res.3", split.by = "Stage", plot_type = "violin", stack = T)
FeatureStatPlot(sr.rabbit, stat.by = c("LMX1A", "MAFB", "MSX1", "RSPO3", "BMP7", "GSC"), 
                group.by = "RNA_snn_res.3", split.by = "Stage", plot_type = "violin", stack = T)


### >>> 4. Preprocess data
# - human seurat object
hs.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.new.final.rds")
hs.embryo$SuperCluster <- str_split_fixed(hs.embryo$CellType, ":", 2)[, 1]
Idents(hs.embryo$SuperCluster)
marker.genes <- MarkerGene(
  sr.obj = hs.embryo, method = "cosg", top.n = 100,
  assay = "RNA", slot = "data", idents = "SuperCluster",
  cols.pal = "Blues", cols.rev = F,
  plt.min = 0, plt.max = 0.5, file.name = "file.name",
  res.out = file.path(res.out, "Marker_gene_of_human_clusters")
)
# - human and pig homologous genes
hs.ra <- read.table("/home/yhw/document/ensembl/Homologous_Genes/release112_GRCh38_OryCun2.0_mart_export.txt", 
                    header = T, stringsAsFactors = F, sep = "\t")
hs.ra <- subset(hs.ra, Rabbit.homology.type == "ortholog_one2one")
# 
hs.ra <- unique(hs.ra[,c(1,8)]) %>%
  dplyr::filter(Gene.name %in% rownames(hs.embryo)) %>%
  dplyr::filter(Rabbit.gene.name %in% rownames(sr.rabbit))
sr.predict <- list()
sr.predict$hs <- hs.embryo[rownames(hs.embryo) %in% hs.ra$Gene.name, ]
sr.predict$ra <- sr.rabbit[rownames(sr.rabbit) %in% hs.ra$Rabbit.gene.name, ]


### >>> 5. Annotating cell types by Seurat
anchors <- FindTransferAnchors(reference = sr.predict$hs, query = sr.predict$ra, dims = 1:30, reference.reduction = "pca")
sr.predict$map <- MapQuery(anchorset = anchors, reference = sr.predict$hs, query = sr.predict$ra,
                           refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
sr.predict$map$CellType.Seurat <- sr.predict$map$predicted.celltype
sr.predict$map$SuperCluster.Seurat <- str_split_fixed(sr.predict$map$predicted.celltype, ":", 2)[,1]
CellDimPlot(srt = sr.predict$map, group.by = c("SuperCluster.Seurat"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank")


### >>> 6. Annotating cell types by SingleR
library("SingleR")
# - transform Seurat object into SingleCellExperiment object
hs.sce <- as.SingleCellExperiment(sr.predict$hs)
ra.sce <- as.SingleCellExperiment(sr.predict$ra)
# - cell type annotation
predictions.sgr <- SingleR(test = ra.sce, assay.type.test = 1, ref = hs.sce, labels = hs.sce$CellType, num.threads = 12)
# - visualization
if (all(colnames(sr.predict$map) == rownames(predictions.sgr))) {
  sr.predict$map$CellType.SingleR <- predictions.sgr$labels
  sr.predict$map$SuperCluster.SingleR <- str_split_fixed(sr.predict$map$CellType.SingleR, ":", 2)[, 1]
}


### >>> 7. Annotating cell types by SciBet (R)
library("scibet")
# - prepare data
ref.set <- GetAssayData(sr.predict$hs, slot = "data", assay = "RNA") %>% as.data.frame()
ref.set <- t(ref.set)
ref.set$label <- sr.predict$hs$CellType
query.set <- GetAssayData(sr.predict$ra, slot = "data", assay = "RNA") %>% as.data.frame()
query.set <- t(query.set)
# - predict cell type
predictions.scibet.4k <- SciBet(ref.set, query.set, k = 4000)
predictions.scibet.3k <- SciBet(ref.set, query.set, k = 3000)
predictions.scibet.2k <- SciBet(ref.set, query.set, k = 2000)
predictions.scibet.1k <- SciBet(ref.set, query.set, k = 1000)
rm(ref.set, query.set)
save.image("/home/laborer/Hs_Embryo_Multiome_Rabbit.RData")
# - visualization
unique(predictions.scibet.1k)
sr.predict$map$CellType.SciBet <- predictions.scibet.1k
sr.predict$map$SuperCluster.SciBet <- str_split_fixed(sr.predict$map$CellType.SciBet, ":", 2)[, 1]


### >>> 8. Compare the annotation between SingleR, Seurat and sciBet
# majority-rule
cell.type <- sr.predict$map@meta.data[,(ncol(sr.predict$map@meta.data)-4):ncol(sr.predict$map@meta.data)]
cell.type$CellType.Stat <- apply(cell.type[,seq(1,4,2)], 1, function(x){length(sort(table(as.character(x)),decreasing=TRUE))})
cell.type$CellType.Final <- apply(cell.type[,seq(1,4,2)], 1, function(x){names(sort(table(as.character(x)),decreasing=TRUE)[1])})
cell.type$CellType.Final[cell.type$CellType.Stat == 2] <- cell.type$CellType.SingleR[cell.type$CellType.Stat == 2]
table(cell.type$CellType.Final) %>% as.data.frame() %>% arrange(Var1)
if (identical(colnames(sr.predict$map), rownames(cell.type)) & identical(colnames(sr.rabbit), rownames(cell.type))) {
  sr.predict$map$CellType.Final <- cell.type$CellType.Final
  sr.predict$map$SuperCluster.Final <- str_split_fixed(sr.predict$map$CellType.Final, ":", 2)[,1]
  sr.rabbit$CellType <- cell.type$CellType.Final
  sr.rabbit$SuperCluster <- str_split_fixed(sr.rabbit$CellType, ":", 2)[,1]
}
pd.col.fix <- pd.col[1:length(unique(sr.predict$hs$SuperCluster))]
names(pd.col.fix) <- sort(unique(sr.predict$hs$SuperCluster))
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_aggregation.pdf"), height = 5, width = 15)
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map, reduction = "ref.umap", group.by = "SuperCluster.Final", label = F, pt.size = 0.15) +
  ggtitle("Rabbit embryos") +
  scale_color_manual(values = pd.col.fix)
p1 + p2
dev.off()
unique(sr.rabbit$SuperCluster) %>% length()

sr.rabbit$SuperCluster  <- as.character(sr.rabbit$SuperCluster)
sr.rabbit$SuperCluster <- factor(sr.rabbit$SuperCluster, 
                                 levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", 
                                            "Craniofacial", "Head mesoderm", "Somite", "Intermediate mesoderm", 
                                            "Limb","Lateral plate mesoderm",
                                            "Heart", "Endothelium", "Blood", "Endoderm","Fibroblast"))
Idents(sr.rabbit) <- sr.rabbit$SuperCluster

marker.genes$markers$names$Endoderm
pd.gene <- c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","STMN2","DCX","TP63","MYH3","ASB5","XIRP2",
             "TUBB3","FOXD3","MPZ","PLP1","S100B","ALX1","ALX3",
             "DLX1","DLX2","SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9","PAX8","WT1","TBX5","PITX1","FOXF1","GATA6",
             "GATA5","NKX2-5","PLVAP","CD53","PLAC8","CD68","CORO1A","FOXA2","CLDN6","PERP","KRT18","VTN","PLG","FGG")
pd.gene <- intersect(pd.gene, rownames(sr.rabbit))
library(SeuratWrappers)
sr.rabbit <- RunALRA(sr.rabbit)
pdf(file.path(res.out, "Dotplot_to_show_rabbit_super_cluster_markers.pdf"), height = 6, width = 14)
DotPlot(sr.rabbit, features = pd.gene, group.by = "SuperCluster", assay = "alra",
        cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 75) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
dev.off()


### >>> CAME analysis
# process homologous genes
homo.genes <- read.table("/home/yhw/document/ensembl/Homologous_Genes/release112_GRCh38_OryCun2.0_mart_export.txt", header = T, sep = "\t")
homo.genes <- homo.genes[homo.genes$Rabbit.homology.type != "", ]
homo.genes <- homo.genes[homo.genes$Gene.name != "", ]
homo.genes <- homo.genes[homo.genes$Rabbit.gene.name != "", ]
write.csv(homo.genes[, c("Gene.name", "Rabbit.gene.name", 
                         "Rabbit.Gene.order.conservation.score",
                         "Rabbit.Whole.genome.alignment.coverage",
                         "Rabbit.orthology.confidence..0.low..1.high.")],
          "/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_human2rabbit_final.csv", row.names = F)
homo.genes.1v1 <- subset(homo.genes, Rabbit.homology.type == "ortholog_one2one")
write.csv(homo.genes.1v1[, c("Gene.name", "Rabbit.gene.name", 
                             "Rabbit.Gene.order.conservation.score",
                             "Rabbit.Whole.genome.alignment.coverage",
                             "Rabbit.orthology.confidence..0.low..1.high.")],
          "/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2rabbit_final.csv", row.names = F)
# output h5 file
library("dior")
rownames(sr.rabbit)
table(sr.rabbit$SuperCluster)
tmp <- sr.rabbit
tmp <- subset(tmp, SuperCluster != "Undefined")
tmp <- CreateSeuratObject(counts = GetAssayData(tmp, slot = "count"), meta.data = tmp@meta.data)
dior::write_h5(tmp, file = paste0("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/CAME_ra_embryo_All.h5"), object.type = 'seurat')
# make scatter plot
came.out <- "/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/came/_temp/('multiome-human-embryo', 'multiome-rabbit-embryo')-(02-21 17.04.43)/"
adata <- dior::read_h5(file=file.path(came.out, "adt_hidden_cell.h5"), 
                       target.object = 'seurat')
table(adata$REF, adata$dataset)
Idents(adata) <- adata$REF
adata$REF <- factor(adata$REF,
                    levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", 
                               "Craniofacial", "Head mesoderm", "Somite", "Intermediate mesoderm", 
                               "Limb","Lateral plate mesoderm",
                               "Heart", "Endothelium", "Blood", "Endoderm"))
supercluster.col <- c("#8c564b", "#d62728", "#ffbb78", "#ff7f0e", "#aec7e8", "#9467bd", "#1f77b4",
                      "#2ca02c", "#ff9896", "#98df8a", "#f7b6d2", "#c49c94", "#e377c2", "#c5b0d5")
dataset.col <- c("#1f77b4", "#98df8a")
names(supercluster.col) <- c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", 
                             "Craniofacial", "Head mesoderm", "Somite", "Intermediate mesoderm", 
                             "Limb","Lateral plate mesoderm", "Heart", "Endothelium", "Blood", "Endoderm")
adata <- subset(adata, downsample = 6000)
p1 <- CellDimPlot(srt = adata, group.by = "REF", palcolor = supercluster.col,
                  pt.size = 0.1, pt.alpha = 1, label_insitu = F, label = F, label_repel = T, 
                  reduction = "UMAP", theme_use = "theme_blank", raster = F)
p2 <- CellDimPlot(srt = adata, group.by = "dataset", palcolor = dataset.col,
                  pt.size = 0.1, pt.alpha = 1, label_insitu = F, label = F, label_repel = T, 
                  reduction = "UMAP", theme_use = "theme_blank", raster = F)
pdf(file.path(came.out, "Scatter_plot_updated.pdf"), height = 6, width = 15)
p1 + p2
dev.off()
# 
came.out <- "/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/came/_temp/('multiome-human-embryo', 'multiome-mouse-embryo')-(02-18 13.09.41)"
adata <- dior::read_h5(file=file.path(came.out, "adt_hidden_cell.h5"), target.object = 'seurat')
adata <- subset(adata, dataset == "multiome-human-embryo")
common.cell <- intersect(colnames(sr.predict$hs), adata@meta.data$original_name)
adata.umap <- Embeddings(adata, reduction = "umap")
rownames(adata.umap) <- adata@meta.data$original_name
adata.umap <- adata.umap[common.cell, ]
colnames(adata.umap) <- c("UMAP_1", "UMAP_2")
table(colnames(sr.predict$hs) %in% adata@meta.data$original_name)
sr.predict$hs.new <- sr.predict$hs[, match(common.cell, colnames(sr.predict$hs))]
sr.predict$hs.new <- RunUMAP(sr.predict$hs.new, dims = 1:20, reduction = "pca", return.model = T)
if (all(colnames(sr.predict$hs.new) == rownames(adata.umap))) {
  sr.predict$hs.new[["umap"]] <- CreateDimReducObject(embeddings = adata.umap, key = "UMAP_", 
                                                      assay = DefaultAssay(sr.predict$hs.new))
}
sr.predict$ra <- NormalizeData(sr.predict$ra, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)
sr.predict$ra <- RunUMAP(sr.predict$ra, dims = 1:20, reduction = "pca", return.model = T)
anchors <- FindTransferAnchors(reference = sr.predict$hs.new, query = sr.predict$ra, 
                               dims = 1:20, reference.reduction = "pca")
sr.predict$map2 <- MapQuery(anchorset = anchors, reference = sr.predict$hs.new, query = sr.predict$ra,
                            refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
DimPlot(sr.predict$map2, reduction = "ref.umap", group.by = "SuperCluster.Final", label = F, pt.size = 0.15) +
  ggtitle("Rabbit embryos") +
  scale_color_manual(values = pd.col.fix)



# =================================
# 4th part: Sub-population analysis ----
# =================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Rabbit_Embryos/Sub-population")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Human embryo analysis
hs.embryo.pop <- subset(sr.predict$hs, SuperCluster %in% c("Neuron", "Neural progenitor"))
colnames(hs.embryo.pop@meta.data)
hs.embryo.pop@meta.data <- hs.embryo.pop@meta.data[, c(-13:-15)]
hs.embryo.pop$Batch <- hs.embryo.pop$Sample
hs.embryo.pop <- NormalizeData(hs.embryo.pop) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 1000)
hs.embryo.pop <- ScaleData(hs.embryo.pop, verbose = FALSE, features = rownames(hs.embryo.pop))
hs.embryo.pop <- RunPCA(hs.embryo.pop, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.embryo.pop))
ElbowPlot(hs.embryo.pop, ndims = 30)
dim.n <- 12
hs.embryo.pop <- RunHarmony(hs.embryo.pop, group.by.vars = "Batch", dims = 1:dim.n)
hs.embryo.pop <- RunUMAP(hs.embryo.pop, reduction = "harmony", dims = 1:dim.n)
DimPlot(hs.embryo.pop, group.by = "SuperCluster")
hs.embryo.pop@meta.data <- hs.embryo.pop@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
                                                       "Batch", "CellType", "SuperCluster")]
hs.embryo.pop@meta.data$Species <- "Human"
# somite
hs.embryo.sm <- subset(sr.predict$hs, SuperCluster %in% c("Somite"))
colnames(hs.embryo.sm@meta.data)
hs.embryo.sm@meta.data <- hs.embryo.sm@meta.data[, c(-13:-15)]
hs.embryo.sm@meta.data <- hs.embryo.sm@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
                                                     "Batch", "CellType", "SuperCluster")]
hs.embryo.sm@meta.data$Species <- "Human"


### >>> 3. Rabbit embryo analysis
ra.embryo.pop <- subset(sr.predict$map, SuperCluster.Final %in% c("Neuron", "Neural progenitor"))
colnames(ra.embryo.pop@meta.data)
ra.embryo.pop@meta.data <- ra.embryo.pop@meta.data[, c(-7:-25)]
ra.embryo.pop$Batch <- ra.embryo.pop$Group
ra.embryo.pop <- NormalizeData(ra.embryo.pop) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 1000)
ra.embryo.pop <- ScaleData(ra.embryo.pop, verbose = FALSE, features = rownames(ra.embryo.pop))
ra.embryo.pop <- RunPCA(ra.embryo.pop, npcs = 30, verbose = FALSE, features = VariableFeatures(ra.embryo.pop))
ElbowPlot(ra.embryo.pop, ndims = 30)
dim.n <- 12
ra.embryo.pop <- RunHarmony(ra.embryo.pop, group.by.vars = "Group", dims = 1:dim.n)
ra.embryo.pop <- RunUMAP(ra.embryo.pop, reduction = "harmony", dims = 1:dim.n)
DimPlot(ra.embryo.pop, group.by = "SuperCluster.Final")
ra.embryo.pop$CellType <- ra.embryo.pop$CellType.Final
ra.embryo.pop$SuperCluster <- ra.embryo.pop$SuperCluster.Final
ra.embryo.pop@meta.data <- ra.embryo.pop@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
                                                       "Batch", "CellType", "SuperCluster")]
ra.embryo.pop@meta.data$Species <- "Rabbit"
colnames(ra.embryo.pop@meta.data)
# somite
ra.embryo.sm <- subset(sr.predict$map, SuperCluster.Final %in% c("Somite"))
colnames(ra.embryo.sm@meta.data)
ra.embryo.sm@meta.data <- ra.embryo.sm@meta.data[, c(-7:-25)]
ra.embryo.sm$Batch <- ra.embryo.sm$Group
ra.embryo.sm$CellType <- ra.embryo.sm$CellType.Final
ra.embryo.sm$SuperCluster <- ra.embryo.sm$SuperCluster.Final
ra.embryo.sm@meta.data <- ra.embryo.sm@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
                                                       "Batch", "CellType", "SuperCluster")]
ra.embryo.sm@meta.data$Species <- "Rabbit"
colnames(ra.embryo.sm@meta.data)


### >>> 4. Mouse embryo analysis
mm.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/mm.embryo.final.version.rds")
mm.embryo$SuperCluster <- str_split_fixed(mm.embryo$CellType.Final,":",2)[,1]
mm.embryo$CellType <- mm.embryo$CellType.Final
mm.embryo$Batch <- mm.embryo$development_stage
mm.embryo.pop <- subset(mm.embryo, SuperCluster %in% c("Neuron", "Neural progenitor"))
colnames(mm.embryo.pop@meta.data)
mm.embryo.pop@meta.data <- mm.embryo.pop@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Batch",
                                                       "percent.mt", "CellType", "SuperCluster")]
count <- GetAssayData(mm.embryo.pop, slot = "count") %>% as.data.frame()
meta <- mm.embryo.pop@meta.data
homo.gene <- read.csv("/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2mouse_final.csv")
rownames(homo.gene) <- homo.gene$gene.name
keep.gene <- intersect(rownames(homo.gene), rownames(count))
homo.gene <- homo.gene[keep.gene, ]
count <- count[keep.gene, ]
if (all(rownames(count) == rownames(homo.gene))) {
  rownames(count) <- homo.gene$human.gene.name
}
mm.embryo.pop <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, homo.gene, keep.gene)
mm.embryo.pop <- NormalizeData(mm.embryo.pop) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mm.embryo.pop <- ScaleData(mm.embryo.pop, verbose = FALSE, features = rownames(mm.embryo.pop))
mm.embryo.pop <- RunPCA(mm.embryo.pop, npcs = 30, verbose = FALSE, features = VariableFeatures(mm.embryo.pop))
ElbowPlot(mm.embryo.pop, ndims = 30)
dim.n <- 12
mm.embryo.pop <- RunHarmony(mm.embryo.pop, group.by.vars = "Batch", dims = 1:dim.n)
mm.embryo.pop <- RunUMAP(mm.embryo.pop, reduction = "harmony", dims = 1:dim.n)
DimPlot(mm.embryo.pop, group.by = "SuperCluster")
mm.embryo.pop@meta.data$Species <- "Mouse"
colnames(mm.embryo.pop@meta.data)
mm.embryo.pop <- mm.embryo.pop[rownames(mm.embryo.pop) %in% rownames(ra.embryo.pop), ]
# somite
mm.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/mm.embryo.final.version.rds")
mm.embryo$SuperCluster <- str_split_fixed(mm.embryo$CellType.Final,":",2)[,1]
mm.embryo$CellType <- mm.embryo$CellType.Final
mm.embryo$Batch <- mm.embryo$development_stage
mm.embryo.sm <- subset(mm.embryo, SuperCluster %in% c("Somite"))
colnames(mm.embryo.sm@meta.data)
mm.embryo.sm@meta.data <- mm.embryo.sm@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Batch",
                                                     "percent.mt", "CellType", "SuperCluster")]
count <- GetAssayData(mm.embryo.sm, slot = "count") %>% as.data.frame()
meta <- mm.embryo.sm@meta.data
homo.gene <- read.csv("/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2mouse_final.csv")
rownames(homo.gene) <- homo.gene$gene.name
keep.gene <- intersect(rownames(homo.gene), rownames(count))
homo.gene <- homo.gene[keep.gene, ]
count <- count[keep.gene, ]
if (all(rownames(count) == rownames(homo.gene))) {
  rownames(count) <- homo.gene$human.gene.name
}
mm.embryo.sm <- CreateSeuratObject(counts = count, meta.data = meta)
mm.embryo.sm <- mm.embryo.sm[rownames(mm.embryo.sm) %in% rownames(ra.embryo.sm), ]
mm.embryo.sm@meta.data$Species <- "Mouse"
rm(count, homo.gene, keep.gene, meta, mm.embryo)


### >>> 5. Pig embryo analysis
ss.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/ss.embryo.final.version.rds")
ss.embryo$SuperCluster <- str_split_fixed(ss.embryo$CellType,":",2)[,1]
ss.embryo$CellType <- ss.embryo$CellType
ss.embryo$Batch <- ss.embryo$Stage
ss.embryo.pop <- subset(ss.embryo, SuperCluster %in% c("Neuron", "Neural progenitor"))
colnames(ss.embryo.pop@meta.data)
ss.embryo.pop@meta.data <- ss.embryo.pop@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Batch",
                                                       "percent.mt", "CellType", "SuperCluster")]
count <- GetAssayData(ss.embryo.pop, slot = "count") %>% as.data.frame()
meta <- ss.embryo.pop@meta.data
homo.gene <- read.csv("/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2pig_final.csv")
homo.gene <- homo.gene[!duplicated(homo.gene$gene.name), ]
rownames(homo.gene) <- homo.gene$gene.name
keep.gene <- intersect(rownames(homo.gene), rownames(count))
homo.gene <- homo.gene[keep.gene, ]
count <- count[keep.gene, ]
if (all(rownames(count) == rownames(homo.gene))) {
  rownames(count) <- homo.gene$human.gene.name
}
ss.embryo.pop <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, homo.gene, keep.gene)
ss.embryo.pop <- NormalizeData(ss.embryo.pop) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
ss.embryo.pop <- ScaleData(ss.embryo.pop, verbose = FALSE, features = rownames(ss.embryo.pop))
ss.embryo.pop <- RunPCA(ss.embryo.pop, npcs = 30, verbose = FALSE, features = VariableFeatures(ss.embryo.pop))
ElbowPlot(ss.embryo.pop, ndims = 30)
dim.n <- 12
ss.embryo.pop <- RunHarmony(ss.embryo.pop, group.by.vars = "Batch", dims = 1:dim.n)
ss.embryo.pop <- RunUMAP(ss.embryo.pop, reduction = "harmony", dims = 1:dim.n)
DimPlot(ss.embryo.pop, group.by = "SuperCluster")
ss.embryo.pop@meta.data$Species <- "Pig"
colnames(ss.embryo.pop@meta.data)
ss.embryo.pop <- ss.embryo.pop[rownames(ss.embryo.pop) %in% rownames(ra.embryo.pop), ]
# somite
ss.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/ss.embryo.final.version.rds")
ss.embryo$SuperCluster <- str_split_fixed(ss.embryo$CellType,":",2)[,1]
ss.embryo$CellType <- ss.embryo$CellType
ss.embryo$Batch <- ss.embryo$Stage
ss.embryo.sm <- subset(ss.embryo, SuperCluster %in% c("Somite"))
colnames(ss.embryo.sm@meta.data)
ss.embryo.sm@meta.data <- ss.embryo.sm@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Batch",
                                                     "percent.mt", "CellType", "SuperCluster")]
count <- GetAssayData(ss.embryo.sm, slot = "count") %>% as.data.frame()
meta <- ss.embryo.sm@meta.data
homo.gene <- read.csv("/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2pig_final.csv")
homo.gene <- homo.gene[!duplicated(homo.gene$gene.name), ]
rownames(homo.gene) <- homo.gene$gene.name
keep.gene <- intersect(rownames(homo.gene), rownames(count))
homo.gene <- homo.gene[keep.gene, ]
count <- count[keep.gene, ]
if (all(rownames(count) == rownames(homo.gene))) {
  rownames(count) <- homo.gene$human.gene.name
}
ss.embryo.sm <- CreateSeuratObject(counts = count, meta.data = meta)
ss.embryo.sm <- ss.embryo.sm[rownames(ss.embryo.sm) %in% rownames(ra.embryo.pop), ]
ss.embryo.sm@meta.data$Species <- "Pig"
rm(count, homo.gene, keep.gene, meta, ss.embryo)


### >>> 6. Integrate all species
ra.embryo.pop[['prediction.score.celltype']] <- "NULL"
sr.subpop <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)), 
                    list(hs.embryo.pop, mm.embryo.pop, ss.embryo.pop, ra.embryo.pop))
count <- GetAssayData(sr.subpop, slot = "count")
meta <- sr.subpop@meta.data
sr.subpop <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, meta)
table(sr.subpop$Species)
sr.subpop <- NormalizeData(sr.subpop) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(verbose = FALSE)
sr.subpop <- RunPCA(sr.subpop, npcs = 30, verbose = FALSE)
ElbowPlot(sr.subpop, ndims = 30)
dim.n <- 7
sr.subpop <- RunUMAP(sr.subpop, reduction = "pca", dims = 1:dim.n)
DimPlot(sr.subpop, group.by = "Species")
sr.list <- SplitObject(sr.subpop, split.by = "Species")
for (i in 1:length(sr.list)) {
  sr.list[[i]] <- NormalizeData(sr.list[[i]], verbose = FALSE)
  sr.list[[i]] <- FindVariableFeatures(sr.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sr.list, dims = 1:20)
sr.subpop <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(sr.subpop) <- "integrated"
sr.subpop <- ScaleData(sr.subpop, verbose = FALSE)
sr.subpop <- RunPCA(sr.subpop, npcs = 30, verbose = FALSE)
ElbowPlot(sr.subpop, ndims = 30)
sr.subpop <- RunUMAP(sr.subpop, reduction = "pca", dims = 1:10)
sr.subpop <- FindNeighbors(sr.subpop, dims = 1:10) %>%
  FindClusters(resolution = seq(0, 4, 0.25))
# plot clustering
pdf(file.path(res.out, "Integrative_clustering_of_progenitor_and_neuron_after_removing_batch_effect.pdf"), height = 6, width = 12)
CellDimPlot(srt = sr.subpop, group.by = c("SuperCluster", "Species"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "Feature_to_show_FP_markers_in_integrative_clustering_of_progenitor_and_neuron.pdf"), height = 10, width = 10)
FeatureDimPlot(srt = sr.subpop, features = c("FOXA2", "SPON1", "ARX", "SHH"), reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "Feature_to_show_RP_markers_in_integrative_clustering_of_progenitor_and_neuron.pdf"), height = 10, width = 10)
FeatureDimPlot(srt = sr.subpop, features = c("LMX1A", "RSPO3", "BMP6"), reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "Feature_to_show_SHH_markers_in_integrative_clustering_of_progenitor_and_neuron.pdf"), height = 10, width = 10)
FeatureDimPlot(srt = sr.subpop, features = c("PTCH1", "GLI2", "GLI3"), reduction = "UMAP", theme_use = "theme_blank")
dev.off()
FeatureDimPlot(srt = sr.subpop, features = c("ROBO3"), reduction = "UMAP", theme_use = "theme_blank")
pdf(file.path(res.out, "Integrative_clustering_res.3_of_progenitor_and_neuron_after_removing_batch_effect.pdf"), height = 6, width = 12)
CellDimPlot(srt = sr.subpop, group.by = c("integrated_snn_res.3", "Species"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()
# find markers
Idents(sr.subpop) <- sr.subpop$integrated_snn_res.3
marker.genes <- MarkerGene(sr.obj = sr.subpop, method = "cosg", top.n = 100,
                           assay = "integrated", slot = "data", idents = "integrated_snn_res.3", 
                           cols.pal = "Purples", cols.rev = F, 
                           plt.min = 0, plt.max = 0.3,
                           file.name = "Integrative_clustering", 
                           res.out = file.path(res.out, "Marker_gene_of_clusters"))
pdf(file.path(res.out, "Integrative_clustering_res.3_labelled_of_progenitor_and_neuron_after_removing_batch_effect.pdf"), height = 8, width = 8)
CellDimPlot(srt = sr.subpop, group.by = c("integrated_snn_res.3"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label_insitu = T, label = T)
dev.off()
pdf(file.path(res.out, "Stacked_barplots_to_show_cell_component_integrative_clustering_res.3_of_progenitor_and_neuron.pdf"), height = 5, width = 12)
CellStatPlot(sr.subpop, stat.by = "Species", group.by = "integrated_snn_res.3", stat_type = "percent", plot_type = "trend")
dev.off()
FeatureDimPlot(srt = sr.subpop, features = c("LGR5"), reduction = "UMAP", theme_use = "theme_blank")
GroupHeatmap(sr.subpop, assay = "integrated",slot = "data", 
             features = c("RORB", "LHX9", "LHX8", "LHX4", "LHX1", "LHX5"), 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
FeatureDimPlot(srt = sr.subpop, features = c("RORB", "LHX9", "LHX8", "LHX4", "LHX1", "LHX5"), reduction = "UMAP", theme_use = "theme_blank")
GroupHeatmap(sr.subpop, assay = "integrated",slot = "data", 
             features = marker.genes$markers$names$`23`[1:30], 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
# load genesets
mgi.genes <- list()
mgi.genes$Abnormal.Gait <- read.table("/home/yhw/document/genesets/MGI/MP0001406_abnormal_gait.txt", sep = "\t") %>% 
  select(V1) %>% as.vector() %>% unlist() %>% unique()
mgi.genes$Impaired.Limb.Coordination <- read.table("/home/yhw/document/genesets/MGI/MP0001524_impaired_limb_coordination.txt", sep = "\t") %>% 
  select(V1) %>% as.vector() %>% unlist() %>% unique()
mgi.genes$RORB.net <- c("GRIA2", "ANKRD44", "CHAT", "DENND1B", "GAL3ST4", "CHST8", "SLC17A6", "RBFOX3", "KCNJ12", 
                        "ENOX1", "GAS2", "MEIS1", "ZFP90", "DCC", "ZBTB38", "SORBS2", "GABBR2", "ACSL4", "EPB41L1", 
                        "JAK1", "C2CD4C", "ZFHX3", "KCTD4", "NR2F1", "CACNA1D", "ONECUT2", "KDM4B")
DefaultAssay(sr.subpop) <- "RNA"
sr.subpop <- ScoreGeneset(sr.obj = sr.subpop, gene.set = mgi.genes)
sr.subpop$RORB <- GetAssayData(sr.subpop, slot = "data", assay = "integrated")["RORB", ]
sr.subpop$SHH <- GetAssayData(sr.subpop, slot = "data", assay = "RNA")["SHH", ]
sr.subpop$PTCH1 <- GetAssayData(sr.subpop, slot = "data", assay = "RNA")["PTCH1", ]
sr.subpop$Abnormal.Gait <- GetAssayData(sr.subpop, slot = "data", assay = "UCell")["Abnormal.Gait", ]
sr.subpop$Impaired.Limb.Coordination <- GetAssayData(sr.subpop, slot = "data", assay = "UCell")["Impaired.Limb.Coordination", ]
sr.subpop$RORB.net <- GetAssayData(sr.subpop, slot = "data", assay = "UCell")["RORB.net", ]
pdf(file.path(res.out, "Heatmap_to_show_gene_expression_and_pathway_activity_split_by_species.pdf"), height = 4, width = 8)
GroupHeatmap(sr.subpop, assay = "ModuleScore", slot = "data", 
             features = c("SHH", "PTCH1", "RORB", names(mgi.genes)), group.by = c("Species"), split.by = "SuperCluster", 
             show_row_names = T, cluster_rows = F)
dev.off()
GroupHeatmap(sr.subpop, assay = "AUCell", slot = "data", 
             features = c("SHH", "PTCH1", "RORB", names(mgi.genes)), group.by = c("Species"), split.by = "SuperCluster", 
             show_row_names = T, cluster_rows = F)
GroupHeatmap(sr.subpop, assay = "UCell", slot = "data", 
             features = c("SHH", "PTCH1", "RORB", names(mgi.genes)), group.by = c("Species"), split.by = "SuperCluster", 
             show_row_names = T, cluster_rows = F)
GroupHeatmap(sr.subpop, assay = "JAS.likelihood", slot = "data", 
             features = c("SHH", "PTCH1", "RORB", names(mgi.genes)), group.by = c("Species"), split.by = "SuperCluster", 
             show_row_names = T, cluster_rows = F)
GroupHeatmap(sr.subpop, assay = "JAS.oddsratio", slot = "data", 
             features = c("SHH", "PTCH1", "RORB", names(mgi.genes)), group.by = c("Species"), split.by = "SuperCluster", 
             show_row_names = T, cluster_rows = F)
sr.subpop@meta.data %>% 
  filter(SuperCluster == "Neuron") %>% 
  ggplot(aes(x = RORB, y = RORB.net)) +
  geom_point(aes(group = Species)) +
  geom_smooth(aes(group = Species), method = "lm", se = TRUE) +
  facet_wrap(.~SuperCluster)
  theme_bw()
# RORB+ cells
sr.subpop@meta.data <- sr.subpop@meta.data %>% 
  mutate(RORB.type = case_when(RORB > 0.25 ~ "RORB.pos",
                               RORB <= 0.25 ~ "RORB.neg"))
pdf(file.path(res.out, "Stacked_barplot_to_show_cell_groups_composition_RORB_split_by_species.pdf"), height = 4, width = 8)
CellStatPlot(sr.subpop, stat.by = "RORB.type", group.by = "Species", split.by = "SuperCluster", 
             label = TRUE, plot_type = "trend", stat_type = "percent")
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_cell_ratio_without_replicates.pdf"), height = 6, width = 10)
pd <- sr.subpop@meta.data %>% 
  filter(CellType.sub %in% c("roof plate", "floor plate", 
                             "pA1", "pA2", "pA3", "pB1",
                             "dA1", "dA2", "dA3", "dB3", "dB4",
                             "p0.rhombomere", "p1.rhombomere", "p2.rhombomere", 
                             "V0.rhombomere", "V1.rhombomere", "V2b.rhombomere",
                             "dp1","dp2","dp3","dp4","dp5","dp6", 
                             "dI1","dI2","dI4","dI5","dI6", 
                             "p0","p1","p2","pMN","p3", 
                             "v0","v1","v2a","v2b","v3"))
p <- RatioPlot(meta.data = subset(sr.subpop@meta.data, Species == "Rabbit"), group = "RORB.type", rep = NULL,
               cell = "CellType.sub", lfc = log2(1.5), pt.col = NULL, pt.size = 5)
p$plot
dev.off()
# somite
ra.embryo.sm[['prediction.score.celltype']] <- "NULL"
sr.subpop.sm <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)), 
                       list(hs.embryo.sm, mm.embryo.sm, ss.embryo.sm, ra.embryo.sm))
count <- GetAssayData(sr.subpop.sm, slot = "count")
meta <- sr.subpop.sm@meta.data
sr.subpop.sm <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, meta)
table(sr.subpop.sm$Species, sr.subpop.sm$CellType)


### >>> 6. Integrate all species (spinal cord)
sr.subpop$CellType.sub <- str_split_fixed(sr.subpop$CellType,":",2)[, 2]
sr.subpop.sc <- subset(sr.subpop, CellType.sub %in% c("dp3", "dp4", "dp5", 
                                                      "dI1", "dI2", "dI4", "dILA/dBLa",
                                                      "p3", 
                                                      "v1", "intermediate V2 precursor", "v2a", "v2b", 
                                                      "lateral motor columns", "medial motor columns"))
DefaultAssay(sr.subpop.sc) <- "RNA"
sr.list <- SplitObject(sr.subpop.sc, split.by = "Species")
for (i in 1:length(sr.list)) {
  sr.list[[i]] <- NormalizeData(sr.list[[i]], verbose = FALSE)
  sr.list[[i]] <- FindVariableFeatures(sr.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sr.list, dims = 1:20)
sr.subpop.sc <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(sr.subpop.sc) <- "integrated"
sr.subpop.sc <- ScaleData(sr.subpop.sc, verbose = FALSE)
sr.subpop.sc <- RunPCA(sr.subpop.sc, npcs = 30, verbose = FALSE)
ElbowPlot(sr.subpop.sc, ndims = 30)
sr.subpop.sc <- RunUMAP(sr.subpop.sc, reduction = "pca", dims = 1:10)
sr.subpop.sc <- FindNeighbors(sr.subpop.sc, dims = 1:10) %>%
  FindClusters(resolution = seq(0, 4, 0.25))
table(sr.subpop.sc$Species)
# plot clustering
pdf(file.path(res.out, "Integrative_clustering_of_spinal_cord_progenitor_and_neuron_after_removing_batch_effect.pdf"), height = 6, width = 12)
CellDimPlot(srt = sr.subpop.sc, group.by = c("SuperCluster", "Species"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "Integrative_clustering_res.3_of_spinal_cord_progenitor_and_neuron_after_removing_batch_effect.pdf"), height = 6, width = 12)
CellDimPlot(srt = sr.subpop.sc, group.by = c("integrated_snn_res.3", "Species"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label = T, label_insitu = T)
dev.off()
pdf(file.path(res.out, "Integrative_clustering_res.3_of_spinal_cord_progenitor_and_neuron_after_annotation.pdf"), height = 6, width = 12)
CellDimPlot(srt = sr.subpop.sc, group.by = c("integrated_snn_res.3", "CellType.sub"), pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label = T, label_insitu = T, label_repel = T)
dev.off()
# find markers
Idents(sr.subpop.sc) <- sr.subpop.sc$integrated_snn_res.3
marker.genes.sc <- MarkerGene(sr.obj = sr.subpop.sc, method = "cosg", top.n = 100,
                              assay = "integrated", slot = "data", idents = "integrated_snn_res.3", 
                              cols.pal = "Purples", cols.rev = F, 
                              plt.min = 0, plt.max = 0.3,
                              file.name = "Integrative_clustering", 
                              res.out = file.path(res.out, "Marker_gene_of_clusters_in_spinal_cord"))
pdf(file.path(res.out, "Stacked_barplots_to_show_cell_component_integrative_clustering_res.3_of_spinal_cord_progenitor_and_neuron.pdf"), height = 5, width = 12)
CellStatPlot(sr.subpop.sc, stat.by = "Species", group.by = "integrated_snn_res.3", stat_type = "percent", plot_type = "trend")
dev.off()
FeatureDimPlot(srt = sr.subpop.sc, features = c("PTCH1", "GLI2", "GLI3"), reduction = "UMAP", theme_use = "theme_blank")
FeatureDimPlot(srt = sr.subpop.sc, features = c("FOXA2", "SPON1", "ARX", "SHH"), reduction = "UMAP", theme_use = "theme_blank", assay = "RNA")
FeatureDimPlot(srt = sr.subpop.sc, features = c("SHH"), split.by = "Species", reduction = "UMAP", theme_use = "theme_blank", assay = "RNA")
GroupHeatmap(sr.subpop.sc, assay = "integrated",slot = "data", 
             features = c("RORB", "LHX9", "LHX8", "LHX4", "LHX1", "LHX5", "SHH", "PTCH1", "CNTN1", "ETV4"), 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
FeatureDimPlot(srt = sr.subpop.sc, features = c("RORB", "LHX9", "LHX8", "LHX4", "LHX1", "LHX5"), assay = "RNA", reduction = "UMAP", theme_use = "theme_blank")
GroupHeatmap(sr.subpop.sc, assay = "integrated",slot = "data", 
             features = marker.genes.sc$markers$names$`3`[1:50], 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
# load genesets
mgi.genes <- list()
mgi.genes$Abnormal.Gait <- read.table("/home/yhw/document/genesets/MGI/MP0001406_abnormal_gait.txt", sep = "\t") %>% 
  select(V1) %>% as.vector() %>% unlist() %>% unique()
mgi.genes$Impaired.Limb.Coordination <- read.table("/home/yhw/document/genesets/MGI/MP0001524_impaired_limb_coordination.txt", sep = "\t") %>% 
  select(V1) %>% as.vector() %>% unlist() %>% unique()
GroupHeatmap(sr.subpop.sc, assay = "integrated",slot = "data", 
             features = mgi.genes$Abnormal.Gait$V1, 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
# 
DefaultAssay(sr.subpop.sc) <- "RNA"
sr.subpop.sc <- ScoreGeneset(sr.obj = sr.subpop.sc, gene.set = mgi.genes)
GroupHeatmap(sr.subpop.sc, assay = "ModuleScore", slot = "data", 
             features = names(mgi.genes), 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
FeatureStatPlot(sr.subpop.sc, stat.by = names(mgi.genes)[1], assay = "AUCell", slot = "data", 
                group.by = "integrated_snn_res.3", split.by = "Species", plot_type = "box")
FeatureStatPlot(sr.subpop.sc, stat.by = "RORB", assay = "integrated", slot = "data", 
                group.by = "Species", plot_type = "box")
for (i in names(marker.genes.sc$markers$names)) {
  print(i)
  print(intersect(marker.genes.sc$markers$names[[i]], mgi.genes$Abnormal.Gait))
}
for (i in names(marker.genes.sc$markers$names)) {
  print(i)
  print(intersect(marker.genes.sc$markers$names[[i]], mgi.genes$Impaired.Limb.Coordination))
}
GroupHeatmap(sr.subpop.sc, assay = "integrated", slot = "data", 
             features = unique(mgi.genes$Abnormal.Gait), 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
GroupHeatmap(sr.subpop.sc, assay = "integrated", slot = "data", 
             features = unique(c(intersect(marker.genes.sc$markers$names$`3`, mgi.genes$Abnormal.Gait),
                                 intersect(marker.genes.sc$markers$names$`4`, mgi.genes$Abnormal.Gait),
                                 intersect(marker.genes.sc$markers$names$`9`, mgi.genes$Abnormal.Gait),
                                 intersect(marker.genes.sc$markers$names$`13`, mgi.genes$Abnormal.Gait),
                                 intersect(marker.genes.sc$markers$names$`32`, mgi.genes$Abnormal.Gait),
                                 intersect(marker.genes.sc$markers$names$`42`, mgi.genes$Abnormal.Gait), "RORB")), 
             group.by = c("integrated_snn_res.3"), split.by = "Species", 
             show_row_names = T, cluster_rows = T)
# RORB+ cells
sr.subpop.sc$RORB <- GetAssayData(sr.subpop.sc, slot = "count", assay = "RNA")["RORB", ]
sr.subpop.sc@meta.data <- sr.subpop.sc@meta.data %>% 
  mutate(RORB.type = case_when(RORB > 0 ~ "RORB.pos",
                               RORB <= 0 ~ "RORB.neg"))
table(sr.subpop.sc$RORB.type, sr.subpop.sc$Species)
pdf(file.path(outdir, "Hind_brain_all_cell_groups_composition_RORB_split_by_species.pdf"), height = 5, width = 4)
CellStatPlot(sr.subpop.sc, stat.by = "RORB.type", group.by = "Species", label = TRUE, plot_type = "trend")
dev.off()
pdf(file.path(outdir, "Hind_brain_all_cell_groups_composition_RORB_split_by_none.pdf"), height = 5, width = 5)
CellStatPlot(sr.subpop.sc, stat.by = "RORB.type", group.by = "integrated_snn_res.3", label = TRUE, plot_type = "trend", split.by = "Species")
dev.off()
pdf(file.path(outdir, "Hind_brain_all_cell_groups_composition_RORB_split_by_supercluster.pdf"), height = 5, width = 8)
CellStatPlot(pd, stat.by = "RORB.type", group.by = "Species", label = TRUE, plot_type = "trend", split.by = "SuperCluster")
dev.off()
pdf(file.path(outdir, "Hind_brain_all_cell_groups_composition_RORB_split_by_stage.pdf"), height = 5, width = 10)
CellStatPlot(pd, stat.by = "RORB.type", group.by = "CellType.sub2", label = TRUE, plot_type = "trend", split.by = "Species")
dev.off()



# ==============================
# 5th part: Spinal cord analysis ----
# ==============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Rabbit_Embryos/Spinal_Cord")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Human spinal cord analysis
library("dior")
hs.spinal <- dior::read_h5(file = '/home/yhw/document/public_data/Spinal_Cord/wxq_data/all_rna.h5', target.object = 'seurat')
hs.spinal <- subset(hs.spinal, cell_type_final_sub %in% c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                                          "dp1", "dp2", "dp3-6",
                                                          "p0-p2", "p3", "pMN", 
                                                          "v0", "v1", "v2a", "v2b", "v3", "MN"))
table(hs.spinal$sample_fromcell_type_final_sub)
Idents(hs.spinal) <- hs.spinal$cell_type_final_sub
hs.spinal <- subset(hs.spinal, downsample = 1000)
hs.spinal <- hs.spinal[rownames(hs.spinal) %in% rownames(sr.subpop), ]
hs.spinal$Species <- "Human.WXQ"
hs.spinal$CellType <- hs.spinal$cell_type_final_sub


### >>> 3. Annotate cell types by Seurat
hs.spinal <- NormalizeData(hs.spinal, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)
ElbowPlot(hs.spinal, ndims = 30)
dims <- 20
hs.spinal <- RunUMAP(hs.spinal, dims = 1:dims, reduction = "pca", return.model = T)
anchors <- FindTransferAnchors(reference = hs.spinal, query = sr.subpop, dims = 1:25, reference.reduction = "pca")
map <- MapQuery(anchorset = anchors, reference = hs.spinal, query = sr.subpop.sc,
                refdata = list(celltype = "cell_type_final_sub"), reference.reduction = "pca", reduction.model = "umap")
map.list <- Annotate.Cell(ref = hs.spinal, query = sr.subpop, gene.set = NULL, 
                          group.by = "CellType", dim.n = 25, weight.n = 50)
map.list$seurat <- map.list$seurat@meta.data
map.list$scpred <- map.list$scpred@meta.data
map.list$singler <- map.list$singler@meta.data
table(map.list$scpred$Batch, map.list$scpred$scpred_prediction)
table(map.list$singler$labels)
table(map.list$seurat$predicted.id)
if (all(colnames(map) == rownames(map.list$scpred))) {
  # add annotation
  map$scpred_prediction <- map.list$scpred$scpred_prediction
  # filter cell
  meta <- map.list$scpred %>%
    rownames_to_column(var = "Barcode") %>%
    dplyr::filter(scpred_prediction != "unassigned") %>% 
    dplyr::group_by(Batch, scpred_prediction) %>%
    dplyr::top_n(n = 500, wt = scpred_max)
  table(meta$Batch, meta$scpred_prediction)
  # subset cells
  map.sub <- map[, colnames(map) %in% meta$Barcode]
}
table(map.sub$Batch, map.sub$scpred_prediction)
map.sub@meta.data <- map.sub@meta.data %>% 
  mutate(Stage = case_when(Batch %in% c("Hs.5W.1", "Hs.5W.2") ~ "Hs-5W",
                           Batch %in% c("Hs.6W.1") ~ "Hs-6W",
                           Batch %in% c("E9.5") ~ "Mm-E9.5",
                           Batch %in% c("E10.5") ~ "Mm-E10.5",
                           Batch %in% c("E11.5") ~ "Mm-E11.5",
                           Batch %in% c("E12.5") ~ "Mm-E12.5",
                           Batch %in% c("Day18") ~ "Pig-D18",
                           Batch %in% c("Day22") ~ "Pig-D22",
                           Batch %in% c("E10-1", "E10-2") ~ "Ra-E10",
                           Batch %in% c("E12-1", "E12-2") ~ "Ra-E12",
                           Batch %in% c("E14-1", "E14-2") ~ "Ra-E14"),
         SuperCluster = case_when(scpred_prediction %in% c("dp1", "dp2", "dp3-6", "p0-p2", "pMN", "p3") ~ "Progenitor",
                                  scpred_prediction %in% c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6", 
                                                           "v0", "v1", "v2a", "v2b", "MN", "v3") ~ "Neuron"))
map.sub$scpred_prediction <- factor(map.sub$scpred_prediction, 
                                    levels = c("dp1", "dp2", "dp3-6", "dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                               "p0-p2", "pMN", "p3", "v0", "v1", "v2a", "v2b", "MN", "v3"))
map.sub$Stage <- factor(map.sub$Stage, levels = c("Hs-5W","Hs-6W","Mm-E9.5","Mm-E10.5","Mm-E11.5","Mm-E12.5",
                                                  "Pig-D18","Pig-D22","Ra-E10","Ra-E12","Ra-E14"))
# plot cell type
pdf(file.path(res.out, "Spinal_cord_cell_clustering_after_integration.pdf"), height = 6, width = 7)
CellDimPlot(srt = map.sub, group.by = "scpred_prediction", pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label = T, label_insitu = T, label_repel = T)
dev.off()
# plot cell composition
pdf(file.path(res.out, "Spinal_cord_cell_groups_composition_after_annotation.pdf"), height = 8, width = 14)
CellStatPlot(map.sub, stat.by = "scpred_prediction", group.by = "Stage", split.by = "SuperCluster",
             label = TRUE, plot_type = "trend", stat_type = "percent")
dev.off()
table(map$Species, map$Batch, map$SuperCluster)
table(map.sub$Species, map.sub$scpred_prediction)


### >>> 4. Integration again
table(map.sub$Species, map.sub$SuperCluster)
sr.list <- SplitObject(map.sub, split.by = "Species")
for (i in 1:length(sr.list)) {
  sr.list[[i]] <- NormalizeData(sr.list[[i]], verbose = FALSE)
  sr.list[[i]] <- FindVariableFeatures(sr.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sr.list, dims = 1:20)
map.inte <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(map.inte) <- "integrated"
map.inte <- ScaleData(map.inte, verbose = FALSE)
map.inte <- RunPCA(map.inte, npcs = 30, verbose = FALSE)
ElbowPlot(map.inte, ndims = 30)
map.inte <- RunUMAP(map.inte, reduction = "pca", dims = 1:10)
map.inte$scpred_prediction <- factor(map.inte$scpred_prediction, 
                                     levels = c("dp1", "dp2", "dp3-6", "dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                                "p0-p2", "pMN", "p3", "v0", "v1", "v2a", "v2b", "MN", "v3"))
map.inte$Stage <- factor(map.inte$Stage, levels = c("Hs-5W","Hs-6W","Mm-E9.5","Mm-E10.5","Mm-E11.5","Mm-E12.5",
                                                    "Pig-D18","Pig-D22","Ra-E10","Ra-E12","Ra-E14"))
pdf(file.path(res.out, "Spinal_cord_cell_clustering_after_integration.pdf"), height = 6, width = 7)
CellDimPlot(srt = map.inte, group.by = "scpred_prediction", pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label = T, label_insitu = T, label_repel = T)
dev.off()
pdf(file.path(res.out, "Spinal_cord_cell_groups_composition_after_integration.pdf"), height = 8, width = 14)
CellStatPlot(map.sub, stat.by = "scpred_prediction", group.by = "Stage", split.by = "SuperCluster",
             label = TRUE, plot_type = "trend", stat_type = "percent")
dev.off()
# 
table(map.sub$Stage, map.sub$predicted.celltype) %>% 
  as.data.frame() %>% 
  spread(key = "Var1", value = "Freq") %>% 
  column_to_rownames(var = "Var2") -> pd
pd[pd == 0] <- 1
library("robCompositions")
pd.dist <- matrix(NA, nrow = ncol(pd), ncol = ncol(pd))
rownames(pd.dist) <- colnames(pd)
colnames(pd.dist) <- colnames(pd)
for (i in colnames(pd)) {
  for (j in colnames(pd)) {
    i.index <- grep(i, colnames(pd))
    j.index <- grep(j, colnames(pd))
    pd.dist[i.index, j.index] <- aDist(pd[, i], pd[, j])
  }
}

pd.dist <- max(pd.dist) - pd.dist
pd.dist[pd.dist > 8] <- 9
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
pdf(file.path(res.out, "Aitchison_distance_of_cell_composition_1.pdf"), height = 7, width = 7)
corrplot::corrplot(pd.dist, 
                   is.corr = F, col = cols(100),
                   col.lim = c(min(pd.dist), max(pd.dist)),
                   tl.pos = "lt", method = "circle", tl.col = "#000000", type = "full",
                   order = "hclust", addrect = 3, tl.cex = 0.75, pch.cex = 0.9)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()
pdf(file.path(res.out, "Aitchison_distance_of_cell_composition_2.pdf"), height = 7, width = 7)
corrplot::corrplot.mixed(pd.dist, 
                         is.corr = F, upper.col = cols(100), lower.col = cols(100),
                         col.lim = c(min(pd.dist), max(pd.dist)),
                         tl.pos = "lt", tl.col = "#000000", lower = 'number', upper = 'circle',
                         order = "hclust", addrect = 4, tl.cex = 0.75, pch.cex = 0.9)
dev.off()
# RORB+ cells
map.inte$RORB <- GetAssayData(map.inte, slot = "data", assay = "integrated")["RORB", ]
map.inte@meta.data <- map.inte@meta.data %>% 
  mutate(RORB.type = case_when(RORB > 0.5 ~ "RORB.pos",
                               RORB <= 0.5 ~ "RORB.neg"))
table(map.inte$RORB.type, map.inte$Species)
pdf(file.path(res.out, "Spinal_cord_cell_groups_composition_RORB_split_by_species.pdf"), height = 5, width = 4)
CellStatPlot(map.inte, stat.by = "RORB.type", group.by = "Species", label = TRUE, plot_type = "trend")
dev.off()
# Aitchison distance
rm(i, j, anchors, pd, pd.list, i.index, j.index)


### >>> 5. Molecular clock
table(map.inte$Species, map.inte$SuperCluster)
# - integrated assay
mfuzz.inte <- list()
DefaultAssay(map.inte) <- "integrated"
for (i in c("Neuron", "Progenitor")) {
  pd <- subset(map.inte, SuperCluster == i)
  pd.gene <- VariableFeatures(pd)
  pd <- AverageExpression(pd, assays = "integrated", group.by = "Stage") %>% as.data.frame()
  pd.meta <- data.frame(Group = gsub("integrated.", "", colnames(pd)),
                        Species = gsub("\\..*", "", gsub("integrated.", "", colnames(pd))),
                        row.names = colnames(pd))
  pd.sva <- CorrectBatch(expr = pd, meta = pd.meta, batch = "Species", low.expr = 0,
                         method = "sva", sva.mode = "data")
  for (j in seq(500, 2000, 250)) {
    p1 <- PCAplot(data = pd[intersect(names(rowSums(pd[HVG.Topn(pd, j), ])>0), rownames(pd)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (Seurat) with Top", j, " Genes"), pt.size = 4)
    p2 <- PCAplot(data = pd.sva[intersect(names(rowSums(pd.sva[HVG.Topn(pd.sva, j), ])>0), rownames(pd.sva)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (Seurat + SVA) with Top", j, " Genes"), pt.size = 4)
    p3 <- PCAplot(data = pd[intersect(names(rowSums(pd[pd.gene, ])>0), rownames(pd)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (Seurat) with Top HVGs"), pt.size = 4)
    p4 <- PCAplot(data = pd.sva[intersect(names(rowSums(pd.sva[pd.gene, ])>0), rownames(pd.sva)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (Seurat + SVA) with Top HVGs"), pt.size = 4)
    write.csv(do.call(rbind, list(p1$data, p2$data, p3$data, p4$data)),
              file.path(res.out, paste0("Integrated_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes.csv")))
    write.csv(pd.sva,
              file.path(res.out, paste0("Integrated_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes_expression_RNA_SVA.csv")))
    pdf(file.path(res.out, paste0("Integrated_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 14, width = 12)
    print(plot_grid(p1$plot + p2$plot + p3$plot + p4$plot, nrow = 2))
    dev.off()
    # Seurat
    mfuzz.inte[[paste0("Seurat_", i, "_", j)]] <- Pipe.Mfuzz(expr.mtx = as.matrix(pd[HVG.Topn(pd, j), grep("Mm", colnames(pd))]),
                                                             res.out = file.path(res.out, "Mfuzz/mouse"),
                                                             prefix = paste0(i, "_Seurat_with_Top_", j, "_Genes"),
                                                             cluster.num = 10, min.acore = 0.5)
    p.list <- list()
    for (k in 1:10) {
      p.list[[k]] <- pd %>%
        rownames_to_column(var = "Genes") %>%
        gather(key = "Group", value = "Value", -Genes) %>%
        mutate(Value = log1p(Value),
               Group = factor(Group, levels = colnames(pd))) %>%
        dplyr::filter(Genes %in% subset(mfuzz.rna[[paste0("RNA_", i, "_", j)]]$cluster, Cluster == k)$SYMBOL) %>%
        ggplot(aes(x = Group, y = Value)) +
        geom_boxplot(aes(fill = Group)) +
        scale_fill_brewer(palette = "Paired") +
        labs(x = paste0("Cluster", k)) +
        theme_bw() +
        theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
              axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
              axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              panel.grid = element_blank(), legend.position = "none")
    }
    pdf(file.path(res.out, paste0("Seurat_BoxPlot_Cluster_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 9, width = 18)
    print(plot_grid(plotlist = p.list, nrow = 2, ncol = 5))
    dev.off()
    # Seurat + SVA
    mfuzz.inte[[paste0("Seurat_SVA_", i, "_", j)]] <- Pipe.Mfuzz(expr.mtx = as.matrix(pd.sva[HVG.Topn(pd.sva, j), grep("Mm", colnames(pd.sva))]),
                                                                 res.out = file.path(res.out, "Mfuzz/mouse"),
                                                                 prefix = paste0(i, "_Seurat_SVA_with_Top_", j, "_Genes"),
                                                                 cluster.num = 10, min.acore = 0.5)
    p.list <- list()
    for (k in 1:10) {
      p.list[[k]] <- pd.sva %>%
        rownames_to_column(var = "Genes") %>%
        gather(key = "Group", value = "Value", -Genes) %>%
        mutate(Value = log1p(Value),
               Group = factor(Group, levels = colnames(pd.sva))) %>%
        dplyr::filter(Genes %in% subset(mfuzz.rna[[paste0("RNA_", i, "_", j)]]$cluster, Cluster == k)$SYMBOL) %>%
        ggplot(aes(x = Group, y = Value)) +
        geom_boxplot(aes(fill = Group)) +
        scale_fill_brewer(palette = "Paired") +
        labs(x = paste0("Cluster", k)) +
        theme_bw() +
        theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
              axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
              axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              panel.grid = element_blank(), legend.position = "none")
    }
    pdf(file.path(res.out, paste0("Seurat_SVA_BoxPlot_Cluster_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 9, width = 18)
    print(plot_grid(plotlist = p.list, nrow = 2, ncol = 5))
    dev.off()
  }
}
# - RNA assay
mfuzz.rna <- list()
DefaultAssay(map.inte) <- "RNA"
for (i in c("Neuron", "Progenitor")) {
  pd <- subset(map.inte, SuperCluster == i)
  pd <- NormalizeData(pd, verbose = FALSE) %>% 
    FindVariableFeatures(selection.method = "vst", 
                         nfeatures = 2000, verbose = FALSE)
  pd.gene <- VariableFeatures(pd)
  pd <- AverageExpression(pd, assays = "RNA", group.by = "Stage") %>% as.data.frame()
  pd.meta <- data.frame(Group = gsub("RNA.", "", colnames(pd)),
                        Species = gsub("\\..*", "", gsub("RNA.", "", colnames(pd))),
                        row.names = colnames(pd))
  pd.sva <- CorrectBatch(expr = pd, meta = pd.meta, batch = "Species", low.expr = 0,
                         method = "sva", sva.mode = "data")
  for (j in seq(500, 2000, 250)) {
    p1 <- PCAplot(data = pd[intersect(names(rowSums(pd[HVG.Topn(pd, j), ])>0), rownames(pd)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (RNA) with Top", j, " Genes"), pt.size = 4)
    p1$data$Title <- paste0("2D-Plot PCA of ", i, " (RNA) with Top", j, " Genes")
    p2 <- PCAplot(data = pd.sva[intersect(names(rowSums(pd.sva[HVG.Topn(pd.sva, j), ])>0), rownames(pd.sva)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (RNA + SVA) with Top", j, " Genes"), pt.size = 4)
    p2$data$Title <- paste0("2D-Plot PCA of ", i, " (RNA + SVA) with Top", j, " Genes")
    p3 <- PCAplot(data = pd[intersect(names(rowSums(pd[pd.gene, ])>0), rownames(pd)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (RNA) with Top HVGs"), pt.size = 4)
    p3$data$Title <- paste0("2D-Plot PCA of ", i, " (RNA) with Top HVGs")
    p4 <- PCAplot(data = pd.sva[intersect(names(rowSums(pd.sva[pd.gene, ])>0), rownames(pd.sva)), ], 
                  label = pd.meta$Group, color.by = pd.meta$Species, 
                  title = paste0("2D-Plot PCA of ", i, " (RNA + SVA) with Top HVGs"), pt.size = 4)
    p4$data$Title <- paste0("2D-Plot PCA of ", i, " (RNA + SVA) with Top HVGs")
    pdf(file.path(res.out, paste0("RNA_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 14, width = 12)
    print(plot_grid(p1$plot + p2$plot + p3$plot + p4$plot, nrow = 2))
    dev.off()
    write.csv(do.call(rbind, list(p1$data, p2$data, p3$data, p4$data)),
              file.path(res.out, paste0("RNA_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes.csv")))
    write.csv(pd.sva,
              file.path(res.out, paste0("RNA_2D-Plot_PCA_of_", i, "_with_Top_", j, "_Genes_expression_RNA_SVA.csv")))
    # RNA
    mfuzz.rna[[paste0("RNA_", i, "_", j)]] <- Pipe.Mfuzz(expr.mtx = as.matrix(pd[HVG.Topn(pd, j), grep("Mm", colnames(pd))]),
                                                         res.out = file.path(res.out, "Mfuzz/mouse"),
                                                         prefix = paste0(i, "_RNA_with_Top_", j, "_Genes"),
                                                         cluster.num = 10, min.acore = 0.5)
    p.list <- list()
    for (k in 1:10) {
      p.list[[k]] <- pd %>%
        rownames_to_column(var = "Genes") %>%
        gather(key = "Group", value = "Value", -Genes) %>%
        mutate(Value = log1p(Value),
               Group = factor(Group, levels = colnames(pd))) %>%
        dplyr::filter(Genes %in% subset(mfuzz.rna[[paste0("RNA_", i, "_", j)]]$cluster, Cluster == k)$SYMBOL) %>%
        ggplot(aes(x = Group, y = Value)) +
        geom_boxplot(aes(fill = Group)) +
        scale_fill_brewer(palette = "Paired") +
        labs(x = paste0("Cluster", k)) +
        theme_bw() +
        theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
              axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
              axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              panel.grid = element_blank(), legend.position = "none")
    }
    pdf(file.path(res.out, paste0("RNA_BoxPlot_Cluster_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 9, width = 18)
    print(plot_grid(plotlist = p.list, nrow = 2, ncol = 5))
    dev.off()
    # RNA + SVA
    mfuzz.rna[[paste0("RNA_SVA_", i, "_", j)]] <- Pipe.Mfuzz(expr.mtx = as.matrix(pd.sva[HVG.Topn(pd.sva, j), grep("Mm", colnames(pd.sva))]),
                                                             res.out = file.path(res.out, "Mfuzz/mouse"),
                                                             prefix = paste0(i, "_RNA_SVA_with_Top_", j, "_Genes"),
                                                             cluster.num = 10, min.acore = 0.5)
    p.list <- list()
    for (k in 1:10) {
      p.list[[k]] <- pd.sva %>%
        rownames_to_column(var = "Genes") %>%
        gather(key = "Group", value = "Value", -Genes) %>%
        mutate(Value = log1p(Value),
               Group = factor(Group, levels = colnames(pd.sva))) %>%
        dplyr::filter(Genes %in% subset(mfuzz.rna[[paste0("RNA_", i, "_", j)]]$cluster, Cluster == k)$SYMBOL) %>%
        ggplot(aes(x = Group, y = Value)) +
        geom_boxplot(aes(fill = Group)) +
        scale_fill_brewer(palette = "Paired") +
        labs(x = paste0("Cluster", k)) +
        theme_bw() +
        theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
              axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
              axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              panel.grid = element_blank(), legend.position = "none")
    }
    pdf(file.path(res.out, paste0("RNA_SVA_BoxPlot_Cluster_of_", i, "_with_Top_", j, "_Genes.pdf")), height = 9, width = 18)
    print(plot_grid(plotlist = p.list, nrow = 2, ncol = 5))
    dev.off()
  }
}
rm(i, j, k, p1, p2, p3, p4, p.list, pd, pd.sva, pd.gene, pd.meta, pd.dist)
# Pseudotime analysis
library("ggridges")
library("ggplot2")
pd.files <- list.files("./R/Graphs/Rabbit_Embryos/Spinal_Cord/", pattern = "*RNA_SVA.csv", full.names = T)
stage.col <- c("#9ECAE1", "#6BAED6",
               "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D",
               "#FC9272", "#FB6A4A",
               "#BCBDDC", "#9E9AC8", "#807DBA")
names(stage.col) <- c("Hs.5W", "Hs.6W", "Mm.E9.5", "Mm.E10.5", "Mm.E11.5", "Mm.E12.5", "Pig.D18", "Pig.D22", "Ra.E10", "Ra.E12", "Ra.E14")
for (file in pd.files) {
  pd.sva <- read.csv(file, row.names = 1)
  results <- list()
  for (i in 1:50) {
    pd <- Pipe.CBPA(expr = pd.sva[sample(1:nrow(pd.sva), i*20), ],
                    growth.function = "linear", zscore = T, cor.method = "spearman")
    results[[i]] <- pd$sample_pseudotime
  }
  pdf(gsub(".csv", "_Pseudotime.pdf", file), height = 4, width = 7)
  print(Reduce(rbind, results) %>% 
          as.data.frame() %>% 
          gather(key = "Stage", value = "Pseudotime") %>% 
          mutate(Stage = gsub("RNA.", "", gsub("integrated.", "", Stage))) %>%  
          mutate(Stage = factor(Stage, levels = names(stage.col)),
                 Species = gsub("\\..*$", "", gsub("RNA.", "", gsub("integrated.", "", Stage)))) %>% 
          ggplot(aes(x = Pseudotime, y = Species, fill = Stage)) +
          geom_density_ridges(scale = 1, rel_min_height = 0) +
          scale_fill_manual(values = stage.col) +
          labs(title = "Weighted Correlation-Based Pseudotime Analysis") +
          xlab("Pseudotime (low ---> high)") +
          theme_bw() +
          theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
                axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
                panel.grid = element_blank(), legend.position = "right")
  )
  dev.off()
}
rm(i, file, pd.sva, results)
# 
pd.gene <- c("NKX6-1","FOXA2","SHH","SPON1", # Ventral.FP
             "LHX1","LHX5","EN1","OTP","LHX3","VSX2","SOX14",# Ventral neuron
             "LBX1","PAX2","LHX2","LHX9","BARHL2","BARHL1","POU4F1", # Dorsal neuron
             "NKX6-1","NKX6-2","NKX2-2","NKX2-8","OLIG1","OLIG2","PRDM8","IRX3","PAX6","DBX2" # Ventral progenitor
)
GroupHeatmap(srt = map.inte2, flip = F, limits = c(-1, 2),
             exp_method = "zscore", assay = "RNA", slot = "data",
             cluster_columns = T, cluster_rows = T,
             features = intersect(unique(pd.gene), rownames(map.inte)),
             group.by = c("Species"), 
             show_row_names = TRUE, row_names_side = "left",
             add_dot = F, add_reticle = F, dot_size = unit(5, "mm"))


### >>> 6. New pig data
sr.pig <- readRDS("/home/yhw/document/public_data/GSE206914_pig_somites_E16_18_21_28/R/CodeData/GSE206914_pig_somites_E16_18_21_28.rds")
table(sr.pig$CellType)
sr.pig <- subset(sr.pig, CellType %in% c("Floor plate", "Roof plate", "Neural progenitor", "Neuron"))
sr.pig <- sr.pig[rownames(sr.pig) %in% rownames(sr.subpop), ]
sr.pig <- NormalizeData(sr.pig, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)
table(sr.pig@meta.data$Stage)
map.list.pig <- Annotate.Cell(ref = hs.spinal, query = sr.pig, gene.set = NULL, 
                              group.by = "CellType", dim.n = 25, weight.n = 50)
map.list.pig$seurat <- map.list.pig$seurat@meta.data
map.list.pig$scpred <- map.list.pig$scpred@meta.data
table(map.list.pig$scpred$scpred_prediction)
if (all(rownames(map.list.pig$scpred) == colnames(sr.pig))) {
  sr.pig$CellType.old <- sr.pig$CellType
  sr.pig$CellType <- map.list.pig$scpred$scpred_prediction
  sr.pig$Stage <- paste0("Pig-", sr.pig$Stage)
}
table(sr.pig$Stage)
sr.pig$Species <- "Pig2"
table(sr.pig$CellType)
sr.pig <- subset(sr.pig, CellType != "unassigned")


### >>> 7. Integration again and agin
map.sub2 <- merge(map.sub, sr.pig)
table(map.sub2$Species)
sr.list <- SplitObject(map.sub2, split.by = "Species")
for (i in 1:length(sr.list)) {
  sr.list[[i]] <- NormalizeData(sr.list[[i]], verbose = FALSE)
  sr.list[[i]] <- FindVariableFeatures(sr.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sr.list, dims = 1:20)
map.inte2 <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(map.inte2) <- "integrated"
map.inte2 <- ScaleData(map.inte2, verbose = FALSE)
map.inte2 <- RunPCA(map.inte2, npcs = 30, verbose = FALSE)
ElbowPlot(map.inte2, ndims = 30)
map.inte2 <- RunUMAP(map.inte2, reduction = "pca", dims = 1:10)
map.inte2@meta.data <- map.inte2@meta.data %>% 
  mutate(scpred_prediction = case_when(Species %in% c("Human", "Mouse", "Pig", "Rabbit") ~ scpred_prediction,
                                       Species %in% c("Pig2") ~ CellType),
         SuperCluster = case_when(scpred_prediction %in% c("dp1", "dp2", "dp3-6", "p0-p2", "pMN", "p3") ~ "Progenitor",
                                  scpred_prediction %in% c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6", 
                                                           "v0", "v1", "v2a", "v2b", "MN", "v3") ~ "Neuron"))
table(map.inte2$Species, map.inte2$scpred_prediction)
map.inte2$scpred_prediction <- factor(map.inte2$scpred_prediction, 
                                     levels = c("dp1", "dp2", "dp3-6", "dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                                "p0-p2", "pMN", "p3", "v0", "v1", "v2a", "v2b", "MN", "v3"))
map.inte2$Stage <- factor(map.inte2$Stage, levels = c("Hs-5W","Hs-6W","Mm-E9.5","Mm-E10.5","Mm-E11.5","Mm-E12.5",
                                                    "Pig-D18","Pig-D22","Pig-E16","Pig-E18","Pig-E21","Pig-E28",
                                                    "Ra-E10","Ra-E12","Ra-E14"))
pdf(file.path(res.out, "Spinal_cord_cell_clustering_after_integration_2.pdf"), height = 6, width = 7)
CellDimPlot(srt = map.inte2, group.by = "scpred_prediction", pt.size = 0.1, 
            reduction = "UMAP", theme_use = "theme_blank", label = T, label_insitu = T, label_repel = T)
dev.off()
pdf(file.path(res.out, "Spinal_cord_cell_groups_composition_after_integration_2.pdf"), height = 8, width = 18)
CellStatPlot(map.inte2, stat.by = "scpred_prediction", group.by = "Stage", split.by = "SuperCluster",
             label = TRUE, plot_type = "trend", stat_type = "percent")
dev.off()
tmp.1 <- subset(map.inte2, Stage %in% c("Hs-6W","Mm-E11.5","Pig-E21","Ra-E12"))
tmp.1 <- subset(tmp.1, scpred_prediction %in% c("dp1", "dp2", "dp3-6", "p0-p2", "pMN", "p3"))
tmp.1$Stage <- factor(tmp.1$Stage, levels = c("Hs-6W","Mm-E11.5","Pig-E21","Ra-E12"))
tmp.1$scpred_prediction <- factor(tmp.1$scpred_prediction, levels = c("dp1", "dp2", "dp3-6", "p0-p2", "pMN", "p3"))
p1 <- CellStatPlot(tmp.1, stat.by = "scpred_prediction", group.by = "Stage", split.by = "SuperCluster",
                   label = TRUE, plot_type = "trend", stat_type = "percent")
tmp.2 <- subset(map.inte2, Stage %in% c("Hs-6W","Mm-E11.5","Pig-E21","Ra-E12"))
tmp.2 <- subset(tmp.2, scpred_prediction %in% c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                                "v0", "v1", "v2a", "v2b", "MN", "v3"))
tmp.2$Stage <- factor(tmp.2$Stage, levels = c("Hs-6W","Mm-E11.5","Pig-E21","Ra-E12"))
tmp.2$scpred_prediction <- factor(tmp.2$scpred_prediction, 
                                  levels = c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6",
                                             "v0", "v1", "v2a", "v2b", "MN", "v3"))
p2 <- CellStatPlot(tmp.2, stat.by = "scpred_prediction", group.by = "Stage", split.by = "SuperCluster",
                   label = TRUE, plot_type = "trend", stat_type = "percent")
pdf(file.path(res.out, "Spinal_cord_cell_groups_composition_after_integration_3.pdf"), height = 8, width = 10)
p1 + p2
dev.off()
# 
table(map.inte2$Stage, map.inte2$predicted.celltype) %>% 
  as.data.frame() %>% 
  spread(key = "Var1", value = "Freq") %>% 
  column_to_rownames(var = "Var2") -> pd
pd[pd == 0] <- 1
library("robCompositions")
pd.dist <- matrix(NA, nrow = ncol(pd), ncol = ncol(pd))
rownames(pd.dist) <- colnames(pd)
colnames(pd.dist) <- colnames(pd)
for (i in colnames(pd)) {
  for (j in colnames(pd)) {
    i.index <- grep(i, colnames(pd))
    j.index <- grep(j, colnames(pd))
    pd.dist[i.index, j.index] <- aDist(pd[, i], pd[, j])
  }
}
pd.dist <- max(pd.dist) - pd.dist
pd.dist[pd.dist > 8] <- 9
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
pdf(file.path(res.out, "Aitchison_distance_of_cell_composition_1_more.pdf"), height = 7, width = 7)
corrplot::corrplot(pd.dist, 
                   is.corr = F, col = cols(100),
                   col.lim = c(min(pd.dist), max(pd.dist)),
                   tl.pos = "lt", method = "circle", tl.col = "#000000", type = "full",
                   order = "hclust", addrect = 3, tl.cex = 0.75, pch.cex = 0.9)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()
pdf(file.path(res.out, "Aitchison_distance_of_cell_composition_2_more.pdf"), height = 7, width = 7)
corrplot::corrplot.mixed(pd.dist, 
                         is.corr = F, upper.col = cols(100), lower.col = cols(100),
                         col.lim = c(min(pd.dist), max(pd.dist)),
                         tl.pos = "lt", tl.col = "#000000", lower = 'number', upper = 'circle',
                         order = "hclust", addrect = 4, tl.cex = 0.75, pch.cex = 0.9)
dev.off()


### >>> 8. Differential expression analysis
map.inte@meta.data <- map.inte@meta.data %>% 
  mutate(Tissue = case_when(scpred_prediction %in% c("dp1", "dp2", "dp3-6") ~ "Dorsal-Progenitor",
                            scpred_prediction %in% c("p0-p2", "pMN", "p3") ~ "Ventral-Progenitor",
                            scpred_prediction %in% c("dl1", "dl2", "dl3", "dl4", "dl5", "dl6") ~ "Dorsal-Neuron",
                            scpred_prediction %in% c("v0", "v1", "v2a", "v2b", "MN", "v3") ~ "Ventral-Neuron")) %>% 
  mutate(DEG.group = paste0(Stage, ".", Tissue))
table(map.inte@meta.data$DEG.group)
table(sr.subpop.sm$Species, sr.subpop.sm$SuperCluster)
sr.subpop.sm$DEG.group <- paste0(sr.subpop.sm$Species, ".", sr.subpop.sm$SuperCluster)
table(sr.subpop.sm$DEG.group)
# DEA for comparison between cell types
deg <- list()
deg.meta <- data.frame(g1 = c("Mm-E9.5.Ventral-Neuron", "Mm-E9.5.Ventral-Progenitor", 
                              "Mm-E10.5.Ventral-Neuron", "Mm-E10.5.Ventral-Progenitor",
                              "Mm-E11.5.Ventral-Neuron", "Mm-E11.5.Ventral-Progenitor",
                              "Mm-E12.5.Ventral-Neuron", "Mm-E12.5.Ventral-Progenitor",
                              "Mm-E9.5.Ventral-Neuron", "Mm-E9.5.Ventral-Progenitor", 
                              "Mm-E10.5.Ventral-Neuron", "Mm-E10.5.Ventral-Progenitor",
                              "Mm-E11.5.Ventral-Neuron", "Mm-E11.5.Ventral-Progenitor",
                              "Mm-E12.5.Ventral-Neuron", "Mm-E12.5.Ventral-Progenitor",
                              "Mm-E9.5.Ventral-Neuron", "Mm-E9.5.Ventral-Progenitor", 
                              "Mm-E10.5.Ventral-Neuron", "Mm-E10.5.Ventral-Progenitor",
                              "Mm-E11.5.Ventral-Neuron", "Mm-E11.5.Ventral-Progenitor",
                              "Mm-E12.5.Ventral-Neuron", "Mm-E12.5.Ventral-Progenitor",
                              "Mm-E9.5.Ventral-Neuron", "Mm-E9.5.Ventral-Progenitor", 
                              "Mm-E10.5.Ventral-Neuron", "Mm-E10.5.Ventral-Progenitor",
                              "Mm-E11.5.Ventral-Neuron", "Mm-E11.5.Ventral-Progenitor",
                              "Mm-E12.5.Ventral-Neuron", "Mm-E12.5.Ventral-Progenitor"), 
                       g2 = c("Pig-D18.Ventral-Neuron", "Pig-D18.Ventral-Progenitor", 
                              "Pig-D18.Ventral-Neuron", "Pig-D18.Ventral-Progenitor",
                              "Pig-D18.Ventral-Neuron", "Pig-D18.Ventral-Progenitor",
                              "Pig-D18.Ventral-Neuron", "Pig-D18.Ventral-Progenitor",
                              "Pig-D22.Ventral-Neuron", "Pig-D22.Ventral-Progenitor", 
                              "Pig-D22.Ventral-Neuron", "Pig-D22.Ventral-Progenitor",
                              "Pig-D22.Ventral-Neuron", "Pig-D22.Ventral-Progenitor",
                              "Pig-D22.Ventral-Neuron", "Pig-D22.Ventral-Progenitor",
                              "Ra-E10.Ventral-Neuron", "Ra-E10.Ventral-Progenitor",
                              "Ra-E10.Ventral-Neuron", "Ra-E10.Ventral-Progenitor",
                              "Ra-E10.Ventral-Neuron", "Ra-E10.Ventral-Progenitor",
                              "Ra-E10.Ventral-Neuron", "Ra-E10.Ventral-Progenitor",
                              "Ra-E12.Ventral-Neuron", "Ra-E12.Ventral-Progenitor",
                              "Ra-E12.Ventral-Neuron", "Ra-E12.Ventral-Progenitor",
                              "Ra-E12.Ventral-Neuron", "Ra-E12.Ventral-Progenitor",
                              "Ra-E12.Ventral-Neuron", "Ra-E12.Ventral-Progenitor"))
deg.meta <- data.frame(g1 = c("Ra-E10.Ventral-Neuron", "Ra-E10.Ventral-Progenitor",
                              "Ra-E12.Ventral-Neuron", "Ra-E12.Ventral-Progenitor"), 
                       g2 = c("Pig-D18.Ventral-Neuron", "Pig-D18.Ventral-Progenitor", 
                              "Pig-D22.Ventral-Neuron", "Pig-D22.Ventral-Progenitor"))
for (i in 1:nrow(deg.meta)) {
  tmp.sr <- subset(map.inte, DEG.group %in% deg.meta[i,])
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "DEG.group", 
                       comparison = c(deg.meta[i, "g1"], deg.meta[i, "g2"]), 
                       sample.n = 300)
  tmp.name <- paste0(deg.meta[i, "g2"], "_vs_", deg.meta[i, "g1"])
  deg[[tmp.name]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                 sample.n = NULL, group.by = "CellType", 
                                 g1 = paste0("1_", deg.meta[i, "g1"]), 
                                 g2 = paste0("2_", deg.meta[i, "g2"]), 
                                 lfc = log2(1.5), sig = 0.05,
                                 res.out = file.path(res.out, paste0("DEG/", tmp.name)))
}
deg.meta <- data.frame(g1 = c("Human.Somite", "Mouse.Somite", "Rabbit.Somite", "Mouse.Somite"), 
                       g2 = c("Pig.Somite", "Pig.Somite", "Pig.Somite", "Rabbit.Somite"))
for (i in 1:nrow(deg.meta)) {
  tmp.sr <- subset(sr.subpop.sm, DEG.group %in% deg.meta[i,])
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "DEG.group", 
                       comparison = c(deg.meta[i, "g1"], deg.meta[i, "g2"]), 
                       sample.n = 300)
  tmp.name <- paste0(deg.meta[i, "g2"], "_vs_", deg.meta[i, "g1"])
  deg[[tmp.name]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                 sample.n = NULL, group.by = "CellType", 
                                 g1 = paste0("1_", deg.meta[i, "g1"]), 
                                 g2 = paste0("2_", deg.meta[i, "g2"]), 
                                 lfc = log2(1.5), sig = 0.05,
                                 res.out = file.path(res.out, paste0("DEG/", tmp.name)))
}
rm(i, tmp.sr, tmp)
hs.mm <- read.csv("/home/yhw/document/ensembl/Homologous_Genes/CAME_gene_matches_1v1_human2mouse.csv")
top.pos <- c("Gli1","Hhip","Ptch2","Ptch1","Gm5127","Foxd1","Foxf2","Id4",
             "Epha3","Frem1","1700018A04Rik","Adamts15","Lama2","Lypd6",
             "Trp53i11","St8sia2","Arhgap20","Gm13476","Ntn1","Ckb","Fez1",
             "Gpr30","Foxc2","H2-Ke6","Ephb2","Adam5","Efnb1","Rfc3","Fgfr2","Nek7")
top.pos <- subset(hs.mm, gene.name %in% top.pos)$human.gene.name
top.neg <- c("Edn1","Gprc5a","Aoc3","Stambpl1","Slc38a4","St6galnac2","Gas1",
             "Grem1","Bnc1","Adamts9","Lhfp","Mylk","Cacnb2","Prkg1","Cadm1",
             "Nrp2","B4galnt1","Pamr1","Npnt","Ak4","Pard6b","Aldoc","Slc1a4",
             "Nfkbiz","Mast4","Maml2","Syde2","Lpar4","Lamb1")
top.neg <- subset(hs.mm, gene.name %in% top.neg)$human.gene.name
top.neg <- c(top.neg, "IFT122", "TTC21B", "GPR161")
ventral.markers <- c("GAD1", "GAD2", "BCL11B", 
                     "NKX2-1", "LHX8", "CHAT", "OLIG2", "ISL1", "SLC18A3",
                     "NKX2-2", "NKX6-1", "DLX1", "DLX2", "DLX5",
                     "MSN", "ADORA2A", "CALB1", "ARPP21", "GPR6", "GPR88", "ASCL1",
                     "CHRM4", "DRD1", "DRD2", "PENK", "PPP1R1B", "GSX2", "FOXP1")
for (i in names(deg)) {
  # # top10
  # pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_Genes.pdf")), height = 6, width = 8)
  # print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = NULL, 
  #                         p.col = "PValue", q.col = "FDR", lfc.col = "logFC", 
  #                         sig = 0.05, lfc = log2(1.5),
  #                         curve.plot = F, nosig.remove = FALSE,
  #                         pt.size = 3, pt.shape = 19, pt.alpha = 1,
  #                         l.width = 0.8, top.n = NULL, gene.size = 3,
  #                         label.gene = c("SHH", "IFT122", "TTC21B", "GPR161", "PTCH1", "PTCH2", "FOXA2"),
  #                         high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  # dev.off()
  # # top.pos
  # pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_Genes_top.pos.pdf")), height = 6, width = 8)
  # print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = NULL, 
  #                         p.col = "PValue", q.col = "FDR", lfc.col = "logFC", 
  #                         sig = 0.05, lfc = log2(1.5),
  #                         curve.plot = F, nosig.remove = FALSE,
  #                         pt.size = 3, pt.shape = 19, pt.alpha = 1,
  #                         l.width = 0.8, top.n = NULL, gene.size = 3,
  #                         label.gene = c("SHH", top.pos),
  #                         high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  # dev.off()
  # # top.neg
  # pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_Genes_top.neg.pdf")), height = 6, width = 8)
  # print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = NULL, 
  #                         p.col = "PValue", q.col = "FDR", lfc.col = "logFC", 
  #                         sig = 0.05, lfc = log2(1.5),
  #                         curve.plot = F, nosig.remove = FALSE,
  #                         pt.size = 3, pt.shape = 19, pt.alpha = 1,
  #                         l.width = 0.8, top.n = NULL, gene.size = 3,
  #                         label.gene = c("SHH", top.neg),
  #                         high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  # dev.off()
  # # tfs
  # pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_TFs.pdf")), height = 6, width = 8)
  # print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = hs.tfs$Symbol, 
  #                         p.col = "PValue", q.col = "FDR", lfc.col = "logFC", 
  #                         sig = 0.05, lfc = log2(1.5),
  #                         curve.plot = F, nosig.remove = FALSE,
  #                         pt.size = 3, pt.shape = 19, pt.alpha = 1,
  #                         l.width = 0.8, top.n = 10, label.gene = NULL, gene.size = 3,
  #                         high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  # dev.off()
  # # custom
  # pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_custom.pdf")), height = 6, width = 8)
  # print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = c("SHH", top.pos, top.neg), 
  #                         p.col = "PValue", q.col = "FDR", lfc.col = "logFC", 
  #                         sig = 0.05, lfc = log2(1.5),
  #                         curve.plot = F, nosig.remove = FALSE,
  #                         pt.size = 3, pt.shape = 19, pt.alpha = 1,
  #                         l.width = 0.8, top.n = 10, label.gene = NULL, gene.size = 3,
  #                         high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  # dev.off()
  # custom
  pdf(file.path(res.out, paste0("DEG/", gsub("\\/", "-", i), "_volcano_ventral_markers.pdf")), height = 6, width = 8)
  print(VisDEG.volcano.v2(pd = deg[[i]]$all, geneset = ventral.markers,
                          p.col = "PValue", q.col = "FDR", lfc.col = "logFC",
                          sig = 0.05, lfc = log2(1.5),
                          curve.plot = F, nosig.remove = FALSE,
                          pt.size = 3, pt.shape = 19, pt.alpha = 1,
                          l.width = 0.8, top.n = 10, label.gene = NULL, gene.size = 3,
                          high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9"))
  dev.off()
}
rm(i)
subset(deg$`Pig-D18.Ventral-Progenitor_vs_Mm-E10.5.Ventral-Progenitor`$all, SYMBOL %in% ventral.markers)
subset(deg$`Pig-D18.Ventral-Neuron_vs_Mm-E10.5.Ventral-Neuron`$all, SYMBOL %in% ventral.markers) %>% arrange(desc(logFC))
subset(deg$`Pig-D22.Ventral-Progenitor_vs_Mm-E11.5.Ventral-Progenitor`$all, SYMBOL %in% ventral.markers)
subset(deg$`Pig-D22.Ventral-Neuron_vs_Mm-E11.5.Ventral-Neuron`$all, SYMBOL %in% ventral.markers) %>% arrange(desc(logFC))
subset(deg$`Ra-E10.Ventral-Progenitor_vs_Mm-E10.5.Ventral-Progenitor`$all, SYMBOL %in% ventral.markers)
subset(deg$`Ra-E12.Ventral-Progenitor_vs_Mm-E11.5.Ventral-Progenitor`$all, SYMBOL %in% ventral.markers)
# ventral markers
pd <- list()
for (i in names(deg)) {
  pd[[gsub(" ", "_", gsub("-", "_", gsub("\\/", "-", i)))]] <- deg[[i]]$all
}
names(pd)
pd <- pd[c(3:4, 13:14, 19:20, 29:30, 33:37)]
pdf(file.path(res.out, "DEG_volcano_ventral_markers.pdf"), height = 10, width = 16)
VisDEG.volcano.multi(deg.list = pd, geneset = c(ventral.markers, grep("^SEMA", rownames(sr.subpop.sc), value = T)), 
                     lfc = log2(1.25), sig = 0.05,
                     top.n = 10, pt.size = 2, jitter.width = 0.35, label.gene = NULL,
                     high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#000000")
dev.off()
# 
pd <- list()
for (i in names(deg)) {
  pd[[gsub(" ", "_", gsub("-", "_", gsub("\\/", "-", i)))]] <- deg[[i]]$all
}
names(pd)
pd <- pd[c(3:4, 13:14, 19:20, 29:30, 33:37)]
pdf(file.path(res.out, "DEG_volcano_impaired_limb_coordination.pdf"), height = 10, width = 16)
VisDEG.volcano.multi(deg.list = pd, geneset = c(mgi.genes$Impaired.Limb.Coordination), 
                     lfc = log2(1.25), sig = 0.05,
                     top.n = 10, pt.size = 2, jitter.width = 0.35, label.gene = NULL,
                     high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#000000")
dev.off()
# SEMA ligands
pd <- list()
for (i in names(deg)) {
  pd[[gsub(" ", "_", gsub("-", "_", gsub("\\/", "-", i)))]] <- deg[[i]]$all
}
names(pd)
pd <- pd[c(37:40)]
pdf(file.path(res.out, "DEG_volcano_SEMAs.pdf"), height = 10, width = 8)
VisDEG.volcano.multi(deg.list = pd, geneset = c(grep("^SEMA", rownames(sr.subpop.sc), value = T)), 
                     lfc = log2(1.25), sig = 0.05,
                     top.n = 10, pt.size = 3, jitter.width = 0.35, label.gene = NULL,
                     high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#000000")
dev.off()
# GO
go <- list()
for (fc in c(1.25, 1.5, 2)) {
  for (i in names(deg)) {
    # for the name of cell type with "/"
    tmp <- gsub("\\/", "-", i)
    go[[paste0(tmp, ".up.", fc)]] <- Pipe.GO(species = "human", 
                                             genelist = subset(deg[[i]]$all, PValue <= 0.05 & logFC >= log2(fc))$SYMBOL,
                                             basename = paste0(tmp, "_up_FC", fc), 
                                             genetype = "SYMBOL",
                                             res.out = file.path(res.out, paste0("GO_Pvalue0.05_FC", fc, "/", tmp)))
    go[[paste0(tmp, ".down.", fc)]] <- Pipe.GO(species = "human", 
                                               genelist = subset(deg[[i]]$all, PValue <= 0.05 & logFC <= -log2(fc))$SYMBOL,
                                               basename = paste0(tmp, "_down_FC", fc),
                                               genetype = "SYMBOL",
                                               res.out = file.path(res.out, paste0("GO_Pvalue0.05_FC", fc, "/", tmp)))
  }
}
rm(fc, i, tmp)
# GSEA
gsea <- list()
for (i in names(deg)) {
  gsea[[paste0(i, "_fc2_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg[[i]]$all, deg.type = "edger",
                                                    lfc = log2(2), sig = 0.05,
                                                    reversed = FALSE, species = "human",
                                                    basename = paste0(i, "_fc2_pvalue0.05"),
                                                    genetype = "SYMBOL", gene.col = "SYMBOL",
                                                    outdir = file.path(res.out, paste0("GSEA_fc2_pvalue0.05/", i)))
  gsea[[paste0(i, "_fc1.5_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg[[i]]$all, deg.type = "edger",
                                                      lfc = log2(1.5), sig = 0.05,
                                                      reversed = FALSE, species = "human",
                                                      basename = paste0(i, "_fc1.5_pvalue0.05"),
                                                      genetype = "SYMBOL", gene.col = "SYMBOL",
                                                      outdir = file.path(res.out, paste0("GSEA_fc1.5_pvalue0.05/", i)))
  gsea[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg[[i]]$all, deg.type = "edger", 
                                         lfc = 0, sig = 1, 
                                         reversed = FALSE, species = "human", 
                                         basename = paste0(i, "_all"), 
                                         genetype = "SYMBOL", gene.col = "SYMBOL", 
                                         outdir = file.path(res.out, paste0("GSEA_all/", i)))
}
rm(i)
# plot multiple GO
# hsa04340	Hedgehog signaling pathway
# R-HSA-5632684	Hedgehog 'on' state
# R-HSA-5358351	Signaling by Hedgehog
# R-HSA-5610787	Hedgehog 'off' state
for (fc in c(1.25, 1.5, 2)) {
  go.data <- go[grep(paste0("\\.", fc, "$"), names(go))]
  go.data <- go.data[grep("Progenitor", names(go.data))]
  go.terms <- c("hsa04340", "R-HSA-5632684", "R-HSA-5358351", "R-HSA-5610787")
  pd <- GO.replot(go.data = go.data, go.terms = go.terms,
                  color.by = "pvalue", size.by = "GeneRatio", cols.pal = brewer.pal(9, "Reds"),
                  multi.go = T, multi.cluster = T, seq.x = names(go.data), seq.y = NULL,
                  col.min = -log10(0.05), col.max = 3, pd.title = "GO analysis from DEGs")
  dev.off()
  pdf(file.path(res.out, paste0("GO_Replot_in_up_and_down_genes_across_all_cell_types_fc", fc, ".pdf")), 
      height = 10, width = 15)
  print(pd$bubble)
  dev.off()
}
# GSEA: SHH
gsea.shh <- Pipe.GSEA.v2(deg = deg, lfc.col = "logFC", gene.col = "SYMBOL", genetype = "SYMBOL", 
                         species = "human", msigdbr.cat = "C2",
                         pathway.name = "HEDGEHOG", outdir = file.path(res.out, "GSEA_v2"))
gsea.shh$`Pig-D22.Ventral-Progenitor_vs_Mm-E11.5.Ventral-Progenitor`[, c("Description", "core_enrichment_name")]
# GSEA: sema
gsea.sema <- Pipe.GSEA.v2(deg = deg, lfc.col = "logFC", gene.col = "SYMBOL", genetype = "SYMBOL", 
                          species = "human", msigdbr.cat = "C2",
                          pathway.name = "SEMA", outdir = file.path(res.out, "GSEA_sema"))
# 
shh.targets <- read.table("./scenic_plus/", header = T)
colnames(shh.targets)



### >>> 9. Differential abundance analysis
# miloR
table(map.inte2$Stage)
tmp <- map.inte2
DefaultAssay(tmp) <- "RNA"
for (i in Seurat::Assays(tmp)[-1:-2]) {
  tmp[[i]] <- NULL
}
Idents(tmp) <- tmp$Species
tmp <- subset(tmp, Stage %in% c("Mm-E10.5", "Mm-E11.5", "Pig-E18", "Pig-E21"))
table(tmp$Species)
tmp$Species <- factor(tmp$Species, levels = c("Mouse", "Pig2"))
tmp$Stage <- factor(tmp$Stage, levels = c("Mm-E10.5", "Mm-E11.5", "Pig-E18", "Pig-E21"))
Idents(tmp) <- tmp$Stage
table(tmp$Stage)
tmp <- subset(tmp, downsample = 3000)
table(tmp$Stage)
saveRDS(tmp, file.path(res.out, "miloR/downsample_1000_analysis_seurat.rds"))
milor <- Pipe.miloR(sr.obj = tmp, down.sample = NULL, reduction = "PCA",
                   group.by = "Species", batch.by = "Stage", cell.by = "scpred_prediction", 
                   outdir = file.path(res.out, "miloR_all"))
Pipe.miloR.replot(data = milor, mode = "beeswarm", pd.seq = c(levels(map.inte$scpred_prediction), "Mixed"), pt.size = 2,
                  p.height = 6, p.width = 9, 
                  low = "#0fa3b1", mid = "#e0e0e0", high = "#c52233",
                  outdir = file.path(res.out, "miloR_all"))
Pipe.miloR.replot(data = milor, mode = "hood", pd.seq = c(levels(map.inte$scpred_prediction), "Mixed"), pt.size = 2,
                  p.height = 5, p.width = 12, 
                  low = "#0fa3b1", mid = "#e0e0e0", high = "#c52233",
                  outdir = file.path(res.out, "miloR_all"))
# output to h5ad
tmp <- readRDS(file.path(res.out, "miloR/downsample_1000_analysis_seurat.rds"))
SeuratToH5ad(sr.obj = tmp, method = "SeuratDisk",
             group.by = c("scpred_prediction", "Species", "Stage"), res.out = file.path(res.out, "miloR"), 
             prefix = "downsample_1000_analysis_seurat",
             fileter.by = NULL, filter.term = NULL, gene.set = NULL,
             de.slot = "count", de.assays = "RNA")


### >>> 10. CAME analysis
# output h5 file
library("dior")
rownames(hs.spinal)
table(hs.spinal$CellType)
tmp <- hs.spinal
tmp <- CreateSeuratObject(counts = GetAssayData(tmp, slot = "count"), meta.data = tmp@meta.data)
dior::write_h5(tmp, file = paste0("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/CAME_hs_spinal_wxq.h5"), object.type = 'seurat')
# 
rownames(map.inte)
table(map.inte$scpred_prediction)
tmp <- map.inte
tmp <- CreateSeuratObject(counts = GetAssayData(tmp, slot = "count", assay = "RNA"), meta.data = tmp@meta.data)
dior::write_h5(tmp, file = paste0("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/CAME_hs_spinal_ours.h5"), object.type = 'seurat')



# ============================
# 6th part: AfterChat analysis ----
# ============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Rabbit_Embryos/AfterChat")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Human spinal cord analysis
sr.rabbit$CellType.sub <- str_split_fixed(sr.rabbit$CellType, ":", 2)[, 2]
rabbit.ns <- sr.rabbit
rabbit.ns@meta.data <- rabbit.ns@meta.data %>%
  mutate(CellType.sub2 = case_when(CellType.sub %in% c("roof plate") ~ "Dorsal.RP",
                                   CellType.sub %in% c("dp1","dp2","dp3","dp4","dp5","dp6") ~ "Dorsal.Progenitor",
                                   CellType.sub %in% c("dI1","dI2","dI4","dI5","dI6") ~ "Dorsal.Neuron",
                                   CellType.sub %in% c("floor plate") ~ "Ventral.FP",
                                   CellType.sub %in% c("p0","p1","p2","pMN","p3") ~ "Ventral.Progenitor",
                                   CellType.sub %in% c("v0","v1","v2a","v2b","v3") ~ "Ventral.Neuron",
                                   CellType.sub %in% c("aPSM", "pPSM") ~ "PSM",
                                   CellType.sub %in% c("sclerotome.early", "sclerotome.late") ~ "Sclerotome",
                                   CellType.sub %in% c("myotome") ~ "Myotome",
                                   CellType.sub %in% c("syndetome") ~ "Syndetome",
                                   CellType.sub %in% c("dermomyotome") ~ "Dermomyotome"))
rabbit.ns <- rabbit.ns[, !is.na(rabbit.ns$CellType.sub2)]
table(rabbit.ns$Stage, rabbit.ns$CellType.sub2)
library("AfterChat")
library("SCopeLoomR")
db.cellchat.hs <- CellChatDB.human
ct.rabbit <- Pipe.CellChat(sc.obj = rabbit.ns, cell.db = db.cellchat.hs,
                           group.by = "CellType.sub2", split.by = NULL)
ct.rabbit <- computeCommunProb(ct.rabbit$Sample, type = "truncatedMean", trim = 0.1) %>%
  computeCommunProbPathway() %>%
  aggregateNet()
AfterChat::PathCentrality(ct.obj = ct.rabbit, outdir = file.path(res.out, "rabbit_embryo_ns"), file.prefix = "rabbit_embryo_ns")
AfterChat::PathInteracion(ct.obj = ct.rabbit, outdir = file.path(res.out, "rabbit_embryo_ns"), file.prefix = "rabbit_embryo_ns")
AfterChat::LRsContribution(ct.obj = ct.rabbit, outdir = file.path(res.out, "rabbit_embryo_ns"), file.prefix = "rabbit_embryo_ns")
AfterChat::LRsInteraction(ct.obj = ct.rabbit, outdir = file.path(res.out, "rabbit_embryo_ns"), file.prefix = "rabbit_embryo_ns",
                          cell.source = c("PSM", "Sclerotome", "Syndetome", "Dermomyotome", "Myotome"), 
                          cell.target = c("Dorsal.Neuron", "Ventral.Neuron"))
AfterChat::PathClustering(ct.obj = ct.rabbit, outdir = file.path(res.out, "rabbit_embryo_ns"), file.prefix = "rabbit_embryo_ns")



# =================
# Last part: Exit R ----
# =================
qsave(sr.rabbit, "sr.rabbit.qs")
save.image("Rabbit_Embryo_Data_Analysis.RData")
rm(list = ls())
