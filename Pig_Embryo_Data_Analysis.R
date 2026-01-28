### ========================
### 1st part: Global setting ----
### ========================
setwd("/home/yhw/bioinfo/project-mine/MultiOmics")
dir.create("R/CodeData", recursive = T)
dir.create("R/Graphs", recursive = T)
dir.create("R/Table", recursive = T)
# - input parameter: gene annotation GTF file
gtf <- "/home/yhw/refgenome/embl/rs107/sus_scrofa/Sus_scrofa.Sscrofa11.1.107_CellRanger_for_Signac.gtf"
# - species
species <- "pig"
# - column name in cell metadata to show group info
gp.name <- "Sample"
# - library path
pkg.lib <- "/home/yhw/software/anaconda3/envs/rstudio-4.2.1/lib/R/library"
.libPaths(pkg.lib)
# - human seurat object
hs.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.new.final.rds")
hs.embryo$SuperCluster <- str_split_fixed(hs.embryo$CellType, ":", 2)[,1]
hs.embryo.sub <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.merge2.sub1.rds")
# - human and pig homologous genes
hs.sus.hm <- read.table("/home/cmq/document/ensembl/homologous/human_pig/human_pig_homologous_genes_ensembl.txt", header = T, stringsAsFactors = F, sep = "\t")
hs.sus.hm <- hs.sus.hm[,c(1,2,14,11:13,5:10)]
colnames(hs.sus.hm) <- c("Hs.id", "Hs.id.version", "Hs.chr", "Hs.start", "Hs.end", "Hs.symbol", "Sus.id", "Sus.symbol", "Sus.chr", "Sus.start", "Sus.end", "homology.type")
hs.sus.hm <- subset(hs.sus.hm, homology.type=="ortholog_one2one")
# - human marker genes
marker <- list()
marker$super.cluster <- c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3","TP63","FOXD3","MPZ","PLP1","S100B","ALX1","ALX3",
                          "DLX1","DLX2","SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9","PAX8","WT1","TBX5","PITX1","FOXF1","GATA6",
                          "GATA5","NKX2-5","PLVAP","CD53","PLAC8","CD68","CORO1A","FOXA2","CLDN6","PERP","KRT18",
                          "COL1A1","COL1A2","COL3A1","POSTN","DCN","PAX2","OSR1")
# - dotplot theme setting
theme_dp <- theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
                  axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
                  axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
                  axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))



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



### ==================================================
### 3rd step: Annotation of cell type on super cluster ----
### ==================================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Sus_embryo")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Prepare various data
# - load the RNA and ATAC data
h5.file <- list.files(path = "/home/yhw/bioinfo/project-mine/MultiOmics", pattern = "filtered_feature_bc_matrix.h5", recursive = T, full.names = T)
h5.file <- h5.file[c(2,4)]
names(h5.file) <- c("Day18", "Day22")
counts <- list()
for (h5 in names(h5.file)) {
  counts[[h5]] <- Read10X_h5(h5.file[h5])
}
fragpath <- grep("Day", list.files(getwd(), "atac_fragments.tsv.gz$", recursive = T, full.names = T), value = T)
fragpath <- fragpath[c(2,4)]
names(fragpath) <- names(h5.file)
# - get gene annotations for hg38
library(ensembldb)
gtf.ens <- ensDbFromGtf(gtf, organism = "Sus_scrofa", genomeVersion = "Sscrofa11", version = "107")
gtf.ens <- EnsDb(gtf.ens)
annotation <- GetGRangesFromEnsDb(ensdb = gtf.ens)
seqlevels(annotation) <- gsub("chrMT", "chrM", paste("chr", seqlevels(annotation), sep = ""))
saveRDS(annotation, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/ss11.rs107.gene.annotation.rds")
# - load human TFs
hs.tfs <- read.table("/home/yhw/document/AnimalTFDBs/Homo_sapiens_TF.txt", header = T, stringsAsFactors = F, sep = "\t")
# - create a Seurat object containing the RNA adata
ss.embryo <- list()
for (sample in names(counts)) {
  ss.embryo[[sample]] <- CreateSeuratObject(counts = counts[[sample]]$`Gene Expression`, assay = "RNA")
}
# - create ATAC assay and add it to the object
for (sample in names(counts)) {
  ss.embryo[[sample]][["ATAC"]] <- CreateChromatinAssay(counts = counts[[sample]]$Peaks,
                                                     sep = c(":", "-"),
                                                     fragments = fragpath[sample],
                                                     annotation = annotation)
}
# - merge datasets
ss.embryo <- merge(ss.embryo$Day18, y = ss.embryo$Day22,
                add.cell.ids = c("Ss.Day18", "Ss.Day22"), project = "Ss.multiome")
# - delete useless variables
rm(counts)


### >>> 3. Quality control
# - nucleosome_signal: Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
# - compute per-cell quality control metrics
DefaultAssay(ss.embryo) <- "ATAC"
ss.embryo <- NucleosomeSignal(ss.embryo)
ss.embryo <- TSSEnrichment(ss.embryo, fast = F)
DefaultAssay(ss.embryo) <- "RNA"
ss.mt.gene <- annotation %>%
  plyranges::filter(seqnames == 'chrM') %>%
  plyranges::select(gene_name)
ss.embryo[["percent.mt"]] <- PercentageFeatureSet(ss.embryo, features = ss.mt.gene$gene_name)
ss.embryo$Stage <- gsub("Ss.", "", str_split_fixed(rownames(ss.embryo@meta.data), "_", 2)[,1])
pdf(file.path(res.out, "Quality_control_before_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = ss.embryo,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 7, split.by = "Stage", pt.size = 0, group.by = "Stage")
dev.off()
# filter out low quality cells
ss.embryo <- subset(
  x = ss.embryo,
  subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000 &
    nucleosome_signal <= 2 &
    TSS.enrichment >= 3 &
    percent.mt <= 30
)
pdf(file.path(res.out, "Quality_control_after_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = ss.embryo,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 7, split.by = "Stage", pt.size = 0, group.by = "Stage")
dev.off()
DefaultAssay(ss.embryo) <- "ATAC"
pdf(file.path(res.out, "tss_enrichment_and_fragment_histogram_after_filtering.pdf"), height = 5, width = 12)
p1 <- TSSPlot(ss.embryo, group.by = "orig.ident") + NoLegend() + scale_color_manual(values = c("#0984e3"))
p2 <- FragmentHistogram(object = ss.embryo, group.by = "orig.ident") + scale_fill_manual(values = c("#0984e3"))
p1+p2
dev.off()


### >>> 4. Dimension reduction and clustering
# - post-processing
DefaultAssay(ss.embryo) <- "RNA"
ss.embryo <- NormalizeData(ss.embryo) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
ss.embryo <- ScaleData(ss.embryo, verbose = FALSE, features = rownames(ss.embryo))
ss.embryo <- RunPCA(ss.embryo, npcs = 50, verbose = FALSE, features = VariableFeatures(ss.embryo))
ElbowPlot(ss.embryo, ndims = 30)
dim.n <- 20
# - dimension reduction and clustering
ss.embryo <- FindNeighbors(ss.embryo, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo <- RunUMAP(ss.embryo, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
# - visualize clustering
pdf(file.path(res.out, "Scatter_plot_to_show_uncorrected_clustering.pdf"), height = 5, width = 12)
DimPlot(ss.embryo, reduction = "umap", group.by = c("Stage", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# - remove batch effect
library(harmony)
ss.embryo <- RunHarmony(ss.embryo, group.by.vars = "Stage")
ss.embryo <- RunUMAP(ss.embryo, reduction = "harmony", dims = 1:dim.n)
ss.embryo <- FindNeighbors(ss.embryo, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_corrected_clustering.pdf"), height = 5, width = 12)
DimPlot(ss.embryo, reduction = "umap", group.by = c("Stage", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.25)
dev.off()
# - find marker genes
ss.embryo.markers <- FindAllMarkers(ss.embryo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


### >>> 5. Annotating cell types by Seurat (R)
# - filtering non-homologous genes
hs.sus.hm.tar <- unique(hs.sus.hm[,c(6,7,8)]) %>%
  dplyr::mutate(Sus.symbol=case_when(Sus.symbol=="" ~ Sus.id, TRUE ~ Sus.symbol)) %>%
  dplyr::filter(Hs.symbol %in% rownames(hs.embryo)) %>%
  dplyr::filter(Sus.symbol %in% rownames(ss.embryo)) %>%
  dplyr::filter(! Sus.symbol %in% "AMELY")
sr.predict <- list()
sr.predict$hs <- hs.embryo[rownames(hs.embryo)%in%hs.sus.hm.tar$Hs.symbol,]
sr.predict$sus <- ss.embryo[rownames(ss.embryo)%in%hs.sus.hm.tar$Sus.symbol,]
# - annotate cells using human dataset as reference (Seurat 4.0.6)
anchors <- FindTransferAnchors(reference = sr.predict$hs, query = sr.predict$sus, dims = 1:30, reference.reduction = "pca")
sr.predict$map <- MapQuery(anchorset = anchors, reference = sr.predict$hs, query = sr.predict$sus,
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
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_seurat.pdf"), height = 5, width = 15)
p1+p2
dev.off()
# - delete useless varibles
rm(p1, p2, anchors)


### >>> 6. Annotating cell types by SingleR (R)
# - transform Seurat object into SingleCellExperiment object
hs.sce <- as.SingleCellExperiment(sr.predict$hs)
ss.sce <- as.SingleCellExperiment(sr.predict$sus)
# - cell type annotation
predictions.sgr <- SingleR(test=ss.sce, assay.type.test=1, ref=hs.sce, labels=hs.sce$CellType, num.threads=16)
# - visualization
if (all(colnames(sr.predict$sus) == rownames(predictions.sgr)) & all(colnames(sr.predict$map) == rownames(predictions.sgr))) {
  # pig object
  sr.predict$sus$CellType.SingleR <- predictions.sgr$labels
  sr.predict$sus$SuperCluster.SingleR <- str_split_fixed(sr.predict$sus$CellType.SingleR, ":", 2)[,1]
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
rm(p1, p2, hs.sce, ss.sce)


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
query.set <- GetAssayData(ss.embryo, slot = "data") %>% t() %>% as.data.frame()
# - predict cell type
predictions.scibet.5k <- SciBet(ref.set, query.set, k = 5000)
predictions.scibet.4k <- SciBet(ref.set, query.set, k = 4000)
predictions.scibet.3k <- SciBet(ref.set, query.set, k = 3000)
predictions.scibet.2k <- SciBet(ref.set, query.set, k = 2000)
predictions.scibet.1k <- SciBet(ref.set, query.set, k = 1000)
# -
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
if (identical(colnames(sr.predict$map), rownames(cell.type)) & identical(colnames(ss.embryo), rownames(cell.type))) {
  sr.predict$map$CellType.Final <- cell.type$CellType.Final
  sr.predict$map$SuperCluster.Final <- str_split_fixed(sr.predict$map$CellType.Final, ":", 2)[,1]
  ss.embryo$CellType <- cell.type$CellType.Final
  ss.embryo$SuperCluster <- str_split_fixed(ss.embryo$CellType, ":", 2)[,1]
}
p1 <- DimPlot(sr.predict$hs, reduction = "umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict$map, reduction = "ref.umap", group.by = "SuperCluster.Final", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_with_predicted_cell_types_by_aggregation.pdf"), height = 5, width = 15)
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
# pig
pd <- ss.embryo@meta.data[,c("SuperCluster","CellType")]
library(webr)
pdf(file.path(res.out, "PieDonut_plot_to_show_pig_cell_type_proportion.pdf"), height = 6, width = 6)
PieDonut(pd, aes(pies = SuperCluster), addDonutLabel = TRUE, showRatioDonut = T,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01), color = "#ffffff",
         ratioByGroup = F, labelposition = 1, pieLabelSize = 2, donutLabelSize = 2, r0 = 0.5) +
  scale_fill_manual(values = pd.col.fix)
dev.off()
pd <- table(pd$SuperCluster) %>% as.data.frame() %>% mutate(Species = "human")
pd.col.fix <- pd.col[1:length(unique(sr.predict$hs$SuperCluster))]
names(pd.col.fix) <- sort(unique(sr.predict$hs$SuperCluster))
pdf(file.path(res.out, "Barplot_plot_to_show_pig_cell_type_proportion.pdf"), height = 10, width = 5)
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
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
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
# - sankey plot
library(ggsankey)
pdf(file.path(res.out, "comparing_of_cell_type_predicted_by_Seurat_and_SingleR_sankey_plot.pdf"), width = 13.5, height = 10)
pd <- sr.predict$map@meta.data[,c("SuperCluster.Seurat", "SuperCluster.SingleR", "SuperCluster.SciBet")]
pd %>%
  make_long(SuperCluster.Seurat, SuperCluster.SingleR, SuperCluster.SciBet) -> pd
ggplot(pd, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  scale_fill_manual(values = pd.col.fix) +
  theme_sankey(base_size = 16)
dev.off()



### =====================================================
### 5th step: Re-clustering to check cell type annotation
### =====================================================

### >>> 1. re-clustering
res.out <- file.path(getwd(), "R/Graphs/Sus_embryo/cell_annotation_final")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
ss.embryo.sub <- list()
pd.col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
            "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
            "#00FBFF","#00BCBF","#D277FD","#AE00FF","#7789FD","#203EFF",
            "#FC8ACD","#FF20A4","#8AFCE9","#1CFFD9","#97FF99","#1CFF20",
            "#7B7B7B","#BFBFBF")
# - Craniofacial
ss.embryo.sub$Craniofacial <- ss.embryo[,grep("Craniofacial", ss.embryo$CellType)]
ss.embryo.sub$Craniofacial <- ss.embryo.sub$Craniofacial %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Craniofacial)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Craniofacial <- RunHarmony(ss.embryo.sub$Craniofacial, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Craniofacial, reduction = "harmony", ndims = 30)
dim.n <- 10
ss.embryo.sub$Craniofacial <- RunUMAP(ss.embryo.sub$Craniofacial, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Craniofacial$CellType.sub <- str_split_fixed(ss.embryo.sub$Craniofacial$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Craniofacial) <- ss.embryo.sub$Craniofacial$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Craniofacial_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Craniofacial, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Head mesoderm
ss.embryo.sub$HeadMesoderm <- ss.embryo[,grep("Head mesoderm", ss.embryo$CellType)]
ss.embryo.sub$HeadMesoderm <- ss.embryo.sub$HeadMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$HeadMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$HeadMesoderm <- RunHarmony(ss.embryo.sub$HeadMesoderm, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$HeadMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
ss.embryo.sub$HeadMesoderm <- RunUMAP(ss.embryo.sub$HeadMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$HeadMesoderm$CellType.sub <- str_split_fixed(ss.embryo.sub$HeadMesoderm$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$HeadMesoderm) <- ss.embryo.sub$HeadMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_HeadMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$HeadMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Lateral plate mesoderm
ss.embryo.sub$LateralPlateMesoderm <- ss.embryo[,grep("Lateral plate mesoderm", ss.embryo$CellType)]
table(ss.embryo.sub$LateralPlateMesoderm$CellType)
ss.embryo.sub$LateralPlateMesoderm <- ss.embryo.sub$LateralPlateMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$LateralPlateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$LateralPlateMesoderm <- RunHarmony(ss.embryo.sub$LateralPlateMesoderm, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$LateralPlateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
ss.embryo.sub$LateralPlateMesoderm <- RunUMAP(ss.embryo.sub$LateralPlateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$LateralPlateMesoderm$CellType.sub <- str_split_fixed(ss.embryo.sub$LateralPlateMesoderm$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$LateralPlateMesoderm) <- ss.embryo.sub$LateralPlateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_LateralPlateMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$LateralPlateMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Limb
ss.embryo.sub$Limb <- ss.embryo[,grep("Limb", ss.embryo$CellType)]
ss.embryo.sub$Limb <- ss.embryo.sub$Limb %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Limb)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Limb <- RunHarmony(ss.embryo.sub$Limb, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Limb, reduction = "harmony", ndims = 30)
dim.n <- 8
ss.embryo.sub$Limb <- RunUMAP(ss.embryo.sub$Limb, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Limb$CellType.sub <- str_split_fixed(ss.embryo.sub$Limb$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Limb) <- ss.embryo.sub$Limb$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Limb_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Limb, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Somite
ss.embryo.sub$Somite <- ss.embryo[,grep("Somite", ss.embryo$CellType)]
ss.embryo.sub$Somite <- ss.embryo.sub$Somite %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Somite)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Somite <- RunHarmony(ss.embryo.sub$Somite, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Somite, reduction = "harmony", ndims = 30)
dim.n <- 10
ss.embryo.sub$Somite <- RunUMAP(ss.embryo.sub$Somite, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Somite$CellType.sub <- str_split_fixed(ss.embryo.sub$Somite$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Somite) <- ss.embryo.sub$Somite$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Somite_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Somite, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Endothelium
ss.embryo.sub$Endothelium <- ss.embryo[,grep("Endothelium", ss.embryo$CellType)]
ss.embryo.sub$Endothelium <- ss.embryo.sub$Endothelium %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Endothelium)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Endothelium <- RunHarmony(ss.embryo.sub$Endothelium, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Endothelium, reduction = "harmony", ndims = 30)
dim.n <- 8
ss.embryo.sub$Endothelium <- RunUMAP(ss.embryo.sub$Endothelium, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Endothelium$CellType.sub <- str_split_fixed(ss.embryo.sub$Endothelium$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Endothelium) <- ss.embryo.sub$Endothelium$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endothelium_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Endothelium, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Blood
library(harmony)
ss.embryo.sub$Blood <- ss.embryo[,grep("Blood", ss.embryo$CellType)]
ss.embryo.sub$Blood <- ss.embryo.sub$Blood %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Blood)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Blood <- RunHarmony(ss.embryo.sub$Blood, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Blood, reduction = "harmony", ndims = 30)
dim.n <- 6
ss.embryo.sub$Blood <- RunUMAP(ss.embryo.sub$Blood, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Blood$CellType.sub <- str_split_fixed(ss.embryo.sub$Blood$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Blood) <- ss.embryo.sub$Blood$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Blood_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Blood, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Intermediate mesoderm
ss.embryo.sub$IntermediateMesoderm <- ss.embryo[,grep("Intermediate mesoderm", ss.embryo$CellType)]
ss.embryo.sub$IntermediateMesoderm <- ss.embryo.sub$IntermediateMesoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$IntermediateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$IntermediateMesoderm <- RunHarmony(ss.embryo.sub$IntermediateMesoderm, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$IntermediateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 5
ss.embryo.sub$IntermediateMesoderm <- RunUMAP(ss.embryo.sub$IntermediateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$IntermediateMesoderm$CellType.sub <- str_split_fixed(ss.embryo.sub$IntermediateMesoderm$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$IntermediateMesoderm) <- ss.embryo.sub$IntermediateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_IntermediateMesoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$IntermediateMesoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Endoderm
ss.embryo.sub$Endoderm <- ss.embryo[,grep("Endoderm", ss.embryo$CellType)]
ss.embryo.sub$Endoderm <- ss.embryo.sub$Endoderm %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Endoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Endoderm <- RunHarmony(ss.embryo.sub$Endoderm, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Endoderm, reduction = "harmony", ndims = 30)
dim.n <- 5
ss.embryo.sub$Endoderm <- RunUMAP(ss.embryo.sub$Endoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Endoderm$CellType.sub <- str_split_fixed(ss.embryo.sub$Endoderm$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Endoderm) <- ss.embryo.sub$Endoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endoderm_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Endoderm, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Neuron
ss.embryo.sub$Neuron <- ss.embryo[,grep("Neuron", ss.embryo$CellType)]
ss.embryo.sub$Neuron <- ss.embryo.sub$Neuron %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Neuron)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Neuron <- RunHarmony(ss.embryo.sub$Neuron, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Neuron, reduction = "harmony", ndims = 30)
dim.n <- 8
ss.embryo.sub$Neuron <- RunUMAP(ss.embryo.sub$Neuron, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Neuron$CellType.sub <- str_split_fixed(ss.embryo.sub$Neuron$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Neuron) <- ss.embryo.sub$Neuron$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Neuron_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Neuron, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Neural progenitor
ss.embryo.sub$NeuralProgenitor <- ss.embryo[,grep("Neural progenitor", ss.embryo$CellType)]
ss.embryo.sub$NeuralProgenitor <- ss.embryo.sub$NeuralProgenitor %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$NeuralProgenitor)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$NeuralProgenitor <- RunHarmony(ss.embryo.sub$NeuralProgenitor, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$NeuralProgenitor, reduction = "harmony", ndims = 30)
dim.n <- 8
ss.embryo.sub$NeuralProgenitor <- RunUMAP(ss.embryo.sub$NeuralProgenitor, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$NeuralProgenitor$CellType.sub <- str_split_fixed(ss.embryo.sub$NeuralProgenitor$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$NeuralProgenitor) <- ss.embryo.sub$NeuralProgenitor$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_NeuralProgenitor_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$NeuralProgenitor, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Schwann cells
ss.embryo.sub$Schwann <- ss.embryo[,grep("Schwann cells", ss.embryo$CellType)]
ss.embryo.sub$Schwann <- ss.embryo.sub$Schwann %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Schwann)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Schwann <- RunHarmony(ss.embryo.sub$Schwann, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Schwann, reduction = "harmony", ndims = 30)
dim.n <- 6
ss.embryo.sub$Schwann <- RunUMAP(ss.embryo.sub$Schwann, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Schwann$CellType.sub <- str_split_fixed(ss.embryo.sub$Schwann$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Schwann) <- ss.embryo.sub$Schwann$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Schwann_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Schwann, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Undefined
ss.embryo.sub$Undefined <- ss.embryo[,grep("undefined", ss.embryo$CellType)]
ss.embryo.sub$Undefined <- ss.embryo.sub$Undefined %>% FindVariableFeatures() %>%
  ScaleData(rownames(ss.embryo.sub$Undefined)) %>% RunPCA(verbose = FALSE)
set.seed(100)
ss.embryo.sub$Undefined <- RunHarmony(ss.embryo.sub$Undefined, group.by.vars = "Stage")
ElbowPlot(ss.embryo.sub$Undefined, reduction = "harmony", ndims = 30)
dim.n <- 10
ss.embryo.sub$Undefined <- RunUMAP(ss.embryo.sub$Undefined, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
ss.embryo.sub$Undefined$CellType.sub <- str_split_fixed(ss.embryo.sub$Undefined$CellType, ":", "2")[,2]
Idents(ss.embryo.sub$Undefined) <- ss.embryo.sub$Undefined$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Undefined_with_predicted_cell_types.pdf"), height = 5, width = 13)
DimPlot(ss.embryo.sub$Undefined, reduction = "umap", group.by = c("Stage","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 1, label.box = T, cols = pd.col)
dev.off()
# - Correlation analysis based on specificity score matrix (Genes + TFs)
spe.score.seurat <- list()
for (cluster in names(hs.embryo.sub)) {
  expr.hs = AverageExpression(hs.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
  expr.ss = AverageExpression(ss.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
  deg.hs = VariableFeatures(hs.embryo.sub[[cluster]])
  deg.ss = VariableFeatures(ss.embryo.sub[[cluster]])
  spe.score.seurat[[cluster]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                 ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                                 Species1 = "human", Species2 = "pig",
                                                 filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", cluster)),
                                                 Height = 12, Width = 12)
}
# all cell types
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.ss = AverageExpression(ss.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.ss = VariableFeatures(ss.embryo)
spe.score.seurat[["all"]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                             ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                             Species1 = "human", Species2 = "pig",
                                             filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")))
CorrComparePlot.again(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                      ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                      objects.to.return = spe.score.seurat[["all"]],
                      Species1 = "human", Species2 = "pig",
                      filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")),
                      Height = 70, Width = 70)
rm(cluster, expr.hs, expr.ss, deg.hs, deg.ss)
# all clusters
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.ss = AverageExpression(ss.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.ss = VariableFeatures(ss.embryo)
spe.score.seurat[["cluster"]] <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                 ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                                 Species1 = "human", Species2 = "pig",
                                                 filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", "All_Super_Cluster")))
# - Correlation analysis based on specificity score matrix (Genes + LRs)
res.out <- "/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/Sus_embryo/cell_annotation_final/LRs"
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
spe.score.seurat.LR <- list()
for (cluster in names(hs.embryo.sub)) {
  expr.hs = AverageExpression(hs.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
  expr.ss = AverageExpression(ss.embryo.sub[[cluster]], assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
  deg.hs = VariableFeatures(hs.embryo.sub[[cluster]])
  deg.ss = VariableFeatures(ss.embryo.sub[[cluster]])
  spe.score.seurat.LR[[cluster]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                      ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                                      Species1 = "human", Species2 = "pig",
                                                      filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", cluster)),
                                                      Height = 12, Width = 12)
}
# all cell types
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.ss = AverageExpression(ss.embryo, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.ss = VariableFeatures(ss.embryo)
spe.score.seurat.LR[["all"]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                  ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                                  Species1 = "human", Species2 = "pig",
                                                  filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")))
CorrComparePlotLR.again(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                        ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                        objects.to.return = spe.score.seurat[["all"]],
                        Species1 = "human", Species2 = "pig",
                        filename = file.path(res.out,paste0("Specificity_score_correlation_plot_of_", "All_cell_types")),
                        Height = 70, Width = 70)
# all clusters
expr.hs = AverageExpression(hs.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.ss = AverageExpression(ss.embryo, assays = "RNA", group.by = "SuperCluster", slot = "data") %>% as.data.frame()
colnames(expr.ss) = paste("Ss", gsub("RNA.", "", colnames(expr.ss)), sep = "_")
deg.hs = VariableFeatures(hs.embryo)
deg.ss = VariableFeatures(ss.embryo)
spe.score.seurat.LR[["cluster"]] <- CorrComparePlotLR(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                                      ExpressionTableSpecies2 = expr.ss, DEgenesSpecies2 = deg.ss,
                                                      Species1 = "human", Species2 = "pig",
                                                      filename = file.path(res.out, paste0("Specificity_score_correlation_plot_of_", "All_Super_Cluster")))
rm(cluster, expr.hs, expr.ss, deg.hs, deg.ss)



### =======================================================
### 6th step: scATAC-seq: integrated clustering (all cells)
### =======================================================

### >>> 1. calculation of gene activity via snaptools
# - create GRange object of gtf
# human
hs.ge.body.tss.u2k <- read.table("/home/yhw/bioinfo/project-mine/MultiOmics/R/Table/gene_annotation/Homo_sapiens_GRCh38_104_gene_body_and_tss_up_2kb.bed", header = F, stringsAsFactors = F, sep = "\t")
colnames(hs.ge.body.tss.u2k) <- c("chr", "start", "end", "strand", "ensembl", "symbol", "type")
hs.ge.body.tss.u2k$symbol <- gsub("_", "-", hs.ge.body.tss.u2k$symbol)
hs.ge.body.tss.u2k <- subset(hs.ge.body.tss.u2k, symbol%in%rownames(hs.embryo))
hs.ge.body.tss.u2k <- GRanges(hs.ge.body.tss.u2k[,1], IRanges(hs.ge.body.tss.u2k[,2], hs.ge.body.tss.u2k[,3]), name=hs.ge.body.tss.u2k[,6], strand = hs.ge.body.tss.u2k[,4])
# pig
ss.ge.body.tss.u2k <- read.table("/home/yhw/bioinfo/project-mine/MultiOmics/R/Table/gene_annotation/Sus_scrofa.Sscrofa11.1.104_CellRanger_For_Signac_gene_body_and_tss_up_2kb.bed", header = F, stringsAsFactors = F, sep = "\t")
colnames(ss.ge.body.tss.u2k) <- c("chr", "start", "end", "strand", "ensembl", "symbol", "type")
ss.ge.body.tss.u2k$symbol <- gsub("_", "-", ss.ge.body.tss.u2k$symbol)
ss.ge.body.tss.u2k <- subset(ss.ge.body.tss.u2k, symbol%in%rownames(ss.embryo))
ss.ge.body.tss.u2k <- GRanges(ss.ge.body.tss.u2k[,1], IRanges(ss.ge.body.tss.u2k[,2], ss.ge.body.tss.u2k[,3]), name=ss.ge.body.tss.u2k[,6], strand = ss.ge.body.tss.u2k[,4])
# - load snap file and quality matrix
snap.path <- list.files(getwd(), "snap$", recursive = T, full.names = T)
names(snap.path) <- c("Day18", "Day22", "Hs.5W.1", "Hs.5W.2", "Hs.4W.1", "Hs.6W.1")
qc.path <- list.files(getwd(), "per_barcode_metrics.csv$", recursive = T, full.names = T)
qc.path <- qc.path[-5]
names(qc.path) <- names(snap.path)
# - create snap object
snap <- list()
for (i in seq(1, length(snap.path))) {
  snap[[i]] <- createSnap(
    file=snap.path[i], sample=names(snap.path)[i], num.cores=1)
  names(snap)[i] <- names(snap.path)[i]
}; rm(i)
# - filter cells based on RNA data and add cell metadata
# pig
for (i in 1:2) {
  tmp <- subset(ss.embryo, Stage==names(snap)[i])
  snap[[i]] <- snap[[i]][snap[[i]]@barcode %in% str_split_fixed(colnames(tmp), "_", 2)[,2],]
  tmp <- read.csv(qc.path[i], header = T)
  rownames(tmp) <- tmp$gex_barcode
  tmp <- tmp[snap[[i]]@barcode, -1]
  if (all(snap[[i]]@barcode==rownames(tmp))) {
    snap[[i]]@metaData <- tmp
  }
}; rm(i, tmp)
# human
for (i in 3:6) {
  tmp <- subset(hs.embryo, Sample==names(snap)[i])
  snap[[i]] <- snap[[i]][snap[[i]]@barcode %in% str_split_fixed(colnames(tmp), "_", 2)[,2],]
  tmp <- read.csv(qc.path[i], header = T)
  rownames(tmp) <- tmp$gex_barcode
  tmp <- tmp[snap[[i]]@barcode, -1]
  if (all(snap[[i]]@barcode==rownames(tmp))) {
    snap[[i]]@metaData <- tmp
  }
}; rm(i, tmp)
# - create raw cell-by-gene matrix
# add bin matrix
for (i in 1:6) {
  snap[[i]] <- addBmatToSnap(snap[[i]])
}; rm(i)
# generation
for (i in 1:2) {
  snap[[i]] <- createGmatFromMat(obj=snap[[i]], input.mat="bmat", genes=ss.ge.body.tss.u2k, do.par=TRUE, num.cores=16)
}
for (i in 3:6) {
  snap[[i]] <- createGmatFromMat(obj=snap[[i]], input.mat="bmat", genes=hs.ge.body.tss.u2k, do.par=TRUE, num.cores=16)
}
# saving
ge.act <- list()
ge.act$raw <- list()
for (i in 1:6) {
  ge.act$raw[[i]] <- snap[[i]]@gmat
  names(ge.act$raw)[i] <- names(snap)[i]
  if (i%in%1:2) {
    rownames(ge.act$raw[[i]]) <- paste("Ss.", names(ge.act$raw)[i], "_", rownames(ge.act$raw[[i]]), sep = "")
  } else if (i%in%3:6) {
    rownames(ge.act$raw[[i]]) <- paste(names(ge.act$raw)[i], "_", rownames(ge.act$raw[[i]]), sep = "")
  }
}; rm(i)
# - normalize the cell-by-gene matrix
# normalization
ge.act$RPM <- list()
for (i in seq(1, length(snap))) {
  snap[[i]] <- scaleCountMatrix(obj=snap[[i]], cov=snap[[i]]@metaData$atac_fragments + 1, mat="gmat", method = "RPM")
  ge.act$RPM[[i]] <- snap[[i]]@gmat
  names(ge.act$RPM)[i] <- names(snap)[i]
  if (i%in%1:2) {
    rownames(ge.act$RPM[[i]]) <- paste("Ss.", names(ge.act$RPM)[i], "_", rownames(ge.act$RPM[[i]]), sep = "")
  } else if (i%in%3:6) {
    rownames(ge.act$RPM[[i]]) <- paste(names(ge.act$RPM)[i], "_", rownames(ge.act$RPM[[i]]), sep = "")
  }
}; rm(i)
saveRDS(ge.act, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.ss.gene.body.tss.u2k.activity.snapatac.rds")

### >>> 2. calculate gene activity via signac
ge.act.sn <- list()
DefaultAssay(hs.embryo) <- "ATAC"
ge.act.sn$hs <- GeneActivity(hs.embryo, assay = "ATAC", extend.upstream = 2000)
DefaultAssay(ss.embryo) <- "ATAC"
ge.act.sn$ss <- GeneActivity(ss.embryo, assay = "ATAC", extend.upstream = 2000)
saveRDS(ge.act.sn, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.ss.gene.body.tss.u2k.activity.signac.rds")

### >>> 3. human: dimension reduction based on ATAC-seq signal on genes
hs.embryo[['ACT']] <- CreateAssayObject(counts = ge.act.sn$hs)
hs.embryo <- NormalizeData(object = hs.embryo, assay = 'ACT', normalization.method = 'LogNormalize', scale.factor = median(hs.embryo$nCount_ACT))
# - standard workflow
DefaultAssay(hs.embryo) <- "ACT"
hs.embryo <- RunTFIDF(hs.embryo)
hs.embryo <- FindTopFeatures(hs.embryo, min.cutoff = 'q0')
hs.embryo <- RunSVD(object = hs.embryo)
DepthCor(hs.embryo)
hs.embryo <- RunUMAP(object = hs.embryo, reduction = 'lsi', dims = 2:30, reduction.name = "atac.umap", return.model = T, min.dist = 0.2)
p1 <- DimPlot(object = hs.embryo, label = TRUE, reduction = "atac.umap", group.by = c("SuperCluster"), cols = pd.col.fix)
p2 <- DimPlot(object = hs.embryo, label = TRUE, reduction = "atac.umap", group.by = c("Stage"), cols = pd.col)
p1+p2
set.seed(123)
hs.embryo <- RunHarmony(hs.embryo, group.by.vars = "Sample", reduction = "lsi", reduction.name="lsi.harmony", project.dim = F)
hs.embryo <- RunUMAP(object = hs.embryo, reduction = 'harmony', dims = 2:30, reduction.name = "atac.umap", return.model = T)
p1 <- DimPlot(object = hs.embryo, label = TRUE, reduction = "atac.umap", group.by = c("SuperCluster"), cols = pd.col.fix)
p2 <- DimPlot(object = hs.embryo, label = TRUE, reduction = "atac.umap", group.by = c("Stage"), cols = pd.col)
p1+p2

### >>> 4. mapping pig data to human data
ss.embryo[['ACT']] <- CreateAssayObject(counts = ge.act.sn$ss)
ss.embryo <- NormalizeData(object = ss.embryo, assay = 'ACT', normalization.method = 'LogNormalize', scale.factor = median(ss.embryo$nCount_ACT))
# - standard workflow
DefaultAssay(ss.embryo) <- "ACT"
ss.embryo <- RunTFIDF(ss.embryo)
ss.embryo <- FindTopFeatures(ss.embryo, min.cutoff = 'q0')
# - subset seurat object in homologous genes
sr.predict.atac <- list()
sr.predict.atac$hs <- hs.embryo[hs.sus.hm.tar$Hs.symbol,]
sr.predict.atac$ss <- ss.embryo[hs.sus.hm.tar$Sus.symbol,]
# - map pig dataset onto the human dataset
# find transfer anchors
anchors <- FindTransferAnchors(reference = sr.predict.atac$hs, query = sr.predict.atac$ss, reference.reduction = "lsi",
                               reduction = "lsiproject", dims = 2:30)
# map query onto the reference dataset
sr.predict.atac$map <- MapQuery(anchorset = anchors, reference = sr.predict.atac$hs, query = sr.predict.atac$ss, refdata = sr.predict.atac$hs$CellType,
                                reference.reduction = "lsi", new.reduction.name = "ref.lsi", reduction.model = 'atac.umap')
p1 <- DimPlot(hs.embryo, reduction = "atac.umap", group.by = "SuperCluster", label = F, pt.size = 0.15) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.predict.atac$map, reduction = "ref.umap", group.by = "SuperCluster.Seurat", label = F, pt.size = 0.15) +
  ggtitle("Pig embryos") +
  scale_color_manual(values = pd.col.fix)
pdf(file.path(res.out, "Integrated_clustering_of_gene_ATAC_signal_with_predicted_cell_types_by_seurat.pdf"), height = 5, width = 15)
p1+p2
dev.off(); rm(p1, p2)



### =====================================================
### 7th step: visualization of super cluster marker genes
### =====================================================

### >>> 1. human
# - RNA assay
Idents(hs.embryo) <- hs.embryo$SuperCluster
pd <- subset(hs.embryo, SuperCluster %in% unique(setdiff(hs.embryo$SuperCluster, c("Undefined", "Fibroblast"))))
DefaultAssay(pd) <- "RNA"
pd <- NormalizeData(pd) %>% ScaleData()
pd@meta.data %>%
  mutate(SuperCluster=factor(SuperCluster,
                             levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", "Craniofacial",
                                        "Head mesoderm", "Somite", "Intermediate mesoderm", "Limb","Lateral plate mesoderm",
                                        "Heart", "Endothelium", "Blood", "Endoderm"))) -> pd@meta.data
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/hs_MPR/Dotplot_to_show_human_super_cluster_markers_grouped_by_cell_type.pdf",
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39],
        assay = "RNA", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()
pd <- SeuratWrappers::RunALRA(pd)
pd <- ScaleData(pd, features = rownames(pd))
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/hs_MPR/Dotplot_to_show_human_super_cluster_markers_grouped_by_cell_type_alra.pdf",
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39],
        assay = "alra", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()


### >>> 2. pig
# - RNA assay
# predicted cell types by Seurat
Idents(ss.embryo) <- ss.embryo$SuperCluster
pd <- subset(ss.embryo, SuperCluster %in% setdiff(ss.embryo$SuperCluster, c("Undefined", "Fibroblast")))
DefaultAssay(pd) <- "RNA"
pd <- NormalizeData(pd) %>% ScaleData()
pd@meta.data %>%
  mutate(SuperCluster=factor(SuperCluster,
                             levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", "Craniofacial",
                                        "Head mesoderm", "Somite", "Intermediate mesoderm", "Limb","Lateral plate mesoderm",
                                        "Heart", "Endothelium", "Blood", "Endoderm"))) -> pd@meta.data
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/Sus_embryo/Dotplot_to_show_pig_super_cluster_markers_grouped_by_cell_type_RNA.pdf",
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39],
        assay = "RNA", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()
pd <- SeuratWrappers::RunALRA(pd)
pd <- ScaleData(pd, features = rownames(pd))
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/Sus_embryo/Dotplot_to_show_pig_super_cluster_markers_grouped_by_cell_type_alra.pdf",
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster", marker$super.cluster[1:39],
        assay = "alra", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()


# predicted cell types by SingleR
ss.embryo@meta.data <- ss.embryo@meta.data %>%
  mutate(SuperCluster.SingleR=str_split_fixed(CellType.SingleR, ":", 2)[,1])
Idents(ss.embryo) <- ss.embryo$SuperCluster.SingleR
pd <- subset(ss.embryo, SuperCluster.SingleR %in% setdiff(ss.embryo$SuperCluster.SingleR, c("Undefined", "Fibroblast")))
DefaultAssay(pd) <- "RNA"
pd <- NormalizeData(pd) %>% ScaleData()
pd@meta.data %>%
  mutate(SuperCluster.SingleR=factor(SuperCluster.SingleR,
                                     levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", "Craniofacial",
                                                "Head mesoderm", "Somite", "Intermediate mesoderm", "Limb","Lateral plate mesoderm",
                                                "Heart", "Endothelium", "Blood", "Endoderm"))) -> pd@meta.data
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/Sus_embryo/Dotplot_to_show_pig_super_cluster_markers_grouped_by_cell_type_SingleR.pdf",
    height = 4.5, width = 12)
DotPlot(pd, group.by = "SuperCluster.SingleR", marker$super.cluster[1:39],
        assay = "RNA", cols = c("#eaeaea", "#fc0330"), scale = T) +
  RotatedAxis() + labs(x = "Hs Super Clsuter Marker Genes", y = "Super Cluster")
dev.off()



### ====================
### Last part: Save Data ----
### ====================
saveRDS(ss.embryo, "ss.embryo.final.version.rds")
save.image("Ss_Embryo_Multiome.RData")
rm(list = ls())
