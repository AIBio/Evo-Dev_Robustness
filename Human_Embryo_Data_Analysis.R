### ========================
### 1st part: Global setting ----
### ========================
setwd("/home/yhw/bioinfo/project-mine/MultiOmics")
dir.create("R/CodeData", recursive = T)
dir.create("R/Graphs", recursive = T)
dir.create("R/Table", recursive = T)
# - input parameter: gene annotation GTF file
gtf <- "/home/yhw/refgenome/gencode/hg38_v38/gencode.v38.primary_assembly.annotation_CellRanger_filter.gtf"
# - species
species <- "human"
# - column name in cell metadata to show group info
gp.name <- "Sample"
# - library path
pkg.lib <- "/home/yhw/software/anaconda3/envs/rstudio-4.2.1/lib/R/library"
.libPaths(pkg.lib)
# - dotplot theme setting
theme_dp <- theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
                  axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
                  axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
                  axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))



### =================
### 2nd part: Library ----
### =================
library(devtools)
library(BiocManager)
### >>> 1. CRAN packages
cran.pks <- c("SoupX", "Seurat", "harmony", "rliger", "stringr", "circlize", "ggplot2", "ggalluvial", "RCurl",
              "ggrepel", "forcats", "dplyr", "tidyr", "tidyverse", "stringr", "Signac", "SnapATAC")
for (pks in cran.pks) {
  #if (!require(pks, character.only = TRUE)) {
  #  install.packages(pks, lib = pkg.lib)
  #  library(pks, character.only = T)
  #} else {
    library(pks, character.only = T)
  #}
}


### >>> 2. Bioconductor packages
bioc.pks <- c("BSgenome.Hsapiens.UCSC.hg38", "SingleR")
for (pks in bioc.pks) {
  if (!require(pks, character.only = TRUE)) {
    BiocManager::install(pks, lib = pkg.lib)
    library(pks, character.only = T)
  } else {
    library(pks, character.only = T)
  }
}


### >>> 3. GitHub packages
if (!require("SeuratDisk", character.only = TRUE)) {
  devtools::install_github('mojaveazure/seurat-disk')
  library("SeuratDisk", character.only = T)
} else {
  library("SeuratDisk", character.only = T)
}
if (!require("CellChat", character.only = TRUE)) {
  devtools::install_github("sqjin/CellChat")
  library("CellChat", character.only = T)
} else {
  library("CellChat", character.only = T)
}
if (!require("SCENIC", character.only = TRUE)) {
  devtools::install_github("aertslab/SCENIC")
  library("SCENIC", character.only = T)
} else {
  library("SCENIC", character.only = T)
}
if (!require("velocyto.R", character.only = TRUE)) {
  devtools::install_github("velocyto-team/velocyto.R")
  library("velocyto.R", character.only = T)
} else {
  library("velocyto.R", character.only = T)
}
if (!require("pagoda2", character.only = TRUE)) {
  devtools::install_github("kharchenkolab/pagoda2")
  library("pagoda2", character.only = T)
} else {
  library("pagoda2", character.only = T)
}
if (!require("leidenbase", character.only = TRUE)) {
  devtools::install_github("cole-trapnell-lab/leidenbase")
  library("leidenbase", character.only = T)
} else {
  library("leidenbase", character.only = T)
}
if (!require("SCopeLoomR", character.only = TRUE)) {
  devtools::install_github("aertslab/SCopeLoomR")
  library("SCopeLoomR", character.only = T)
} else {
  library("SCopeLoomR", character.only = T)
}



### ==================================================
### 3rd step: Annotation of cell type on super cluster ----
### ==================================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/hs_MPR")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Prepare various data
# - load the RNA and ATAC data
h5.file <- list.files(getwd(), "filtered_feature_bc_matrix.h5", recursive = T, full.names = T)
names(h5.file) <- c("Hs.5W.2", "Hs.5W.3", "Hs.4W.1", "Hs.6W.1")
counts <- list()
for (h5 in names(h5.file)) {
  counts[[h5]] <- Read10X_h5(h5.file[h5])
}
fragpath <- list.files(getwd(), "atac_fragments.tsv.gz$", recursive = T, full.names = T)
names(fragpath) <- names(h5.file)
# - get gene annotations for hg38
#gtf.ens <- ensDbFromGtf(gtf, organism = "Homo_sapiens", genomeVersion = "GRCh38", version = "104")
#gtf.ens <- EnsDb(gtf.ens)
annotation <- GetGRangesFromEnsDb(ensdb = gtf.ens)
seqlevels(annotation) <- gsub("chrMT", "chrM", paste("chr", seqlevels(annotation), sep = ""))
saveRDS(annotation, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hg38.gene.annotation.rds")
# - load human TFs
hs.tfs <- read.table("/home/yhw/document/AnimalTFDBs/Homo_sapiens_TF.txt", header = T, stringsAsFactors = F, sep = "\t")
# - create a Seurat object containing the RNA adata
hs.mpr <- list()
for (sample in names(counts)) {
  hs.mpr[[sample]] <- CreateSeuratObject(counts = counts[[sample]]$`Gene Expression`, assay = "RNA")
}
# - create ATAC assay and add it to the object
for (sample in names(counts)) {
  hs.mpr[[sample]][["ATAC"]] <- CreateChromatinAssay(counts = counts[[sample]]$Peaks,
                                                     sep = c(":", "-"),
                                                     fragments = fragpath[sample],
                                                     annotation = annotation)
}
# - merge datasets
hs.mpr <- merge(hs.mpr$Hs.4W.1, y = c(hs.mpr$Hs.5W.2, hs.mpr$Hs.5W.3, hs.mpr$Hs.6W.1),
                add.cell.ids = c("Hs.4W.1", "Hs.5W.1", "Hs.5W.2", "Hs.6W.1"), project = "Hs.multiome")
# - load marker genes (classified by tissue)
marker <- list()
marker$super.cluster <- c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3","PERP","KRT18","FOXD3","MPZ","PLP1","S100B","ALX1","ALX3",
                          "DLX1","DLX2","SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9","PAX8","WT1","TBX5","PITX1","FOXF1","GATA6",
                          "GATA5","NKX2-5","PLVAP","CD53","PLAC8","CD68","CORO1A","FOXA2","CLDN6","EMX2","PAX2","OSR1")
marker$nv <- c("BARHL1","BARHL2","ZIC1","LHX5","PAX6", "FGF8","FGF17","SIX3", "DMRTA1","DMRTA2","EMX2","PAX6","DMRT3",
               "NKX2-1","SOX2", "WNT4","WNT2B", "SHH","WNT5A", "MEIS2","SP5","TAL2","MAB21L2", "FGF8","FOXG1",
               "GATA3","DLX5","LMX1A","OC90", "RAX","SIX3","PAX6","TBX3","BMP4", "MITF","PAX2","DCT",
               "VSX2", "NKX2-1","DLX1","DLX2","DLX5", "LMX1A", "MSX1","OLIG3","ATOH1", "MSX1","OLIG3",
               "MSX1","OLIG3","ASCL1","PHOX2B", "ASCL1","PTF1A", "PAX3","PHOX2B", "PAX3","DBX2", "PAX3","DBX2","PRDM13",
               "PAX6","DBX2", "PRDM8","NKX6-2", "PRDM8","NKX6-1","ASCL1", "OLIG2","NKX6-1", "NKX2-2","PHOX2B",
               "SHH","FOXA2", "LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3", "MSX1","OLIG3","ATOH1","HOXB6",
               "MSX1","OLIG3","HOXB6", "MSX1","OLIG3","ASCL1","HOXB6", "ASCL1","PAX3","HOXB6",
               "PAX3","DBX2","HOXB6", "PAX3","DBX2","PRDM13","HOXB6", "PAX6","DBX2","HOXB6",
               "PRDM8","NKX6-2","HOXB6", "PAX6","IRX3","NKX6-1","HOXB6", "NKX6-1","OLIG1","OLIG2","HOXB6",
               "NKX6-1","NKX2-2","NKX2-8","HOXB6", "NKX6-1","FOXA2","SHH","SPON1","HOXB6", "EN1","PAX8","FGF17","CNPY1",
               "ASCL1","BARHL1","BARHL2","DLX1","DLX2","DLX5","DLX6-AS1","DMRT3","EN1","EOMES","EVX1","FOXA2","FOXP1","GAD2","GATA2","GATA3",
               "HOXB6","ISL1","ISL2","LBX1","LHX1","LHX2","LHX3","LHX5","LHX9","NKX2-1","NKX2-2","NKX6-1","OTP","PAX2","PHOX2A","PHOX2B",
               "PITX2","POU4F1","SIM1","SOX14","SOX21","TAL2","TBR1","TLX1","TLX3","VSX1","VSX2","SIX1","FOXD3","FOXP2","BHLHE22","DRGX")
marker$sw <- c("ASCL1","DCX","ELAVL4","FOXD3","HAND2","ISL1","MITF","MPZ","PAX3","PHOX2A","PHOX2B","PLP1","PMEL","PMP22","SOX10","TFAP2A","TFAP2B","TUBB3")
marker$epidermis <- c("FGF8","PITX2","HESX1","SIX6")
marker$craniofacial <- c("ALX1","ALX3","PITX2","TCF21","LHX2","MSC","DLX2","DLX3","DLX4","DLX5","DLX6","DLX1","DLX2","DLX1","DLX2","DLX5","DLX6","HAND2",
                         "DLX3","DLX4","DLX5","DLX6","HOXA2","HOXA3","DLX1","DLX3","HOXA3","HOXB2","DLX1","DLX2","DLX5","DLX6","HOXA2","HOXA3")
marker$hm <- c("PITX2","MYF5","MSX2","PRRX1","CYP26C1","CSPG4","KCNJ8","ABCC9")
marker$somite <- c("PAX1","PAX9","PAX1","PAX9","NKX3-2","SOX9","SOX9","SCX","ETV4","PAX3","DMRT2","MYF5","MYF6","MYOD1",
                   "LBX1","PAX3","SIX1","FOXC1","FOXC2","SHISA2","WNT5A","CYP26A1","SOX2","MESP1")
marker$im <- c("PAX8","EMX2","LHX1")
marker$somatic.lpm <- c("IRX3","GATA6")
marker$limb <- c("TBX5","PITX1")
marker$splanchnic.lpm <- c("TBX1","WNT2","NKX6-1","LHX2","HLX","FOXF1","BARX1","ISL1","NKX2-5","TBX5","NPPA","MYH6","NKX2-5",
                           "MYL2","MYH7","SLC8A1","GATA6","NKX2-5","HAND2","HEY2","SHOX2","TBX18","TBX18","WT1","TNNI3","TNNT1",
                           "TCF21","SNAI2","TCF21","SOX9","THY1","HAND2","NPPA","GATA4","TGFBI","HEY2","TWIST1")
marker$endothelium <- c("CDH5","CD34","PECAM1","CDH5","CD34","PECAM1","NOTCH1","NOTCH4","DLL4",
                        "HEY1","HEY2","TMEM100","GJA5","GJA4","ELN","XBP1","STAB2","CD36")
marker$blood <- c("CLEC10A","CD1C","PLAC8","HDC","LMO4","CLC","CD34","LTB","IL7R","CD52","JCHAIN","CD68","CD14","CD163","CSF1R","ITGAM","ID3",
                  "NR1H3","CD9","PECAM1","ITGA2B","GP1BA","PPBP","S100A8","S100A9","LYZ","PRTN3","AZU1","MPO","HBA1","GYPB")
marker$endoderm <- c("PAX1","HOXA3","GLI3","SOX2","SOX2","KLF5","NKX2-1","SOX2","NKX2-1","SOX9","ETV5","MYCN","IRX1","IRX2",
                     "IRX3","HHEX","KLF5","ONECUT2","FOXA1","PDX1","PTF1A","NKX6-1","CDX2","GATA4","XBP1","CEBPA","ALB")
marker$pgc <- c("POU5F1","NANOS3","DPPA3","TFAP2C")
marker$fibroblast <- c("COL1A2","ALB")
marker$neuron <- c("ASCL1","ATOH1","BARHL1","BARHL2","BHLHE22","DLX1","DLX2","DLX5","DLX6-AS1","DMRT3","DRGX","EN1","EOMES","EVX1","FOXA2","FOXD3",
                   "FOXP1","FOXP2","GAD2","GATA2","GATA3","HOXB6","ISL1","ISL2","LBX1","LHX1","LHX2","LHX3","LHX5","LHX9","LMX1B","NKX2-1","NKX2-2",
                   "NKX6-1","OTP","PAX2","PHOX2A","PHOX2B","PITX2","POU4F1","SIM1","SIX1","SOX14","SOX21","TAL2","TBR1","TLX1","TLX3","VSX1","VSX2")
marker$np <- c("ASCL1","ATOH1","BARHL1","BARHL2","BMP4","CNPY1","DBX2","DCT","DLX1","DLX2","DLX5","DMRT3","DMRTA1","DMRTA2","EMX2","EN1","FGF17","FGF8",
               "FOXA2","FOXG1","GATA3","HOXB6","IRX3","LHX5","LMX1A","MAB21L2","MAFB","MEIS2","MITF","MSX1","MSX2","NKX2-1","NKX2-2","NKX2-8","NKX6-1",
               "NKX6-2","OC90","OLIG1","OLIG2","OLIG3","PAX2","PAX3","PAX6","PAX8","PHOX2B","PRDM13","PRDM8","PTF1A","RAX","RSPO1","RSPO3","SHH","SIX3",
               "SOX2","SP5","SPON1","TAL2","TBX3","VSX2","WNT2B","WNT4","WNT5A","ZIC1")
# - load marker genes (classified by cell type)
ct.mk <- list()
# cerebrum
ct.mk$cerebrum <- data.frame(CellType=c(rep("Excitatory neurons", 3), rep("Inhibitory neurons", 3), rep("Oligodendrocytes", 4),
                                        rep("Microglia", 3), rep("Astrocytes", 3), rep("Vascular endothelial cells", 2),
                                        rep("Limbic system neurons", 3), rep("SKOR2_NPSR1 positive cells", 2), rep("Megakaryocytes", 2)),
                             Marker=c("SLC17A7", "NEUROD6", "MAB21L1", "GAD1", "RELN", "CALB1", "MOG", "MAG", "SOX10", "OLIG1",
                                      "CX3CR1", "P2RY12", "TMEM119", "ALDH1L1", "SLC1A3", "AQP4", "CLDN5", "SLCO1C1", "TCF7L2", "LHX9", "NECAB2",
                                      "CALB1", "TFAP2A", "PPBP", "PF4"))
# cerebellum
ct.mk$cerebellum <- data.frame(CellType=c(rep("Inhibitory interneurons", 3), rep("Oligodendrocytes", 2), rep("Astrocytes", 3),
                                        rep("Granule neurons", 1), rep("Purkinje neurons", 4), rep("Microglia", 3),
                                        rep("Vascular endothelial cells", 2), rep("Unipolar brush cells", 1), rep("SLC24A4_PEX5L positive cells", 3)),
                               Marker=c("GAD1", "RELN", "PAX2", "MOG", "MAG", "ALDH1L1", "SLC1A3", "AQP4", "PAX6", "PCP4", "GAD1", "GAD2", "NECAB2",
                                        "CX3CR1", "P2RY12", "TMEM119", "CDH5", "KDR", "EOMES", "SLC24A4", "PEX5L", "L3MBTL4"))
# lung
ct.mk$lung <- data.frame(CellType=c(rep("Bronchiolar and alveolar epithelial cells", 6),
                                    rep("Stromal cells", 10), rep("Ciliated epithelial cells", 4),
                                    rep("Neuroendocrine cells", 1), rep("Squamous epithelial cells", 2), rep("Visceral neurons", 3),
                                    rep("Myeloid cells", 14), rep("Lymphoid cells", 8),
                                    rep("Megakaryocytes", 4), rep("Vascular endothelial cells", 9),
                                    rep("Lymphatic endothelial cells", 2), rep("Mesothelial cells", 2), rep("CSH1_CSH2 positive cells", 2)),
                          Marker=c("EPCAM","CDH1","AGER","SFTPB","SFTPA1","SFTPC","PDGFRB","COL1A1",
                                   "DCN","ACTA2","TAGLN","MYH11","COL3A1","SNAI2",
                                   "TCF21","LUM","TPPP3","FOXJ1","CD24","SCGB1A1","ASCL1",
                                   "TMPRSS11B","KRT13","PHOX2B","PRPH","CHRNA3","H2-AA","CD24",
                                   "CLEC10A","CST3","IRF8","CCL2","FSCN1","ITGAX","CD68",
                                   "C1QA","C1QC","C1QB","CD68","CD52","CD3D","CD3G","LEF1","CD3E",
                                   "GZMA","CCL5","NKG7","KLRB1","PPBP","PF4","ITGA2B",
                                   "ITGB3","PECAM1","CD34","FLI1","CAV1","CD93","KDR","TMEM100","CLDN5",
                                   "VWF","LYVE1","PDPN","COL8A1","MSLN","CSH1","CSH2"))
# heart
ct.mk$heart <- data.frame(CellType=c(rep("Cardiomyocytes", 5), rep("Stromal cells", 5),
                                     rep("Vascular endothelial cells", 5), rep("Smooth muscle cells", 3),
                                     rep("Endocardial cells", 1), rep("Epicardial fat cells", 1),
                                     rep("Lymphoid cells", 9), rep("Myeloid cells", 11),
                                     rep("Lymphatic endothelial cells", 2), rep("Visceral neurons", 3),
                                     rep("Schwann cells", 2), rep("SATB2_LRRC7 positive cells", 2),
                                     rep("ELF3_AGBL2 positive cells", 8), rep("CLC_IL5RA positive cells", 4),
                                     rep("Megakaryocytes", 4), rep("Erythroblasts", 3)),
                          Marker=c("NPPA","MYL7","TNNI3","TNNT2","ACTN2","TCF21","COLLA2","COL1A1","DCN",
                                   "SOX9","FABP4","CDH5","CAV1","CDH5","KDR","ACTA2","TAGLN",
                                   "MYH11","NPR3","PRG4","CD79A","CD79B","MS4A1","GZMA","CCL5","NKG7",
                                   "CD3D","CD3G","LEF1","C1QA","C1QC","C1QB","IRF8","CCL2","FSCN1",
                                   "ITGAX","CCL3","CCL9","S100A8","S100A9","LYVE1","PDPN","PHOX2B","PRPH",
                                   "CHRNA3","MPZ","PLP1","SATB2","DISP3","SCGB3A2","SFTPB",
                                   "ZMYND10","AGR2","ELF3","IL20RA","CFAP57","AGBL2","CLC","IL5RA","PRG2",
                                   "MS4A3","PPBP","PF4","ITGA2B","ITGB3","GYPA","HBA","HBB"))
# placenta
ct.mk$placenta <- data.frame(CellType=c(rep("Syncytiotrophoblasts and villous cytotrophoblasts", 4),
                                        rep("Stromal cells", 4), rep("Vascular endothelial cells", 2),
                                        rep("Trophoblast giant cells", 3), rep("Myeloid cells", 11),
                                        rep("Extravillous trophoblasts", 5), rep("IGFBP1_DKK1 positive cells", 2),
                                        rep("Lymphoid cells", 9),rep("Smooth muscle cells", 3),
                                        rep("PAEP_MECOM positive cells", 2), rep("AFP_ALB positive cells", 3),
                                        rep("Megakaryocytes", 4)),
                             Marker=c("SLC22A11","PAGE4","PARP1","PERP","COL3A1","LUM","COL1A1","DCN",
                                      "CDH5","KDR","PSG4","PSG9","PSG6","C1QA","C1QC","C1QB","IRF8","CCL2",
                                      "FSCN1","ITGAX","CCL3","CCL9","S100A8","S100A9","HLA-G","DIO2",
                                      "HTRA4","MFAP5","HPGD","IGFBP1","DKK1","CD79A","CD79B","MS4A1","GZMA",
                                      "CCL5","NKG7","CD3D","CD3G","LEF1","ACTA2","TAGLN","MYH11","PAEP",
                                      "MECOM","AFP","ALB","APOA2","PPBP","PF4","ITGA2B","ITGB3"))
# adrenal
ct.mk$adrenal <- data.frame(CellType=c(rep("Stromal cells", 4), rep("Vascular endothelial cells", 2),
                                       rep("Lymphoid cells", 9), rep("Myeloid cells", 11),
                                       rep("Megakaryocytes", 4), rep("Erythroblasts", 3),
                                       rep("CSH1_CSH2 positive cells", 4), rep("SLC26A4_PAEP positive cells", 4),
                                       rep("Adrenocortical cells", 2), rep("Sympathoblasts", 3),
                                       rep("Schwann cells", 2), rep("Chromaffin cells", 3)),
                            Marker=c("COL3A1","LUM","COL1A1","DCN","CDH5","KDR","CD79A","CD79B",
                                     "MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G","LEF1","C1QA","C1QC","C1QB",
                                     "IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9","S100A8","S100A9",
                                     "PPBP","PF4","ITGA2B","ITGB3","GYPA","HBA","HBB","CSH1","CSH2","PAPPA2",
                                     "CGA","SLC26A4","PAEP","VMO1","CXCL5","APOE","APOA1","CARTPT",
                                     "PNMT","DBH","MPZ","PLP1","TH","PHOX2B","CHGB"))
# liver
ct.mk$liver <- data.frame(CellType=c(rep("Vascular endothelial cells", 2), rep("Lymphoid cells", 9),
                                     rep("Myeloid cells", 11), rep("Megakaryocytes", 4),
                                     rep("Erythroblasts", 3), rep("Stellate cells", 6),
                                     rep("Hepatoblasts", 11), rep("Mesothelial cells", 2),
                                     rep("Hematopoietic stem cells", 1)),
                          Marker=c("CDH5","KDR","CD79A","CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D",
                                   "CD3G","LEF1","C1QA","C1QC","C1QB","IRF8","CCL2","FSCN1",
                                   "ITGAX","CCL3","CCL9","S100A8","S100A9","PPBP","PF4","ITGA2B",
                                   "ITGB3","GYPA","HBA","HBB","ACTA2","COL1A1","COL1A2","TAGLN",
                                   "COL3A1","SPARC","ALB","TTR","APOA1","SERPINA1C","KRT19","EPCAM",
                                   "KRT8","KRT18","AFP","GPC3","DLK1","COL8A1","MSLN","CD34"))
# muscle
ct.mk$muscle <- data.frame(CellType=c(rep("Vascular endothelial cells", 2), rep("Lymphoid cells", 9),
                                      rep("Myeloid cells", 11), rep("Smooth muscle cells", 3),
                                      rep("Lymphatic endothelial cells", 2), rep("Schwann cells", 2),
                                      rep("Megakaryocytes", 2), rep("Erythroblasts", 3),
                                      rep("Stromal cells", 4), rep("Skeletal muscle cells", 2), rep("Satellite cells", 3)),
                           Marker=c("CDH5","KDR","CD79A","CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G",
                                    "LEF1","C1QA","C1QC","C1QB","IRF8","CCL2","FSCN1","ITGAX",
                                    "CCL3","CCL9","S100A8","S100A9","ACTA2","TAGLN","MYH11","LYVE1","PDPN",
                                    "MPZ","PLP1","PPBP","PF4","GYPA","HBA","HBB","COL3A1","LUM",
                                    "COL1A1","DCN","MYH3","MYH8","CALCR","MYF5","PAX7"))
# spleen
ct.mk$spleen <- data.frame(CellType=c(rep("Stromal cells", 4), rep("Vascular endothelial cells", 2),
                                      rep("Lymphoid cells", 9), rep("Myeloid cells", 11),
                                      rep("Megakaryocytes", 4), rep("Erythroblasts", 3),
                                      rep("Mesothelial cells", 2), rep("AFP_ALB positive cells", 3),
                                      rep("STC2_TLX1 positive cells", 3)),
                           Marker=c("COL3A1","LUM","COL1A1","DCN","CDH5","KDR","CD79A","CD79B",
                                    "MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G","LEF1","C1QA",
                                    "C1QC","C1QB","IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9",
                                    "S100A8","S100A9","PPBP","PF4","ITGA2B","ITGB3","GYPA","HBA",
                                    "HBB","COL8A1","MSLN","AFP","ALB","APOA2","MARCKS","BGN","EMILIN1"))
# thymus
ct.mk$thymus <- data.frame(CellType=c(rep("Antigen presenting cells", 5), rep("Thymic epithelial cells", 4),
                                      rep("Vascular endothelial cells", 2), rep("Stromal cells", 2), rep("Thymocytes", 5)),
                           Marker=c("HLA-DRB1","HLA-DPA1","C1QA","C1QC","C1QB","KRT5","KRT8","CDH1",
                                    "ASCL1","CDH5","KDR","GPC3","DLK1","IL2RB","RUNX3","PTCRA","RAG2","RAG1"))
# eye
ct.mk$eye <- data.frame(CellType=c(rep("Retinal progenitors and Muller glia", 8),
                                   rep("Ganglion cells", 4), rep("Amacrine cells", 2), rep("Horizontal cells", 4),
                                   rep("Photoreceptor cells", 4), rep("Bipolar cells", 3),
                                   rep("Microglia", 3), rep("Astrocytes", 3),
                                   rep("Vascular endothelial cells", 2), rep("Stromal cells", 4),
                                   rep("Smooth muscle cells", 3), rep("Retinal pigment cells", 2),
                                   rep("Lens fibre cells", 3), rep("Corneal and conjunctival epithelial cells", 4), rep("Skeletal muscle cells", 3),
                                   rep("PDE11A_TAFA2 positive cells", 2)),
                        Marker=c("SOX2","HES1","NAP1L1","DEK","RLBP1","SOX9","CRYM","GATB","GAP43","POU4F2",
                                 "PRPH","FKBP1B","GAD2","LHX9","LHX1","ONECUT1","GRIA3","TMOD1",
                                 "CRX","AIPL1","MATK","MAOA","VSX1","VSX2","TRPM1","CX3CR1","C1QC","GPR34",
                                 "GFAP","ANGPT1","IFI44L","CDH5","KDR","COL3A1","LUM","COL1A1","DCN",
                                 "MYH11","ACTA2","TAGLN","ENPP2","PLD5","CRIM1","CRYBB3","CRYAB","KRT12",
                                 "KRT13","KRT4","AQP3","MYH3","NEB","ACTC1","PDE11A","TAFA2"))
# pancreas
ct.mk$pancreas <- data.frame(CellType=c(rep("Acinar cells", 2), rep("Stromal cells", 4),
                                        rep("Ductal cells", 2), rep("Lymphoid cells", 9),
                                        rep("Vascular endothelial cells", 2), rep("Islet endocrine cells", 11),
                                        rep("ENS glia", 2), rep("ENS neurons", 3),
                                        rep("Erythroblasts", 3), rep("Myeloid cells", 11),
                                        rep("CCL19_CCL21 positive cells", 3), rep("Mesothelial cells", 2),
                                        rep("Lymphatic endothelial cells", 2), rep("Smooth muscle cells", 3)),
                             Marker=c("CPA1","PTF1A","COL3A1","LUM","COL1A1","DCN","KRT19","HNF1B","CD79A",
                                      "CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G","LEF1","CDH5","KDR",
                                      "MAFA","PRSS53","NKX6-1","PDX1","INS","GCG","SST","PDX1","ISL1","PPY",
                                      "ARX","MPZ","PLP1","DSCAM","FGF13","CHRNA3","GYPA","HBA","HBB","C1QA",
                                      "C1QC","C1QB","IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9","S100A8",
                                      "S100A9","CCL19","CCL21","PDPN","COL8A1","MSLN","STAB2","MMRN1","MYH11","ACTA2","TAGLN"))
# kidney
ct.mk$kidney <- data.frame(CellType=c(rep("Vascular endothelial cells", 2), rep("Mesangial cells", 2),
                                      rep("Metanephric cells", 12), rep("Megakaryocytes", 2),
                                      rep("Lymphoid cells", 9), rep("Myeloid cells", 11), rep("Erythroblasts", 3), rep("Stromal cells", 4),
                                      rep("Ureteric bud cells", 5)),
                           Marker=c("CDH5","KDR","VIM","ACTA2","PODXL","PTPRO","WT1","VIL1","RNF24",
                                    "COL19A1","SLC5A12","SLC12A1","KCNE1","SCNN1A","AQP2","SLC12A3","PPBP","PF4",
                                    "CD79A","CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G","LEF1",
                                    "C1QA","C1QC","C1QB","IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9","S100A8",
                                    "S100A9","GYPA","HBA","HBB","COL3A1","LUM","COL1A1","DCN","PSCA","HPGD","DHRS2","UPK1A","S100P"))
ct.mk$kidney <- rbind(ct.mk$kidney,
                      data.frame(CellType=c(rep("ureteric epithelial cells",1), rep("ureteric bud tip", 2), rep("renal stromal cell", 7),
                                            rep("cap mesenchyme", 6), rep("renal vesicle", 5), rep("ureteric stalk", 1)),
                                 Marker=c("HOXB7", "RET", "WNT11", "FAT1", "FGF8", "FOXD1", "PDX1", "TCF21", "RARA", "RARB",
                                          "FGFR1", "FGFR2", "GDNF", "SIX2", "WT1", "WNT4", "FGF8", "PAX8", "WNT4", "WT1", "SIX2", "WNT9B")))
# stomach
ct.mk$stomach <- data.frame(CellType=c(rep("Goblet cells", 4), rep("Parietal and chief cells", 8),
                                       rep("Squamous epithelial cells", 1), rep("Stromal cells", 4),
                                       rep("MUC13_DMBT1 positive cells", 4), rep("Lymphoid cells", 9),
                                       rep("Vascular endothelial cells", 2), rep("PDE1C_ACSM3 positive cells", 3),
                                       rep("Myeloid cells", 11), rep("Erythroblasts", 3), rep("ENS glia", 2), rep("ENS neurons", 3),
                                       rep("Ciliated epithelial cells", 3), rep("Neuroendocrine cells", 5),
                                       rep("Lymphatic endothelial cells", 2), rep("Mesothelial cells", 2)),
                            Marker=c("GKN2","TFF1","MUCL3","TFF2","APOA1","ATP4B","ATP4A","FABP3","LDHB",
                                     "CKB","PGC","GIF","KRT4","COL3A1","LUM","COL1A1","DCN","MUC13","DMBT1",
                                     "MT1G","APOB","CD79A","CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D",
                                     "CD3G","LEF1","CDH5","KDR","PDE1C","ACSM3","TFCP2L1","C1QA","C1QC","C1QB",
                                     "IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9","S100A8","S100A9",
                                     "GYPA","HBA","HBB","MPZ","PLP1","DSCAM","FGF13","CHRNA3","CDHR3","DNAH12","DNAH3",
                                     "GAST","GHRL","CPE","CHGB","CHGA","STAB2","MMRN1","COL8A1","MSLN"))
# intestine
ct.mk$intestine <- data.frame(CellType=c(rep("Intestinal epithelial cells", 7), rep("Stromal cells", 4), rep("ENS glia", 2), rep("ENS neurons", 3),
                                         rep("Myeloid cells", 11), rep("Lymphoid cells", 9),
                                         rep("Vascular endothelial cells", 2), rep("Smooth muscle cells", 3),
                                         rep("Chromaffin cells", 2), rep("Erythroblasts", 3),
                                         rep("Lymphatic endothelial cells", 2), rep("Mesothelial cells", 2)),
                              Marker=c("KRT20","SLC26A3","SPDEF","TFF3","MANF","CCL9","LGR5","COL3A1","LUM","COL1A1",
                                       "DCN","MPZ","PLP1","DSCAM","FGF13","CHRNA3","C1QA","C1QC",
                                       "C1QB","IRF8","CCL2","FSCN1","ITGAX","CCL3","CCL9","S100A8","S100A9","CD79A",
                                       "CD79B","MS4A1","GZMA","CCL5","NKG7","CD3D","CD3G","LEF1",
                                       "CDH5","KDR","ACTA2","TAGLN","MYH11","CHGA","CHGB","GYPA","HBA","HBB","STAB2","MMRN1","COL8A1","MSLN"))


### >>> 3. Quality control
# - nucleosome_signal: Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
# - compute per-cell quality control metrics
DefaultAssay(hs.mpr) <- "ATAC"
hs.mpr <- NucleosomeSignal(hs.mpr)
hs.mpr <- TSSEnrichment(hs.mpr, fast = F)
DefaultAssay(hs.mpr) <- "RNA"
hs.mpr[["percent.mt"]] <- PercentageFeatureSet(hs.mpr, pattern = "^MT-")
hs.mpr$Sample <- str_split_fixed(rownames(hs.mpr@meta.data), "_", 2)[,1]
hs.mpr$Stage <- str_split_fixed(rownames(hs.mpr@meta.data), "\\.", 3)[,2]
pdf(file.path(res.out, "Quality_control_before_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = hs.mpr,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 7, split.by = "Stage", pt.size = 0, group.by = "Stage")
dev.off()
# filter out low quality cells
hs.mpr <- subset(
  x = hs.mpr,
  subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000 &
    nucleosome_signal <= 2 &
    TSS.enrichment >= 3 &
    percent.mt <= 30
)
pdf(file.path(res.out, "Quality_control_after_filtering.pdf"), height = 5, width = 20)
VlnPlot(object = hs.mpr,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 7, pt.size = 0, group.by = "Stage")
dev.off()
DefaultAssay(ss.embryo) <- "ATAC"
pdf(file.path(res.out, "tss_enrichment_and_fragment_histogram_after_filtering.pdf"), height = 5, width = 12)
p1 <- TSSPlot(hs.mpr, group.by = "orig.ident") + NoLegend() + scale_color_manual(values = c("#0984e3"))
p2 <- FragmentHistogram(object = hs.mpr, group.by = "orig.ident") + scale_fill_manual(values = c("#0984e3"))
p1+p2
dev.off()

### >>> 4. Peak calling (No Run)
DefaultAssay(hs.mpr) <- "ATAC"
# - call peaks using MACS2
peaks <- CallPeaks(hs.mpr, macs2.path = "/home/yhw/software/anaconda3/envs/dna.p37/bin/macs2")
# - remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# - quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(hs.mpr),
  features = peaks,
  cells = colnames(hs.mpr)
)
# - create a new assay using the MACS2 peak set and add it to the Seurat object
hs.mpr[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)


### >>> 7. Annotating super-cluster cell types by Seurat tools
# - post-processing
DefaultAssay(hs.mpr) <- "RNA"
hs.mpr <- NormalizeData(hs.mpr) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.mpr <- ScaleData(hs.mpr, verbose = FALSE, features = rownames(hs.mpr))
hs.mpr <- RunPCA(hs.mpr, npcs = 50, verbose = FALSE, features = VariableFeatures(hs.mpr))
ElbowPlot(hs.mpr, ndims = 30)
dim.n <- 25
# - dimension reduction and clustering
hs.mpr <- FindNeighbors(hs.mpr, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr <- RunUMAP(hs.mpr, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
# - visualize clustering
pdf(file.path(res.out, "Scatter_plot_to_show_uncorrected_clustering.pdf"), height = 5, width = 12.5)
DimPlot(hs.mpr, reduction = "umap", group.by = c("Sample", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# - remove batch effect
library(harmony)
hs.mpr <- RunHarmony(hs.mpr, group.by.vars = "Sample")
hs.mpr <- RunUMAP(hs.mpr, reduction = "harmony", dims = 1:dim.n)
hs.mpr <- FindNeighbors(hs.mpr, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_corrected_clustering.pdf"), height = 5, width = 12)
DimPlot(hs.mpr, reduction = "umap", group.by = c("Sample", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.25)
dev.off()
FeaturePlot(hs.mpr, features = c("nCount_RNA", "nFeature_RNA"))
# - find marker genes
hs.mpr.markers <- FindAllMarkers(hs.mpr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# - define 4w cells
pd.gene <- hs.mpr.markers %>% dplyr::filter(cluster == "3") %>% dplyr::top_n(12, avg_log2FC)
FeaturePlot(hs.mpr, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 4)
pd.gene <- hs.mpr.markers %>% dplyr::filter(cluster == "32") %>% dplyr::top_n(12, avg_log2FC)
FeaturePlot(hs.mpr, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 4)
pd.gene <- hs.mpr.markers %>% dplyr::filter(cluster == "22") %>% dplyr::top_n(12, avg_log2FC)
FeaturePlot(hs.mpr, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 4)
FeaturePlot(hs.mpr, features = c("HOXC8","HHEX"),
            slot = "scale.data", reduction = "umap", ncol = 4)
pd.gene <- c("SPINK1","BRI3","NPW","FST","MAP2K2",
             "MT1H","IGF2","AMN","APOM","DLK1","RSPO3","GSTP1","BAAT",
             "AFP","CYP3A7","F2","AHSG","APOA1","LIPC","AGT",
             "ALB","ANGPTL3","CFHR1","LEPR","APOB","GC","VTN","AKR1C1",
             "CYP2E1","SAA1","HP","C9","ADH1B","F9","CES1","TAT","AZGP1","APCS","NNMT")
DotPlot(hs.mpr, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# - define super cluster (Human A single-cell transcriptome atlas of human early embryogenesis)
pd.gene <- c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3","PERP","KRT18","FOXD3","MPZ","PLP1","S100B","ALX1","ALX3",
             "DLX1","DLX2","SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9","PAX8","WT1","TBX5","PITX1","FOXF1","GATA6",
             "GATA5","NKX2-5","PLVAP","CD53","PLAC8","CD68","CORO1A","FOXA2","CLDN6","EMX2","PAX2","OSR1")
pdf(file.path(res.out, "Dotplot_to_show_human_super_cluster_markers.pdf"), height = 10.5, width = 14)
DotPlot(hs.mpr, features = pd.gene, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
png(file.path(res.out, "Featureplot_to_show_human_super_cluster_markers.png"), res = 300, width = 7500, height = 7000)
FeaturePlot(hs.mpr, features = pd.gene,
            slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# nervous system
FeaturePlot(hs.mpr, features = c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("5","11","14","15","28","29","32","34","35") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# schwann
FeaturePlot(hs.mpr, features = c("FOXD3", "MPZ", "PLP1", "S100B"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("25") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# endothelium
FeaturePlot(hs.mpr, features = c("PLVAP"), slot = "scale.data", reduction = "umap", ncol=4, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("10","20","23","27") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# blood
FeaturePlot(hs.mpr, features = c("CD53","PLAC8","CD68","CORO1A"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("13","26","36") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# Epidermis
FeaturePlot(hs.mpr, features = c("PERP","KRT18","CLDN6","PAWR","FOXA2"), slot = "scale.data", reduction = "umap", order = T)
FeaturePlot(hs.mpr, features = c("MET","CDH1"), slot = "scale.data", reduction = "umap")
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("30") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# endoderm
FeaturePlot(hs.mpr, features = c("SPINK1","MT1H","APOM","AFP","APOA1","ALB","AGT","VTN","DCDC2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("22", "31", "37") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# craniofacial
FeaturePlot(hs.mpr, features = c("ALX1","ALX3","DLX1","DLX2","SNAI2","TWIST1"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("0","1","9","12") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# head mesoderm
FeaturePlot(hs.mpr, features = c("SNAI2","TWIST1","FOXC1","FOXC2"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("4") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# somite
FeaturePlot(hs.mpr, features = c("SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("3","6") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# Intermediate mesoderm
FeaturePlot(hs.mpr, features = c("PAX8","EMX2","PAX2","OSR1"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("8") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# Lateral plate mesoderm
FeaturePlot(hs.mpr, features = c("GATA5","GATA6"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("16", "18", "19") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# limb
FeaturePlot(hs.mpr, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol=4, order = T)
# Heart
FeaturePlot(hs.mpr, features = c("TNNI3","TNNT2","ACTN2","TAGLN","MYL4","MYL7","TCF21","SMARCD3"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("24") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# Erythroid
FeaturePlot(hs.mpr, features = c("HBA1","HBA2","HBG1","HBG2"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("33") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# Mesenchyme
FeaturePlot(hs.mpr, features = c("KRT8","KRT18","AHNAK","BMP4","TAGLN","ACTA2"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("7") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("24") ~ as.character(seurat_clusters), TRUE ~ "Non")) %>%
  mutate(Clustering=factor(Clustering, levels = c("Non","24")))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5, order = T)
# visualization of lineage clustering
hs.mpr@meta.data <- hs.mpr@meta.data %>%
  mutate(Lineage=case_when(seurat_clusters%in%c(0:4,6,8:10,12:13,16:21,23:24,26:27,33,36,39) ~ "Mesoderm",
                           seurat_clusters%in%c(5,11,14,15,25,28:29,32,34:35) ~ "Ectoderm",
                           seurat_clusters%in%c(22,31,37) ~ "Endoderm",
                           seurat_clusters%in%c(7,30,38) ~ "Undefined"))
DimPlot(hs.mpr, reduction = "umap", group.by = c("Lineage"), label = TRUE, repel = TRUE, pt.size = 0.25)



### ====================================
### 4th step: Annotation of sub-clusters ----
### ====================================

### >>> 1. Ectoderm1 (neural progenitor + neuron): 1st round clustering ----
hs.mpr.sub1 <- list()
hs.mpr.sub1$ecto1 <- subset(hs.mpr, seurat_clusters%in%c("5","11","14","15","28","29","32","34","35"))
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$ecto1 <- FindVariableFeatures(hs.mpr.sub1$ecto1) %>% ScaleData(rownames(hs.mpr.sub1$ecto1)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$ecto1 <- RunHarmony(hs.mpr.sub1$ecto1, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$ecto1, reduction = "harmony", ndims = 30)
dim.n <- 20
hs.mpr.sub1$ecto1 <- RunUMAP(hs.mpr.sub1$ecto1, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_nervous_system_corrected_clustering.pdf"), height = 5, width = 17)
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# - visualization of specific cluster
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("22") ~ "22", TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - find markers
hs.mpr.sub1.markers <- list()
hs.mpr.sub1.markers$ecto1 <- FindAllMarkers(hs.mpr.sub1$ecto1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# - plot markers
FeaturePlot(hs.mpr.sub1$ecto1, features = c("OTX2","FEZF2","SOX2","PAX6","HES5","ELAVL4","TUBB3"), slot = "data", reduction = "umap", ncol = 4)
pdf(file.path(res.out, "Dotplot_to_show_nervous_system_subcluster_markers.pdf"), height = 10, width = 25.5)
DotPlot(hs.mpr.sub1$ecto1, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(marker$nv))) +
  labs(x = "Hs Nervous System Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub1$ecto1, features = sort(unique(marker$nv)), slot = "data", reduction = "umap", ncol = 6)
dev.off()
# - define cell type
# c0
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 0) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("0") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c1,18
FeaturePlot(hs.mpr.sub1$ecto1, features = c("ASCL1", "PAX3", "HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("1","18") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c2
FeaturePlot(hs.mpr.sub1$ecto1, features = c("LHX1", "LHX5", "PAX2"), slot = "scale.data", reduction = "umap", ncol = 2)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("2") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c3
FeaturePlot(hs.mpr.sub1$ecto1, features = c("DMRTA1","DMRTA2","EMX2","PAX6","DMRT3"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("3") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c4
FeaturePlot(hs.mpr.sub1$ecto1, features = c("PAX3","PRDM13","ASCL1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("4") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c5,7,9
FeaturePlot(hs.mpr.sub1$ecto1, features = c("LHX1","LHX5","POU4F1","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("5","7","9") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c6
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 6) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("6") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c8
FeaturePlot(hs.mpr.sub1$ecto1, features = c("ISL1","ISL2","LHX3","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("8") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
VlnPlot(hs.mpr.sub1$ecto1, features = c("ISL1","ISL2","LHX3","HOXB6"), flip = T, stack = T)
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 8) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
# c10
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 10) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("10") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c11
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 11) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("11") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c12
FeaturePlot(hs.mpr.sub1$ecto1, features = c("PAX6","DBX2"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 12) %>% dplyr::top_n(n = 20, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 5)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("12") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c13
FeaturePlot(hs.mpr.sub1$ecto1, features = c("WNT2B","WNT4"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.sub1$ecto1, features = c("DMRTA1","DMRTA2","EMX2","PAX6","DMRT3"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("13") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c15,17,22
FeaturePlot(hs.mpr.sub1$ecto1, features = c("LHX2","LHX9","BARHL2","BARHL1","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("15","17","22") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c16
FeaturePlot(hs.mpr.sub1$ecto1, features = c("MSX1","OLIG3","ATOH1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("16") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c19
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 19) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("19") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
VlnPlot(hs.mpr.sub1$ecto1, features = c("DLX1","DLX2","DLX5", "GATA2","GATA3"), stack = T, flip = T)
# c20,25
FeaturePlot(hs.mpr.sub1$ecto1, features = c("SIX1","ISL1","TLX3","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.sub1$ecto1, features = c("TLX3","POU4F1","ISL1","PHOX2B"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("20","25") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c21
FeaturePlot(hs.mpr.sub1$ecto1, features = c("LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.sub1$ecto1, features = c("LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3"),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("21") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c23
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 23) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = c("TLX3","PHOX2B","POU4F1","ISL1"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("23") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c24
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 24) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("24") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c26
FeaturePlot(hs.mpr.sub1$ecto1, features = c("LBX1","LHX1","LHX5","PAX2","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("26") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c27
FeaturePlot(hs.mpr.sub1$ecto1, features = c("NKX6-1","NKX2-2","NKX2-8","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("27") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c28
hs.mpr.sub1.markers$ecto1 %>% dplyr::filter(cluster == 28) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto1, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("28") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c29
FeaturePlot(hs.mpr.sub1$ecto1, features = c("NKX6-1","FOXA2","SHH","SPON1","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("29") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c30
FeaturePlot(hs.mpr.sub1$ecto1, features = c("NKX6-1","PHOX2A","PHOX2B","ISL1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("30") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c31
FeaturePlot(hs.mpr.sub1$ecto1, features = c("VSX1","ASCL1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("31") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - defined classified/unclassified clusters
hs.mpr.sub1$ecto1@meta.data <- hs.mpr.sub1$ecto1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(23,30,0,7,10,11,19,20,25) ~ "UnClassified", TRUE ~ "Classified"))
DimPlot(hs.mpr.sub1$ecto1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - marker visualization
FeaturePlot(hs.mpr.sub1$ecto1, features = c("ATOH1","BHLHE22","DRGX","FOXD3","FOXP2","LBX1","LHX1","LHX5","LMX1B","PAX2","PHOX2B","TLX3"))
DotPlot(hs.mpr.sub1$ecto1, features = c("ATOH1","BHLHE22","DRGX","FOXD3","FOXP2","LBX1","LHX1","LHX5","LMX1B","PAX2","PHOX2B","TLX3"),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ecto1, features = c("FOXD3","FOXP2"))
DotPlot(hs.mpr.sub1$ecto1, features = c("FOXD3","FOXP2"),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp


### >>> 2. Ectoderm1: 2nd round clustering ----
# - cluster 2,26
hs.mpr.sub2 <- list()
hs.mpr.sub2$nv1 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("2", "26"))
hs.mpr.sub2$nv1@meta.data <- hs.mpr.sub2$nv1@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv1 <- FindVariableFeatures(hs.mpr.sub2$nv1) %>% ScaleData(rownames(hs.mpr.sub2$nv1)) %>% RunPCA(verbose = FALSE)
set.seed(200)
hs.mpr.sub2$nv1 <- RunHarmony(hs.mpr.sub2$nv1, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv1, reduction = "harmony", ndims = 30)
dim.n <- 5
hs.mpr.sub2$nv1 <- RunUMAP(hs.mpr.sub2$nv1, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv1, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv1, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(marker$nv))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c2c26.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv1, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers <- list()
hs.mpr.sub2.markers$nv1 <- FindAllMarkers(hs.mpr.sub2$nv1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
hs.mpr.sub2$nv1@meta.data <- hs.mpr.sub2$nv1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("3","7") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv1, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv1, features = c("LHX1","LHX5","EN1","OTP","HOXB6","PAX2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv1, features = c("PAX6","IRX3","NKX6-1","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv1, features = c("LBX1","LHX1","LHX5","PAX2","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv1, features = c("GATA2","GATA3","LHX1","LHX5","HOXB6","VSX2", "SOX14", "SOX21"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv1, features = c("NKX2-1","DLX1","DLX2","DLX5","DLX6-AS1","GAD2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
VlnPlot(hs.mpr.sub2$nv1, features = c("LBX1","LHX1","LHX5","PAX2","HOXB6"), stack = T, flip = T)
VlnPlot(hs.mpr.sub2$nv1, features = c("NKX2-1","DLX1","DLX2","DLX5","DLX6-AS1","GAD2"), stack = T, flip = T)
DotPlot(hs.mpr.sub2$nv2, features = c("B3GAT1","AP-2A", "TWIST1", "SOX9", "MYC", "ETS1", "DLX1", "DLX2", "CRABP1", "EPHA2", "ITGB1"),
        cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv1, features = c("B3GAT1"))
DotPlot(hs.mpr.sub2$nv1, features = c("LBX1","PAX2","LHX1","LHX5","BHLHE22"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# - cluster 15,17,22
hs.mpr.sub2$nv2 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("15", "17", "22"))
hs.mpr.sub2$nv2@meta.data <- hs.mpr.sub2$nv2@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv2 <- FindVariableFeatures(hs.mpr.sub2$nv2) %>% ScaleData(rownames(hs.mpr.sub2$nv2)) %>% RunPCA(verbose = FALSE)
set.seed(210)
hs.mpr.sub2$nv2 <- RunHarmony(hs.mpr.sub2$nv2, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv2, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.sub2$nv2 <- RunUMAP(hs.mpr.sub2$nv2, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv2, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv2, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c15c17c22.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv2, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv2 <- FindAllMarkers(hs.mpr.sub2$nv2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
hs.mpr.sub2$nv2@meta.data <- hs.mpr.sub2$nv2@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(5) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv2, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv2, features = c("LHX2","LHX9","BARHL2","BARHL1","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv2, features = c("LHX1","LHX5","POU4F1","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv2, features = c("LBX1","TLX1","TLX3","PHOX2A","HOXB6","PHOX2B","POU4F1","ISL1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv2, features = c("MSX1","OLIG3","ATOH1","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
VlnPlot(hs.mpr.sub2$nv2, features = c("NKX2-1","DLX1","DLX2","DLX5","DLX6-AS1","GAD2"), stack = T, flip = T)
hs.mpr.sub2.markers$nv2 %>% dplyr::filter(cluster == 4) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv2, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = T)
DotPlot(hs.mpr.sub2$nv2, features = c("B3GAT1"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv2, features = c("B3GAT1"))
# - cluster 23,30
hs.mpr.sub2$nv3 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("23", "30"))
hs.mpr.sub2$nv3@meta.data <- hs.mpr.sub2$nv3@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv3 <- FindVariableFeatures(hs.mpr.sub2$nv3) %>% ScaleData(rownames(hs.mpr.sub2$nv3)) %>% RunPCA(verbose = FALSE)
set.seed(220)
hs.mpr.sub2$nv3 <- RunHarmony(hs.mpr.sub2$nv3, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv3, reduction = "harmony", ndims = 30)
dim.n <- 5
hs.mpr.sub2$nv3 <- RunUMAP(hs.mpr.sub2$nv3, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv3, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv3, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c23c30.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv3, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv3 <- FindAllMarkers(hs.mpr.sub2$nv3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
hs.mpr.sub2.markers$nv3 %>% dplyr::filter(cluster == 2) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv3, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$nv3, features = c("NKX6-1","PHOX2B","ISL1","PHOX2A"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv3, features = c("GATA2","GATA3","PHOX2A","PHOX2B","SPON1","TBX3"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv3, features = c("TLX3","PHOX2B","POU4F1","ISL1"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F, pt.size = 0.5)
hs.mpr.sub2$nv3@meta.data <- hs.mpr.sub2$nv3@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(0,1,3,5) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv3, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - cluster 8,24
hs.mpr.sub2$nv4 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("8", "24"))
hs.mpr.sub2$nv4@meta.data <- hs.mpr.sub2$nv4@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv4 <- FindVariableFeatures(hs.mpr.sub2$nv4) %>% ScaleData(rownames(hs.mpr.sub2$nv4)) %>% RunPCA(verbose = FALSE)
set.seed(230)
hs.mpr.sub2$nv4 <- RunHarmony(hs.mpr.sub2$nv4, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv4, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$nv4 <- RunUMAP(hs.mpr.sub2$nv4, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2.4)
DimPlot(hs.mpr.sub2$nv4, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv4, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c8c24.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv4, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv4 <- FindAllMarkers(hs.mpr.sub2$nv4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
hs.mpr.sub2.markers$nv4 %>% dplyr::filter(cluster == 0) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv4, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.sub2$nv4, features = c("NKX6-1","PHOX2B","ISL1","PHOX2A"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv4, features = c("ISL1","ISL2","FOXP1","LHX3","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv4, features = c("LHX1","LHX5","POU4F1","HOXB6","GATA2","GATA3","PAX2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv4, features = c("NKX6-1","ISL1","PHOX2B","PHOX2A"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv4, features = marker$sw, slot = "scale.data", reduction = "umap", ncol = 5, order = F, pt.size = 0.5)
hs.mpr.sub2$nv4@meta.data <- hs.mpr.sub2$nv4@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(2,6,7,9,10) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv4, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv4, features = c("SNAI2","TWIST1","TBX5","PITX1","PAX3","PMP22","TFAP2B","MITF","PLP1"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# - cluster 14,16
hs.mpr.sub2$nv5 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("14", "16"))
hs.mpr.sub2$nv5@meta.data <- hs.mpr.sub2$nv5@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv5 <- FindVariableFeatures(hs.mpr.sub2$nv5) %>% ScaleData(rownames(hs.mpr.sub2$nv5)) %>% RunPCA(verbose = FALSE)
set.seed(240)
hs.mpr.sub2$nv5 <- RunHarmony(hs.mpr.sub2$nv5, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv5, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.sub2$nv5 <- RunUMAP(hs.mpr.sub2$nv5, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv5, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv5, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c14c16.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv5, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv5 <- FindAllMarkers(hs.mpr.sub2$nv5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
hs.mpr.sub2.markers$nv5 %>% dplyr::filter(cluster == 7) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv5, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$nv5, features = c("MSX1","OLIG3","ATOH1","PHOX2B","ASCL1"), slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv5, features = c("PRDM13","PAX3","DBX2"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
hs.mpr.sub2$nv5@meta.data <- hs.mpr.sub2$nv5@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(3,5) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv5, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv5, features = c("ASCL1","PTF1A"), slot = "scale.data", reduction = "umap", ncol = 3, order = T, pt.size = 0.5)
# - cluster 4
hs.mpr.sub2$nv6 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("4"))
hs.mpr.sub2$nv6@meta.data <- hs.mpr.sub2$nv6@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv6 <- FindVariableFeatures(hs.mpr.sub2$nv6) %>% ScaleData(rownames(hs.mpr.sub2$nv6)) %>% RunPCA(verbose = FALSE)
set.seed(250)
hs.mpr.sub2$nv6 <- RunHarmony(hs.mpr.sub2$nv6, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv6, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.sub2$nv6 <- RunUMAP(hs.mpr.sub2$nv6, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 3.8)
DimPlot(hs.mpr.sub2$nv6, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv6, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c4.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv6, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv6 <- FindAllMarkers(hs.mpr.sub2$nv6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
FeaturePlot(hs.mpr.sub2$nv6, features = c("MSX1","OLIG3","ATOH1","PHOX2B","ASCL1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("FOXA2","PITX2","ISL2","ISL1","LHX3","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("LHX3","VSX2","SOX14","SOX21","VSX1","HOXB6", "PTF1A","ASCL1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("LBX1","LHX1","LHX5","PAX2","EN1","OTP","HOXB6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("NKX2-1","DLX1","DLX2","DLX5"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("SOX14","SOX21","VSX2","PRDM13","SIM1","NKX6-1","IRX3","OLIG3","LHX9"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
hs.mpr.sub2.markers$nv6 %>% dplyr::filter(cluster == 5) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv6, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.sub2$nv6@meta.data <- hs.mpr.sub2$nv6@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(5) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv6, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv6, features = c("PAX3","DBX2","PRDM13"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
# - cluster 3,13,28
hs.mpr.sub2$nv7 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("3", "28", "13"))
hs.mpr.sub2$nv7@meta.data <- hs.mpr.sub2$nv7@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv7 <- FindVariableFeatures(hs.mpr.sub2$nv7) %>% ScaleData(rownames(hs.mpr.sub2$nv7)) %>% RunPCA(verbose = FALSE)
set.seed(260)
hs.mpr.sub2$nv7 <- RunHarmony(hs.mpr.sub2$nv7, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv7, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.sub2$nv7 <- RunUMAP(hs.mpr.sub2$nv7, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv7, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv7, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c3c13c28.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv7, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv7 <- FindAllMarkers(hs.mpr.sub2$nv7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
FeaturePlot(hs.mpr.sub2$nv7, features = c("MITF","PAX2","DCT"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv7, features = c("LMX1A","SHH","FOXA2","WNT5A"), slot = "data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv7, features = c("WNT4","WNT2B","DMRTA1","DMRTA2","EMX2","PAX6","DMRT3"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv7, features = c("GYPA","GYPB","HBE1","HBZ","HBG2","HBA2","HBG1","HBA1","SLC25A21","ANK1","SLC4A1","SLC25A37"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv7, features = c("SIX3","FGF8"), slot = "scale.data", reduction = "umap", ncol = 5, order = F, pt.size = 0.5)
hs.mpr.sub2.markers$nv7 %>% dplyr::filter(cluster == 2) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv7, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
DotPlot(hs.mpr.sub2$nv7, features = c("GYPA","GYPB","HBE1","HBZ","HBG2","HBA2","HBG1","HBA1","SLC25A21","ANK1","SLC4A1","SLC25A37"),
        cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("GYPA","GYPB","HBE1","HBZ","HBG2","HBA2","HBG1","HBA1","SLC25A21","ANK1","SLC4A1","SLC25A37")) +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
hs.mpr.sub2$nv7@meta.data <- hs.mpr.sub2$nv7@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(1,4,5) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv7, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - cluster 1,6,12,18,27
hs.mpr.sub2$nv8 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("1", "6", "12", "18", "27"))
hs.mpr.sub2$nv8@meta.data <- hs.mpr.sub2$nv8@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv8 <- FindVariableFeatures(hs.mpr.sub2$nv8) %>% ScaleData(rownames(hs.mpr.sub2$nv8)) %>% RunPCA(verbose = FALSE)
set.seed(270)
hs.mpr.sub2$nv8 <- RunHarmony(hs.mpr.sub2$nv8, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv8, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$nv8 <- RunUMAP(hs.mpr.sub2$nv8, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv8, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv8, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c1c6c12c18c27.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv8, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv8 <- FindAllMarkers(hs.mpr.sub2$nv8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# visualization
hs.mpr.sub2.markers$nv8 %>% dplyr::filter(cluster == 11) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv8, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.sub2$nv8, features = c("PAX6","DBX2"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv8, features = c("ASCL1","PAX3","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv8, features = c("NKX6-1","NKX2-2","NKX2-8","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv8, features = c("PRDM8","NKX6-1","ASCL1"), slot = "scale.data", reduction = "umap", ncol = 5, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv8, features = c("MITF","MSX1","MSX2","PRDM13","ZIC1"), slot = "scale.data", reduction = "umap", ncol = 5, order = F, pt.size = 0.5)
hs.mpr.sub2$nv8@meta.data <- hs.mpr.sub2$nv8@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(0,2:5,7:9,11:12,14) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv8, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - cluster 0,7,10,11,19
hs.mpr.sub2$nv9 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("0", "7", "10", "11", "19"))
hs.mpr.sub2$nv9@meta.data <- hs.mpr.sub2$nv9@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv9 <- FindVariableFeatures(hs.mpr.sub2$nv9) %>% ScaleData(rownames(hs.mpr.sub2$nv9)) %>% RunPCA(verbose = FALSE)
set.seed(280)
hs.mpr.sub2$nv9 <- RunHarmony(hs.mpr.sub2$nv9, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv9, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$nv9 <- RunUMAP(hs.mpr.sub2$nv9, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 3)
DimPlot(hs.mpr.sub2$nv9, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv9, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c0c7c10c11c19.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv9, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# visualize specific cluster
hs.mpr.sub2$nv9@meta.data <- hs.mpr.sub2$nv9@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(6,11,15,20) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv9, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$nv9 <- FindAllMarkers(hs.mpr.sub2$nv9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$nv9 %>% filter(cluster == "2") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$nv9, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot markers
DotPlot(hs.mpr.sub2$nv9, features = unique(marker$super.cluster), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.sub2$nv9, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv9, features = c("WNT2B","WNT4"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv9, features = c("ASCL1","DBX2","HOXB6","IRX3","NKX2-2","NKX2-8","NKX6-1","NKX6-2","PAX3","PAX6",
                                          "PRDM13","PRDM8","MSX1","MSX2","OLIG3","PHOX2B","SOX2","SOX21","TBX3","WNT5A"),
            slot = "scale.data", reduction = "umap", ncol = 6, order = T, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv9, features = c("ASCL1","DBX2","HOXB6","IRX3","NKX2-2","NKX2-8","NKX6-1","NKX6-2","PAX3","PAX6",
                                      "PRDM13","PRDM8","MSX1","MSX2","OLIG3","PHOX2B","SOX2","SOX21","TBX3","WNT5A"),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.sub2$nv9, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.sub2$nv9, features = unique(c("CCNA1","HES2","DOK5","VSX1","ELFN1","GBE1","ASCL1","MYBL1","MT1F","NPTX2",
                                             "MT3","CDH13","FOXN4","VSX1","KIF19","C8ORF46","MYBL1","DOK5","DGKB","NEUROD1")), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv9, features = c("BARHL1","BARHL2","BHLHE22","DMRT3","EVX1","FOXD3","FOXP2","HOXB6","LBX1",
                                          "LHX1","LHX2","LHX5","LHX9","LMX1B","PAX2","PHOX2A","POU4F1","DRGX","TLX1","TLX3"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv9, features = c("BARHL1","BARHL2","BHLHE22","DMRT3","EVX1","FOXD3","FOXP2","HOXB6","LBX1",
                                      "LHX1","LHX2","LHX5","LHX9","LMX1B","PAX2","PHOX2A","POU4F1","DRGX","TLX1","TLX3"),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv9, features = c("ISL2","ISL1","FOXP1","HOXB6",
                                          "GATA2","GATA3",
                                          "SIX1","ISL1","TLX3",
                                          "EN1"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv9, features = unique(c("NTN5","SUSD2","USP18","ZNF385C","CLPTM1L","EYA2","PPP1R17","SDCCAG8","ASS1","FAM110A",
                                             "C2ORF91","TBC1D3B","HMX1","TRRHDE","CORO2A","THSD7B","CACNA1A","IQSEC1","NEDD1","PPP1R17",
                                             "PPP1R1C","NTRK1","CASP3","SYNPR","NTRK3","ADAM11","KCNJ12","BCAS3","NPR2","ZNF385C",
                                             "GKAP1","D21S2088E","S100A1","SST","TLX2","SELM","RPRM","CIB2","SIX1","THBS3")),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv9, features = c("GSX2"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv9, features = unique(marker$fibroblast), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100, col.min = -1) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv9, features = unique(marker$hm),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T, pt.size = 0.5)
# - cluster 20,25
hs.mpr.sub2$nv10 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("20", "25"))
hs.mpr.sub2$nv10@meta.data <- hs.mpr.sub2$nv10@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv10 <- FindVariableFeatures(hs.mpr.sub2$nv10) %>% ScaleData(rownames(hs.mpr.sub2$nv10)) %>% RunPCA(verbose = FALSE)
set.seed(290)
hs.mpr.sub2$nv10 <- RunHarmony(hs.mpr.sub2$nv10, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv10, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.sub2$nv10 <- RunUMAP(hs.mpr.sub2$nv10, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$nv10, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv10, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(pd.gene))) +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c20c25.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv10, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv10 <- FindAllMarkers(hs.mpr.sub2$nv10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
FeaturePlot(hs.mpr.sub2$nv10, features = c("SIX1","ISL1","TLX3","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
FeaturePlot(hs.mpr.sub1$ecto1, features = c("SIX1","ISL1","TLX3","HOXB6"), slot = "scale.data", reduction = "umap", ncol = 2, order = F, pt.size = 0.5)
VlnPlot(hs.mpr.sub2$nv10, features = c("SIX1","ISL1","TLX3","HOXB6"), stack = T, flip = T)
FeaturePlot(hs.mpr.sub2$nv10, features = c("PAX3","SOX2","ISL2"), slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
hs.mpr.sub2.markers$nv10 %>% dplyr::filter(cluster == 4) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv10, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
hs.mpr.sub2$nv10@meta.data <- hs.mpr.sub2$nv10@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(0,2:5,7:9,11:12,14) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv10, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - cluster 5,9
hs.mpr.sub2$nv11 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("5", "9"))
hs.mpr.sub2$nv11@meta.data <- hs.mpr.sub2$nv11@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv11 <- FindVariableFeatures(hs.mpr.sub2$nv11) %>% ScaleData(rownames(hs.mpr.sub2$nv11)) %>% RunPCA(verbose = FALSE)
set.seed(300)
hs.mpr.sub2$nv11 <- RunHarmony(hs.mpr.sub2$nv11, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv11, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$nv11 <- RunUMAP(hs.mpr.sub2$nv11, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 3.6)
DimPlot(hs.mpr.sub2$nv11, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv11, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  theme_dp
pdf(file.path(res.out, "Featureplot_to_show_nervous_system_subcluster_markers_in_sub_c5c9.pdf"), height = 90, width = 40)
FeaturePlot(hs.mpr.sub2$nv11, features = sort(unique(marker$nv)), slot = "scale.data", reduction = "umap", ncol = 6)
dev.off()
# find markers
hs.mpr.sub2.markers$nv11 <- FindAllMarkers(hs.mpr.sub2$nv11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# visualization
FeaturePlot(hs.mpr.sub2$nv11, features = c("ATOH1","BARHL1","BARHL2","BHLHE22","EN1","EVX1","FOXD3","FOXP2","HOXB6","LBX1","LHX1",
                                           "LHX2","LHX5","LHX9","LMX1B","OTP","PAX2","PHOX2A","PHOX2B","POU4F1","DRGX","TLX1","TLX3"),
            slot = "scale.data", reduction = "umap", ncol = 6, order = T, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv11, features = c("FOXD3","FOXP2"),
            slot = "scale.data", reduction = "umap", ncol = 6, order = T, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv11, features = c("LHX1","LHX5"),
            slot = "scale.data", reduction = "umap", ncol = 6, order = T, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv11, features = c("FOXD3","FOXP2"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv11, features = c("LHX1","LHX5","EN1","OTP","HOXB6","LBX1","PAX2","POU4F1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F, pt.size = 0.5)
hs.mpr.sub2.markers$nv11 %>% dplyr::filter(cluster == 2) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub2$nv11, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4, order = F)
hs.mpr.sub2$nv11@meta.data <- hs.mpr.sub2$nv11@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c(0,2:5,7:9,11:12,14) ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub2$nv11, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# - cluster 21
hs.mpr.sub2$nv12 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("21"))
hs.mpr.sub2$nv12@meta.data <- hs.mpr.sub2$nv12@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv12 <- FindVariableFeatures(hs.mpr.sub2$nv12) %>% ScaleData(rownames(hs.mpr.sub2$nv12)) %>% RunPCA(verbose = FALSE)
set.seed(400)
hs.mpr.sub2$nv12 <- RunHarmony(hs.mpr.sub2$nv12, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv12, reduction = "harmony", ndims = 30)
dim.n <- 6
hs.mpr.sub2$nv12 <- RunUMAP(hs.mpr.sub2$nv12, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub2$nv12, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv12, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# find markers
hs.mpr.sub2.markers$nv12 <- FindAllMarkers(hs.mpr.sub2$nv12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# visualization
FeaturePlot(hs.mpr.sub2$nv12, features = c("LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3"), slot = "scale.data", reduction = "umap", ncol = 6)
DotPlot(hs.mpr.sub2$nv12, features = c("LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$nv12, features = c("ATOH1","BARHL1","BARHL2","BHLHE22","EN1","EVX1","FOXD3","FOXP2","HOXB6","LBX1","LHX1",
                                           "LHX2","LHX5","LHX9","LMX1B","OTP","PAX2","PHOX2A","PHOX2B","POU4F1","DRGX","TLX1","TLX3"),
            slot = "scale.data", reduction = "umap", ncol = 6, order = T, pt.size = 0.5)
FeaturePlot(hs.mpr.sub2$nv12, features = c("WNT2B"),
            slot = "scale.data", reduction = "umap", ncol = 1, order = T, pt.size = 0.5)
pd.gene <- subset(hs.mpr.sub2.markers$nv12, cluster == "2") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$nv12, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 6)
DotPlot(hs.mpr.sub2$nv12, features = sort(unique(marker$splanchnic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# - cluster 29
hs.mpr.sub2$nv13 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("29"))
hs.mpr.sub2$nv13@meta.data <- hs.mpr.sub2$nv13@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv13 <- FindVariableFeatures(hs.mpr.sub2$nv13) %>% ScaleData(rownames(hs.mpr.sub2$nv13)) %>% RunPCA(verbose = FALSE)
set.seed(500)
hs.mpr.sub2$nv13 <- RunHarmony(hs.mpr.sub2$nv13, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv13, reduction = "harmony", ndims = 30)
dim.n <- 4
hs.mpr.sub2$nv13 <- RunUMAP(hs.mpr.sub2$nv13, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub2$nv13, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv13, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# find markers
hs.mpr.sub2.markers$nv13 <- FindAllMarkers(hs.mpr.sub2$nv13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# visualization
FeaturePlot(hs.mpr.sub2$nv13, features = c("NKX6-1","FOXA2","SHH","SPON1","HOXB6"), slot = "data", reduction = "umap", ncol = 6)
FeaturePlot(hs.mpr.sub2$nv13, features = c("MITF","MEIS2","HOXB6"), slot = "data", reduction = "umap", ncol = 6)
FeaturePlot(hs.mpr.sub2$nv13, features = c("DMRTA1","DMRTA2","EMX2","PAX6","DMRT3"), slot = "data", reduction = "umap", ncol = 6)
FeaturePlot(hs.mpr.sub2$nv7, features = c("NKX6-1","FOXA2","SHH","SPON1","HOXB6"), slot = "data", reduction = "umap", ncol = 6)
FeaturePlot(hs.mpr.sub2$nv7, features = c("MITF","MEIS2","HOXB6"), slot = "data", reduction = "umap", ncol = 6)
VlnPlot(hs.mpr.sub2$nv7, features = c("MITF","MEIS2","HOXB6"), stack = T, flip = T, slot = "data")
VlnPlot(hs.mpr.sub2$nv13, features = c("MITF","MEIS2","HOXB6"), stack = T, flip = T, slot = "data")
# - cluster 31
hs.mpr.sub2$nv14 <- subset(hs.mpr.sub1$ecto1, seurat_clusters%in%c("31"))
hs.mpr.sub2$nv14@meta.data <- hs.mpr.sub2$nv14@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$nv14 <- FindVariableFeatures(hs.mpr.sub2$nv14) %>% ScaleData(rownames(hs.mpr.sub2$nv14)) %>% RunPCA(verbose = FALSE)
set.seed(600)
hs.mpr.sub2$nv14 <- RunHarmony(hs.mpr.sub2$nv14, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$nv14, reduction = "harmony", ndims = 30)
dim.n <- 4
hs.mpr.sub2$nv14 <- RunUMAP(hs.mpr.sub2$nv14, reduction = "harmony", dims = 1:dim.n, n.neighbors = 26) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub2$nv14, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
DotPlot(hs.mpr.sub2$nv14, features = sort(unique(marker$nv)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# find markers
hs.mpr.sub2.markers$nv14 <- FindAllMarkers(hs.mpr.sub2$nv14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# visualization
FeaturePlot(hs.mpr.sub2$nv14, features = c("PRDM8", "NKX6-1", "ASCL1", "VSX1", "PTF1A"), slot = "data", reduction = "umap", ncol = 6)
FeaturePlot(hs.mpr.sub2$nv14, features = c("SLC8A1","GATA6","NKX2-5","HAND2","HEY2"), slot = "data", reduction = "umap", ncol = 6)
pd.gene <- subset(hs.mpr.sub2.markers$nv14, cluster == "1") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$nv14, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 6)
DotPlot(hs.mpr.sub2$nv14, features = sort(unique(marker$splanchnic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp


### >>> 3. Ectoderm2 (schwann cell): 1st round clustering ----
hs.mpr.sub1$ecto2 <- subset(hs.mpr, seurat_clusters%in%c("25"))
hs.mpr.sub1$ecto2@meta.data <- hs.mpr.sub1$ecto2@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$ecto2 <- FindVariableFeatures(hs.mpr.sub1$ecto2) %>% ScaleData(rownames(hs.mpr.sub1$ecto2)) %>% RunPCA(verbose = FALSE)
set.seed(110)
hs.mpr.sub1$ecto2 <- RunHarmony(hs.mpr.sub1$ecto2, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$ecto2, reduction = "harmony", ndims = 30)
dim.n <- 20
hs.mpr.sub1$ecto2 <- RunUMAP(hs.mpr.sub1$ecto2, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_schwann_cell_corrected_clustering.pdf"), height = 5, width = 17)
DimPlot(hs.mpr.sub1$ecto2, reduction = "umap", group.by = c("Sample", "ident", "raw_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# find markers
hs.mpr.sub1.markers$ecto2 <- FindAllMarkers(hs.mpr.sub1$ecto2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# plot markers
FeaturePlot(hs.mpr.sub1$ecto2, features = c("FOXD3", "MPZ", "PLP1", "S100B"), slot = "data", reduction = "umap", ncol = 4)
pdf(file.path(res.out, "Dotplot_to_show_schwann_cell_subcluster_markers.pdf"), height = 4, width = 9)
DotPlot(hs.mpr.sub1$ecto2, features = sort(unique(marker$sw)), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 1.5,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(unique(marker$sw))) +
  labs(x = "Hs Schwann Cell Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pdf(file.path(res.out, "Featureplot_to_show_schwann_cell_subcluster_markers.pdf"), height = 26, width = 35)
FeaturePlot(hs.mpr.sub1$ecto2, features = sort(unique(marker$sw)), slot = "data", reduction = "umap", ncol = 5)
dev.off()
# c0,4,6
FeaturePlot(hs.mpr.sub1$ecto2, features = c("SOX10","MPZ","PLP1","FOXD3","PMP22"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
hs.mpr.sub1$ecto2@meta.data <- hs.mpr.sub1$ecto2@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("0","4","6") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto2, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c1
hs.mpr.sub1.markers$ecto2 %>% dplyr::filter(cluster == 1) %>% dplyr::top_n(n = 10, wt = avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC)) -> pd
FeaturePlot(hs.mpr.sub1$ecto2, features = pd$gene, slot = "scale.data", reduction = "umap", ncol = 4)
# c2
FeaturePlot(hs.mpr.sub1$ecto2, features = c("HAND2","ASCL1","PHOX2A","PHOX2B","DCX","TUBB3","ELAVL4"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
hs.mpr.sub1$ecto2@meta.data <- hs.mpr.sub1$ecto2@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("2") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto2, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)
# c7
FeaturePlot(hs.mpr.sub1$ecto2, features = c("MITF","TFAP2A","TFAP2B","PAX3","PMEL"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
hs.mpr.sub1$ecto2@meta.data <- hs.mpr.sub1$ecto2@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("7") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$ecto2, reduction = "umap", group.by = c("Clustering"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5)


### >>> 4. Mesoderm1 (LPM+HM+Somite+Limb+Craniofacial): 1st round clustering ----
hs.mpr.sub1$meso1 <- subset(hs.mpr, seurat_clusters%in%c("0","1","2","3","4","6","8","9","12","16","17","18","19","21","39"))
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$meso1 <- FindVariableFeatures(hs.mpr.sub1$meso1) %>%
  ScaleData(rownames(hs.mpr.sub1$meso1)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$meso1 <- RunHarmony(hs.mpr.sub1$meso1, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$meso1, reduction = "harmony", ndims = 30)
dim.n <- 20
hs.mpr.sub1$meso1 <- RunUMAP(hs.mpr.sub1$meso1, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_meso1_corrected_clustering.pdf"), height = 5, width = 17)
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
# find markers
hs.mpr.sub1.markers$meso1 <- FindAllMarkers(hs.mpr.sub1$meso1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$meso1 %>% filter(cluster == "20") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$meso1, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub1$meso1, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1")
DotPlot(hs.mpr.sub1$meso1, features = pd.gene, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ALX1","ALX3","PITX2","TCF21","LHX2","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","HAND2")
DotPlot(hs.mpr.sub1$meso1, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
# head mesoderm
FeaturePlot(hs.mpr.sub1$meso1, features = c("SNAI2","TWIST1","FOXC1","FOXC2"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.sub1$meso1, features = c("PITX2","MYF5","MSX2","PRRX1","CYP26C1","CSPG4","KCNJ8","ABCC9"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("1","11","23","30") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# somite
FeaturePlot(hs.mpr.sub1$meso1, features = c("SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub1$meso1, features = c("PAX1","PAX9","NKX3-2","SOX9","SCX","ETV4", "PAX3","DMRT2","MYF5","MYF6","MYOD1","LBX1","SIX1",
                                            "FOXC1","FOXC2","SHISA2","WNT5A","CYP26A1","SOX2","MESP1"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("2","4","5","17","19","29","33","34") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# Lateral plate mesoderm
FeaturePlot(hs.mpr.sub1$meso1, features = c("GATA5","GATA6"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
FeaturePlot(hs.mpr.sub1$meso1, features = c("IRX3","GATA6","TBX1","WNT2","NKX6-1","LHX2","HLX","FOXF1","BARX1","ISL1","NKX2-5","TBX5",
                                            "NPPA","MYH6","MYL2","MYH7","SLC8A1","HAND2","HEY2","SHOX2","TBX18","WT1","TNNI3","TNNT1",
                                            "TCF21","SNAI2","SOX9","THY1","NPPA","GATA4","TGFBI","TWIST1"),
            slot = "scale.data", reduction = "umap", ncol=6, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("10","12","22","25","26","27","35") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# Intermediate mesoderm
FeaturePlot(hs.mpr.sub1$meso1, features = c("PAX8","WT1"), slot = "scale.data", reduction = "umap", ncol=2, order = T)
FeaturePlot(hs.mpr.sub1$meso1, features = c("PAX8","EMX2","LHX1","OSR1"), slot = "scale.data", reduction = "umap", ncol=4, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("32") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# limb
FeaturePlot(hs.mpr.sub1$meso1, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("9","18") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# craniofacial
FeaturePlot(hs.mpr.sub1$meso1, features = c("ALX1","ALX3","DLX1","DLX2","SNAI2","TWIST1"), slot = "scale.data", reduction = "umap", ncol=3, order = T)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("6","7","14","15","16","20","21","21","24","31") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
FeaturePlot(hs.mpr.sub1$meso1, features = c("ALX1","ALX3","PITX2","TCF21","LHX2","MSC","HOXA2","HOXA3","HOXB2",
                                            "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","HAND2"), slot = "scale.data", reduction = "umap", ncol=4, order = T)
# neuron: lateral motor columns
FeaturePlot(hs.mpr.sub1$meso1, features = c("RBFOX1","CTNNA2","NRXN1","NRXN3"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("0") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
# undefined
hs.mpr.sub1$meso1@meta.data <- hs.mpr.sub1$meso1@meta.data %>%
  mutate(Clustering=case_when(seurat_clusters%in%c("3","8","13") ~ as.character(seurat_clusters), TRUE ~ "Non"))
DimPlot(hs.mpr.sub1$meso1, reduction = "umap", group.by = c("Clustering"), ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)


### >>> 5. Mesoderm1: 2nd round clustering ----
# - head mesoderm
hs.mpr.sub2$hm <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("1","11","23","30"))
hs.mpr.sub2$hm@meta.data <- hs.mpr.sub2$hm@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$hm <- FindVariableFeatures(hs.mpr.sub2$hm) %>% ScaleData(rownames(hs.mpr.sub2$hm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$hm <- RunHarmony(hs.mpr.sub2$hm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$hm, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$hm <- RunUMAP(hs.mpr.sub2$hm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub2$hm, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"), ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$hm <- FindAllMarkers(hs.mpr.sub2$hm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$hm %>% filter(cluster == "0") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$hm, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$hm, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub2$hm, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("PITX2","MYF5","MSX2","PRRX1","CYP26C1","CSPG4","KCNJ8","ABCC9")
DotPlot(hs.mpr.sub2$hm, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
FeaturePlot(hs.mpr.sub2$hm, features = c("SNAI2","TWIST1","FOXC1","FOXC2"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$hm, features = c("PITX2","MYF5","MSX2","PRRX1","CYP26C1","CSPG4","KCNJ8","ABCC9"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$hm, features = c("ALX1","ALX3","PAX1","PAX9","PITX1","SOX9"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$hm, features = c("DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","COL6A6"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# - somite
hs.mpr.sub2$somite <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("2","4","5","17","19","29","33","34"))
hs.mpr.sub2$somite@meta.data <- hs.mpr.sub2$somite@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$somite <- FindVariableFeatures(hs.mpr.sub2$somite) %>% ScaleData(rownames(hs.mpr.sub2$somite)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$somite <- RunHarmony(hs.mpr.sub2$somite, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$somite, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$somite <- RunUMAP(hs.mpr.sub2$somite, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub2$somite, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$somite <- FindAllMarkers(hs.mpr.sub2$somite, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$somite %>% filter(cluster == "2") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$somite, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$somite, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub2$somite, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("PAX1","PAX9","NKX3-2","SOX9","SCX","ETV4","PAX3","DMRT2","MYF5","MYF6","MYOD1",
             "LBX1","SIX1","FOXC1","FOXC2","SHISA2","WNT5A","CYP26A1","SOX2","MESP1")
DotPlot(hs.mpr.sub2$somite, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
FeaturePlot(hs.mpr.sub2$somite, features = c("PAX1","PAX9","NKX3-2","SOX9"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$somite, features = c("FOXC1","FOXC2","SHISA2","WNT5A","CYP26A1","PAX3","SOX2","ETV4",
                                             "SCX","DMRT2","LBX1","SIX1","MYF5","MYF6","MYOD1"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.sub2$somite, features = c("FOXF1","GATA6","GATA5","TBX1","MSX2","PRRX1","TGFBI","TWIST1"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# - Lateral plate mesoderm
hs.mpr.sub2$lpm <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("10","12","22","25","26","27","35"))
hs.mpr.sub2$lpm@meta.data <- hs.mpr.sub2$lpm@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$lpm <- FindVariableFeatures(hs.mpr.sub2$lpm) %>% ScaleData(rownames(hs.mpr.sub2$lpm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$lpm <- RunHarmony(hs.mpr.sub2$lpm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$lpm, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.sub2$lpm <- RunUMAP(hs.mpr.sub2$lpm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$lpm, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$lpm <- FindAllMarkers(hs.mpr.sub2$lpm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$lpm %>% filter(cluster == "16") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$lpm, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$lpm, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub2$lpm, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
DotPlot(hs.mpr.sub2$lpm, features = unique(marker$endoderm), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$endoderm)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
# atrioventricular canal
FeaturePlot(hs.mpr.sub2$lpm, features = c("SLC8A1","GATA6","HAND2","HEY2"), slot = "scale.data", reduction = "umap", ncol = 4, order = F)
# epicardium
FeaturePlot(hs.mpr.sub2$lpm, features = c("TBX18","WT1","TNNI3","TNNT1"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# hepatic stellate cell
FeaturePlot(hs.mpr.sub2$lpm, features = c("LHX2","HLX","FOXF1","AFP"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# epicardial derived cell
FeaturePlot(hs.mpr.sub2$lpm, features = c("SNAI2","TCF21","SOX9","THY1"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# somite
FeaturePlot(hs.mpr.sub2$lpm, features = c("PAX1","PAX9","SOX9","FOXC1","FOXC2","SHISA2"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# limb
FeaturePlot(hs.mpr.sub2$lpm, features = c("ALX1","ALX3","TWIST1","TBX5","PITX1"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# - craniofacial
hs.mpr.sub2$craniofacial <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("6","7","14","15","16","20","21","24","31"))
hs.mpr.sub2$craniofacial@meta.data <- hs.mpr.sub2$craniofacial@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$craniofacial <- FindVariableFeatures(hs.mpr.sub2$craniofacial) %>% ScaleData(rownames(hs.mpr.sub2$craniofacial)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$craniofacial <- RunHarmony(hs.mpr.sub2$craniofacial, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$craniofacial, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$craniofacial <- RunUMAP(hs.mpr.sub2$craniofacial, reduction = "harmony", dims = 1:dim.n, n.neighbors = 30, min.dist = 0.1) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$craniofacial, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$craniofacial <- FindAllMarkers(hs.mpr.sub2$craniofacial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$craniofacial %>% filter(cluster == "5") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$craniofacial, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$craniofacial, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub2$craniofacial, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ALX1","ALX3","PITX2","TCF21","LHX2","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","HAND2")
DotPlot(hs.mpr.sub2$craniofacial, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
DotPlot(hs.mpr.sub2$craniofacial, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$splanchnic.lpm)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
DotPlot(hs.mpr.sub2$craniofacial, features = unique(marker$somite), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$somite)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
# craniofacial
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("ALX1","ALX3","DLX1","DLX2","SNAI2","TWIST1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# frontonasal mesenchyme
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("ALX1","ALX3"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# PA1/2 hand2
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("HAND2"), slot = "scale.data", reduction = "umap", ncol = 1, order = T)
# PA1/2 proximal
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("DLX1","DLX2"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# PA1/2 middle
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("DLX1","DLX2","DLX5","DLX6"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# PA3/4 middle
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("DLX1","DLX2","DLX5","DLX6","HOXA2","HOXA3"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# Head mesoderm: undefined
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("CYP26A1","CYP26C1"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# atrioventricular canal.2
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("TBX5","GATA6","ERBB4","TSHZ2"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F)
# atrioventricular canal.2
FeaturePlot(hs.mpr.sub2$craniofacial, features = c("PAX1","PAX3"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = F)
# - mixed
hs.mpr.sub2$mixed <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("0","3","8","9","13","18","32"))
hs.mpr.sub2$mixed@meta.data <- hs.mpr.sub2$mixed@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$mixed <- FindVariableFeatures(hs.mpr.sub2$mixed) %>% ScaleData(rownames(hs.mpr.sub2$mixed)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$mixed <- RunHarmony(hs.mpr.sub2$mixed, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$mixed, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$mixed <- RunUMAP(hs.mpr.sub2$mixed, reduction = "harmony", dims = 1:dim.n, n.neighbors = 30, min.dist = 0.1) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$mixed, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$mixed <- FindAllMarkers(hs.mpr.sub2$mixed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$mixed %>% filter(cluster == "4") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$mixed, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$mixed, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub2$mixed, features = sort(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = sort(pd.gene)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 11, angle = 45),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 11, angle = 0))
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$craniofacial), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$craniofacial)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$splanchnic.lpm)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
# limb
FeaturePlot(hs.mpr.sub2$mixed, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
# craniofacial
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$craniofacial), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$mixed, features = c("SNAI2","TCF21","SOX9","THY1"), slot = "scale.data", reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.sub2$mixed, features = c("ALX1","ALX3","HOXA2","HOXA3","HOXB2","HAND2","DLX1","DLX2","DLX3","DLX4","DLX5","DLX6"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.sub2$mixed, features = c("ALX1","ALX3","DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","HAND2",
                                            "HOXA2","HOXA3","HOXB2","LHX2","MSC","PITX2","TCF21"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# intermediate mesoderm
FeaturePlot(hs.mpr.sub2$mixed, features = c("WT1","PAX8","HOXB7"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$im), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.sub2$mixed, features = unique(c("KLK6","WT1","ALDH1A2","APOBEC3C","SERINC5","PAX8","BST2","PODXL","PLAGL1","BCAM",
                                               "MGAT4C","SHISA3","PLAGL1","GMDS","TSHZ2","HOXC6","IL11RA","BST2","CXCL12","HOXA7",
                                               "NRK","MGAT4C","MIR503HG","PLAGL1","KLK6","PHACTR2","WDR86","FAM89A","ALDH1A2","APOBEC3C")),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# somatic LPM (anterior somatic LPM)
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$somatic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# somite
DotPlot(hs.mpr.sub2$mixed, features = unique(c(marker$somite)),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$mixed, features = c("SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# splanchnic LPM
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$mixed, features = c("HLX","LHX2","FOXF1","AFP"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# limb: forelimb + hindlimb
FeaturePlot(hs.mpr.sub2$mixed, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
# endothelium: liver sinusoidal endothelial cell
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$mixed, features = c("STAB2","CD36","XBP1","CD34"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
# nervous system
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$nv), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.sub2$mixed, features = c("SOX2","PAX6","EMX2","OTX2","HES5","HOXB6","SIX3","TUBB3","ELAVL3","LHX5","POU4F1","ISL2",
                                        "HESX1","FGF8","TFAP2A","TFAP2B","ISL1","SOX10","MPZ","ALX3"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# schwann cells
DotPlot(hs.mpr.sub2$mixed, features = unique(marker$sw), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub2$mixed, features = c("MITF","TFAP2A","TFAP2B","PAX3","PMEL"), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.sub2$mixed, features = c("ASCL1","DCX","ELAVL4","HAND2","PHOX2A","PHOX2B","TUBB3"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# - limb
hs.mpr.sub2$limb <- subset(hs.mpr.sub1$meso1, seurat_clusters %in% c("28"))
hs.mpr.sub2$limb@meta.data <- hs.mpr.sub2$limb@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$limb <- FindVariableFeatures(hs.mpr.sub2$limb) %>% ScaleData(rownames(hs.mpr.sub2$limb)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$limb <- RunHarmony(hs.mpr.sub2$limb, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$limb, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$limb <- RunUMAP(hs.mpr.sub2$limb, reduction = "harmony", dims = 1:dim.n, n.neighbors = 30, min.dist = 0.1) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$limb, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$limb <- FindAllMarkers(hs.mpr.sub2$limb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$limb %>% filter(cluster == "8") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$limb, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$limb, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
# limb: forelimb + hindlimb
FeaturePlot(hs.mpr.sub2$limb, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
# Neuron
FeaturePlot(hs.mpr.sub2$limb, features = c("TUBB3","OSR1"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
# - Erythroid
hs.mpr.sub2$erythroid <- subset(hs.mpr, seurat_clusters %in% c("33"))
hs.mpr.sub2$erythroid@meta.data <- hs.mpr.sub2$erythroid@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub2$erythroid <- FindVariableFeatures(hs.mpr.sub2$erythroid) %>% ScaleData(rownames(hs.mpr.sub2$erythroid)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub2$erythroid <- RunHarmony(hs.mpr.sub2$erythroid, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub2$erythroid, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub2$erythroid <- RunUMAP(hs.mpr.sub2$erythroid, reduction = "harmony", dims = 1:dim.n, n.neighbors = 30, min.dist = 0.1) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub2$erythroid, reduction = "umap", group.by = c("Sample", "raw_clusters", "seurat_clusters"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub2.markers$erythroid <- FindAllMarkers(hs.mpr.sub2$erythroid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub2.markers$erythroid %>% filter(cluster == "0") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub2$erythroid, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub2$erythroid, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
DotPlot(hs.mpr.sub2$erythroid, features = marker$blood, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$blood) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
DotPlot(hs.mpr.sub2$erythroid, features = c("HBA1","HBA2","HBG1","HBG2"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("HBA1","HBA2","HBG1","HBG2")) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
# erythroid
FeaturePlot(hs.mpr.sub2$erythroid, features = c("HBA1","HBA2","HBG1","HBG2"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.sub2$erythroid, features = c("CD34","SOX2"), slot = "scale.data", reduction = "umap", ncol = 1, order = T)


### >>> 6. Mesoderm2 (endothelium): 1st round clustering ----
hs.mpr.sub1$meso2 <- subset(hs.mpr, seurat_clusters%in%c("10","20","23","27"))
hs.mpr.sub1$meso2@meta.data <- hs.mpr.sub1$meso2@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$meso2 <- FindVariableFeatures(hs.mpr.sub1$meso2) %>%
  ScaleData(rownames(hs.mpr.sub1$meso2)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$meso2 <- RunHarmony(hs.mpr.sub1$meso2, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$meso2, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub1$meso2 <- RunUMAP(hs.mpr.sub1$meso2, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
DimPlot(hs.mpr.sub1$meso2, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub1.markers$meso2 <- FindAllMarkers(hs.mpr.sub1$meso2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$meso2 %>% filter(cluster == "0") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$meso2, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub1$meso2, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
pd.gene <- c("CDH5","CD34","PECAM1","NOTCH1","NOTCH4","DLL4","HEY1","HEY2",
             "TMEM100","GJA5","GJA4","ELN","XBP1","STAB2","CD36")
DotPlot(hs.mpr.sub1$meso2, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub1$meso2, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
# vascular endothelium
FeaturePlot(hs.mpr.sub1$meso2, features = c("CDH5","CD34","PECAM1"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# arterial endothelium
FeaturePlot(hs.mpr.sub1$meso2, features = c("NOTCH1","NOTCH4","DLL4","HEY1","HEY2","TMEM100","GJA5","GJA4","ELN"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F)
# liver sinusoidal endothelial cell
FeaturePlot(hs.mpr.sub1$meso2, features = c("XBP1","STAB2","CD36"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# endocardial derived cell
FeaturePlot(hs.mpr.sub1$meso2, features = c("GATA4","TGFBI","HEY2","TWIST1","HAND2"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
DotPlot(hs.mpr.sub1$meso2, features = unique(marker$nv), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$nv)) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
# fibroblast?
FeaturePlot(hs.mpr.sub1$meso2, features = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.sub1$meso2, features = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1","ALB"),
        cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1","ALB")) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
DotPlot(hs.mpr.sub1$meso2, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$endothelium)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp


### >>> 7. Mesoderm3 (blood): 1st round clustering ----
hs.mpr.sub1$meso3 <- subset(hs.mpr, seurat_clusters%in%c("13","26","36"))
hs.mpr.sub1$meso3@meta.data <- hs.mpr.sub1$meso3@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$meso3 <- FindVariableFeatures(hs.mpr.sub1$meso3) %>%
  ScaleData(rownames(hs.mpr.sub1$meso3)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$meso3 <- RunHarmony(hs.mpr.sub1$meso3, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$meso3, reduction = "harmony", ndims = 30)
dim.n <- 6
hs.mpr.sub1$meso3 <- RunUMAP(hs.mpr.sub1$meso3, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub1$meso3, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub1.markers$meso3 <- FindAllMarkers(hs.mpr.sub1$meso3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$meso3 %>% filter(cluster == "8") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$meso3, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# plot marker genes
DotPlot(hs.mpr.sub1$meso3, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("CLEC10A","CD1C","PLAC8","HDC","LMO4","CLC","CD34","LTB","IL7R","CD52","JCHAIN",
             "CD68","CD14","CD163","CSF1R","ITGAM","ID3","NR1H3","CD9","PECAM1","ITGA2B",
             "GP1BA","PPBP","S100A8","S100A9","LYZ","PRTN3","AZU1","MPO","HBA1","GYPB")
DotPlot(hs.mpr.sub1$meso3, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("ABCC9","BARX1","CSPG4","CYP26A1","CYP26C1","DMRT2","EMX2","ETV4","FOXC1","FOXC2","FOXF1","GATA4","GATA6",
             "HAND2","HEY2","HLX","IRX3","ISL1","KCNJ8","LBX1","LHX1","LHX2","MESP1","MSX2","MYF5","MYF6","MYH6","MYH7",
             "MYL2","MYOD1","NKX2-5","NKX3-2","NKX6-1","NPPA","PAX1","PAX3","PAX8","PAX9","PITX1","PITX2","PRRX1","SCX",
             "SHISA2","SHOX2","SIX1","SLC8A1","SNAI2","SOX2","SOX9","TBX1","TBX18","TBX5","TCF21","TGFBI","THY1","TNNI3",
             "TNNT1","TWIST1","WNT2","WNT5A","WT1","ALX1","ALX3","MSC","HOXA2","HOXA3","HOXB2",
             "DLX1","DLX2","DLX3","DLX4","DLX5","DLX6")
DotPlot(hs.mpr.sub1$meso3, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
# dendritic cell
FeaturePlot(hs.mpr.sub1$meso3, features = c("CLEC10A","CD1C","PLAC8"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# eosino/basophil/mast cell progenitor
FeaturePlot(hs.mpr.sub1$meso3, features = c("HDC","LMO4","CLC"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# lymphocyte
FeaturePlot(hs.mpr.sub1$meso3, features = c("CD34","LTB","IL7R","CD52","JCHAIN"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# megakaryocyte
FeaturePlot(hs.mpr.sub1$meso3, features = c("CD9","PECAM1","ITGA2B","GP1BA","PPBP"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# macrophage
FeaturePlot(hs.mpr.sub1$meso3, features = c("CD68","CD14","CD163","CSF1R","ITGAM"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# neutrophil
FeaturePlot(hs.mpr.sub1$meso3, features = c("S100A8","S100A9","LYZ"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# mixed: 2,8,10,11
FeaturePlot(hs.mpr.sub1$meso3, features = c("ID3","NR1H3","PITX1"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# fibroblast?
FeaturePlot(hs.mpr.sub1$meso3, features = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.sub1$meso3, features = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1","ALB"),
        cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("TNNT2","ACTN2","COL1A1","COL1A2","COL3A1","POSTN","DCN","PDGFRA","SOX9","WT1","ALB")) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
# atrioventricular canal
DotPlot(hs.mpr.sub1$meso3, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$splanchnic.lpm)) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
FeaturePlot(hs.mpr.sub1$meso3, features = c("SLC8A1","GATA6","NKX2-5","HAND2","HEY2"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)


### >>> 8. Mesoderm4 (Heart): 1st round clustering ----
hs.mpr.sub1$meso4 <- subset(hs.mpr, seurat_clusters%in%c("24"))
hs.mpr.sub1$meso4@meta.data <- hs.mpr.sub1$meso4@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$meso4 <- FindVariableFeatures(hs.mpr.sub1$meso4) %>%
  ScaleData(rownames(hs.mpr.sub1$meso4)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$meso4 <- RunHarmony(hs.mpr.sub1$meso4, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$meso4, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.sub1$meso4 <- RunUMAP(hs.mpr.sub1$meso4, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub1$meso4, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub1.markers$meso4 <- FindAllMarkers(hs.mpr.sub1$meso4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$meso4 %>% filter(cluster == "8") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$meso4, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.sub1$meso4, features = c("PITX2","MYF5","MYOD1","MYF6"), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# plot marker genes
# super cluster
DotPlot(hs.mpr.sub1$meso4, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") + theme_dp
# cardiomyocytes marker in mouse embryo
FeaturePlot(hs.mpr.sub1$meso4, features = c("SMARCD3","ACTA2","TAGLN","MYL7","MYL4","TNNT2","NKX2-5"),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.sub1$meso4, features = c("SMARCD3","ACTA2","TAGLN","MYL7","MYL4","TNNT2","NKX2-5"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("SMARCD3","ACTA2","TAGLN","MYL7","MYL4","TNNT2","NKX2-5")) +
  labs(x = "Mm Cardiomyocytes Marker Gene", y = "Cluster") +
  theme_dp
# hepatocyte
FeaturePlot(hs.mpr.sub1$meso4, features = c("XBP1","CEBPA","ALB"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# myotome
DotPlot(hs.mpr.sub1$meso4, features = unique(c("MYH3","MYLPF","TNNC2","IL17B","MYL1","MYOG","KLHL41","MAP3K7CL",
                                               "MYF5","RAPSN","PRL","CDH15","MEOX1","NEB","MYF6","SGCA")),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$meso4, features = unique(c("MYH3","MYLPF","TNNC2","IL17B","MYL1","MYOG","KLHL41","MAP3K7CL")),
            slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# head muscle
DotPlot(hs.mpr.sub1$meso4, features = unique(c("MYF5","ZNF385D","SPX","SHISA2","CDH15","PITX2","FLRT3","MXRA8","PGF","NEB")),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# migrating hypaxial muscle
DotPlot(hs.mpr.sub1$meso4, features = unique(c("MET","VMO1","CCDC140","NAV3","PDGFC","PDGFA","SIX1","LBX1","PAX3","PCDH17",
                                               "ITIH5","PDGFC","PART1","SIX2","VMO1","PLK2","NAV3","LBX1","SIX1","NEB")),
        cols = c("#eaeaea", "#fc0330"), scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# limb AER
FeaturePlot(hs.mpr.sub1$meso4, features = c("TP63"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)


### >>> 9. Endoderm: 1st round clustering ----
hs.mpr.sub1$endo <- subset(hs.mpr, seurat_clusters%in%c("22","31","37"))
hs.mpr.sub1$endo@meta.data <- hs.mpr.sub1$endo@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$endo <- FindVariableFeatures(hs.mpr.sub1$endo) %>%
  ScaleData(rownames(hs.mpr.sub1$endo)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$endo <- RunHarmony(hs.mpr.sub1$endo, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$endo, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.sub1$endo <- RunUMAP(hs.mpr.sub1$endo, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub1$endo, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# find markers
hs.mpr.sub1.markers$endo <- FindAllMarkers(hs.mpr.sub1$endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$endo %>% filter(cluster == "10") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$endo, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
# plot marker genes
DotPlot(hs.mpr.sub1$endo, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker$super.cluster) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("PAX1","HOXA3","GLI3","SOX2","KLF5","NKX2-1","SOX9","ETV5","MYCN","IRX1","IRX2","IRX3",
             "HHEX","ONECUT2","FOXA1","PDX1","PTF1A","NKX6-1","CDX2","GATA4","XBP1","CEBPA","ALB")
DotPlot(hs.mpr.sub1$endo, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = pd.gene) +
  labs(x = "Hs Marker Gene", y = "Cluster") +
  theme_dp
pd.gene <- c("SPINK1","BRI3","NPW","FST","MAP2K2",
             "MT1H","IGF2","AMN","APOM","DLK1","RSPO3","GSTP1","BAAT",
             "AFP","CYP3A7","F2","AHSG","APOA1","LIPC","AGT",
             "ALB","ANGPTL3","CFHR1","LEPR","APOB","GC","VTN","AKR1C1",
             "CYP2E1","SAA1","HP","C9","ADH1B","F9","CES1","TAT","AZGP1","APCS","NNMT")
DotPlot(hs.mpr.sub1$endo, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"), col.min = 0, col.max = 3,
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# thymus
FeaturePlot(hs.mpr.sub1$endo, features = c("PAX1","HOXA3","GLI3","SOX2"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
# hepatocyte
FeaturePlot(hs.mpr.sub1$endo, features = c("XBP1","CEBPA","ALB"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)
# lung distal epithelium
FeaturePlot(hs.mpr.sub1$endo, features = c("NKX2-1","SOX9","ETV5","MYCN","IRX1","IRX2","IRX3"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F)
# fibroblast?
FeaturePlot(hs.mpr.sub1$endo, features = c("COL3A1","COL2A1"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
FeaturePlot(hs.mpr.sub1$endo, features = c("IRX3","GATA6"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.sub1$endo, features = c("WNT2","NKX6-1"), slot = "scale.data", reduction = "umap", ncol = 2, order = F)
DotPlot(hs.mpr.sub1$endo, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$splanchnic.lpm)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# stomach
pd.gene <- c("HHEX","KLF5","ONECUT2","FOXA1","FXYD3","DCDC2","STARD10","BEX5","ELF3","CTSH")
FeaturePlot(hs.mpr.sub1$endo, features = pd.gene, slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.sub1$endo, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# duodenum
pd.gene <- c("GATA4","CDX2","LGALS3","HABP2","STK35","ARSE","GPR160","ANXA4","PDE11A","SAYSD1","TJP2")
FeaturePlot(hs.mpr.sub1$endo, features = pd.gene, slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.sub1$endo, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# pancreas
pd.gene <- c("PDX1","PTF1A","NKX6-1","CPA2","SLC4A4","MARVELD3","TM4SF4","HNF1B","STAT1","DPEP1","CHST9","LYPD6B")
FeaturePlot(hs.mpr.sub1$endo, features = pd.gene, slot = "scale.data", reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.sub1$endo, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# renal epithelium
FeaturePlot(hs.mpr.sub1$endo, features = c("PAX8","WT1","EMX2","LHX1"), slot = "scale.data", reduction = "umap", ncol = 4, order = F)
DotPlot(hs.mpr.sub1$endo, features = c("PAX8","WT1","EMX2","LHX1"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = c("PAX8","WT1","EMX2","LHX1")) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
# liver sinusoidal endothelial cell
DotPlot(hs.mpr.sub1$endo, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = unique(marker$endothelium)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
FeaturePlot(hs.mpr.sub1$endo, features = c("XBP1","STAB2","CD36"), slot = "scale.data", reduction = "umap", ncol = 3, order = F)


### >>> 10. Undefined cluster: 1st round clustering ----
hs.mpr.sub1$ud <- subset(hs.mpr, seurat_clusters%in%c("7","30","38"))
hs.mpr.sub1$ud@meta.data <- hs.mpr.sub1$ud@meta.data %>% mutate(raw_clusters=seurat_clusters)
hs.mpr.sub1$ud <- FindVariableFeatures(hs.mpr.sub1$ud) %>%
  ScaleData(rownames(hs.mpr.sub1$ud)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.sub1$ud <- RunHarmony(hs.mpr.sub1$ud, group.by.vars = "Sample")
ElbowPlot(hs.mpr.sub1$ud, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.sub1$ud <- RunUMAP(hs.mpr.sub1$ud, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 1)
DimPlot(hs.mpr.sub1$ud, reduction = "umap", group.by = c("Sample", "raw_clusters", "ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
# - find markers
hs.mpr.sub1.markers$ud <- FindAllMarkers(hs.mpr.sub1$ud, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.sub1.markers$ud %>% filter(cluster == "1") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.sub1$ud, features = pd.gene$gene, slot = "scale.data", reduction = "umap", ncol = 5, order = T)
# - plot marker genes
# super cluster
DotPlot(hs.mpr.sub1$ud, features = marker$super.cluster, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# head mesoderm: pericyte
FeaturePlot(hs.mpr.sub1$ud, features = c("SNAI2","TWIST1","FOXC1","FOXC2"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
DotPlot(hs.mpr.sub1$ud, features = marker$hm, cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ud, features = c("CSPG4","KCNJ8","ABCC9"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# somite: neuromesodermal progenitors
FeaturePlot(hs.mpr.sub1$ud, features = c("SNAI2","TWIST1","FOXC1","FOXC2","PAX1","PAX9"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
DotPlot(hs.mpr.sub1$ud, features = unique(marker$somite), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ud, features = c("PAX1","PAX9","SOX2","NKX3-2","SOX9"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
VlnPlot(hs.mpr.sub1$ud, features = c("PAX1","PAX9","SOX2","NKX3-2","SOX9"), slot = "data", flip = T, stack = T)
# craniofacial: undefined
FeaturePlot(hs.mpr.sub1$ud, features = c("ALX1","ALX3","DLX1","DLX2"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
DotPlot(hs.mpr.sub1$ud, features = unique(marker$craniofacial), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
# endothelium: liver sinusoidal endothelial cell
FeaturePlot(hs.mpr.sub1$ud, features = c("KRT18","PLVAP"), slot = "scale.data", reduction = "umap", ncol = 4, order = F)
DotPlot(hs.mpr.sub1$ud, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ud, features = c("XBP1","STAB2","CD36"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# epidermis
FeaturePlot(hs.mpr.sub1$ud, features = c("PERP","KRT18","CLDN6"), slot = "scale.data", reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.sub1$ud, features = unique(marker$epidermis), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ud, features = c(marker$epidermis,"TP63"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
# placenta marker in human fetal
DotPlot(hs.mpr.sub1$ud, features = unique(ct.mk$placenta$Marker), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.sub1$ud, features = c("HLA-G", "DIO2", "HTRA4", "MFAP5", "HPGD"), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.sub1$ud, features = c("SLC22A11", "PAGE4", "PARP1", "PERP", "PSG4", "PSG9", "PSG6"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
# LPM
FeaturePlot(hs.mpr.sub1$ud, features = c("SNAI2","TWIST1","FOXF1","GATA6","GATA5"), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.sub1$ud, features = unique(marker$splanchnic.lpm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp


### >>> 11. Rename cells ----
# Endoderm
new.cluster.ids <- c("Endothelium:liver sinusoidal endothelial cell","Endothelium:liver sinusoidal endothelial cell",
                     "Endothelium:liver sinusoidal endothelial cell","Gut:duodenum","Liver:hepatocyte",
                     "Intermediate mesoderm:renal epithelium","Endothelium:liver sinusoidal endothelial cell",
                     "Intermediate mesoderm:renal epithelium","Craniofacial:undefined","Liver:hepatocyte","Undefined:undefined")
names(new.cluster.ids) <- levels(hs.mpr.sub1$endo)
hs.mpr.sub1$endo <- RenameIdents(hs.mpr.sub1$endo, new.cluster.ids)
hs.mpr.sub1$endo$CellType <- Idents(hs.mpr.sub1$endo)
# Mesoderm2-Endothelium
new.cluster.ids <- c("Fibroblast:fibroblast","Endothelium:liver sinusoidal endothelium","Endothelium:vascular endothelium",
                     "Endothelium:vascular endothelium","Endothelium:vascular endothelium","Fibroblast:fibroblast",
                     "Endothelium:vascular endothelium","Endothelium:arterial endothelium","Endothelium:vascular endothelium",
                     "Endothelium:vascular endothelium","Endothelium:arterial endothelium","Endothelium:liver sinusoidal endothelium",
                     "Endothelium:liver sinusoidal endothelium","Fibroblast:fibroblast","Fibroblast:fibroblast",
                     "Lateral plate mesoderm:endocardial derived cell","Fibroblast:fibroblast","Endothelium:liver sinusoidal endothelium",
                     "Endothelium:vascular endothelium","Endothelium:liver sinusoidal endothelium","Endothelium:arterial endothelium",
                     "Endothelium:arterial endothelium","Fibroblast:fibroblast")
names(new.cluster.ids) <- levels(hs.mpr.sub1$meso2)
hs.mpr.sub1$meso2 <- RenameIdents(hs.mpr.sub1$meso2, new.cluster.ids)
hs.mpr.sub1$meso2$CellType <- Idents(hs.mpr.sub1$meso2)
# Mesoderm3-Blood
new.cluster.ids <- c("Blood:macrophage","Blood:macrophage","Fibroblast:fibroblast","Blood:macrophage","Blood:dendritic cell",
                     "Blood:macrophage","Blood:granulocytes cell progenitor","Blood:lymphocyte","Fibroblast:fibroblast","Blood:macrophage",
                     "Lateral plate mesoderm:atrioventricular canal","Neural progenitor:undefined","Blood:megakaryocyte","Blood:neutrophil")
names(new.cluster.ids) <- levels(hs.mpr.sub1$meso3)
hs.mpr.sub1$meso3 <- RenameIdents(hs.mpr.sub1$meso3, new.cluster.ids)
hs.mpr.sub1$meso3$CellType <- Idents(hs.mpr.sub1$meso3)
# Mesoderm4-Muscle
new.cluster.ids <- c("Somite:migrating hypaxial muscle","Head mesoderm:head muscle","Head mesoderm:head muscle","Liver:hepatocyte",
                     "Head mesoderm:head muscle","Somite:myotome","Heart:cardiomyocytes","Somite:myotome",
                     "Peripheral surface ectoderm:limb AER","Head mesoderm:head muscle")
names(new.cluster.ids) <- levels(hs.mpr.sub1$meso4)
hs.mpr.sub1$meso4 <- RenameIdents(hs.mpr.sub1$meso4, new.cluster.ids)
hs.mpr.sub1$meso4$CellType <- Idents(hs.mpr.sub1$meso4)
# Mesoderm-Head mesoderm
new.cluster.ids <- c("Craniofacial:frontonasal mesenchyme","Head mesoderm:undefined","Head mesoderm:cranium","Head mesoderm:undefined",
                     "Head mesoderm:pericyte","Somite:sclerotome.early","Head mesoderm:cranium","Head mesoderm:cranium",
                     "Head mesoderm:undefined","Head mesoderm:pericyte")
names(new.cluster.ids) <- levels(hs.mpr.sub2$hm)
hs.mpr.sub2$hm <- RenameIdents(hs.mpr.sub2$hm, new.cluster.ids)
hs.mpr.sub2$hm$CellType <- Idents(hs.mpr.sub2$hm)
# Mesoderm-Somite
new.cluster.ids <- c("Somite:sclerotome.early","Somite:posterior presomitic mesoderm","Somite:posterior presomitic mesoderm",
                     "Somite:sclerotome.late","Somite:chondrogenic progenitors","Somite:syndetome","Somite:posterior presomitic mesoderm",
                     "Somite:migrating hypaxial muscle","Somite:posterior presomitic mesoderm","Somite:chondrogenic progenitors",
                     "Somite:sclerotome.late","Somite:posterior presomitic mesoderm","Somite:dermomyotome")
names(new.cluster.ids) <- levels(hs.mpr.sub2$somite)
hs.mpr.sub2$somite <- RenameIdents(hs.mpr.sub2$somite, new.cluster.ids)
hs.mpr.sub2$somite$CellType <- Idents(hs.mpr.sub2$somite)
# Mesoderm-Lateral plate mesoderm
new.cluster.ids <- c("Lateral plate mesoderm:hepatic stellate cell","Lateral plate mesoderm:epicardial derived cell",
                     "Lateral plate mesoderm:atrioventricular canal","Somite:undefined","Somite:undefined","Limb:forelimb",
                     "Limb:forelimb","Lateral plate mesoderm:hepatic stellate cell","Lateral plate mesoderm:atrioventricular canal",
                     "Limb:forelimb","Limb:forelimb","Limb:forelimb","Lateral plate mesoderm:hepatic stellate cell",
                     "Lateral plate mesoderm:atrioventricular canal","Somite:undefined","Lateral plate mesoderm:undefined",
                     "Lateral plate mesoderm:undefined","Lateral plate mesoderm:epicardium","Lateral plate mesoderm:hepatic stellate cell",
                     "Lateral plate mesoderm:epicardial derived cell","Lateral plate mesoderm:hepatic stellate cell",
                     "Lateral plate mesoderm:respiratory mesoderm")
names(new.cluster.ids) <- levels(hs.mpr.sub2$lpm)
hs.mpr.sub2$lpm <- RenameIdents(hs.mpr.sub2$lpm, new.cluster.ids)
hs.mpr.sub2$lpm$CellType <- Idents(hs.mpr.sub2$lpm)
# Mesoderm-Craniofacial
new.cluster.ids <- c("Head mesoderm:undefined","Craniofacial:PA1/2 middle","Craniofacial:PA1/2 distal","Craniofacial:PA3/4 middle",
                     "Craniofacial:PA3/4 distal","Lateral plate mesoderm:atrioventricular canal","Somite:undefined","Craniofacial:PA1/2 hand2",
                     "Craniofacial:PA1/2 proximal","Craniofacial:PA3/4 distal","Craniofacial:PA1/2 hand2","Somite:undefined",
                     "Craniofacial:frontonasal mesenchyme","Craniofacial:PA3/4 middle","Craniofacial:PA1/2 hand2","Head mesoderm:undefined",
                     "Craniofacial:PA1/2 middle","Craniofacial:frontonasal mesenchyme","Craniofacial:PA3/4 distal","Craniofacial:PA3/4 middle")
names(new.cluster.ids) <- levels(hs.mpr.sub2$craniofacial)
hs.mpr.sub2$craniofacial <- RenameIdents(hs.mpr.sub2$craniofacial, new.cluster.ids)
hs.mpr.sub2$craniofacial$CellType <- Idents(hs.mpr.sub2$craniofacial)
# Mesoderm-Mixed
new.cluster.ids <- c("Neuron:undefined","Neuron:undefined","Craniofacial:frontonasal mesenchyme","Craniofacial:PA3/4 middle",
                     "Craniofacial:undefined","Craniofacial:PA1/2 middle","Craniofacial:PA1/2 middle","Limb:hindlimb",
                     "Craniofacial:undefined","Limb:forelimb","Neuron:undefined","Craniofacial:PA1/2 hand2","Neuron:undefined",
                     "Intermediate mesoderm:undefined","Undefined:undefined","Craniofacial:PA1/2 hand2",
                     "Endothelium:liver sinusoidal endothelial cell","Lateral plate mesoderm:hepatic stellate cell")
names(new.cluster.ids) <- levels(hs.mpr.sub2$mixed)
hs.mpr.sub2$mixed <- RenameIdents(hs.mpr.sub2$mixed, new.cluster.ids)
hs.mpr.sub2$mixed$CellType <- Idents(hs.mpr.sub2$mixed)
# Mesoderm-Limb
new.cluster.ids <- c("Limb:hindlimb","Limb:hindlimb","Limb:hindlimb","Limb:hindlimb","Limb:hindlimb","Limb:hindlimb",
                     "Limb:hindlimb","Limb:hindlimb","Undefined:undefined")
names(new.cluster.ids) <- levels(hs.mpr.sub2$limb)
hs.mpr.sub2$limb <- RenameIdents(hs.mpr.sub2$limb, new.cluster.ids)
hs.mpr.sub2$limb$CellType <- Idents(hs.mpr.sub2$limb)
# Mesoderm-Erythroid
new.cluster.ids <- c("Blood:erythroid","Blood:erythroid","Blood:erythroid","Blood:erythroid",
                     "Blood:erythroid","Blood:erythroid","Blood:erythroid",
                     "Blood:erythroid","Blood:erythroid")
names(new.cluster.ids) <- levels(hs.mpr.sub2$erythroid)
hs.mpr.sub2$erythroid <- RenameIdents(hs.mpr.sub2$erythroid, new.cluster.ids)
hs.mpr.sub2$erythroid$CellType <- Idents(hs.mpr.sub2$erythroid)
# Epidermis
new.cluster.ids <- c("Liver:hepatocyte","Liver:hepatocyte","Craniofacial:undefined",
                     "Liver:hepatocyte","Somite:neuromesodermal progenitors",
                     "Head mesoderm:pericyte","Liver:hepatocyte","Epidermis:epidermis",
                     "Endothelium:liver sinusoidal endothelial cell",
                     "Epidermis:epidermis","Liver:hepatocyte","Placenta:extravillous trophoblasts","Undefined:undefined")
names(new.cluster.ids) <- levels(hs.mpr.sub1$ud)
hs.mpr.sub1$ud <- RenameIdents(hs.mpr.sub1$ud, new.cluster.ids)
hs.mpr.sub1$ud$CellType <- Idents(hs.mpr.sub1$ud)
# Ectoderm1
new.cluster.ids <- c("Neuron:v1","Neuron:GABAergic neuron",
                     "Neuron:V2b.rhombomere","Neuron:V0.rhombomere","Neuron:dB4",
                     "Neuron:dl4","Neuron:GABAergic neuron","Neuron:V0.rhombomere",
                     "Neuron:V1.rhombomere","Neuron:dB4","Neuron:V2b.rhombomere")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv1)
hs.mpr.sub2$nv1 <- RenameIdents(hs.mpr.sub2$nv1, new.cluster.ids)
hs.mpr.sub2$nv1$CellType <- Idents(hs.mpr.sub2$nv1)
# Ectoderm2
new.cluster.ids <- c("Neuron:dA1","Neuron:dl1","Neuron:dl1","Neuron:dl1","Neuron:dA1","Neuron:dB3","Neural progenitor:pA1",
                     "Neuron:dl1","Neuron:dl1","Neural progenitor:pA1","Neuron:dl2","Neuron:dA1")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv2)
hs.mpr.sub2$nv2 <- RenameIdents(hs.mpr.sub2$nv2, new.cluster.ids)
hs.mpr.sub2$nv2$CellType <- Idents(hs.mpr.sub2$nv2)
# Ectoderm3
new.cluster.ids <- c("Neuron:dA3","Neuron:dA3","Neuron:V2b.rhombomere","Neuron:dA3","Neuron:MNv","Neuron:dA3")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv3)
hs.mpr.sub2$nv3 <- RenameIdents(hs.mpr.sub2$nv3, new.cluster.ids)
hs.mpr.sub2$nv3$CellType <- Idents(hs.mpr.sub2$nv3)
# Ectoderm4
new.cluster.ids <- c("Schwann cells:melanocyte","Schwann cells:schwann progenitor",
                     "Neuron:lateral motor columns","Neuron:v2b",
                     "Neuron:dA2","Neuron:dl2","Neuron:lateral motor columns",
                     "Neuron:medial motor columns","Neuron:MNv",
                     "Neuron:lateral motor columns","Neuron:lateral motor columns","Neuron:V0.rhombomere")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv4)
hs.mpr.sub2$nv4 <- RenameIdents(hs.mpr.sub2$nv4, new.cluster.ids)
hs.mpr.sub2$nv4$CellType <- Idents(hs.mpr.sub2$nv4)
# Ectoderm5
new.cluster.ids <- c("Neural progenitor:pA1","Neural progenitor:pA1",
                     "Neural progenitor:pA2","Neural progenitor:pA3",
                     "Neural progenitor:pA2","Neural progenitor:pA3",
                     "Neural progenitor:pA1","Neural progenitor:pB4","Neural progenitor:pB4")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv5)
hs.mpr.sub2$nv5 <- RenameIdents(hs.mpr.sub2$nv5, new.cluster.ids)
hs.mpr.sub2$nv5$CellType <- Idents(hs.mpr.sub2$nv5)
# Ectoderm6
new.cluster.ids <- c("Neural progenitor:pB1","Neuron:dA2","Neuron:dB4","Neuron:V0.rhombomere","Neural progenitor:pB4",
                     "Neural progenitor:GABAergic neuron precursor","Neural progenitor:pB4","Neuron:midbrain Pitx2 neurons",
                     "Neural progenitor:pB1","Neuron:v2a","Neuron:medial motor columns","Neural progenitor:pB1",
                     "Undefined:undefined","Neural progenitor:pB1")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv6)
hs.mpr.sub2$nv6 <- RenameIdents(hs.mpr.sub2$nv6, new.cluster.ids)
hs.mpr.sub2$nv6$CellType <- Idents(hs.mpr.sub2$nv6)
# Ectoderm7
new.cluster.ids <- c("Neural progenitor:mesencephalon","Blood:erythroid",
                     "Undefined:undefined","Undefined:undefined","Blood:erythroid",
                     "Blood:erythroid","Neural progenitor:ventral diencephalon and ZLI",
                     "Neural progenitor:dorsal telencephalon", "Neural progenitor:dorsal diencephalon",
                     "Neural progenitor:dorsal telencephalon","Neural progenitor:roof plate.rhombomere",
                     "Neural progenitor:retinal pigment epithelium")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv7)
hs.mpr.sub2$nv7 <- RenameIdents(hs.mpr.sub2$nv7, new.cluster.ids)
hs.mpr.sub2$nv7$CellType <- Idents(hs.mpr.sub2$nv7)
# Ectoderm8
new.cluster.ids <- c("Neural progenitor:dp4","Neural progenitor:p0.rhombomere",
                     "Neural progenitor:dp4","Neural progenitor:dp4","Neural progenitor:dp4",
                     "Neural progenitor:dp4","Neural progenitor:p2.rhombomere",
                     "Neural progenitor:dp4","Neural progenitor:dp4","Neural progenitor:dp4",
                     "Neural progenitor:p3","Neural progenitor:dp4","Neural progenitor:dp4",
                     "Neural progenitor:p0.rhombomere","Neural progenitor:dp4")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv8)
hs.mpr.sub2$nv8 <- RenameIdents(hs.mpr.sub2$nv8, new.cluster.ids)
hs.mpr.sub2$nv8$CellType <- Idents(hs.mpr.sub2$nv8)
# Ectoderm9
new.cluster.ids <- c("Undefined:undefined","Undefined:undefined","Undefined:undefined",
                     "Undefined:undefined","Neural progenitor:dp3",
                     "Neural progenitor:dorsal diencephalon","Undefined:undefined",
                     "Neuron:dl1","Neural progenitor:dp3","Undefined:undefined",
                     "Neural progenitor:dp5","Undefined:undefined","Neuron:dILA/dBLa",
                     "Undefined:undefined","Undefined:undefined",
                     "Undefined:undefined","Neuron:dI2","Neuron:dA4",
                     "Undefined:undefined","Undefined:undefined","Undefined:undefined")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv9)
hs.mpr.sub2$nv9 <- RenameIdents(hs.mpr.sub2$nv9, new.cluster.ids)
hs.mpr.sub2$nv9$CellType <- Idents(hs.mpr.sub2$nv9)
# Ectoderm10
new.cluster.ids <- c("Neuron:trigeminal ganglia","Neuron:dorsal root ganglia",
                     "Neuron:dorsal root ganglia","Neuron:dorsal root ganglia",
                     "Neuron:dorsal root ganglia","Neuron:dorsal root ganglia",
                     "Neuron:dorsal root ganglia","Neuron:trigeminal ganglia",
                     "Neuron:dorsal root ganglia")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv10)
hs.mpr.sub2$nv10 <- RenameIdents(hs.mpr.sub2$nv10, new.cluster.ids)
hs.mpr.sub2$nv10$CellType <- Idents(hs.mpr.sub2$nv10)
# Ectoderm11
new.cluster.ids <- c("Neuron:dA2","Neuron:dA2","Neuron:dA2","Neuron:dA2","Neuron:dA2",
                     "Neuron:dA2","Neuron:dA2","Neuron:dA2","Neuron:dA2",
                     "Neuron:dA2","Neuron:dA2","Neuron:dA2","Neuron:V0.rhombomere","Neuron:dA2","Neuron:dA2")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv11)
hs.mpr.sub2$nv11 <- RenameIdents(hs.mpr.sub2$nv11, new.cluster.ids)
hs.mpr.sub2$nv11$CellType <- Idents(hs.mpr.sub2$nv11)
# Ectoderm12
new.cluster.ids <- c("Neural progenitor:roof plate","Neural progenitor:roof plate","Neural progenitor:roof plate",
                     "Neural progenitor:roof plate","Neural progenitor:roof plate.rhombomere")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv12)
hs.mpr.sub2$nv12 <- RenameIdents(hs.mpr.sub2$nv12, new.cluster.ids)
hs.mpr.sub2$nv12$CellType <- Idents(hs.mpr.sub2$nv12)
# Ectoderm13
new.cluster.ids <- c("Neural progenitor:floor plate","Neural progenitor:floor plate")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv13)
hs.mpr.sub2$nv13 <- RenameIdents(hs.mpr.sub2$nv13, new.cluster.ids)
hs.mpr.sub2$nv13$CellType <- Idents(hs.mpr.sub2$nv13)
# Ectoderm14
new.cluster.ids <- c("Neuron:intermediate V2 precursor","Neuron:intermediate V2 precursor")
names(new.cluster.ids) <- levels(hs.mpr.sub2$nv14)
hs.mpr.sub2$nv14 <- RenameIdents(hs.mpr.sub2$nv14, new.cluster.ids)
hs.mpr.sub2$nv14$CellType <- Idents(hs.mpr.sub2$nv14)
# Ectoderm-Schwann cells
new.cluster.ids <- c("Schwann cells:schwann progenitor","Schwann cells:undefined","Schwann cells:sympathetic neuron",
                     "Schwann cells:undefined","Schwann cells:schwann progenitor","Schwann cells:undefined",
                     "Schwann cells:schwann progenitor","Schwann cells:melanocyte","Schwann cells:undefined",
                     "Schwann cells:undefined","Schwann cells:undefined")
names(new.cluster.ids) <- levels(hs.mpr.sub1$ecto2)
hs.mpr.sub1$ecto2 <- RenameIdents(hs.mpr.sub1$ecto2, new.cluster.ids)
hs.mpr.sub1$ecto2$CellType <- Idents(hs.mpr.sub1$ecto2)
#
rm(new.cluster.ids)


### >>> 12. Merge cells in subclusters ----
tmp.list <- list(hs.mpr.sub1$meso2,hs.mpr.sub1$meso3,hs.mpr.sub1$meso4,hs.mpr.sub1$ud,
                 hs.mpr.sub2$hm,hs.mpr.sub2$somite,hs.mpr.sub2$lpm,
                 hs.mpr.sub2$craniofacial,hs.mpr.sub2$mixed,hs.mpr.sub2$limb,hs.mpr.sub2$erythroid,
                 hs.mpr.sub2$nv1,hs.mpr.sub2$nv2,hs.mpr.sub2$nv3,hs.mpr.sub2$nv4,
                 hs.mpr.sub2$nv5,hs.mpr.sub2$nv6,hs.mpr.sub2$nv7,hs.mpr.sub2$nv8,
                 hs.mpr.sub2$nv9,hs.mpr.sub2$nv10,hs.mpr.sub2$nv11,hs.mpr.sub2$nv12,
                 hs.mpr.sub2$nv13,hs.mpr.sub2$nv14,hs.mpr.sub1$ecto2)
hs.mpr.merge <- merge(hs.mpr.sub1$endo, y = tmp.list)
table(colnames(hs.mpr.merge) %in% colnames(hs.mpr))
length(hs.mpr.merge$CellType)-length(na.omit(hs.mpr.merge$CellType))
hs.mpr.merge@meta.data <- hs.mpr.merge@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC",
                                                    "nucleosome_signal","nucleosome_percentile","TSS.enrichment","TSS.percentile",
                                                    "percent.mt","Sample","Stage","raw_clusters","CellType")]


### >>> 13. Clustering again ----
hs.mpr.merge.sub1 <- list()
hs.mpr.merge.markers <- list()
# - Craniofacial
hs.mpr.merge.sub1$Craniofacial <- hs.mpr.merge[,grep("Craniofacial", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Craniofacial <- FindVariableFeatures(hs.mpr.merge.sub1$Craniofacial) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Craniofacial)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Craniofacial <- RunHarmony(hs.mpr.merge.sub1$Craniofacial, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Craniofacial, reduction = "harmony", ndims = 30)
dim.n <- 11
hs.mpr.merge.sub1$Craniofacial <- RunUMAP(hs.mpr.merge.sub1$Craniofacial, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Craniofacial$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Craniofacial$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Craniofacial, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Craniofacial <- FindAllMarkers(hs.mpr.merge.sub1$Craniofacial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Craniofacial %>% dplyr::filter(cluster == "9") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Craniofacial, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Craniofacial, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$craniofacial), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$craniofacial), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$endothelium), slot = "scale.data", reduction = "umap", ncol = 5, order = T)
DotPlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Craniofacial, features = unique(marker$limb), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Craniofacial, features = c("TBX5","PITX1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.merge.sub1$Craniofacial@meta.data %>%
  mutate(CellType = case_when(seurat_clusters == 1 | seurat_clusters == 9 ~ "Limb:forelimb",
                              seurat_clusters == 2 | seurat_clusters == 18 ~ "Undefined:undefined",
                              seurat_clusters == 17 ~ "Lateral plate mesoderm:sinoatrial node",
                              seurat_clusters == 20 ~ "Craniofacial:PA3/4 proximal",
                              seurat_clusters == 22 ~ "Endothelium:liver sinusoidal endothelial cell",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$Craniofacial@meta.data
hs.mpr.merge.sub1$Craniofacial$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Craniofacial$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Craniofacial, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Head mesoderm
hs.mpr.merge.sub1$HeadMesoderm <- hs.mpr.merge[,grep("Head mesoderm", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$HeadMesoderm <- FindVariableFeatures(hs.mpr.merge.sub1$HeadMesoderm) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$HeadMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$HeadMesoderm <- RunHarmony(hs.mpr.merge.sub1$HeadMesoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$HeadMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge.sub1$HeadMesoderm <- RunUMAP(hs.mpr.merge.sub1$HeadMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$HeadMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$HeadMesoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$HeadMesoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$HeadMesoderm <- FindAllMarkers(hs.mpr.merge.sub1$HeadMesoderm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$HeadMesoderm %>% dplyr::filter(cluster == "0") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$HeadMesoderm, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$HeadMesoderm, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$HeadMesoderm, features = unique(marker$hm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$HeadMesoderm, features = unique(marker$hm), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
hs.mpr.merge.sub1$HeadMesoderm@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(0,5,6,10,11,12,13,15,19) ~ "Head mesoderm:pericyte",
                              seurat_clusters %in% c(3,4,7,8,14,16) ~ "Head mesoderm:cranium",
                              seurat_clusters %in% c(9,17,20) ~ "Head mesoderm:head muscle",
                              seurat_clusters %in% c(1,18) ~ "Head mesoderm:head mesoderm",
                              seurat_clusters == 2 ~ "Head mesoderm:undefined",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$HeadMesoderm@meta.data
hs.mpr.merge.sub1$HeadMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$HeadMesoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$HeadMesoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Lateral plate mesoderm
hs.mpr.merge.sub1$LateralPlateMesoderm <- hs.mpr.merge[,grep("Lateral plate mesoderm", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$LateralPlateMesoderm <- FindVariableFeatures(hs.mpr.merge.sub1$LateralPlateMesoderm) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$LateralPlateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$LateralPlateMesoderm <- RunHarmony(hs.mpr.merge.sub1$LateralPlateMesoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.merge.sub1$LateralPlateMesoderm <- RunUMAP(hs.mpr.merge.sub1$LateralPlateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$LateralPlateMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$LateralPlateMesoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$LateralPlateMesoderm <- FindAllMarkers(hs.mpr.merge.sub1$LateralPlateMesoderm, only.pos = TRUE,
                                                            min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$LateralPlateMesoderm %>% dplyr::filter(cluster == "10") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = unique(c(marker$splanchnic.lpm,marker$somatic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("GATA4","TGFBI","HEY2","TWIST1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("LHX2","HLX","FOXF1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("TBX18","WT1","TNNI3","TNNT1","TCF21"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("SNAI2","TCF21","SOX9","THY1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("WNT2","NKX6-1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("SLC8A1","GATA6","NKX2-5","HAND2","HEY2"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("PITX1","TBX5"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("CDH5","CD34","PECAM1","XBP1","STAB2","CD36"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
DotPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = unique(marker$blood), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("CD68","CD14","CD163","CSF1R","ITGAM"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
DotPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = unique(marker$craniofacial), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$LateralPlateMesoderm, features = c("ALX1","ALX3"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
hs.mpr.merge.sub1$LateralPlateMesoderm@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(0,3,5,16,17) ~ "Lateral plate mesoderm:hepatic stellate cell",
                              seurat_clusters %in% c(1,12,19) ~ "Limb:undefined",
                              seurat_clusters %in% c(4,8,9) ~ "Lateral plate mesoderm:atrioventricular canal",
                              seurat_clusters %in% c(6,7,15) ~ "Lateral plate mesoderm:epicardial derived cell",
                              seurat_clusters == 2 ~ "Lateral plate mesoderm:respiratory mesoderm",
                              seurat_clusters == 13 ~ "Blood:macrophage",
                              seurat_clusters == 11 ~ "Endothelium:vascular endothelium",
                              seurat_clusters == 18 ~ "Endothelium:liver sinusoidal endothelium",
                              seurat_clusters == 14 ~ "Lateral plate mesoderm:epicardium",
                              seurat_clusters == 10 ~ "Undefined:undefined")) -> hs.mpr.merge.sub1$LateralPlateMesoderm@meta.data
hs.mpr.merge.sub1$LateralPlateMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$LateralPlateMesoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$LateralPlateMesoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Limb
hs.mpr.merge.sub1$Limb <- hs.mpr.merge[,grep("Limb", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Limb <- FindVariableFeatures(hs.mpr.merge.sub1$Limb) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Limb)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Limb <- RunHarmony(hs.mpr.merge.sub1$Limb, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Limb, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge.sub1$Limb <- RunUMAP(hs.mpr.merge.sub1$Limb, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Limb$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Limb$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Limb, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Limb <- FindAllMarkers(hs.mpr.merge.sub1$Limb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Limb %>% dplyr::filter(cluster == "10") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Limb, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Limb, features = unique(c(marker$splanchnic.lpm,marker$somatic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("GATA4","TGFBI","HEY2","TWIST1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("SHOX2","TBX18"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("SLC8A1","GATA6","NKX2-5","HAND2","HEY2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("TBX18","WT1","TNNI3","TNNT1","TCF21","GATA6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("SOX6","EZR"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("GATA6","IRX3"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.merge.sub1$Limb, features = unique(marker$limb), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Limb, features = c("PITX1","TBX5"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.merge.sub1$Limb@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(3,4,9,11,12,13) ~ "Limb:hindlimb",
                              seurat_clusters %in% c(0,6,8) ~ "Limb:forelimb",
                              seurat_clusters %in% c(2,5) ~ "Lateral plate mesoderm:sinoatrial node",
                              seurat_clusters %in% c(10,14) ~ "Lateral plate mesoderm:epicardium",
                              seurat_clusters == 1 ~ "Lateral plate mesoderm:anterior somatic LPM",
                              seurat_clusters == 7 ~ "Lateral plate mesoderm:endocardial derived cell")) -> hs.mpr.merge.sub1$Limb@meta.data
hs.mpr.merge.sub1$Limb$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Limb$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Limb, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Somite
hs.mpr.merge.sub1$Somite <- hs.mpr.merge[,grep("Somite", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Somite <- FindVariableFeatures(hs.mpr.merge.sub1$Somite) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Somite)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Somite <- RunHarmony(hs.mpr.merge.sub1$Somite, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Somite, reduction = "harmony", ndims = 30)
dim.n <- 15
hs.mpr.merge.sub1$Somite <- RunUMAP(hs.mpr.merge.sub1$Somite, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Somite$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Somite$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Somite, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$Somite, reduction = "umap", group.by = c("CellType.sub","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Somite <- FindAllMarkers(hs.mpr.merge.sub1$Somite, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Somite %>% dplyr::filter(cluster == "18") %>% dplyr::top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Somite, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Somite, features = unique(marker$somite), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("PAX1","PAX9","NKX3-2","SOX9"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("MYF5","MYF6","MYOD1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("PAX3","DMRT2"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("LBX1","PAX3","SIX1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("SCX","ETV4"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("FOXC1","FOXC2","SHISA2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("WNT5A","CYP26A1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.merge.sub1$Somite, features = unique(c(marker$splanchnic.lpm,marker$somatic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("SLC8A1","GATA6","NKX2-5","HAND2","HEY2"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("TBX18","WT1","TNNI3","TNNT1","TCF21","GATA6"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("SHOX2","TBX18"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
DotPlot(hs.mpr.merge.sub1$Somite, features = unique(marker$limb), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Somite, features = c("TBX5","PITX1"),
            slot = "scale.data", reduction = "umap", ncol = 2, order = T)
hs.mpr.merge.sub1$Somite@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(6,7,14,17) ~ "Somite:aPSM",
                              seurat_clusters %in% c(3,5) ~ "Somite:pPSM",
                              seurat_clusters %in% c(0,8,9) ~ "Somite:syndetome",
                              seurat_clusters %in% c(1,15) ~ "Somite:sclerotome.late",
                              seurat_clusters %in% c(13,22) ~ "Somite:sclerotome.early",
                              seurat_clusters %in% c(4,18) ~ "Somite:migrating hypaxial muscle",
                              seurat_clusters %in% c(2,19) ~ "Somite:chondrogenic progenitors",
                              seurat_clusters == 10 ~ "Lateral plate mesoderm:atrioventricular canal",
                              seurat_clusters == 11 ~ "Somite:undefined",
                              seurat_clusters == 12 ~ "Limb:undefined",
                              seurat_clusters == 16 ~ "Lateral plate mesoderm:epicardium",
                              seurat_clusters == 20 ~ "Somite:NMP",
                              seurat_clusters == 21 ~ "Somite:dermomyotome",
                              seurat_clusters == 23 ~ "Somite:myotome")) -> hs.mpr.merge.sub1$Somite@meta.data
hs.mpr.merge.sub1$Somite$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Somite$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Somite, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Endothelium
hs.mpr.merge.sub1$Endothelium <- hs.mpr.merge[,grep("Endothelium", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Endothelium <- FindVariableFeatures(hs.mpr.merge.sub1$Endothelium) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Endothelium)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Endothelium <- RunHarmony(hs.mpr.merge.sub1$Endothelium, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Endothelium, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.merge.sub1$Endothelium <- RunUMAP(hs.mpr.merge.sub1$Endothelium, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Endothelium$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Endothelium$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Endothelium, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Endothelium <- FindAllMarkers(hs.mpr.merge.sub1$Endothelium, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Endothelium %>% dplyr::filter(cluster == "0") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Endothelium, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Endothelium, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Endothelium, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Endothelium, features = unique(marker$endothelium)[7:12],
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Endothelium, features = c("XBP1","STAB2","CD36"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Endothelium, features = c("CDH5","CD34","PECAM1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
hs.mpr.merge.sub1$Endothelium@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(1,3,5,7,8,12) ~ "Endothelium:vascular endothelium",
                              seurat_clusters %in% c(4,6,15,20,21) ~ "Endothelium:liver sinusoidal endothelium",
                              seurat_clusters %in% c(10,11,17,18,19) ~ "Endothelium:arterial endothelium",
                              seurat_clusters %in% c(0,2,9,13,14,16) ~ "Endoderm:hepatocyte")) -> hs.mpr.merge.sub1$Endothelium@meta.data
hs.mpr.merge.sub1$Endothelium$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Endothelium$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Endothelium, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Endoderm
hs.mpr.merge.sub1$Endoderm <- hs.mpr.merge[,c(grep("Gut", hs.mpr.merge$CellType),
                                              grep("Liver", hs.mpr.merge$CellType),
                                              grep("Intermediate mesoderm", hs.mpr.merge$CellType),
                                              grep("Epidermis", hs.mpr.merge$CellType),
                                              grep("Placenta", hs.mpr.merge$CellType))]
hs.mpr.merge.sub1$Endoderm <- FindVariableFeatures(hs.mpr.merge.sub1$Endoderm) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Endoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Endoderm <- RunHarmony(hs.mpr.merge.sub1$Endoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Endoderm, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.merge.sub1$Endoderm <- RunUMAP(hs.mpr.merge.sub1$Endoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Endoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Endoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Endoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Endoderm <- FindAllMarkers(hs.mpr.merge.sub1$Endoderm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Endoderm %>% dplyr::filter(cluster == "14") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Endoderm, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Endoderm, features = unique(marker$endoderm), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = unique(marker$Endoderm),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("PAX8","EMX2","LHX1"), slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("PAX8","WT1"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("SOX2","KLF5"), slot = "scale.data", reduction = "umap", ncol = 2, order = T)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("SOX2","NKX2-1","SOX9","ETV5","MYCN","IRX1","IRX2","IRX3"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("HHEX","KLF5","ONECUT2","FOXA1","CDX2","GATA4"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("XBP1","CEBPA","ALB","TTR","APOA2","AFP","APOA1"),
            slot = "scale.data", reduction = "umap", ncol = 4, order = F)
DotPlot(hs.mpr.merge.sub1$Endoderm, features = unique(marker$im), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Endoderm, features = unique(marker$epidermis), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Endoderm, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Endoderm, features = c("HLA-G","DIO2","HTRA4","MFAP5","HPGD"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Endoderm, features = c("HLA-G","DIO2","HTRA4","MFAP5","HPGD"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = T)
hs.mpr.merge.sub1$Endoderm@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(0,1,2,3,4,7,15) ~ "Endoderm:hepatocyte",
                              seurat_clusters %in% c(5,12) ~ "Endoderm:foregut",
                              seurat_clusters %in% c(6,9) ~ "Intermediate mesoderm:undefined",
                              seurat_clusters %in% c(10,13) ~ "Intermediate mesoderm:renal epithelium",
                              seurat_clusters %in% c(8,11) ~ "Endoderm:stomach/duodenum",
                              seurat_clusters %in% c(14) ~ "Placenta:extravillous trophoblasts")) -> hs.mpr.merge.sub1$Endoderm@meta.data
hs.mpr.merge.sub1$Endoderm$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Endoderm$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Endoderm, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Blood
hs.mpr.merge.sub1$Blood <- hs.mpr.merge[,grep("Blood", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Blood <- FindVariableFeatures(hs.mpr.merge.sub1$Blood) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Blood)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Blood <- RunHarmony(hs.mpr.merge.sub1$Blood, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Blood, reduction = "harmony", ndims = 30)
dim.n <- 8
hs.mpr.merge.sub1$Blood <- RunUMAP(hs.mpr.merge.sub1$Blood, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Blood$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Blood$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Blood, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Blood <- FindAllMarkers(hs.mpr.merge.sub1$Blood, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Blood %>% dplyr::filter(cluster == "11") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Blood, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Blood, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Blood, features = c(unique(marker$blood),"HBE1","HBZ","HBG2","HBA2","HBG1"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
hs.mpr.merge.sub1$Blood@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(2,4,11,15,16,19) ~ "Blood:erythroid",
                              seurat_clusters %in% c(0,1,3,5,9,12,18) ~ "Blood:macrophage",
                              seurat_clusters %in% c(13) ~ "Blood:megakaryocyte",
                              seurat_clusters %in% c(14) ~ "Blood:neutrophil",
                              seurat_clusters %in% c(8,17) ~ "Blood:lymphocyte",
                              seurat_clusters %in% c(10) ~ "Blood:granulocyte progenitor",
                              seurat_clusters %in% c(6,7) ~ "Blood:dendritic cell")) -> hs.mpr.merge.sub1$Blood@meta.data
hs.mpr.merge.sub1$Blood$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Blood$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Blood, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Neuron
table(hs.mpr.merge$CellType)
hs.mpr.merge.sub1$Neuron <- hs.mpr.merge[,grep("Neuron:", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Neuron <- FindVariableFeatures(hs.mpr.merge.sub1$Neuron) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Neuron)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Neuron <- RunHarmony(hs.mpr.merge.sub1$Neuron, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Neuron, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.merge.sub1$Neuron <- RunUMAP(hs.mpr.merge.sub1$Neuron, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Neuron$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Neuron$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Neuron, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$Neuron, reduction = "umap", group.by = c("CellType.sub"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$Neuron, reduction = "umap", group.by = c("CellType.sub", "ident"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
hs.mpr.merge.markers$Neuron <- FindAllMarkers(hs.mpr.merge.sub1$Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Neuron %>% dplyr::filter(cluster == "6") %>% top_n(10, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("NKX6-1","PHOX2B","ISL1","PHOX2A"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
VlnPlot(hs.mpr.merge.sub1$Neuron, features = c("NKX6-1","PHOX2B","ISL1","PHOX2A"), flip = T, stack = T)
VlnPlot(hs.mpr.merge.sub1$Neuron, features = c("LHX2", "LHX9", "BARHL1","BARHL2", "HOXB6"), flip = T, stack = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("LHX2", "LHX9", "BARHL1","BARHL2", "HOXB6","EVX1"), slot = "scale.data",
            reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("FOXD3","FOXP2","LHX1","LHX5"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("GATA2","GATA3","PHOX2A","PHOX2B"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("LHX3","VSX2","SOX14","SOX21","VSX1","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("FOXD3","FOXP2", "LHX1", "LHX5"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("LHX2","LHX9","LHX1","LHX5","TLX3","PHOX2B","POU4F1","ISL1","BARHL1","BARHL2","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 6, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("ISL2","ISL1","FOXP1","HOXB6","LHX3","GATA2","GATA3","LHX1","LHX5"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("LHX2","LHX9","BARHL2","BARHL1","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("LHX1","LHX5","POU4F1","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("FOXD3","FOXP2"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.merge.sub1$Neuron, features = c("SOX9","FOXC1","FOXC2","SHISA2","WNT5A","TFPI","PRRX1"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)

DotPlot(hs.mpr.merge.sub1$Neuron, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Neuron, features = c(unique(marker$neuron)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Neuron, features = c(unique(marker$pgc)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp

subset(hs.mpr.merge.sub1$Neuron@meta.data, seurat_clusters=="20")$CellType.sub %>% table()
subset(hs.mpr.merge.sub1$Neuron@meta.data, CellType.sub=="dA4")$seurat_clusters %>% table()

hs.mpr.merge.sub1$Neuron@meta.data %>%
  mutate(Plot = case_when(seurat_clusters %in% c(13) ~ as.character(seurat_clusters),
                          TRUE ~ "Non")) -> hs.mpr.merge.sub1$Neuron@meta.data
DimPlot(hs.mpr.merge.sub1$Neuron, reduction = "umap", group.by = c("Plot"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5)
hs.mpr.merge.sub1$Neuron@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(0,1,2,3) ~ "Undefined:undefined",
                              seurat_clusters %in% c(11) ~ "Neuron:trigeminal ganglia",
                              seurat_clusters %in% c(15) ~ "Neuron:dorsal root ganglia",
                              seurat_clusters %in% c(21) ~ "Neuron:MNv",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$Neuron@meta.data
hs.mpr.merge.sub1$Neuron@meta.data %>%
  mutate(CellType = case_when(CellType %in% c("Neuron:dl1") ~ "Neuron:dI1",
                              CellType %in% c("Neuron:dl2") ~ "Neuron:dI2",
                              CellType %in% c("Neuron:dl4") ~ "Neuron:dI4",
                              CellType %in% c("Neuron:dl5") ~ "Neuron:dI5",
                              CellType %in% c("Neuron:dl6") ~ "Neuron:dI6",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$Neuron@meta.data
hs.mpr.merge.sub1$Neuron$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Neuron$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Neuron, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Neural Progenitor
table(hs.mpr.merge$CellType)
hs.mpr.merge.sub1$NeuralProgenitor <- hs.mpr.merge[,grep("Neural progenitor", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$NeuralProgenitor <- FindVariableFeatures(hs.mpr.merge.sub1$NeuralProgenitor) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$NeuralProgenitor)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$NeuralProgenitor <- RunHarmony(hs.mpr.merge.sub1$NeuralProgenitor, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "harmony", ndims = 30)
dim.n <- 15
hs.mpr.merge.sub1$NeuralProgenitor <- RunUMAP(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$NeuralProgenitor$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$NeuralProgenitor$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "umap", group.by = c("CellType.sub","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "umap", group.by = c("CellType.sub"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "umap", group.by = c("ident"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$NeuralProgenitor <- FindAllMarkers(hs.mpr.merge.sub1$NeuralProgenitor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$NeuralProgenitor %>% dplyr::filter(cluster == "0") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("CD68","CD14","CD163","CSF1R","ITGAM"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("LMX1A","MAFB","MSX1","MSX2","RSPO1","RSPO3","WNT4","WNT2B"), slot = "scale.data",
            reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("ASCL1","PAX3","HOXB6","PTF1A"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("PAX3","DBX2","PRDM13","PAX6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("PRDM8","NKX6-1","ASCL1","PAX6","DBX2","NKX2-2","NKX2-8","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("MSX1","OLIG3","ATOH1","PHOX2B","PAX3","DBX2","PRDM13","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("MSX1","OLIG3","ASCL1","HOXB6","PAX3","DBX2"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("MSX1","OLIG3","ATOH1","PHOX2B","ASCL1","PHOX2B"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("PRDM8","NKX6-1","ASCL1","NKX6-2"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("HHIP","RNF220"), slot = "scale.data", reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("SHH","WNT5A","FOXA2","NKX6-1","SPON1","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("MEIS2","SP5","TAL2","MAB21L2"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("MITF","PAX2","DCT"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("ASCL1","PAX3","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("PAX1","PAX9","FOXC1","FOXC2","SHISA2","TFPI","FN1","WNT4","WNT2B"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = unique(marker$craniofacial), slot = "scale.data",
            reduction = "umap", ncol = 6, order = T)
FeaturePlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c("TUBB3","ELAVL4"), slot = "scale.data",
            reduction = "umap", ncol = 2, order = F)

DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$np)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$neuron)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$blood)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$somite)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$hm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$craniofacial)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$endothelium)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$NeuralProgenitor, features = c(unique(marker$epidermis)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp

subset(hs.mpr.merge.sub1$NeuralProgenitor@meta.data, seurat_clusters=="18")$CellType.sub %>% table()
subset(hs.mpr.merge.sub1$NeuralProgenitor@meta.data, CellType.sub=="dp5")$seurat_clusters %>% table()

hs.mpr.merge.sub1$NeuralProgenitor@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(6,13) ~ "Craniofacial:undefined",
                              seurat_clusters %in% c(11) ~ "Neuron:dA2",
                              seurat_clusters %in% c(15) ~ "Somite:aPSM",
                              seurat_clusters %in% c(16) ~ "Neuron:dA1",
                              seurat_clusters %in% c(17) ~ "Blood:macrophage",
                              seurat_clusters %in% c(19) ~ "NeuralProgenitor:p1. rhombomere",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$NeuralProgenitor@meta.data
hs.mpr.merge.sub1$NeuralProgenitor$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$NeuralProgenitor$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$NeuralProgenitor, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# Schwann cells
table(hs.mpr.merge$CellType)
hs.mpr.merge.sub1$Schwann <- hs.mpr.merge[,grep("Schwann cells:", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Schwann <- FindVariableFeatures(hs.mpr.merge.sub1$Schwann) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Schwann)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Schwann <- RunHarmony(hs.mpr.merge.sub1$Schwann, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Schwann, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.merge.sub1$Schwann <- RunUMAP(hs.mpr.merge.sub1$Schwann, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Schwann$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Schwann$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Schwann, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$Schwann, reduction = "umap", group.by = c("CellType.sub"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
DimPlot(hs.mpr.merge.sub1$Schwann, reduction = "umap", group.by = c("CellType.sub", "ident"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.5)
hs.mpr.merge.markers$Schwann <- FindAllMarkers(hs.mpr.merge.sub1$Schwann, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Schwann %>% dplyr::filter(cluster == "7") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = unique(marker$sw), slot = "scale.data",
            reduction = "umap", ncol = 5, order = T)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("SOX10","FOXD3","MPZ","PLP1","PMP22"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("ASCL1","PHOX2B","ISL1","HAND2"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("MITF","TFAP2A","TFAP2B","PAX3","PMEL"), slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("HAND2","ASCL1","PHOX2A","PHOX2B","DCX","TUBB3","ELAVL4"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = F)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("ISL2","ISL1","FOXP1","HOXB6"), slot = "scale.data",
            reduction = "umap", ncol = 4, order = T)
FeaturePlot(hs.mpr.merge.sub1$Schwann, features = c("SLC8A1"), slot = "scale.data",
            reduction = "umap", ncol = 1, order = T)
DotPlot(hs.mpr.merge.sub1$Schwann, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c(unique(marker$sw)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c(unique(marker$neuron)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c(unique(marker$somite)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c("SLC1A3","CLDN11","KIAA1755","LYPD1","CYP1B1","LAMC3","AGPS","CHST2","LUM","ITPR2"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c("TAGLN","WFDC1","GMDS","FAM198B","CXCL12","CSRP1","FOXC2","ITGA1","VCL","FRMD4A"), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
DotPlot(hs.mpr.merge.sub1$Schwann, features = c(unique(marker$splanchnic.lpm)), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
subset(hs.mpr.merge.sub1$Schwann@meta.data, seurat_clusters=="8")$CellType.sub %>% table()
subset(hs.mpr.merge.sub1$Schwann@meta.data, CellType.sub=="sympathetic neuron")$seurat_clusters %>% table()
hs.mpr.merge.sub1$Schwann@meta.data %>%
  mutate(CellType = case_when(seurat_clusters %in% c(0) ~ "Somite:undefined",
                              seurat_clusters %in% c(1,2,3,4,9,12) ~ "Schwann cells:Schwann progenitor",
                              seurat_clusters %in% c(6,10,13) ~ "Neuron:LMC",
                              seurat_clusters %in% c(5) ~ "Schwann cells:sympathetic neuron",
                              seurat_clusters %in% c(7,8) ~ "Undefined:undefined",
                              seurat_clusters %in% c(11) ~ "Schwann cells:melanocyte",
                              TRUE ~ CellType)) -> hs.mpr.merge.sub1$Schwann@meta.data
hs.mpr.merge.sub1$Schwann$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Schwann$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Schwann, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# - Undefined
hs.mpr.merge.sub1$Undefined <- hs.mpr.merge[,grep("Undefined:", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Undefined <- FindVariableFeatures(hs.mpr.merge.sub1$Undefined) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Undefined)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Undefined <- RunHarmony(hs.mpr.merge.sub1$Undefined, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Undefined, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge.sub1$Undefined <- RunUMAP(hs.mpr.merge.sub1$Undefined, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Undefined$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Undefined$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Undefined, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Undefined <- FindAllMarkers(hs.mpr.merge.sub1$Undefined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Undefined %>% dplyr::filter(cluster == "7") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Undefined, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Undefined, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
hs.mpr.merge.sub1$Undefined@meta.data %>%
  mutate(CellType = case_when(TRUE ~ CellType)) -> hs.mpr.merge.sub1$Undefined@meta.data
hs.mpr.merge.sub1$Undefined$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Undefined$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Undefined, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# Fibroblast
hs.mpr.merge.sub1$Fibroblast <- hs.mpr.merge[,grep("Fibroblast:", hs.mpr.merge$CellType)]
hs.mpr.merge.sub1$Fibroblast <- FindVariableFeatures(hs.mpr.merge.sub1$Fibroblast) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Fibroblast)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Fibroblast <- RunHarmony(hs.mpr.merge.sub1$Fibroblast, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Fibroblast, reduction = "harmony", ndims = 30)
dim.n <- 7
hs.mpr.merge.sub1$Fibroblast <- RunUMAP(hs.mpr.merge.sub1$Fibroblast, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Fibroblast$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Fibroblast$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Fibroblast, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Fibroblast <- FindAllMarkers(hs.mpr.merge.sub1$Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Fibroblast %>% dplyr::filter(cluster == "7") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Fibroblast, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Fibroblast, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Fibroblast, features = c("COL1A1","COL1A2","COL3A1","POSTN","DCN"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
DotPlot(hs.mpr.merge.sub1$Fibroblast, features = unique(marker$endothelium), cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
FeaturePlot(hs.mpr.merge.sub1$Fibroblast, features = unique(marker$endothelium)[7:12],
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Fibroblast, features = c("XBP1","STAB2","CD36"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
FeaturePlot(hs.mpr.merge.sub1$Fibroblast, features = c("CDH5","CD34","PECAM1"),
            slot = "scale.data", reduction = "umap", ncol = 3, order = F)
hs.mpr.merge.sub1$Fibroblast@meta.data %>%
  mutate(CellType = case_when(TRUE ~ "Fibroblast:fibroblast")) -> hs.mpr.merge.sub1$Fibroblast@meta.data
hs.mpr.merge.sub1$Fibroblast$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Fibroblast$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Fibroblast, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
# Mixed
hs.mpr.merge.sub1$Mixed <- hs.mpr.merge[,c(grep("Heart:cardiomyocytes", hs.mpr.merge$CellType),
                                           grep("Peripheral surface ectoderm:limb AER", hs.mpr.merge$CellType))]
hs.mpr.merge.sub1$Mixed <- FindVariableFeatures(hs.mpr.merge.sub1$Mixed) %>%
  ScaleData(rownames(hs.mpr.merge.sub1$Mixed)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge.sub1$Mixed <- RunHarmony(hs.mpr.merge.sub1$Mixed, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge.sub1$Mixed, reduction = "harmony", ndims = 30)
dim.n <- 5
hs.mpr.merge.sub1$Mixed <- RunUMAP(hs.mpr.merge.sub1$Mixed, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge.sub1$Mixed$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Mixed$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Mixed, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)
hs.mpr.merge.markers$Mixed <- FindAllMarkers(hs.mpr.merge.sub1$Mixed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pd.gene <- hs.mpr.merge.markers$Mixed %>% dplyr::filter(cluster == "7") %>% top_n(20, avg_log2FC)
FeaturePlot(hs.mpr.merge.sub1$Mixed, features = pd.gene$gene, slot = "scale.data",
            reduction = "umap", ncol = 5, order = F)
DotPlot(hs.mpr.merge.sub1$Mixed, features = marker$super.cluster[1:38], cols = c("#eaeaea", "#fc0330"),
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() + theme_dp
table(hs.mpr.merge.sub1$Mixed$CellType.sub)
hs.mpr.merge.sub1$Mixed@meta.data %>%
  mutate(CellType = case_when(TRUE ~ CellType)) -> hs.mpr.merge.sub1$Mixed@meta.data
hs.mpr.merge.sub1$Mixed$CellType.sub <- str_split_fixed(hs.mpr.merge.sub1$Mixed$CellType, ":", "2")[,2]
DimPlot(hs.mpr.merge.sub1$Mixed, reduction = "umap", group.by = c("Sample","CellType.sub","ident"),
        ncol = 3, label = TRUE, repel = TRUE, pt.size = 0.5, label.box = T)


### >>> 14. Merge cells in re-clustering ----
tmp.list <- list(hs.mpr.merge.sub1$NeuralProgenitor, hs.mpr.merge.sub1$Schwann,
                 hs.mpr.merge.sub1$Craniofacial, hs.mpr.merge.sub1$HeadMesoderm, hs.mpr.merge.sub1$LateralPlateMesoderm,
                 hs.mpr.merge.sub1$Limb, hs.mpr.merge.sub1$Somite, hs.mpr.merge.sub1$Endothelium, hs.mpr.merge.sub1$Blood,
                 hs.mpr.merge.sub1$Endoderm, hs.mpr.merge.sub1$Undefined, hs.mpr.merge.sub1$Fibroblast, hs.mpr.merge.sub1$Mixed)
hs.mpr.merge2 <- merge(hs.mpr.merge.sub1$Neuron, y = tmp.list); rm(tmp.list)
table(colnames(hs.mpr.merge2) %in% colnames(hs.mpr))
table(colnames(hs.mpr.merge) %in% colnames(hs.mpr.merge2))
table(hs.mpr.merge2$CellType) %>% as.data.frame() -> cell.stat
write.csv(cell.stat, "/home/yhw/bioinfo/project-mine/MultiOmics/R/Table/Hs_embryo_cell_type_stat.csv", row.names = F)
hs.mpr.merge2@meta.data <- hs.mpr.merge2@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC",
                                                      "nucleosome_signal","nucleosome_percentile","TSS.enrichment","TSS.percentile",
                                                      "percent.mt","Sample","Stage","raw_clusters","CellType")]

### >>> 15. Clustering again ----
res.out <- file.path(getwd(), "R/Graphs/hs_MPR/cell_annotation_seurat")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }
hs.mpr.merge2.sub1 <- list()
hs.mpr.merge2.markers <- list()
pd.col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
            "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
            "#00FBFF","#00BCBF","#D277FD","#AE00FF","#7789FD","#203EFF",
            "#FC8ACD","#FF20A4","#8AFCE9","#1CFFD9","#97FF99","#1CFF20",
            "#7B7B7B","#BFBFBF")
# - Craniofacial
hs.mpr.merge2.sub1$Craniofacial <- hs.mpr.merge2[,grep("Craniofacial", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Craniofacial <- FindVariableFeatures(hs.mpr.merge2.sub1$Craniofacial) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Craniofacial)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Craniofacial <- RunHarmony(hs.mpr.merge2.sub1$Craniofacial, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Craniofacial, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Craniofacial <- RunUMAP(hs.mpr.merge2.sub1$Craniofacial, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Craniofacial$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Craniofacial$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Craniofacial) <- hs.mpr.merge2.sub1$Craniofacial$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Craniofacial_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Craniofacial, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Craniofacial <- FindAllMarkers(hs.mpr.merge2.sub1$Craniofacial, only.pos = TRUE,
                                                     min.pct = 0.25, logfc.threshold = 0.25)
# - Head mesoderm
hs.mpr.merge2.sub1$HeadMesoderm <- hs.mpr.merge2[,grep("Head mesoderm", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$HeadMesoderm <- FindVariableFeatures(hs.mpr.merge2.sub1$HeadMesoderm) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$HeadMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$HeadMesoderm <- RunHarmony(hs.mpr.merge2.sub1$HeadMesoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$HeadMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$HeadMesoderm <- RunUMAP(hs.mpr.merge2.sub1$HeadMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$HeadMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$HeadMesoderm$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$HeadMesoderm) <- hs.mpr.merge2.sub1$HeadMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_HeadMesoderm_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$HeadMesoderm, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$HeadMesoderm <- FindAllMarkers(hs.mpr.merge2.sub1$HeadMesoderm, only.pos = TRUE,
                                                     min.pct = 0.25, logfc.threshold = 0.25)
# - Lateral plate mesoderm
hs.mpr.merge2.sub1$LateralPlateMesoderm <- hs.mpr.merge2[,grep("Lateral plate mesoderm", hs.mpr.merge2$CellType)]
table(hs.mpr.merge2.sub1$LateralPlateMesoderm$CellType)
hs.mpr.merge2.sub1$LateralPlateMesoderm <- FindVariableFeatures(hs.mpr.merge2.sub1$LateralPlateMesoderm) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$LateralPlateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$LateralPlateMesoderm <- RunHarmony(hs.mpr.merge2.sub1$LateralPlateMesoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$LateralPlateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$LateralPlateMesoderm <- RunUMAP(hs.mpr.merge2.sub1$LateralPlateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$LateralPlateMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$LateralPlateMesoderm$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$LateralPlateMesoderm) <- hs.mpr.merge2.sub1$LateralPlateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_LateralPlateMesoderm_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$LateralPlateMesoderm, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$LateralPlateMesoderm <- FindAllMarkers(hs.mpr.merge2.sub1$LateralPlateMesoderm, only.pos = TRUE,
                                                             min.pct = 0.25, logfc.threshold = 0.25)
# - Limb
hs.mpr.merge2.sub1$Limb <- hs.mpr.merge2[,grep("Limb", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Limb <- FindVariableFeatures(hs.mpr.merge2.sub1$Limb) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Limb)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Limb <- RunHarmony(hs.mpr.merge2.sub1$Limb, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Limb, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Limb <- RunUMAP(hs.mpr.merge2.sub1$Limb, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Limb$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Limb$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Limb) <- hs.mpr.merge2.sub1$Limb$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Limb_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Limb, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Limb <- FindAllMarkers(hs.mpr.merge2.sub1$Limb, only.pos = TRUE,
                                             min.pct = 0.25, logfc.threshold = 0.25)
# - Somite
hs.mpr.merge2.sub1$Somite <- hs.mpr.merge2[,grep("Somite", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Somite <- FindVariableFeatures(hs.mpr.merge2.sub1$Somite) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Somite)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Somite <- RunHarmony(hs.mpr.merge2.sub1$Somite, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Somite, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Somite <- RunUMAP(hs.mpr.merge2.sub1$Somite, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Somite$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Somite$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Somite) <- hs.mpr.merge2.sub1$Somite$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Somite_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Somite, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Somite <- FindAllMarkers(hs.mpr.merge2.sub1$Somite, only.pos = TRUE,
                                               min.pct = 0.25, logfc.threshold = 0.25)
# - Endothelium
hs.mpr.merge2.sub1$Endothelium <- hs.mpr.merge2[,grep("Endothelium", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Endothelium <- FindVariableFeatures(hs.mpr.merge2.sub1$Endothelium) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Endothelium)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Endothelium <- RunHarmony(hs.mpr.merge2.sub1$Endothelium, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Endothelium, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Endothelium <- RunUMAP(hs.mpr.merge2.sub1$Endothelium, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Endothelium$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Endothelium$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Endothelium) <- hs.mpr.merge2.sub1$Endothelium$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endothelium_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Endothelium, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Endothelium <- FindAllMarkers(hs.mpr.merge2.sub1$Endothelium, only.pos = TRUE,
                                                    min.pct = 0.25, logfc.threshold = 0.25)
# - Blood
hs.mpr.merge2.sub1$Blood <- hs.mpr.merge2[,grep("Blood", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Blood <- FindVariableFeatures(hs.mpr.merge2.sub1$Blood) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Blood)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Blood <- RunHarmony(hs.mpr.merge2.sub1$Blood, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Blood, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Blood <- RunUMAP(hs.mpr.merge2.sub1$Blood, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Blood$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Blood$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Blood) <- hs.mpr.merge2.sub1$Blood$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Blood_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Blood, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Blood <- FindAllMarkers(hs.mpr.merge2.sub1$Blood, only.pos = TRUE,
                                              min.pct = 0.25, logfc.threshold = 0.25)
# - Intermediate mesoderm
hs.mpr.merge2.sub1$IntermediateMesoderm <- hs.mpr.merge2[,grep("Intermediate mesoderm", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$IntermediateMesoderm <- FindVariableFeatures(hs.mpr.merge2.sub1$IntermediateMesoderm) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$IntermediateMesoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$IntermediateMesoderm <- RunHarmony(hs.mpr.merge2.sub1$IntermediateMesoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$IntermediateMesoderm, reduction = "harmony", ndims = 30)
dim.n <- 5
hs.mpr.merge2.sub1$IntermediateMesoderm <- RunUMAP(hs.mpr.merge2.sub1$IntermediateMesoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$IntermediateMesoderm$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$IntermediateMesoderm$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$IntermediateMesoderm) <- hs.mpr.merge2.sub1$IntermediateMesoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_IntermediateMesoderm_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$IntermediateMesoderm, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$IntermediateMesoderm <- FindAllMarkers(hs.mpr.merge2.sub1$IntermediateMesoderm, only.pos = TRUE,
                                                             min.pct = 0.25, logfc.threshold = 0.25)
# - Endoderm
hs.mpr.merge2.sub1$Endoderm <- hs.mpr.merge2[,grep("Endoderm", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Endoderm <- FindVariableFeatures(hs.mpr.merge2.sub1$Endoderm) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Endoderm)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Endoderm <- RunHarmony(hs.mpr.merge2.sub1$Endoderm, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Endoderm, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Endoderm <- RunUMAP(hs.mpr.merge2.sub1$Endoderm, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Endoderm$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Endoderm$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Endoderm) <- hs.mpr.merge2.sub1$Endoderm$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Endoderm_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Endoderm, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Endoderm <- FindAllMarkers(hs.mpr.merge2.sub1$Endoderm, only.pos = TRUE,
                                                 min.pct = 0.25, logfc.threshold = 0.25)
# - Neuron
hs.mpr.merge2.sub1$Neuron <- hs.mpr.merge2[,grep("Neuron", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Neuron <- FindVariableFeatures(hs.mpr.merge2.sub1$Neuron) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Neuron)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Neuron <- RunHarmony(hs.mpr.merge2.sub1$Neuron, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Neuron, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Neuron <- RunUMAP(hs.mpr.merge2.sub1$Neuron, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Neuron$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Neuron$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Neuron) <- hs.mpr.merge2.sub1$Neuron$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Neuron_with_cell_types_by_seurat.pdf"), height = 5, width = 15)
DimPlot(hs.mpr.merge2.sub1$Neuron, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.15, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Neuron <- FindAllMarkers(hs.mpr.merge2.sub1$Neuron, only.pos = TRUE,
                                               min.pct = 0.25, logfc.threshold = 0.25)
# - Neural progenitor
hs.mpr.merge2.sub1$NeuralProgenitor <- hs.mpr.merge2[,grep("Neural progenitor", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$NeuralProgenitor <- FindVariableFeatures(hs.mpr.merge2.sub1$NeuralProgenitor) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$NeuralProgenitor)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$NeuralProgenitor <- RunHarmony(hs.mpr.merge2.sub1$NeuralProgenitor, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$NeuralProgenitor, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$NeuralProgenitor <- RunUMAP(hs.mpr.merge2.sub1$NeuralProgenitor, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$NeuralProgenitor$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$NeuralProgenitor$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$NeuralProgenitor) <- hs.mpr.merge2.sub1$NeuralProgenitor$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_NeuralProgenitor_with_cell_types_by_seurat.pdf"), height = 5, width = 15)
DimPlot(hs.mpr.merge2.sub1$NeuralProgenitor, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.15, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$NeuralProgenitor <- FindAllMarkers(hs.mpr.merge2.sub1$NeuralProgenitor, only.pos = TRUE,
                                                         min.pct = 0.25, logfc.threshold = 0.25)
# - Schwann cells
hs.mpr.merge2.sub1$Schwann <- hs.mpr.merge2[,grep("Schwann cells", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Schwann <- FindVariableFeatures(hs.mpr.merge2.sub1$Schwann) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Schwann)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Schwann <- RunHarmony(hs.mpr.merge2.sub1$Schwann, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Schwann, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Schwann <- RunUMAP(hs.mpr.merge2.sub1$Schwann, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Schwann$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Schwann$CellType, ":", "2")[,2]
Idents(hs.mpr.merge2.sub1$Schwann) <- hs.mpr.merge2.sub1$Schwann$CellType.sub
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Schwann_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Schwann, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Schwann <- FindAllMarkers(hs.mpr.merge2.sub1$Schwann, only.pos = TRUE,
                                                min.pct = 0.25, logfc.threshold = 0.25)
# - Undefined
hs.mpr.merge2.sub1$Undefined <- hs.mpr.merge2[,grep("undefined", hs.mpr.merge2$CellType)]
hs.mpr.merge2.sub1$Undefined <- FindVariableFeatures(hs.mpr.merge2.sub1$Undefined) %>%
  ScaleData(rownames(hs.mpr.merge2.sub1$Undefined)) %>% RunPCA(verbose = FALSE)
set.seed(100)
hs.mpr.merge2.sub1$Undefined <- RunHarmony(hs.mpr.merge2.sub1$Undefined, group.by.vars = "Sample")
ElbowPlot(hs.mpr.merge2.sub1$Undefined, reduction = "harmony", ndims = 30)
dim.n <- 10
hs.mpr.merge2.sub1$Undefined <- RunUMAP(hs.mpr.merge2.sub1$Undefined, reduction = "harmony", dims = 1:dim.n) %>%
  FindNeighbors(reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.mpr.merge2.sub1$Undefined$CellType.sub <- str_split_fixed(hs.mpr.merge2.sub1$Undefined$CellType, ":", "2")[,2]
dev.off()
pdf(file.path(res.out, "Sub-clustering_of_Undefined_with_cell_types_by_seurat.pdf"), height = 5, width = 13)
DimPlot(hs.mpr.merge2.sub1$Undefined, reduction = "umap", group.by = c("Sample","ident"),
        ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.75, label.box = T, cols = pd.col)
dev.off()
hs.mpr.merge2.markers$Undefined <- FindAllMarkers(hs.mpr.merge2.sub1$Undefined, only.pos = TRUE,
                                                  min.pct = 0.25, logfc.threshold = 0.25)


### >>> 16. add cell type information into harmony treated object
if (all(colnames(hs.mpr.merge2)[match(colnames(hs.mpr), colnames(hs.mpr.merge2))] == colnames(hs.mpr))) {
  hs.mpr$CellType <- hs.mpr.merge2$CellType
}
hs.embryo <- subset(hs.mpr, CellType != "Placenta:extravillous trophoblasts")
DefaultAssay(hs.embryo) <- "RNA"
hs.embryo <- NormalizeData(hs.embryo) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.embryo <- ScaleData(hs.embryo, verbose = FALSE, features = rownames(hs.embryo))
hs.embryo <- RunPCA(hs.embryo, npcs = 50, verbose = FALSE, features = VariableFeatures(hs.embryo))
ElbowPlot(hs.embryo, ndims = 30)
dim.n <- 25
hs.embryo <- FindNeighbors(hs.embryo, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.embryo <- RunUMAP(hs.embryo, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
hs.embryo <- RunHarmony(hs.embryo, group.by.vars = "Sample")
hs.embryo <- RunUMAP(hs.embryo, reduction = "harmony", dims = 1:dim.n, return.model = T)
hs.embryo <- FindNeighbors(hs.embryo, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 2)
pdf(file.path(res.out, "Scatter_plot_to_show_corrected_clustering.pdf"), height = 5, width = 12.5)
DimPlot(hs.embryo, reduction = "umap", group.by = c("Sample", "ident"), ncol = 2, label = TRUE, repel = TRUE, pt.size = 0.25)
dev.off()
hs.embryo@meta.data %>%
  mutate(CellType=case_when(CellType%in%"NeuralProgenitor:p1. rhombomere" ~ "Neural progenitor:p1.rhombomere",
                            TRUE ~ CellType)) -> hs.embryo@meta.data
hs.embryo@meta.data %>%
  mutate(CellType=case_when(CellType%in%"Endothelium:liver sinusoidal endothelial cell" ~ "Endothelium:liver sinusoidal endothelium",
                            CellType%in%"Neuron:LMC" ~ "Neuron:lateral motor columns",
                            TRUE ~ CellType)) -> hs.embryo@meta.data
saveRDS(hs.embryo, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.mpr.new.final.rds")



### ===========================================================
### 5th step: Calculation of ATAC-seq signal on repeat elements ----
### ===========================================================

### >>> 1. load repeat annotation file
re.anno <- read.table("/home/cmq/genome/ucsc/human/hg38/repeat/GRCh38_RepeatMasker_repeat_anno_quan.bed", header = F, stringsAsFactors = F)
colnames(re.anno) <- c("chr", "start", "end", "id", "strand")
re.anno <- separate(re.anno, "id", c("class","family","subfamily","locus"), sep = ":", remove = F)
re.anno <- subset(re.anno, class%in%c("LINE","SINE","LTR","Retroposon"))
re.anno.gr <- GRanges(re.anno[,1], IRanges(re.anno[,2], re.anno[,3]), name=re.anno$id, strand = re.anno$strand)

### >>> 2. calculating activity of repeats via snapatac
# - load snap file
snap.path <- list.files(getwd(), "snap$", recursive = T, full.names = T)
names(snap.path) <- c("Day18", "Day22", "Hs.5W.1", "Hs.5W.2", "Hs.4W.1", "Hs.6W.1")
qc.path <- list.files(getwd(), "per_barcode_metrics.csv$", recursive = T, full.names = T)
qc.path <- qc.path[-5]
names(qc.path) <- names(snap.path)
snap.path <- snap.path[-1:-2]
qc.path <- qc.path[-1:-2]
# - create snap object
library(SnapATAC)
snap <- list()
for (i in seq(1, length(snap.path))) {
  snap[[i]] <- createSnap(
    file=snap.path[i], sample=names(snap.path)[i], num.cores=1)
  names(snap)[i] <- names(snap.path)[i]
}; rm(i)
for (i in 1:4) {
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
for (i in 1:4) {
  snap[[i]] <- addBmatToSnap(snap[[i]])
}; rm(i)
# generation
for (i in 1:4) {
  snap[[i]] <- createGmatFromMat(obj=snap[[i]], input.mat="bmat", genes=re.anno.gr, do.par=TRUE, num.cores=16)
}
saveRDS(snap, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.snap.rds")
saveRDS(re.anno.gr, "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.repeat.elements.grange.rds")
# saving
ge.act <- list()
ge.act$raw <- list()
for (i in 1:6) {
  ge.act$raw[[i]] <- snap[[i]]@gmat
  names(ge.act$raw)[i] <- names(snap)[i]
  rownames(ge.act$raw[[i]]) <- paste(names(ge.act$raw)[i], "_", rownames(ge.act$raw[[i]]), sep = "")
}; rm(i)



### ===============================
### 6th step: Integrate all species ----
### ===============================

# human
hs <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.embryo.final.version.rds")
hs[['ATAC']] <- NULL
count <- GetAssayData(hs, slot = "count")
meta <- hs@meta.data
hs <- CreateSeuratObject(counts = count, meta.data = meta, project = "hs")
# mouse
mm <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/mm.embryo.final.version.rds")
count <- GetAssayData(mm, slot = "count")
rownames(count) <- str_to_upper(rownames(count))
meta <- mm@meta.data
mm <- CreateSeuratObject(counts = count, meta.data = meta, project = "mm")
mm$SuperCluster <- mm$SuperCluster.Final
# pig
ss <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/ss.emb.final.with.alra.rds")
ss[['ATAC']] <- NULL
count <- GetAssayData(ss, slot = "count", assay = "RNA")
meta <- ss@meta.data
ss <- CreateSeuratObject(counts = count, meta.data = meta, project = "ss")
rm(count, meta)
# filter genes
common.gene <- Reduce(intersect, list(rownames(hs), rownames(mm), rownames(ss)))
hs <- hs[rownames(hs) %in% common.gene, ]
mm <- mm[rownames(mm) %in% common.gene, ]
ss <- ss[rownames(ss) %in% common.gene, ]
# merge
sr.merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name, y@project.name)), list(hs, mm, ss))
qsave(sr.merge, "/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/RDS/Hs_Embryo_Multiome_all_species.qs")
scConversion(sc.in = sr.merge, sc.out = "/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/all_species.h5ad",
             from = "seurat", to = "anndata")
sr.merge <- qread("/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/all_species.qs")
tmp <- CreateSeuratObject(counts = GetAssayData(sr.merge, slot = "count"), meta.data = sr.merge@meta.data)
dior::write_h5(tmp, file = paste0("/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/RDS/Hs_Embryo_Multiome_all_species.h5"), object.type = 'seurat')
library(dior)
table(sr.merge$Stage, sr.merge$CellType.sub)
# load corrected data
adata <- dior::read_h5(file="/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/came/_temp/('human', 'mouse')-(01-01 19.47.17)/adt_hidden_cell.h5", 
                       target.object = 'seurat')
adata$Species <- gsub("_.*", "", gsub("SeuratProject_", "", adata$original_name))
table(adata$REF, adata$Species)
table(sr.merge$Stage, sr.merge$SuperCluster)
adata$CellType <- paste0(adata$Species, "-", adata$REF)
Idents(adata) <- adata$CellType
# filter cells
DimPlot(adata, reduction = "umap", group.by = c("REF"), raster = F,
        pt.size = 0.25, label.box = T)
tmp <- data.frame(xmin = c(3, 0, -5), xmax = c(9, 4, 7), 
                  ymin = c(-2, 10, 1), ymax = c(4, 15, 7), 
                  group = c("endo", "sch", "np"))
tmp <- ExtractCellByPos(object = adata, object.type = "seurat", dim.name = "umap",
                        group.coor = tmp, group.name = tmp$group, pt.size = 0.25)
del.cells <- c(intersect(tmp$id$endo, 
                         colnames(adata)[adata$REF %in% c("Somite", "Intermediate mesoderm", "Lateral plate mesoderm")]),
               intersect(tmp$id$sch, 
                         colnames(adata)[adata$REF %in% c("Neuron", "Somite", "Lateral plate mesoderm")]),
               intersect(tmp$id$np, 
                         colnames(adata)[adata$REF %in% c("Craniofacial", "Heart", "Limb", "Intermediate mesoderm", 
                                                          "Somite", "Lateral plate mesoderm", "Blood")]))
adata.fl <- adata[, !(colnames(adata) %in% del.cells)]
# subsample
adata.fl <- subset(adata.fl, downsample = 5000)
table(adata.fl$Species)
adata.fl <- subset(adata.fl, REF != "Undefined")
table(adata.fl$Species)
adata.fl$REF <- factor(adata.fl$REF,
                       levels = c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", 
                                  "Craniofacial", "Head mesoderm", "Somite", "Intermediate mesoderm", 
                                  "Limb","Lateral plate mesoderm",
                                  "Heart", "Endothelium", "Blood", "Endoderm"))
dev.off()
library(cols4all)
supercluster.col <- Fix.Colors(sr.obj = adata, group.by = "REF", colors = sample(c4a("classic20", 14)))
supercluster.col <- c("#8c564b", "#d62728", "#ffbb78", "#ff7f0e", "#aec7e8", "#9467bd", "#1f77b4",
                      "#2ca02c", "#ff9896", "#98df8a", "#f7b6d2", "#c49c94", "#e377c2", "#c5b0d5", "#E6E8EA")
names(supercluster.col) <- c("Neural progenitor", "Neuron", "Peripheral surface ectoderm", "Schwann cells", 
                             "Craniofacial", "Head mesoderm", "Somite", "Intermediate mesoderm", 
                             "Limb","Lateral plate mesoderm", "Heart", "Endothelium", "Blood", "Endoderm")
adata.fl <- qread("/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/RDS/Hs_Embryo_Multiome_all_species_plot.qs")
table(adata.fl$Species)
p1 <- CellDimPlot(srt = subset(adata.fl, Species == "hs"), group.by = "REF", palcolor = supercluster.col,
                  pt.size = 0.1, pt.alpha = 1, label_insitu = F, label = F, label_repel = T, 
                  reduction = "UMAP", theme_use = "theme_blank", raster = F)
p2 <- CellDimPlot(srt = subset(adata.fl, Species == "mm"), group.by = "REF", palcolor = supercluster.col,
                  pt.size = 0.1, pt.alpha = 1, label_insitu = F, label = F, label_repel = T, 
                  reduction = "UMAP", theme_use = "theme_blank", raster = F)
p3 <- CellDimPlot(srt = subset(adata.fl, Species == "ss"), group.by = "REF", palcolor = supercluster.col,
                  pt.size = 0.1, pt.alpha = 1, label_insitu = F, label = F, label_repel = T, 
                  reduction = "UMAP", theme_use = "theme_blank", raster = F)
pdf("/home/yhw/bioinfo/project-mine/MultiOmics/R/Graphs/came/Scatter_plot_after_correction.pdf", height = 6, width = 20)
p1 + p2 + p3
dev.off()
qsave(adata.fl, "/home/yhw/bioinfo/project-mine/MultiOmics/CodeData/RDS/Hs_Embryo_Multiome_all_species_plot.qs")



### ====================
### Last part: Save Data ----
### ====================
hs.mpr@meta.data %>%
  mutate(CellType=case_when(CellType%in%"NeuralProgenitor:p1. rhombomere" ~ "Neural progenitor:p1.rhombomere",
                            TRUE ~ CellType)) -> hs.mpr@meta.data
hs.mpr@meta.data %>%
  mutate(CellType=case_when(CellType%in%"Endothelium:liver sinusoidal endothelial cell" ~ "Endothelium:liver sinusoidal endothelium",
                            CellType%in%"Neuron:LMC" ~ "Neuron:lateral motor columns",
                            TRUE ~ CellType)) -> hs.mpr@meta.data
hs.mpr.merge2@meta.data %>%
  mutate(CellType=case_when(CellType%in%"NeuralProgenitor:p1. rhombomere" ~ "Neural progenitor:p1.rhombomere",
                            TRUE ~ CellType)) -> hs.mpr.merge2@meta.data
hs.mpr.merge2@meta.data %>%
  mutate(CellType=case_when(CellType%in%"Endothelium:liver sinusoidal endothelial cell" ~ "Endothelium:liver sinusoidal endothelium",
                            CellType%in%"Neuron:LMC" ~ "Neuron:lateral motor columns",
                            TRUE ~ CellType)) -> hs.mpr.merge2@meta.data

hs.mpr.merge2.sub1$Neuron@meta.data %>%
  mutate(CellType=case_when(CellType%in%"Neuron:LMC" ~ "Neuron:lateral motor columns",
                            TRUE ~ CellType),
         CellType.sub=case_when(CellType.sub%in%"LMC" ~ "lateral motor columns",
                                TRUE ~ CellType.sub)) -> hs.mpr.merge2.sub1$Neuron@meta.data

hs.mpr.merge2.sub1$Endothelium@meta.data %>%
  mutate(CellType=case_when(CellType%in%"Endothelium:liver sinusoidal endothelial cell" ~ "Endothelium:liver sinusoidal endothelium",
                            TRUE ~ CellType),
         CellType.sub=case_when(CellType.sub%in%"liver sinusoidal endothelial cell" ~ "liver sinusoidal endothelium",
                                TRUE ~ CellType.sub)) -> hs.mpr.merge2.sub1$Endothelium@meta.data



# =================
# Last part: Exit R ----
# =================
saveRDS(hs.embryo, "hs.mpr.new.final.rds")
save.image("Hs_Embryo_Multiome.RData")
rm(list = ls())

