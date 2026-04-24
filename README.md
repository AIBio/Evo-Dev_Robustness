# Evo-Dev_Robustness

Analysis Code about Developmental Robustness and Plasticity

## Evo-Chat

Evo-Chat extends CellChat by incorporating receiver-side transcription factor support to recalibrate pathway-specific communication probabilities. This repository provides a minimal example workflow for:
(1) inferring ligand–receptor interactions with CellChat,
(2) estimating downstream TF activity from transcriptome or chromatin data, and
(3) reweighting pathway activity with Evo-Chat.

### Install from local source package

```r
library(devtools)

# Build and install the package from source
devtools::document()
devtools::build()

install.packages(
  "/home/yhw/EvoChat_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

---

### Load example data

```r
library(Seurat)
library(SeuratData)
library(AfterChat)

# Toy dataset
data("pbmc3k")

# Example embryonic hindbrain Seurat object
test <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/Seurat_hs_embryo_hindbrain.rds")
```

---

### Step 1. Infer TF activity from transcriptome data with decoupleR

```r
# Toy example
act.pbmc <- Pipe.decoupleR(
  sr.obj = pbmc3k[, sample(1:ncol(pbmc3k), 100)],
  species = "human"
)

# Hindbrain example
act.hindbrain <- Pipe.decoupleR(
  sr.obj = test,
  species = "human"
)
```

---

### Step 2. Load external TF activity matrices (optional)

#### 2.1 SCENIC regulon activity

```r
auc.mtx <- read.table(
  "/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/SCENIC/hs_hindbrain/pySCENIC_auc_mtx.csv",
  sep = ",",
  header = TRUE
)

colnames(auc.mtx) <- gsub("\\.\\.\\.", "", colnames(auc.mtx))
colnames(auc.mtx) <- gsub("\\.", "-", colnames(auc.mtx))
auc.mtx$Cell <- test$CellType.sub
```

#### 2.2 Signac/chromVAR motif activity

```r
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    all_versions = FALSE
  )
)

hs.embryo <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hs.embryo.final.version.rds")
DefaultAssay(hs.embryo) <- "ATAC"

# Call peaks
peaks <- CallPeaks(
  hs.embryo,
  macs2.path = "/home/yhw/software/anaconda3/envs/dna.p37/bin/macs2"
)

# Remove non-standard chromosomes and blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(peaks, blacklist_hg38_unified, invert = TRUE)

# Quantify peak counts
macs2_counts <- FeatureMatrix(
  fragments = Fragments(hs.embryo),
  features = peaks,
  cells = colnames(hs.embryo)
)

annotation <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/RDS/hg38.gene.annotation.rds")

hs.embryo[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  annotation = annotation
)

DefaultAssay(hs.embryo) <- "peaks"
hs.embryo <- FindTopFeatures(hs.embryo, min.cutoff = 5)
hs.embryo <- RunTFIDF(hs.embryo)
hs.embryo <- RunSVD(hs.embryo)

# Add motifs and run chromVAR
hs.embryo <- AddMotifs(
  object = hs.embryo,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

hs.embryo <- RunChromVAR(
  object = hs.embryo,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(hs.embryo) <- "chromvar"

# Extract motif/TF activity
signac <- GetTfAct(
  sc.obj = hs.embryo,
  tf.format = "jaspar.id",
  filter.data = TRUE,
  filter.by = "CellType.sub",
  filter.terms = c(
    "floor plate", "roof plate",
    "p0.rhombomere", "p1.rhombomere", "p2.rhombomere",
    "pA1", "pA2", "pA3", "pB1"
  )
)
```

---

### Step 3. Run CellChat

```r
library(CellChat)

db <- CellChatDB.human

ct <- Pipe.CellChat(
  sc.obj = test,
  cell.db = db,
  group.by = "CellType.sub",
  split.by = NULL,
  filter.data = TRUE,
  filter.by = "CellType.sub",
  filter.terms = c(
    "floor plate", "roof plate",
    "p0.rhombomere", "p1.rhombomere", "p2.rhombomere",
    "pA1", "pA2", "pA3", "pB1"
  )
)
```

---

### Step 4. Estimate pathway-specific downstream support

Evo-Chat uses receiver-side TF activity to estimate pathway-specific support, which is then used to reweight CellChat interaction probabilities.

#### 4.1 Example: BMP/WNT support from SCENIC TF activity

```r
scale.f1 <- CalNormFactor(
  tf.data = auc.mtx[, -1],
  ave.method = "all",
  scale.data = TRUE,
  scale.method = c("range-01"),
  threshold = 0.5,
  genelist = c(
    "OLIG3", "MSX1", "MSX2", "SP5",
    "LMX1A", "SMAD1", "SMAD4", "SMAD5"
  ),
  group.info = as.character(test$CellType.sub)
)

bmp_wnt_support <- scale.f1$Ensembl.Mean[levels(ct$Sample@idents)]
```

#### 4.2 Example: HH support from SCENIC TF activity

```r
scale.f2 <- CalNormFactor(
  tf.data = auc.mtx[, -1],
  ave.method = "all",
  scale.data = TRUE,
  scale.method = c("range-01"),
  threshold = 0.5,
  genelist = c("DBX1", "DBX2", "GLI1", "GLI2"),
  group.info = as.character(test$CellType.sub)
)

hh_support <- scale.f2$Ensembl.Mean[levels(ct$Sample@idents)]
```

#### 4.3 Alternative: HH support from chromVAR motif activity

```r
scale.f2.chromvar <- CalNormFactor(
  tf.data = signac$TFact,
  ave.method = "all",
  scale.data = TRUE,
  scale.method = c("range-01"),
  threshold = 0.65,
  genelist = c("MA0675.1", "MA0069.1", "MA1990.1"),
  group.info = as.character(signac$GroupInfo)
)

hh_support_chromvar <- scale.f2.chromvar$Ensembl.Mean[levels(ct$Sample@idents)]
```

### 4.4 Build pathway-specific scaling factors

```r
scale.fs <- list(
  BMP = bmp_wnt_support,
  WNT = bmp_wnt_support,
  HH  = hh_support
)
```

---

### Step 5. Reweight CellChat interactions with Evo-Chat

```r
ct$Sample <- AdjustProPval(
  ct.obj = ct$Sample,
  type = "truncatedMean",
  trim = 0.1,
  scale.factor = scale.fs,
  scale.path = names(scale.fs)
)
```

This step recalibrates pathway-specific interaction probabilities and adjusted p values using receiver-side TF support.

---

### Step 6. Downstream analyses and visualization

```r
AfterChat::PathCentrality(
  ct.obj = ct$Sample,
  outdir = "~/test",
  file.prefix = "Sample"
)

AfterChat::PathClustering(
  ct.obj = ct$Sample,
  outdir = "~/test",
  file.prefix = "Sample"
)

AfterChat::PathInteracion(
  ct.obj = ct$Sample,
  outdir = "~/test",
  file.prefix = "Sample"
)

AfterChat::LRsContribution(
  ct.obj = ct$Sample,
  outdir = "~/test",
  file.prefix = "Sample"
)

AfterChat::LRsInteraction(
  ct.obj = ct$Sample,
  outdir = "~/test",
  file.prefix = "Sample",
  cell.source = c("floor plate", "roof plate"),
  cell.target = c(
    "p0.rhombomere", "p1.rhombomere", "p2.rhombomere",
    "pA1", "pA2", "pA3", "pB1"
  )
)
```

---

### Notes

- Replace all local file paths with your own file locations.
- `Pipe.decoupleR()`, `Pipe.CellChat()`, `CalNormFactor()`, and `AdjustProPval()` are core functions in the Evo-Chat workflow.
- Receiver-side TF activity can be estimated from transcriptome-based regulon inference (e.g., SCENIC/decoupleR) or chromatin-based motif activity (e.g., Signac/chromVAR).
- Pathway-specific reweighting is currently illustrated here for `BMP`, `WNT`, and `HH`.

---

### Citation

If you use Evo-Chat in your work, please cite the associated manuscript when available.
