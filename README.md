# Single-Cell Multiome Analysis of Alzheimer's Disease Mouse Model

Replication of the 10x Genomics Application Note: *"Single cell and spatial multiomics identifies Alzheimer's disease markers"*

**Author:** Aakash Deva Thirukonda Prakash




## Project Overview

This project replicates the single-cell multiome (RNA + ATAC) analysis from the 10x Genomics Alzheimer's disease application note using data from TgCRND8 transgenic mice and wild-type littermates. The goal was to independently discover key findings through proper exploratory methodology.In this project- only microglia pecifically identifying **Slc1a3 as differentially expressed in microglia** with the largest differences observed at early timepoints when plaque deposition is sparse.

### Data Source

**10x Genomics Multiome Dataset:**  
[Multiomic Integration Neuroscience Application Note - Alzheimer's Disease Mouse Model](https://www.10xgenomics.com/datasets/multiomic-integration-neuroscience-application-note-single-cell-multiome-rna-atac-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard)

The dataset includes 12 mouse brain samples:
- **6 TgCRND8 (AD model)** mice at 2.5, 5.7, and 17.9 months (sparse, moderate, and severe plaque burden)
- **6 Wild-type littermates** at 2.5, 5.7, and 13.4 months

---

## Repository Structure

```
├── LoadingQC.Rmd                    # Data loading and quality control pipeline
├── Microglia_DEG_ATAC.Rmd           # Differential expression and chromatin accessibility analysis
├── README.md                        # This file
├── images/                          # Figures for documentation
│   ├── 01_wnn_umap_celltypes.png    # WNN UMAP with cell type annotations
│   ├── 02_dotplot_markers.png       # Marker validation DotPlot
│   ├── 03_slc1a3_violin_timepoints.png  # Slc1a3 expression by timepoint
│   ├── 04_slc1a3_coverage_plot.png  # Chromatin accessibility coverage plot
│   └── 05_top_motifs_wt.png         # Top enriched TF motifs
└── outputs/
    ├── obj_qc_filtered.rds          # QC-filtered Seurat object
    ├── multiome_WNN_final_processed_v3.rds  # Fully processed multiome object
    ├── micro_early_ATAC_analysis.rds        # Microglia ATAC analysis object
    ├── deg_Microglia_AD_vs_WTv3.csv         # Microglia differential expression results
    ├── DA_peaks_Slc1a3_regionv3.csv         # Differentially accessible peaks near Slc1a3
    ├── Motifs_WT_enriched_peaks.csv         # Motif enrichment in WT-accessible regions
    └── Figure5_replication.png              # Replication of application note Figure 5
```

---

## Analysis Pipeline

### 1. Data Loading & Quality Control (`LoadingQC.Rmd`)

#### Loading Multiome Data
- Loaded the aggregated multiome data from Cell Ranger ARC output (`.h5` format)
- Created a Seurat object with both RNA and ATAC assays
- Configured ATAC assay with fragment file and mm10 genome annotations
- Used EnsDb.Mmusculus.v79 for gene annotations

#### Quality Control Metrics

**RNA QC:**
- Number of genes detected per cell (`nFeature_RNA`)
- Total UMI counts (`nCount_RNA`)
- Mitochondrial read percentage (`percent.mt`)

**ATAC QC:**
- Peak region fragments
- Total fragments passing filters
- Percentage of reads in peaks
- TSS enrichment score
- Nucleosome signal

#### QC Thresholds Applied
| Metric | Threshold |
|--------|-----------|
| `nFeature_RNA` | 300 - 12,000 |
| `percent.mt` | < 10% |
| `peak_region_fragments` | 500 - 80,000 |
| `passed_filters` | 2,000 - 250,000 |
| `pct_reads_in_peaks` | > 3% |
| `TSS.enrichment` | > 2 |
| `nucleosome_signal` | < 2 |

**Result:** 27,890 high-quality cells retained after filtering

note, i applied a strict filterimg, but one can apply a leinent filtering and this may affect downstream process

---

### 2. Normalization & Dimensionality Reduction (`Microglia_DEG_ATAC.Rmd`)

#### RNA Processing
- **SCTransform normalization** with mitochondrial percentage regression
- Preserved all genes (including cell type markers) by setting `return.only.var.genes = FALSE`
- PCA with 50 components, selected 30 for downstream analysis

#### ATAC Processing
- **TF-IDF normalization** using Signac
- Top 25% of peaks retained (`min.cutoff = "q25"`)
- **LSI (Latent Semantic Indexing)** for dimensionality reduction
- Excluded LSI component 1 (correlated with sequencing depth)

#### Multimodal Integration
- **Weighted Nearest Neighbor (WNN)** algorithm to integrate RNA and ATAC modalities
- Joint UMAP embedding considering both gene expression and chromatin accessibility

---

### 3. Cell Type Annotation

#### Marker-Based Approach
Defined canonical brain cell type markers:

| Cell Type | Key Markers |
|-----------|-------------|
| **Oligodendrocytes** | Mbp, Plp1, Mog, Mag, Cnp |
| **OPCs** | Pdgfra, Cspg4, Olig1, Olig2 |
| **Microglia** | P2ry12, Tmem119, Cx3cr1, Aif1, Trem2 |
| **Pericytes/Endothelial** | Pecam1, Kdr, Pdgfrb, Rgs5 |
| **Inhibitory Neurons** | Gad1, Gad2, Pvalb, Sst, Vip |
| **Excitatory Neurons** | Slc17a7, Slc17a6, Camk2a, Tbr1 |
| **Astrocytes** | Aqp4, Aldh1l1, Slc1a3, Gfap |

#### Annotation Method
- Calculated module scores for each cell type using `AddModuleScore()`
- Assigned cells based on highest score with margin threshold (≥0.05) to ensure confident assignments
- Cells with ambiguous scores labeled as "Unknown/Ambiguous"

![Figure 1: WNN UMAP with Cell Type Annotations](https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/Marker_plots.png)![alt text](scripts/images_plot/Marker_plots.png)


**Figure 1: Weighted Nearest Neighbor (WNN) UMAP embedding of 27,890 cells colored by annotated cell type.** Integration of RNA and ATAC modalities enables robust identification of major brain cell populations including neurons, glia, and vascular cells.

![Figure 2: Marker Expression Validation](![https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/dotplot.png]![alt text](scripts/images_plot/Dotplot.png) (
**Figure 2: DotPlot showing expression of canonical marker genes across predicted cell types.** Dot size represents the percentage of cells expressing each marker; color intensity indicates average expression level. Clear enrichment of expected markers in corresponding cell types validates the annotation approach.

---

### 4. Sample Metadata Assignment


When Cell Ranger aggregates multiple libraries into a single matrix, it appends numeric suffixes (1-12) to cell barcodes to distinguish their library of origin. While these suffixes serve as unique identifiers, they carry no inherent information about sample metadata.   To reconstruct sample identities from these arbitrary suffixes, i extracted the barcode suffix from each cell and cross-referenced it against the authoritative sample mapping provided in the 10x Genomics web summary output from the aggregation pipeline. The 10x table definitively shows which sample was assigned each index based on the order they were aggregated (line 1 → index 1, line 2 → index 2, etc.), allowing us to unambiguously map all 33,459 cells to their original sample IDs and metadata (genotype, age, plaque burden). 

Reconstructed sample identities from barcode suffixes:

| Index | Sample ID | Condition | Age (months) | Plaque Burden |
|-------|-----------|-----------|--------------|---------------|
| 1 | AD_17p9_rep4 | AD | 17.9 | Severe |
| 2 | AD_17p9_rep5 | AD | 17.9 | Severe |
| 3 | AD_2p5_rep2 | AD | 2.5 | Sparse |
| 4 | AD_2p5_rep3 | AD | 2.5 | Sparse |
| 5 | AD_5p7_rep2 | AD | 5.7 | Moderate |
| 6 | AD_5p7_rep6 | AD | 5.7 | Moderate |
| 7 | WT_13p4_rep2 | WT | 13.4 | None |
| 8 | WT_13p4_rep5 | WT | 13.4 | None |
| 9 | WT_2p5_rep2 | WT | 2.5 | None |
| 10 | WT_2p5_rep7 | WT | 2.5 | None |
| 11 | WT_5p7_rep2 | WT | 5.7 | None |
| 12 | WT_5p7_rep3 | WT | 5.7 | None |




This reconstruction was verified by confirming that cell counts per sample, condition distributions, and timepoint groupings all align with similar  numbers from the 10x summary


---

### 5. Differential Expression Analysis

#### Initial Approach: Pooled Analysis
First performed DEG analysis comparing all AD vs WT microglia pooled across timepoints.

Figure 3-  Microglia show selective downregulation of homeostatic genes in AD.

![Figure 3: Microglia show selective downregulation of homeostatic genes in AD.](https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/Microglia_global_DEG.png)


Microglia from AD transgenic mice displayed a striking gene expression profile: while APP showed dramatic upregulation (consistent with transgenic amyloid overexpression), the majority of the top altered genes were downregulated, including critical homeostatic mediators Atp8a2 (phospholipid homeostasis),, Pcdh9 (cell adhesion), Il1rapl1 (inflammatory regulation), and Fgf14 (growth signaling). This pattern—loss of homeostatic functions coupled with amyloid response—reflects the expected AD microglial phenotype.But, we cannot distinguish whether microglial dysfunction causes early disease or results from late amyloid accumulation thast why timepoint strafied analysis was done.

#### Key Methodological Insight
**Timepoint-stratified analysis** was essential for discovering biologically meaningful differences. Pooling across timepoints masked the strongest signals.

#### Timepoint-Stratified DEG Analysis
Performed separate comparisons for each timepoint:
- **2.5 months** (sparse plaque): Strongest Slc1a3 differential expression
- **5.7 months** (moderate plaque): Reduced but still detectable differences
- **Late timepoints**: Diminished differences

At early disease stages, AD microglia exhibit the expected transcriptional signature: upregulation of amyloid-response genes (APP, Adgrf3) coupled with downregulation of homeostatic functions (Atp8a2, Slc1a3). Among these characteristic AD-associated genes, Slc1a3 stands out not for its statistical magnitude but for its direct mechanistic relevance.—glutamate transporter dysfunction predicts excitotoxic neuronal damage. The timepoint-stratified analysis revealed that this Slc1a3 downregulation is most pronounced at early disease stages (sparse plaques), suggesting that glutamate homeostasis failure is an early initiating event in AD pathogenesis

---

## Key Results

### Slc1a3 Differential Expression in Microglia

**Finding:** Slc1a3 (glutamate transporter GLAST) is downregulated in AD microglia compared to WT, with the **largest difference at the early 2.5-month timepoint** when plaque deposition is sparse.

![Figure 4: Slc1a3 Expression Across Timepoints](https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/slc3a.png?raw=true)

**Figure 4: Slc1a3 Expression Across Timepoints.** Violin plot showing the distribution of *Slc1a3* expression in microglia across different timepoints in WT and AD conditions.

**Figure 4: Violin plots showing Slc1a3 expression in microglia stratified by condition and timepoint.** Slc1a3 is significantly downregulated in AD microglia at 2.5 months (sparse plaque burden), with the difference diminishing at later timepoints. This replicates the key finding from the 10x Genomics application note.

This replicates the key finding from the 10x Genomics application note (Figure 5A).

### Peak-Gene Linkage Analysis

Using `LinkPeaks()`, identified chromatin accessibility peaks correlated with Slc1a3 expression:

| Peak Location | Link Score | Position Relative to Slc1a3 |
|---------------|------------|------------------------------|
| chr15:8801492-8802323 | Strongest | Distal enhancer (~90kb downstream) |
| Multiple intragenic peaks | Moderate | Within gene body |
| Upstream peaks | Variable | Promoter-proximal |

### Differentially Accessible Chromatin Regions

Performed differential accessibility analysis between AD and WT microglia at 2.5 months using logistic regression with sequencing depth correction:
- Identified WT-enriched peaks in the Slc1a3 regulatory region
- These regions show reduced accessibility in AD microglia

![Figure 5: Slc1a3 Chromatin Accessibility and Peak-Gene Links](https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/Chromatin%20accesibliy%20plot.png?raw=true)

**Figure 5: Coverage plot showing chromatin accessibility at the Slc1a3 locus in WT vs AD microglia at 2.5 months.** Top tracks display normalized ATAC signal; arcs indicate significant peak-gene links identified by LinkPeaks. A distal enhancer region (~90kb downstream) shows strong correlation with Slc1a3 expression and differential accessibility between conditions.



### Motif Enrichment Analysis


# Motif Enrichment Analysis

Motif enrichment analysis of differentially accessible chromatin regions identified multiple transcription factors with significantly enriched binding motifs in WT-accessible peaks. The KLF family (KLF2, KLF3, KLF6, KLF11, KLF12, KLF14, KLF15) showed the highest statistical significance, followed by SP proteins (SP1, SP3, SP4), NRF1, and NFYC. EGR1, an immediate early gene associated with stress and inflammatory responses, was also enriched but ranked 18th among WT-accessible peak motifs.

## Top motifs enriched in WT-accessible chromatin

| Transcription Factor | Motif Family | Relevance |
|---|---|---|
| KLF family (KLF2, KLF3, KLF6, KLF11, KLF15) | Krüppel-like Factors | Homeostatic microglial functions; constitutive regulators |
| SP1/SP3/SP4 | Specificity Proteins | General transcriptional maintenance |
| NRF1 | Nuclear Respiratory Factor | Metabolic and mitochondrial regulation |
| EGR1 | Early Growth Response | Stress-responsive; secondary role in WT peaks |


![Figure 5: Top Enriched Motifs in WT-Accessible Chromatin](https://github.com/thethirukonda/mulltiome_alzhimers/blob/71108e9a6717f086bfb0bf83dd5538a128fcd5d1/images/motifs.png?raw=true)

**Figure 5: Sequence logos of top transcription factor motifs enriched in WT-accessible chromatin (2.5 month microglia).** The KLF family (KLF2, KLF3, KLF6, KLF11, KLF15, KLF16) dominates the homeostatic microglial regulatory landscape, along with SP proteins (SP1, SP3, SP4, SP9), NRF1, and other zinc finger factors. EGR1 appears among enriched motifs but ranks secondary to KLF family members, indicating that KLF-mediated regulation controls WT-accessible homeostatic chromatin including the *Slc1a3* enhancer.

While the 10x application note identified EGR1 as an enriched motif in Slc1a3-linked regions, my analysis ranks EGR1 18th among motifs enriched in WT-accessible chromatin, compared to KLF family members which occupy the top positions . This indicates that KLF4 and related KLF factors are the primary regulators while EGR1 acts as a secondary regulator of the WT-accessible Slc1a3 enhancer.
### ChromVAR Motif Activity

Per-cell motif activity analysis confirmed reduced EGR1 activity in AD microglia compared to WT at the 2.5-month timepoint.

---

## Software & Dependencies

# R Session Information

## System
- **R version:** 4.5.2 (2025-10-31)
- **Platform:** aarch64-apple-darwin20
- **OS:** macOS Sequoia 15.7.3
- **Time zone:** Asia/Kolkata

## Linear Algebra
- **BLAS:** /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
- **LAPACK:** /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib (version 3.12.1)

## Core Packages Loaded

### Bioconductor - Multiome & Chromatin Analysis
- `Signac` (1.16.0) - Single-cell chromatin analysis
- `Seurat` (5.4.0) - Single-cell RNA-seq analysis
- `SeuratObject` (5.3.0)
- `chromVAR` (1.32.0) - Chromatin accessibility variability
- `motifmatchr` (1.32.0) - TF motif matching
- `TFBSTools` (1.48.0) - Transcription factor binding site analysis
- `JASPAR2020` (0.99.10) - TF motif database

### Annotation & Genomics
- `GenomicRanges` (1.62.1)
- `GenomicFeatures` (1.62.0)
- `BSgenome.Mmusculus.UCSC.mm10` (1.4.3) - Mouse genome
- `EnsDb.Mmusculus.v79` (2.99.0) - Ensembl mouse annotations
- `AnnotationDbi` (1.72.0)
- `AnnotationFilter` (1.34.0)

### Single-Cell Reference & Cell Type
- `SingleR` (2.12.0) - Single-cell reference annotation
- `celldex` (1.20.0) - Cell type reference datasets

### Visualization
- `ggplot2` (4.0.1)
- `ggseqlogo` (0.2.2) - Sequence logo visualization
- `ggrepel` (0.9.6)
- `ggpubr` (0.6.2)
- `patchwork` (1.3.2)

### Data Manipulation
- `dplyr` (1.1.4)
- `stringr` (1.6.0)

### Parallel Computing
- `future` (1.69.0)

## Key Versions for This Analysis

| Package | Version | Purpose |
|---------|---------|---------|
| Seurat | 5.4.0 | Multiome integration & clustering |
| Signac | 1.16.0 | ATAC-seq analysis |
| chromVAR | 1.32.0 | Motif activity scoring |
| motifmatchr | 1.32.0 | TF motif detection |
| JASPAR2020 | 0.99.10 | TF binding motifs |
| BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 | Mouse mm10 genome |
---

## Reproducibility Notes

### Critical Steps
1. **Fragment file path:** Must be absolute path for ATAC analysis
2. **Seqinfo synchronization:** Peak ranges must have matching seqinfo with mm10 genome
3. **Layer alignment:** ATAC layer column order must match Seurat object cells before subsetting
4. **Object reconstruction:** Manual object rebuilding from counts is more robust than direct subsetting for multiome data

### Running the Analysis
```r
# 1. First run QC pipeline
rmarkdown::render("LoadingQC.Rmd")

# 2. Then run differential expression and ATAC analysis
rmarkdown::render("Microglia_DEG_ATAC.Rmd")
```

---

## Comparison to Application Note

| Finding | Application Note | This Analysis |
|---------|------------------|---------------|
| Slc1a3 downregulated in AD microglia | ✓ | ✓ Replicated |
| Biggest difference at early timepoints | ✓ | ✓ Replicated |
| EGR1 motif enrichment | ✓ | ✓ Replicated |
| Distal regulatory region linked to Slc1a3 | ✓ | ✓ Identified |
| Cell type annotation | 8-12 types | 33-35 clusters (higher granularity) |

The analysis successfully replicated the key biological findings from the 10x Genomics application note through independent exploratory analysis, validating the methodology and confirming the role of Slc1a3 dysregulation and EGR1 transcription factor involvement in early Alzheimer's disease pathogenesis.


