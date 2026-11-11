# Code map: manuscript pipeline steps → scripts in this repository

This document maps each step of the analysis pipeline (as deposited in the GenPilot archive form,
sheet *3文章分析流程*, 19 steps) to the corresponding script file in this repository.

**Legend**
- ✅ code available in this repository
- ❌ no script exists for this step (confirmed by the authors)

> The scripts as archived differed from the published Methods in genome build, MACS2 `--keep-dup`
> and BigWig normalisation. **All three have been corrected so that the code matches the
> manuscript** (mm10 / `--keep-dup 1` / RPKM). Details in `CAVEATS.md §0`.

---

## Overview

| Step | Pipeline stage | Software (manuscript) | Script | Status |
|---|---|---|---|---|
| step_01 | Bulk RNA-seq alignment & quantification | STAR v2.7.3a; RSEM v1.3.1 | `01_upstream_seq/alignment/00_STAR_RSEM_align.sh` | ✅ |
| step_02 | Bulk RNA-seq differential expression | DESeq2 v1.38.3 | `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` | ✅ |
| step_03 | Bulk RNA-seq GO enrichment | clusterProfiler v4.7.1.003 | `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` | ✅ |
| step_04 | Bulk RNA-seq co-expression network | WGCNA v1.72-1 | `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` | ✅ |
| step_05 | CUT&Tag / ATAC-seq QC pre-processing | fastp v0.12.4 | `01_upstream_seq/fastp/02_fastp_dedup_batch.sh`, `01_upstream_seq/fastq_download/01_fasterq_dump_and_fastp.sh` | ✅ |
| step_06 | CUT&Tag / ATAC-seq alignment & filtering | Bowtie2 v2.4.4; Samtools v1.9; sambamba v0.8.0 | `01_upstream_seq/alignment/03_ATAC_bowtie2.sh`, `04_CUTTag_bowtie2.sh`, `05_CUTTag_bowtie2_GSE124718.sh`, `06_CUTTag_bowtie2_TE_mode_spikein.sh` | ✅ |
| step_07 | CUT&Tag / ATAC-seq peak calling | MACS2 v2.2.7.1 | `01_upstream_seq/alignment/03_ATAC_bowtie2.sh`, `04_CUTTag_bowtie2.sh` (section *Calling peak*) | ✅ |
| step_08 | Peak annotation & overlap | ChIPseeker v1.34.1; ChIPpeakAnno v3.32.0 | `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` (`Peak_Anno()`) | ✅ |
| step_09 | Motif discovery | HOMER v4.11.1 | `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` (`Compare_Homer_known_motif_results()`) | ✅ |
| step_10 | BigWig generation | deepTools v3.4.3 | `01_upstream_seq/bam_bigwig/08_bamCoverage_RPKM_bigwig.sh`, bamCoverage blocks in `03` / `04` / `06` | ✅ |
| step_11 | scRNA-seq demultiplex & UMI counting | PISA v0.2; STAR v2.7.3a; sambamba v0.7.0 | **无 (no script)** | ❌ |
| step_12 | scRNA-seq QC, doublet / ambient RNA removal | DoubletFinder v2.0.3; SoupX v1.6.1; Seurat v5.0.1 | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` | ✅ |
| step_13 | scRNA-seq normalization, dim. reduction, clustering | Seurat v5.0.1 | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` | ✅ |
| step_14 | scRNA-seq signature scoring | AUCell v1.20.2 | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` | ✅ |
| step_15 | scRNA-seq integration (Harmony) | Harmony v1.2.0 | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` | ✅ |
| step_16 | scRNA-seq differential expression | Seurat `FindMarkers` v5.0.1 | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` | ✅ |
| step_17 | scRNA-seq trajectory inference (PAGA) | scanpy v1.9.3 (Python v3.9.12) | `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` (SCP wrapper) | ✅ |
| step_18 | IP-MS database search & quantification | MaxQuant v1.6.0.1 | **无 (no script)** | ❌ |
| step_19 | IP-MS differential interaction analysis | DEP v1.20.0; limma v3.54.1 | `05_ip_ms/IP_MS_Klf245_venn_analysis.R` | ✅ |

**Summary: 17 of the 19 steps have code in this repository (all ✅, reconciled with the manuscript).
2 steps have no script — step_11 (PISA demultiplexing / UMI counting) and step_18 (MaxQuant database
search), both confirmed by the authors as not existing as scripts.** The resulting data for those
two steps are available from GEO GSE269422 and PRIDE PXD086085 respectively.

---

## Detailed mapping

### step_01 — Bulk RNA-seq alignment & quantification (✅)
**File:** `01_upstream_seq/alignment/00_STAR_RSEM_align.sh`
STAR v2.7.3a against the mm10 RSEM-ready STAR index (`mm10_gencode_M25_rsem_star`), with
`--quantMode TranscriptomeSAM`, followed by `rsem-calculate-expression --bam` for gene/transcript
quantification. Runs 4 samples in parallel (`GNU parallel --jobs 4`); each sequencing run of a
sample is merged beforehand.

STAR parameters (as in the script): `--outFilterType BySJout`, `--outFilterMultimapNmax 20`,
`--outFilterMismatchNmax 999`, `--outFilterMismatchNoverLmax 0.04`, `--alignIntronMin 20`,
`--alignIntronMax 1000000`, `--alignMatesGapMax 1000000`, `--alignSJoverhangMin 8`,
`--alignSJDBoverhangMin 1`, `--sjdbScore 1`, `--genomeLoad LoadAndKeep`.

Quantification: RSEM gene counts and TPM.

### step_02 — Bulk RNA-seq differential expression (✅)
**File:** `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R`
DESeq2 v1.38.3 on R v4.0.2; thresholds p.adj < 0.1, fold change > 2, TPM > 2 in at least one sample.

### step_03 — Bulk RNA-seq GO enrichment (✅)
**File:** `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R`
Uses `clusterProfiler` (library call at the top of the script) plus the custom helper
`Compare_Cluster_GO_results()` / `Plot_GREAT_GO_results()` from `WXY_R_function.R`.

### step_04 — Bulk RNA-seq co-expression network (✅)
**File:** `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R`
WGCNA v1.72-1; `blockwiseModules(minModuleSize=200, power=18, networkType="signed", TOMType="signed")`
after filtering TPM < 10 and outlier genes. Public datasets GSE123691 / GSE110950 / GSE86586.

### step_05 — CUT&Tag / ATAC-seq QC (✅)
**Files:**
- `01_upstream_seq/fastp/02_fastp_dedup_batch.sh` — batch fastp with PCR-deduplication
  (`-D`), parallelised with `xargs -P`, one JSON + HTML report per sample.
- `01_upstream_seq/fastq_download/01_fasterq_dump_and_fastp.sh` — SRA download
  (`prefetch` + `fastq-dump --split-3 --gzip`) followed by fastp.

### step_06 — CUT&Tag / ATAC-seq alignment & filtering (✅)
**Files:**
- `01_upstream_seq/alignment/03_ATAC_bowtie2.sh` — ATAC-seq branch
- `01_upstream_seq/alignment/04_CUTTag_bowtie2.sh` — contains **both** an `###ATAC###` and a
  `###CUTTag###` Bowtie2 block
- `01_upstream_seq/alignment/05_CUTTag_bowtie2_GSE124718.sh` — single-end mode for the public
  GSE124718 dataset
- `01_upstream_seq/alignment/06_CUTTag_bowtie2_TE_mode_spikein.sh` — TE-mode branch with
  E. coli spike-in normalisation, unique alignment (`-k 1`), **no** MAPQ or duplicate filtering,
  resumable (tracks completed samples)

Common post-alignment chain: `samtools view -bS -q 30` → `sambamba markdup` → `samtools sort` →
`samtools index`. Also runs `gatk CollectInsertSizeMetrics` for fragment-size QC and
`bedtools intersect` against the mm10 blacklist (`mm10.blacklist.bed`).

✅ Reconciled: the Bowtie2 index path is now `mm10` (the archived original pointed at an `hg38`
directory, which was a copied variable — every other artefact in those scripts was mouse).

### step_07 — Peak calling (✅)
**Files:** `03_ATAC_bowtie2.sh`, `04_CUTTag_bowtie2.sh` — section *Calling peak*
```
macs2 callpeak -f BAMPE -g mm -B --SPMR --nomodel --keep-dup 1 --extsize 200 -q 1e-05
```
Followed by `macs2 bdgcmp -m logLR -p 0.01` and `bdg2bw.sh` to convert signal to BigWig.

✅ Reconciled: `--keep-dup 1` (the archived original used `--keep-dup all`).

### step_08 — Peak annotation & overlap (✅)
**File:** `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R`
Custom `Peak_Anno(peak, up = 20000, down = 20000)` — annotates peaks via ChIPseeker and retains
only sites within ±20 kb of a TSS. Overlaps between Klf2 / Klf4 / Klf5 peak sets are computed and
plotted with `Venner_plot()`.

### step_09 — Motif discovery (✅)
**File:** `02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R`
`Compare_Homer_known_motif_results()` — parses HOMER `knownResults.txt` output
(`findMotifsGenome.pl`, run against **mm10**). HOMER itself is a command-line tool; only the
downstream comparison / plotting code is in this repository.

### step_10 — BigWig generation (✅)
**Files:** `01_upstream_seq/bam_bigwig/08_bamCoverage_RPKM_bigwig.sh`,
plus `bamCoverage` blocks in `03_ATAC_bowtie2.sh` / `04_CUTTag_bowtie2.sh` /
`06_CUTTag_bowtie2_TE_mode_spikein.sh`
```
bamCoverage --normalizeUsing RPKM --effectiveGenomeSize 2652783500 \
            --minMappingQuality 30 --ignoreDuplicates
```
✅ Reconciled: RPKM normalisation and the mm10 effective genome size (2652783500). The TE-mode
script additionally emits a CPM track, used for visualisation only.

### step_11 — scRNA-seq demultiplex & UMI counting (❌ 无)
PISA v0.2 → STAR v2.7.3a → sambamba v0.7.0, run as command-line calls on the server.
**No script exists for this step.** The resulting count matrix is available from GEO GSE269422.

### step_12 — scRNA-seq QC, doublets, ambient RNA (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
`DoubletFinder`, `SoupX`, Seurat QC. Filters: > 5 % mitochondrial reads removed; cells with
< 2000 or > 12000 genes removed; genes expressed in < 3 cells removed.

### step_13 — Normalization, dimensionality reduction, clustering (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
`NormalizeData` / `ScaleData`; top 2000 HVGs (`FindVariableFeatures`); `RunPCA`;
`ScoreJackStraw`; `FindNeighbors` / `FindClusters`; `RunTSNE` / `RunUMAP`.

### step_14 — Signature scoring (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
`AUCell_buildRankings()` + `AUCell_calcAUC()` on nPSC / TSC / PGC / BMP / Wnt gene sets
(derived from the WGCNA analysis and public datasets).

### step_15 — Integration with public in vivo embryo data (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
`RunHarmony` via the SCP wrapper (`method = HarmonyIntegration`); UMAP + clustering on the merged
object. Public datasets GSE123046 / GSE148032.

### step_16 — Cluster differential expression (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
`FindAllMarkers` / `FindMarkers` (`logfc.threshold = 1`, `only.pos = TRUE`, `min.pct = 0.25`),
followed by custom GO comparison helpers.

### step_17 — Trajectory inference (✅)
**File:** `04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R`
scanpy-based PAGA via SCP (`tl.pca`, `pp.neighbors`, `tl.diffmap`, `tl.leiden`, `tl.paga`);
non-reprogramming branches annotated with `SCP::RunKNNPredict(bulk_ref = "ref_scMCA",
filter_lowfreq = 500)`. Also contains Slingshot, Monocle and `RunDynamicFeatures` blocks used in
the study.

### step_18 — IP-MS database search & quantification (❌ 无)
MaxQuant v1.6.0.1 (label-free quantification + match-between-runs) against the Mouse Fasta
database; Easy-nLC 1200 + Fusion Lumos. **MaxQuant is a GUI / parameter-file tool, so no script
exists for this step.** Output matrices are available from PRIDE **PXD086085**.

### step_19 — IP-MS differential interaction analysis (✅)
**File:** `05_ip_ms/IP_MS_Klf245_venn_analysis.R`
Reads the MaxQuant iBAQ result tables (KLF5 vs IgG, 3 replicates each; adj.P.val < 0.1,
log2FC > 3 after removing contaminants and reverse hits) and generates the Klf2 / Klf4 / Klf5
overlap Venn diagrams.

---

## Additional analyses present in the code but not listed as separate pipeline steps

The scRNA-seq script also contains blocks for **Slingshot** trajectory inference, **Monocle**
pseudotime, **SCENIC** regulon activity, and **`RunDynamicFeatures`** dynamic-expression analysis.
These support the trajectory / regulon figures in the manuscript.