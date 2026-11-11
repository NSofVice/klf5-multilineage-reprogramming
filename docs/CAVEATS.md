# Known caveats and inconsistencies

The scripts in this repository are the working versions used during the study. They were written
for a specific HPC environment and were not refactored for general release. The points below are
documented so that readers are not misled.

---

## 0. Scripts reconciled with the manuscript (resolved)

The archived working scripts differed from the published Methods in three places. All three have
been **corrected in this repository so that the code matches the manuscript**:

| Item | Archived working script | This repository (= manuscript) | Files touched |
|---|---|---|---|
| Genome build | `hg38` Bowtie2 index path, `annotatePeaks.pl ... hg38` | **mm10** (`mm10` Bowtie2 index, HOMER `mm10`) | `03_ATAC_bowtie2.sh`, `04_CUTTag_bowtie2.sh`, `05_CUTTag_bowtie2_GSE124718.sh`, `06_CUTTag_bowtie2_TE_mode_spikein.sh` |
| Effective genome size (deepTools) | `2913022398` (hg38) | **`2652783500`** (mm10) | same four scripts |
| MACS2 duplicate handling | `--keep-dup all` | **`--keep-dup 1`** | `03_ATAC_bowtie2.sh`, `04_CUTTag_bowtie2.sh` |
| BigWig normalisation | `--normalizeUsing RPGC` | **`--normalizeUsing RPKM`** | `06_CUTTag_bowtie2_TE_mode_spikein.sh`, `08_bamCoverage_RPKM_bigwig.sh` |

The `hg38` index path in the archived originals was almost certainly a copied variable from a human
dataset: every downstream artefact (`mm10.blacklist.bed`, `mm10_chromInfo.txt`, MACS2 `-g mm`,
HOMER `mm10`, the `spe=mm10` variable) is mouse. Only the Bowtie2 `-x` argument and the deepTools
effective genome size were human, and those are what has been corrected.

The TE-mode script also emits a CPM BigWig in addition to the RPKM track; CPM was used for
visualisation only and is retained.

---

## 1. Hard-coded absolute paths

All scripts contain absolute paths from the original environment, e.g.:

- `/mnt/1/genome/mm10/...`, `/mnt/6/WXY/...`, `/mnt/11/WXY_classic_data/...`, `/mnt/12/wu_xinyu/...`
- `/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/...`
- `F:/Backup_20240703/秦组/Data processing/...`

These must be replaced before running.

## 2. Missing helper library

The two large R scripts both begin with:

```r
source("/home/data/t020406/R/WXY_R_function.R")     # scRNA-seq script
source("F:/Backup_20240703/秦组/Data processing/WXY_R_function/WXY_R_function.R")  # figure script
```

`WXY_R_function.R` is a custom helper library (peak annotation `Peak_Anno()`, GREAT-GO plotting
`Plot_GREAT_GO_results()`, HOMER parsing `Compare_Homer_known_motif_results()`, heat-map helpers
`Plot_Flodchange_Dot_Heatmap()`, `Venner_plot()`, and many others). **It is not included in this
repository** because the archived copy was not available at the time of release. Without it the R
scripts will not run. See `lib/README.md`.

## 3. Job-scheduler headers

`01_upstream_seq/` scripts start with SGE/PBS headers (`#PBS -l walltime=...`, `#PBS -q cv2`).
They are intended for `qsub` submission on the original cluster; remove the headers to run them as
plain shell scripts.

## 4. Mixed test variables in the scRNA-seq script

`04_scrna_seq/scRNAseq_Seurat_SCP_pipeline.R` is a working script that was extended over the
course of the project. Variable names from tutorial code (`pancreas_sub`, `monocle`, `seurat`) are
reused, and a few typos survive (e.g. `as.charater`). The analysis logic is intact; the naming is
not cleaned up.

## 5. Steps with no retained script

Two steps of the deposited pipeline have **no code** to deposit, as confirmed by the authors:

- **step_11** — scRNA-seq demultiplexing and UMI counting (PISA v0.2). Run as command-line calls on
  the server; no script was archived.
- **step_18** — IP-MS database search (MaxQuant v1.6.0.1). MaxQuant is a GUI / parameter-file tool;
  no script exists for this step.

The resulting data for both steps are available from the public repositories listed in the
manuscript (GEO GSE269422 for the scRNA-seq count matrix; PRIDE PXD086085 for the MaxQuant output).

## 6. Not runnable elsewhere without editing

Given points 1, 2 and 3 above, these scripts are provided as a faithful record of the analysis
performed, not as a turnkey pipeline. Editing is expected.