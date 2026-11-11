# Klf5 enables multilineage reprogramming from mouse epiblast stem cells
### Analysis code accompanying the manuscript (Communications Biology, COMMSBIO-26-0442C)

This repository contains the analysis scripts used in the study. It is provided for
transparency and reproducibility of the published work.

---

## Contents

| Directory | Description | Data type |
|---|---|---|
| `01_upstream_seq/` | Raw-read processing: bulk RNA-seq STAR+RSEM alignment/quantification; FASTQ download, fastp QC, Bowtie2 alignment, deduplication, sorting, RPKM BigWig generation | Bulk RNA-seq, CUT&Tag, ATAC-seq |
| `02_bulk_rna_seq/` | Differential expression, GO enrichment, WGCNA, peak annotation, motif analysis, figure generation | Bulk RNA-seq + CUT&Tag + ATAC-seq |
| `03_cut_and_tag_atac/` | Pointer to the CUT&Tag / ATAC-seq scripts (see note below) | CUT&Tag, ATAC-seq |
| `04_scrna_seq/` | Seurat/SCP pipeline: QC, clustering, AUCell scoring, Harmony integration, DEG, trajectory | scRNA-seq |
| `05_ip_ms/` | IP-MS interaction analysis (Klf2/Klf4/Klf5 overlap) | IP-MS proteomics |
| `docs/` | `CODE_MAP.md` (step-to-script mapping), `CAVEATS.md` (known inconsistencies) | — |
| `lib/` | Third-party / custom helper library required by the R scripts | — |

## Data availability

| Data type | Repository | Accession |
|---|---|---|
| Bulk RNA-seq, ATAC-seq, CUT&Tag | NCBI GEO | **GSE281941** |
| scRNA-seq | NCBI GEO | **GSE269422** |
| IP-MS proteomics | PRIDE / ProteomeXchange | **PXD086085** |

Raw sequencing data are **not** deposited in this repository; they are available from the
accessions above.

## Reference genome

**mm10 (mouse).** Bowtie2 indices, deepTools `--effectiveGenomeSize 2652783500`, HOMER and MACS2
(`-g mm`) are all set to mm10, matching the manuscript.

## Software

See `docs/CODE_MAP.md` for the software and version used at each step. Summary:

- **Alignment / processing**: fastp v0.12.4, Bowtie2 v2.4.4, Samtools v1.9, sambamba v0.8.0, STAR v2.7.3a, RSEM v1.3.1
- **Peak calling / epigenomics**: MACS2 v2.2.7.1, HOMER v4.11.1, deepTools v3.4.3, ChIPseeker v1.34.1, ChIPpeakAnno v3.32.0, bedtools, GATK
- **R / Bioconductor**: R v4.0.2–v4.3.2, DESeq2 v1.38.3, clusterProfiler v4.7.1.003, WGCNA v1.72-1, Seurat v5.0.1, DoubletFinder v2.0.3, SoupX v1.6.1, AUCell v1.20.2, Harmony v1.2.0
- **Python**: Python v3.9.12, scanpy v1.9.3
- **Proteomics**: MaxQuant v1.6.0.1, DEP v1.20.0, limma v3.54.1

## Important notes before running

1. **Absolute paths are hard-coded.** The scripts were written for a specific HPC/cluster
   environment (SGE/PBS `#PBS` headers, `/mnt/...` paths, and `/home/data/t020406/...` paths).
   They must be edited before they can run elsewhere.
2. **The R scripts depend on a custom helper library** (`WXY_R_function.R`) that is **not**
   included in this repository — see `lib/README.md`.
3. **See `docs/CAVEATS.md`.** The archived working scripts differed from the published Methods in
   genome build, MACS2 `--keep-dup` and BigWig normalisation; these have been corrected here so
   that the code matches the manuscript (mm10 / `--keep-dup 1` / RPKM).

## Coverage

Of the 19 pipeline steps deposited in the GenPilot archive (*3文章分析流程*), **17 have code in this
repository**; step_11 (PISA demultiplexing) and step_18 (MaxQuant) have no script to deposit and are
marked as such in `docs/CODE_MAP.md`.

## Citation

If you use this code, please cite:

> Wu X, Zhou H, Gao Y, et al. Klf5 enables multilineage reprogramming from mouse epiblast stem
> cells. *Communications Biology* (in press).

## License

MIT — see `LICENSE`.
