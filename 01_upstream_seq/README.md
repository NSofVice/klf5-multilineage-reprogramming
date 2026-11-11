# 01_upstream_seq/ — raw sequencing data processing

Bulk RNA-seq, CUT&Tag and ATAC-seq raw-read processing. Scripts are numbered in execution order.

**Genome build: mm10 (mouse).** All indices, effective genome sizes and annotation calls are set to
mm10, matching the manuscript.

| # | Script | Purpose |
|---|---|---|
| 0 | `alignment/00_STAR_RSEM_align.sh` | Bulk RNA-seq: STAR v2.7.3a (`--quantMode TranscriptomeSAM`) → RSEM v1.3.1, `parallel --jobs 4` |
| 1 | `fastq_download/01_fasterq_dump_and_fastp.sh` | SRA `prefetch` + `fastq-dump --split-3 --gzip`, then fastp |
| 2 | `fastp/02_fastp_dedup_batch.sh` | Batch fastp with PCR deduplication (`-D`), parallel via `xargs -P` |
| 3 | `alignment/03_ATAC_bowtie2.sh` | ATAC-seq: Bowtie2 → MAPQ filter → sambamba markdup → sort/index → MACS2 → RPKM BigWig |
| 4 | `alignment/04_CUTTag_bowtie2.sh` | Same chain; contains both an ATAC and a CUT&Tag Bowtie2 parameter block |
| 5 | `alignment/05_CUTTag_bowtie2_GSE124718.sh` | Single-end mode for the public GSE124718 dataset |
| 6 | `alignment/06_CUTTag_bowtie2_TE_mode_spikein.sh` | TE-mode: unique alignments (`-k 1`), E. coli spike-in rate, RPKM + CPM BigWig, resumable |
| 7 | `bam_bigwig/07_samtools_sort_index.sh` | Re-filter BAMs (`-q 20 -F 0x04`), sort and index |
| 8 | `bam_bigwig/08_bamCoverage_RPKM_bigwig.sh` | Batch RPKM-normalised BigWig generation |

## Typical run order

```
Bulk RNA-seq:  00 (STAR + RSEM)
CUT&Tag/ATAC:  01 (download + fastp) → 02 (fastp dedup) → 03 or 04 (align + peak) → 08 (BigWig)
```

## Key parameters (as in the manuscript)

```
Bowtie2        -x <mm10 index>   (ATAC: --end-to-end; CUT&Tag: --local; TE-mode: -k 1)
samtools view  -q 30
macs2 callpeak -f BAMPE -g mm --nomodel --keep-dup 1 --extsize 200 -q 1e-05
bamCoverage    --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --ignoreDuplicates
```

## Requirements

`STAR`, `RSEM`, `fastp`, `Bowtie2`, `samtools`, `sambamba`, `macs2`, `bedtools`, `deeptools`
(`bamCoverage`), `GNU parallel`, and GATK (optional, for fragment-size metrics). Job-scheduler
headers at the top of each script are for SGE/PBS (`qsub`).

See `../docs/CAVEATS.md` before running.