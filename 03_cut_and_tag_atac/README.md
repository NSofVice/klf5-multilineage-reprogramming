# 03_cut_and_tag_atac/

The CUT&Tag and ATAC-seq analysis is split across two directories rather than duplicated here:

| Analysis stage | Location |
|---|---|
| Raw-read QC, alignment, deduplication, peak calling, BigWig | `../01_upstream_seq/` |
| Peak annotation, overlap analysis, motif comparison, figures | `../02_bulk_rna_seq/Figure6A_bulk_CUTTag_ATAC_analysis.R` |

The figure-generation script (`Figure6A_...R`) handles bulk RNA-seq, CUT&Tag and ATAC-seq
together because the manuscript figures present them side by side; it is therefore filed under
`02_bulk_rna_seq/`. This directory exists as a pointer only.
