# lib/ — custom helper library

The R scripts in this repository depend on a custom helper library:

```
WXY_R_function.R
```

Referenced in the code as:

```r
source("/home/data/t020406/R/WXY_R_function.R")
source("F:/Backup_20240703/秦组/Data processing/WXY_R_function/WXY_R_function.R")
```

## Status: NOT included

The library file was not available in the archived code directories at the time this repository
was assembled, so it could not be deposited here.

## What it contains (functions called by the released scripts)

| Function | Used for |
|---|---|
| `Peak_Anno()` | ChIPseeker-based peak annotation, ±20 kb TSS window |
| `Plot_GREAT_GO_results()` | Plot GREAT GO enrichment output |
| `Compare_Homer_known_motif_results()` | Parse and compare HOMER `knownResults.txt` |
| `Plot_Flodchange_Dot_Heatmap()` | Fold-change dot / heat-map figures |
| `Btable_insert_into_Atable()` | Table merging helper |
| `Venner_plot()` | Venn diagram figures |
| `Compare_Cluster_GO_results()` | Per-cluster GO comparison |

## Reproducing the analysis without the library

Every function above is a thin wrapper around standard Bioconductor tools. To reproduce the
analysis, substitute the corresponding public API calls:

- `Peak_Anno()` → `ChIPseeker::annotatePeak()` + `ChIPpeakAnno`
- motif parsing → `read.delim("knownResults.txt")` from HOMER
- Venn diagrams → `VennDiagram` / `ggVennDiagram`
- GO enrichment → `clusterProfiler::enrichGO()` / `compareCluster()`

If you need the original library for exact reproduction, please contact the corresponding author.
