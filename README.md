## Low Yield Single-Cell Differential Expression (Lyscde)

Calculate differential expression between sets of single cells. This method is normalization free and will automatically correct for library size differences. The input is raw gene count data.

## Installation

The software is developed for R. Install and import `lyscde` from GitHub.

```R
devtools::install_github("MaayanLab/lyscde/lyscde")
library("lyscde")
```
## Run differential gene expression

Pass two matrices or dataframes with gene counts.

```R
res = lyscde::diffexp(counts1, counts2)
```

## Plotting

Plot gene counts for select gene. Counts are not normalized to library size.

```R
lyscde::plotcounts("MARCKSL1", counts1, counts2, l1="Condition 1", l2="Condiction 2")
```
