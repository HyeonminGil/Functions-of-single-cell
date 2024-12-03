## 1. with gprofiler2

``` R
library(gprofiler2)

gene_name <- gconvert(
    rownames(dds),
    organism = "mmusculus",
    target = "ENTREZGENE",
    filter_na = TRUE)

matched_idx <- match(rownames(dds), gene_name$input)
length(matched_idx)

valid_rows <- !is.na(matched_idx)
table(valid_rows)

dds <- dds[valid_rows, ]
dds

rownames(dds) <- gene_name$target[matched_idx[valid_rows]]
head(rownames(dds))
```
