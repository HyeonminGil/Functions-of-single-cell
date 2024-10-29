```R
library(Matrix)
matrix_gex <- Matrix(gex[["RNA"]]$counts, sparse = T)
writeMM(obj = matrix_gex, file=glue("{dirdata}seurat_GEX_raw_count_matrix.mtx"))
```