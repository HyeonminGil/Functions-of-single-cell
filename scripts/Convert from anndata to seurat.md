### In Python
```python
## Load packages

## Directory path
dirdata =
namedata = IATEAT_APC_for_R.loom
nameUMAP = 

## Copy anndata
adata_tmp = adata_IE_APC.copy()

## Covert to csr matrix
# raw count matrix
adata_tmp.X = adata_tmp.layers["counts"]
adata_tmp.X = scipy.sparse.csr_matrix(adata_tmp.X)
adata_tmp.raw = adata_tmp.copy()

# data count matrix
adata_tmp.X = adata_tmp.layers["log1p_norm"]
adata_tmp.X = scipy.sparse.csr_matrix(adata_tmp.X)

## Delete the layers
del adata_tmp.uns
# del adata_tmp.var
del adata_tmp.varm
del adata_tmp.layers
del adata_tmp.obsp

## Save UMAP embedding
df_umap = pd.DataFrame(adata_tmp.obsm["X_umap"], columns=["x","y"])
df_umap.index = adata_tmp.obs.index
df_umap.to_csv(f'{dirdata}IATEAT_APC_umap_for_R_meta.csv')

## Save anndata as loom
adata_tmp.write_loom(f'{dirdata}IATEAT_APC_for_R.loom')
```

### In R
```R
library(SeuratDisk)
library(Seurat)
library(glue)
library(dplyr)
library(qs)

tmp <- Connect(glue("{dirdata}IATEAT_APC_for_R.loom"), mode = "r")

# Convert to Seurat object
seurat <- as.Seurat(tmp, features = "gene_ids", cells = "barcode")

# UMAP
mat_UMAP <- read.csv(glue("{dirdata}IATEAT_APC_umap_for_R_meta.csv"), row.names = 1) %>%
  rename("UMAP_1" = "x", "UMAP_2" = "y") %>%
  as("matrix")
seurat[["umap"]] <- CreateDimReducObject(embeddings = mat_UMAP, key = "UMAP_", global = TRUE, assay = "RNA")

# Save seurat object as *.qs
qsave(seurat, glue("{dirdata}IATEAT_APC_for_R.qs"))
```



