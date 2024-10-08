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
