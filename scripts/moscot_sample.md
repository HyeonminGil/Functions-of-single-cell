`Python
dirdata = "/home/hmgil/project/ARS_aorta/data/10x_multiome/out/3_annotation/"

adata_RNA_YC12W = sc.read_h5ad(f'{dirdata}02_seurat_EE_no_LEC_pca_Young_NCD_12W.h5ad')
adata_RNA_YC12W.layers['counts'] = adata_RNA_YC12W.raw.X.copy()
adata_RNA_YC12W.layers['log1p_norm'] = adata_RNA_YC12W.X.copy()


`
