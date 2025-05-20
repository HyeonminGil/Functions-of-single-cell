```Python
dirdata = "/home/hmgil/project/ARS_aorta/data/10x_multiome/out/3_annotation/"

adata_RNA_YC12W = sc.read_h5ad(f'{dirdata}02_seurat_EE_no_LEC_pca_Young_NCD_12W.h5ad')
adata_RNA_YC12W.layers['counts'] = adata_RNA_YC12W.raw.X.copy()
adata_RNA_YC12W.layers['log1p_norm'] = adata_RNA_YC12W.X.copy()

mapping_dict = {
    0: "EC.Clu",
    1: "EC.Lrg1",
    2: "EC.Dkk2",
    3: "EC.Cd36",
    4: "EC.Lyve1",
    5: "EC.Smoc2",
    6: "EndMT"
}

adata_RNA_YC12W.obs["ct_level3"] = adata_RNA_YC12W.obs["ct_level3"].replace(mapping_dict)

pca_YC12W = pd.read_csv(f"{dirdata}02_seurat_EE_no_LEC_pca_Young_NCD_12W.csv", index_col=0)
pca_YC12W = pca_YC12W.loc[adata_RNA_YC12W.obs_names]
adata_RNA_YC12W.obsm["X_pca"] = pca_YC12W.to_numpy()

adata_RNA_YC12W.obs["condition"] = pd.to_numeric(adata_RNA_YC12W.obs["condition"])
adata_RNA_YC12W.obs["ct_level3"] = pd.Categorical(adata_RNA_YC12W.obs["ct_level3"])

tp0_YC12W_pca = TemporalProblem(adata_RNA_YC12W)
tp0_YC12W_pca = tp0_YC12W_pca.prepare(time_key="condition", joint_attr="X_pca")

# Compute the cost matrix
import networkx as nx

dfs = {}
batch_column = "condition"
unique_batches = [1, 2]
for i in range(len(unique_batches) - 1):
    batch1 = unique_batches[i]
    batch2 = unique_batches[i + 1]

    indices = np.where(
        (adata_RNA_YC12W.obs[batch_column] == batch1) | (adata_RNA_YC12W.obs[batch_column] == batch2)
    )[0]
    adata_RNA_YC12W_subset = adata_RNA_YC12W[indices]
    sc.pp.neighbors(adata_RNA_YC12W_subset, use_rep="X_pca", n_neighbors=30)
    G = nx.from_numpy_array(adata_RNA_YC12W_subset.obsp["connectivities"].todense())
    assert nx.is_connected(G)

    dfs[(batch1, batch2)] = pd.DataFrame(
        index=adata_RNA_YC12W_subset.obs_names,
        columns=adata_RNA_YC12W_subset.obs_names,
        data=adata_RNA_YC12W_subset.obsp["connectivities"].todense().astype("float"),
    )
tp0_YC12W_pca[1, 2].set_graph_xy((dfs[1, 2]).astype("float"), t=100.0)

# Solve the temporalproblem
tp0_YC12W_pca = tp0_YC12W_pca.solve(epsilon=1e-3, scale_cost="mean", max_iterations=1e7)

# Identifying ancestors and descendants of cells
order_celltypes = [
    "EC.Clu",
    "EC.Lrg1",
    "EC.Cd36",
    "EC.Smoc2",
    "EndMT"
]

ct_desc = tp0_YC12W_pca.cell_transition(
    1,
    2,
    {"ct_level3": order_celltypes},
    {"ct_level3": order_celltypes},
    forward=False,
    key_added="transitions_2_1_large",
)
ct_desc = tp0_YC12W_pca.cell_transition(
    1,
    2,
    {"ct_level3": order_celltypes},
    {"ct_level3": order_celltypes},
    forward=True,
    key_added="transitions_1_2_large",
)

# Create a 1x2 grid of subplots
fig, axes = plt.subplots(ncols=2, figsize=(12, 6))

axes[0] = mpl.cell_transition(
    tp0_YC12W_pca,
    fontsize=7,
    figsize=(5, 5),
    return_fig=True,
    ax=axes[0],
    key="transitions_2_1_large",
)

axes[1] = mpl.cell_transition(
    tp0_YC12W_pca,
    fontsize=7,
    figsize=(5, 5),
    return_fig=True,
    ax=axes[1],
    key="transitions_1_2_large",
)

fig.subplots_adjust(wspace=0.8)
````
