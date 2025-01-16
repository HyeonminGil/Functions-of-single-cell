```R
## Downsamping function
downsample_cells <- function(seurat_obj, group_column, n_cells) {
  
  # create metadata
  meta_data <- seurat_obj@meta.data
  
  # Sampling n_cells for each group
  sampled_cells <- meta_data %>%
    dplyr::mutate(barcode = rownames(.)) %>%
    dplyr::group_by(!!sym(group_column)) %>%
    dplyr::slice_sample(n = n_cells) %>% # random sampling across group
    dplyr::pull(barcode)
  
  # Create a new Seurat object with sampled cells
  downsampled_obj <- subset(seurat_obj, cells = sampled_cells)
  
  return(downsampled_obj)
}

## Run a downsmpling function
group_column <- "ct_level1"
n_cells <- 1500

seurat_down <- downsample_cells(subset(seurat, subset = ct_level1 != "MC"), group_column, n_cells)
seurat_down
```
