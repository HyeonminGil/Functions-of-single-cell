### 1. Parallel processing
```R
options(future.seed = TRUE)
# set memory usage limit
getOption("future.globals.maxSize")  # NULL
options(future.globals.maxSize = 5 * 1024 ^ 3)  # 5 GB RAM
getOption("future.globals.maxSize")
# multicore
plan()  # sequential
plan("multicore", workers = 30)
plan()
```
---
### 2. Functions of DEG analysis within the same cell type across conditions
It can be applied to RNA or ATAC sequencing data.
```R
FindMarkers_CD <- function(InputData, CelltypeList, GroupList, GroupName) {
    # Subset seurat object
    tmp_input <- subset(InputData, subset = condition %in% GroupList)
    print(ncol(tmp_input))

    for (N in 1:length(CelltypeList)) {
        ct <- CelltypeList[N]

        n_cell1 <- nrow(tmp_input[[]] %>% filter(condition == GroupList[1], celltype == ct))
        n_cell2 <- nrow(tmp_input[[]] %>% filter(condition == GroupList[2], celltype == ct))

        if (n_cell1 < 10) {
            print(glue("The number of {ct} in {GroupList[1]} is lower 10."))
            next
        } else if (n_cell2 < 10) {
            print(glue("The number of {ct} in {GroupList[2]} is lower 10."))
            next
        } else if (n_cell1 + n_cell2 < 30) {
            print(glue("The number of {ct} is lower 30."))
            next
        }

        print(glue("The number of {ct} is {n_cell1 + n_cell2}"))
        DEG_tmp <- FindMarkers(
            tmp_input, subset.ident = ct,
            ident.1 = GroupList[1], ident.2 = GroupList[2], group.by = "condition",
            min.pct = 0, logfc.threshold = 0
            )
        DEG_tmp$condition <- GroupName
        DEG_tmp$celltype <- ct

        if (N == 1) {
            DEG_out <- DEG_tmp
        } else(
            DEG_out <- rbind(DEG_out, DEG_tmp)
        )
        }
        return(DEG_out)
    }
```
---
### 3. Application in R
```R
# Run the functions
execution.time <- system.time(
    DEG_KO_diet <- FindMarkers_CD(seurat, celltype_list, c("condition1", "condition1"), "KO_diet")
    )
# Check for running time
execution.time
```
