```R
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(glue)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p <- UMAPPlot(seurat, group.by = 'celltype') + 
    theme_void() + 
    ggtitle(NULL) +
    theme(text = element_text(size = 18)) +
    scale_color_manual(values = col_vector)
ggsave(glue('{dirout}umap_celltype.png'), p, width = 5.5, height = 4)
```

### Reference for color palette
https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
