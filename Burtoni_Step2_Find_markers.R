burtoni.markers <- FindAllMarkers(
  burtoni.snseq.combined.sct,
  only.pos = TRUE,
  min.pct = 0.25, 
  group.by = "seurat_clusters",
  logfc.threshold = 0.25
)


# Top markers
top.markers_burtoni <- burtoni.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

head(top.markers_burtoni, 30)


top.markers_burtoni %>%
  group_by(cluster) %>%
  summarise(top_genes = paste(gene, collapse = ", ")) %>%
  print(n = Inf)

# based on
levels(burtoni.snseq.combined.sct)
length(levels(burtoni.snseq.combined.sct))
top10_burtoni <- burtoni.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%         
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()


burtoni.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_burtoni
