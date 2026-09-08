library(Seurat)
library(SeuratObject)
library(monocle3)
library(SeuratWrappers)
library(SingleCellExperiment)
library(igraph)
library(Matrix)
library(readxl)
library(ggplot2)
library(dplyr)
library(tradeSeq)
library(SummarizedExperiment)
library(patchwork)

i <- 3
nn <- c(1, 2, 3, 4)
extract_genes <- function(sheet) {
  bimodal <- c(sheet[[1]][!is.na(sheet[[1]])], 
               sheet[[3]][!is.na(sheet[[3]])])
}
if (nn[i] == 1) {
  expression_matrix <- read_excel("Fibroblasts_c57.xlsx")
  name1 <- read_excel("gene name.xlsx", sheet = "F c57")
  DATASET_NAME <- "Fibroblasts_C57"
} else if (nn[i] == 2) {
  expression_matrix <- read_excel("Fibroblasts_cast.xlsx")
  name1 <- read_excel("gene name.xlsx", sheet = "F cast")
  DATASET_NAME <- "Fibroblasts_Cast"
} else if (nn[i] == 3) {
  expression_matrix <- read_excel("Embryonic_c57.xlsx")
  name1 <- read_excel("gene name.xlsx", sheet = "E c57")
  DATASET_NAME <- "Embryonic_C57"
} else if (nn[i] == 4) {
  expression_matrix <- read_excel("Embryonic_cast.xlsx")
  name1 <- read_excel("gene name.xlsx", sheet = "E cast")
  DATASET_NAME <- "Embryonic_Cast"
}
bimodal_genes <- extract_genes(name1)
expression_matrix[is.na(expression_matrix)] <- 0
expression_matrix <- as.data.frame(expression_matrix)
rownames(expression_matrix) <- expression_matrix[, 1]
expression_matrix <- expression_matrix[, -1]
expression_matrix <- as.matrix(expression_matrix)

seurat_obj <- CreateSeuratObject(
  counts = expression_matrix,
  project = "Mouse_Pseudotime",
  min.cells = 3,
  min.features = 200
)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
p_qc <- VlnPlot(
  seurat_obj,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ),
  pt.size = 0,
  ncol = 3
)
print(p_qc)
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 & 
    nFeature_RNA < 6000 & 
    percent.mt < 20
)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
ordering_genes <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = ordering_genes)
seurat_obj <- RunPCA(seurat_obj, features = ordering_genes)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj$seurat_clusters <- as.character(Idents(seurat_obj))

cds <- as.cell_data_set(  seurat_obj,  assay = "RNA", reductions = c("pca", "umap"), default.reduction = "umap")
rowData(cds)$gene_short_name <- rownames(cds)
cds <- cluster_cells(cds, reduction_method = "UMAP")
partition_result <- table( partitions(cds))
print(partition_result)
cds <- learn_graph( cds, use_partition = TRUE,close_loop = FALSE)
plot_cells(cds, 
           color_cells_by = "seurat_clusters", 
           label_cell_groups = FALSE,
           label_branch_points = TRUE, 
           label_leaves = TRUE)

if (nn[i] %in% c(1, 2)) {
  early_markers <- c("Cd34",  "Pdgfra",  "Pi16", "Dpp4", "Col14a1" )
  activated_markers <- c(  "Acta2","Tagln",  "Postn","Ctgf")
  cell_type <- "Fibroblast"
} else if (nn[i] %in% c(3, 4)) {
  early_markers <- c( "Pou5f1", "Sox2","Nanog",  "Klf4", "Myc" )
  activated_markers <- c( "Cdx2", "Gata4", "Gata6")
  cell_type <- "Embryonic"
}
normalized_expression <- LayerData( seurat_obj,  assay = "RNA",layer = "data")
calculate_module_zscore <- function(  genes,
                                      expression_matrix,
                                      module_name) { genes_present <- intersect(  genes,  rownames(expression_matrix)  )
                                      genes_missing <- setdiff( genes, genes_present)
                                      expression_sub <- as.matrix(  expression_matrix[ genes_present,
                                                                                       ,
                                                                                       drop = FALSE])
                                      expression_z <- t( scale(  t(expression_sub) ))
                                      expression_z[!is.finite(expression_z) ] <- 0
                                      score <- colMeans( expression_z, na.rm = TRUE)
                                      return(score)
}
early_score <- calculate_module_zscore(genes = early_markers,
                                       expression_matrix = normalized_expression,
                                       module_name = "Early-state module")
activated_score <- calculate_module_zscore( genes = activated_markers,
                                            expression_matrix = normalized_expression,
                                            module_name = "Activated-state module")
seurat_obj$early_score <- early_score[colnames(seurat_obj)]
seurat_obj$activated_score <- activated_score[  colnames(seurat_obj)]
seurat_obj$orientation_score <- seurat_obj$early_score - seurat_obj$activated_score
trajectory_graph <- principal_graph( cds)[["UMAP"]]
graph_nodes <- igraph::V( trajectory_graph)$name
leaf_nodes <- graph_nodes[ igraph::degree(  trajectory_graph ) == 1]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix( closest_vertex)
closest_vertex <- closest_vertex[ colnames(cds),
                                  ,
                                  drop = FALSE]
raw_cell_node <- as.character( closest_vertex[, 1])
if (all(raw_cell_node %in% graph_nodes)) { cell_node <- raw_cell_node} else {
  node_index <- suppressWarnings(   as.integer(raw_cell_node) )
  valid_node_index <-   !is.na(node_index) &   node_index >= 1 & node_index <= length(graph_nodes)
  if (!all(valid_node_index)) { invalid_values <- unique(  raw_cell_node[ !valid_node_index   ] ) }
  cell_node <- graph_nodes[ node_index]
}
names(cell_node) <- colnames(cds)
unmatched_nodes <- setdiff( unique(cell_node),graph_nodes)
node_to_leaf_distance <- igraph::distances( trajectory_graph,
                                            v = graph_nodes, 
                                            to = leaf_nodes, 
                                            weights = NA)

node_row_index <- match( cell_node,rownames(node_to_leaf_distance))
if (anyNA(node_row_index)) { failed_nodes <- unique( cell_node[    is.na(node_row_index)  ] )}
cell_to_leaf_distance <- node_to_leaf_distance[ node_row_index,
                                                ,
                                                drop = FALSE]
rownames(cell_to_leaf_distance) <- names( cell_node)
usable_leaf <- colSums( is.finite(cell_to_leaf_distance)) > 0
leaf_nodes <- colnames( cell_to_leaf_distance)[usable_leaf]
cell_to_leaf_distance <- cell_to_leaf_distance[
  ,
  usable_leaf,
  drop = FALSE
]
reachable_cell <- rowSums(is.finite(cell_to_leaf_distance)) > 0
cell_to_leaf_distance <- cell_to_leaf_distance[ reachable_cell,
                                                ,
                                                drop = FALSE]

nearest_leaf_index <- max.col( -cell_to_leaf_distance, ties.method = "first")
nearest_leaf_distance <- cell_to_leaf_distance[ cbind( seq_len(nrow(cell_to_leaf_distance) ), nearest_leaf_index)]
cell_leaf_assignment <- data.frame( cell = rownames( cell_to_leaf_distance ),
                                    leaf = colnames(  cell_to_leaf_distance )[ nearest_leaf_index],
                                    graph_distance = nearest_leaf_distance,
                                    stringsAsFactors = FALSE)
print( table( cell_leaf_assignment$leaf))
N_ENDPOINT_CELLS <- 15
endpoint_cells <- cell_leaf_assignment |>
  filter( is.finite(graph_distance) ) |>
  group_by(leaf) |>
  arrange( graph_distance,  .by_group = TRUE ) |>
  slice_head(  n = N_ENDPOINT_CELLS) |>
  ungroup()
print(table(  endpoint_cells$leaf))
cell_information <- FetchData(
  seurat_obj,
  vars = c(
    "seurat_clusters",
    "early_score",
    "activated_score",
    "orientation_score",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  )
)

cell_information$cell <- rownames( cell_information)
rownames(cell_information) <- NULL
cell_information <- cell_information |>
  transmute(  cell = as.character(cell),
              seurat_cluster = as.character( seurat_clusters ),
              early_score = as.numeric(  early_score ),
              activated_score = as.numeric(  activated_score ),
              orientation_score = as.numeric(  orientation_score),
              nCount_RNA = as.numeric( nCount_RNA),
              nFeature_RNA = as.numeric( nFeature_RNA),
              percent.mt = as.numeric(  percent.mt)
  )

endpoint_cells$cell <- as.character( endpoint_cells$cell)
n_matched <- sum(  endpoint_cells$cell %in%  cell_information$cell)
if (n_matched != nrow(endpoint_cells)) { unmatched_cells <- setdiff(endpoint_cells$cell,
                                                                    cell_information$cell
)
}
endpoint_cells <- endpoint_cells |>
  left_join( cell_information, by = "cell")
get_majority_cluster <- function(x) { x <- x[   !is.na(x) ]
if (length(x) == 0) { return(   NA_character_ ) }
names( sort( table(x), decreasing = TRUE) )[1]
}
endpoint_summary <- endpoint_cells |>
  group_by(leaf) |>
  summarise( n_cells = n(), n_matched = sum(  !is.na(nCount_RNA)),
             median_early_score = median( early_score,  na.rm = TRUE ),
             median_activated_score = median( activated_score, na.rm = TRUE),
             median_orientation_score = median(  orientation_score, na.rm = TRUE),
             mean_orientation_score = mean(  orientation_score,  na.rm = TRUE),
             main_cluster = get_majority_cluster(  seurat_cluster),
             median_nCount_RNA = median(  nCount_RNA,  na.rm = TRUE),
             median_nFeature_RNA = median(  nFeature_RNA,  na.rm = TRUE),
             median_percent_mt = median(  percent.mt,  na.rm = TRUE),
             .groups = "drop") |>
  arrange( desc(  median_orientation_score))
print( endpoint_summary, width = Inf)

candidate_root_node <- endpoint_summary$leaf[1]
cds_rooted <- order_cells( cds,
                           reduction_method = "UMAP",
                           root_pr_nodes = candidate_root_node
)
p_pseudotime <- plot_cells(
  cds_rooted,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = TRUE,
  label_leaves = TRUE,
  graph_label_size = 3
) +
  ggtitle(
    paste0(
      DATASET_NAME,
      "; root = ",
      candidate_root_node
    )
  )

print( p_pseudotime)
pseudotime_raw <- pseudotime( cds_rooted)
pseudotime_raw <- pseudotime_raw[  colnames(seurat_obj)]
seurat_obj$pseudotime_raw <- pseudotime_raw

valid_cells <- names(pseudotime_raw)[is.finite(pseudotime_raw)]

counts_matrix <- counts(cds_rooted)
valid_cells <- colnames(counts_matrix)[
  colnames(counts_matrix) %in% valid_cells
]

pseudotime_valid <- pseudotime_raw[valid_cells]

pseudotime_tradeSeq <- pseudotime_valid -
  min(pseudotime_valid)

counts_valid <- counts_matrix[
  ,
  valid_cells,
  drop = FALSE
]

detected_cells <- Matrix::rowSums(
  counts_valid > 0
)

detection_fraction <- detected_cells /
  length(valid_cells)

genes_for_test <- names(detected_cells)[
  detected_cells >= max(
    5,
    ceiling(0.05 * length(valid_cells))
  )
]

trajectory_test <- graph_test(
  cds_rooted[genes_for_test, ],
  neighbor_graph = "principal_graph",
  cores = 1
) |>
  as.data.frame()

trajectory_test$gene <- rownames(trajectory_test)
rownames(trajectory_test) <- NULL

trajectory_test$detection_fraction <-
  detection_fraction[trajectory_test$gene]

if ("status" %in% colnames(trajectory_test)) {
  trajectory_test <- trajectory_test |>
    filter(status == "OK")
}

trajectory_genes_high <- trajectory_test |>
  filter(
    !is.na(q_value),
    !is.na(morans_I),
    q_value < 0.05,
    morans_I >= 0.1,
    detection_fraction >= 0.10
  )
bimodal_genes_clean <- unique(
  trimws(as.character(bimodal_genes))
)

bimodal_genes_clean <- bimodal_genes_clean[
  !is.na(bimodal_genes_clean) &
    bimodal_genes_clean != ""
]

trajectory_bimodal_result <- trajectory_genes_high |>
  filter(gene %in% bimodal_genes_clean)

trajectory_bimodal_genes <-
  trajectory_bimodal_result$gene

pseudotime_matrix <- matrix(
  pseudotime_tradeSeq,
  ncol = 1,
  dimnames = list(
    valid_cells,
    "Lineage1"
  )
)

cell_weights <- matrix(
  1,
  nrow = length(valid_cells),
  ncol = 1,
  dimnames = list(
    valid_cells,
    "Lineage1"
  )
)

counts_tradeSeq <- counts_valid[
  genes_for_test,
  ,
  drop = FALSE
]

trajectory_graph <- principal_graph(
  cds_rooted
)[["UMAP"]]

set.seed(1234)

sce_gam <- tradeSeq::fitGAM(
  counts = counts_tradeSeq,
  pseudotime = pseudotime_matrix,
  cellWeights = cell_weights,
  nknots = 5,
  verbose = TRUE,
  parallel = FALSE,
  sce = TRUE
)

saveRDS(
  sce_gam,
  paste0(DATASET_NAME, "_tradeSeq_model.rds")
)


association_result <- tradeSeq::associationTest(
  sce_gam,
  global = TRUE,
  lineages = FALSE
) |>
  as.data.frame()

association_result$gene <-
  rownames(association_result)

association_result$q_value <-
  p.adjust(
    association_result$pvalue,
    method = "BH"
  )

tradeSeq_significant <- association_result |>
  filter(
    !is.na(q_value),
    q_value < 0.05
  )

high_confidence_bimodal_result <-
  trajectory_bimodal_result |>
  inner_join(
    tradeSeq_significant |>
      select(
        gene,
        tradeSeq_waldStat = waldStat,
        tradeSeq_pvalue = pvalue,
        tradeSeq_qvalue = q_value
      ),
    by = "gene"
  )

high_confidence_bimodal_genes <-
  high_confidence_bimodal_result$gene

write.csv(
  high_confidence_bimodal_result,
  paste0(
    DATASET_NAME,
    "_high_confidence_bimodal_genes.csv"
  ),
  row.names = FALSE
)


plot_tradeSeq_gene <- function(gene_name) {
  
  tradeSeq::plotSmoothers(
    sce_gam,
    counts = counts_tradeSeq,
    gene = gene_name
  ) +
    labs(
      title = gene_name,
      x = "Pseudotime",
      y = "Log(expression + 1)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "italic"
      )
    )
}

if (length(high_confidence_bimodal_genes) > 0) {
  print(
    plot_tradeSeq_gene(
      high_confidence_bimodal_genes[21]
    )
  )
}
if (length(high_confidence_bimodal_genes) > 0) {
  
  pdf_file <- paste0(DATASET_NAME, "_high_confidence_genes_trends.pdf")
  pdf(pdf_file, width = 8, height = 6)
  
  for (gene in high_confidence_bimodal_genes) {
    p <- plot_tradeSeq_gene(gene)
    print(p)
  }
  
  dev.off()
  
}
