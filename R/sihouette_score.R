#' @import Seurat cluster ggplot2 dplyr viridis
library(Seurat)
library(cluster)
library(ggplot2)
library(dplyr)
library(viridis)

#' Calculate Silhouette Scores for a Seurat Object
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col A column name in the metadata specifying clusters.
#' @return A data.frame containing the silhouette scores for PCA, UMAP, and t-SNE embeddings.
#' @description This function calculates silhouette scores for each cell in a Seurat object based on the clusters specified in a metadata column. It adds the silhouette scores to the metadata for PCA, UMAP, and t-SNE embeddings.
#' @export
#' @examples
#' library(Seurat)
#' # Assuming seurat_obj is a Seurat object and "Cell_line" is a metadata column
#' scores <- clusterTH_silhouette_score(seurat_obj, "Cell_line")
#' head(scores)
clusterTH_silhouette_score <- function(seurat_obj, cluster_col) {
  # Ensure the column exists in metadata
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in metadata"))
  }

  # Function to calculate silhouette scores
  calculate_silhouette <- function(embedding, clusters) {
    distance_matrix <- dist(embedding)
    silhouette_scores <- silhouette(as.numeric(clusters), dist = distance_matrix)
    return(silhouette_scores[, 3])
  }

  # Get clusters from metadata
  clusters <- as.factor(seurat_obj@meta.data[[cluster_col]])
  levels(clusters) <- 1:length(levels(clusters))

  # Calculate silhouette scores for PCA, UMAP, and t-SNE, and convert them to heterogeneity scores (the higher, the more heterogeneity)
  clusterTH_pca <- -1*calculate_silhouette(Embeddings(seurat_obj[['pca']]), clusters)
  clusterTH_umap <- -1*calculate_silhouette(Embeddings(seurat_obj[['umap']]), clusters)
  clusterTH_tsne <- -1*calculate_silhouette(Embeddings(seurat_obj[['tsne']]), clusters)

  # Create a data.frame with the silhouette scores
  scores <- data.frame(
    barcode = rownames(seurat_obj@meta.data),
    clusterTH_pca = clusterTH_pca,
    clusterTH_umap = clusterTH_umap,
    clusterTH_tsne = clusterTH_tsne
  )

  # Add the silhouette scores to the Seurat object metadata
  seurat_obj@meta.data$clusterTH_pca <- clusterTH_pca
  seurat_obj@meta.data$clusterTH_umap <- clusterTH_umap
  seurat_obj@meta.data$clusterTH_tsne <- clusterTH_tsne

  return(scores)
}

#' Plot Silhouette Scores
#'
#' @param seurat_obj A Seurat object with silhouette scores in metadata.
#' @param cluster_col A column name in the metadata specifying clusters.
#' @description This function generates plots for silhouette scores.
#' @export
#' @examples
#' library(Seurat)
#' # Assuming seurat_obj is a Seurat object and silhouette scores are calculated
#' plot_silhouette_scores(seurat_obj, "Cell_line")
plot_silhouette_scores <- function(seurat_obj, cluster_col) {
  # Ensure the column exists in metadata
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in metadata"))
  }

  # Plotting for PCA
  mean_clusterTH_pca <- mean(seurat_obj@meta.data$clusterTH_pca)
  df_pca <- seurat_obj@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(!!sym(cluster_col), -clusterTH_pca) %>%
    mutate(barcode = factor(barcode, levels = barcode))

  p1 <- ggplot(df_pca) +
    geom_col(aes(barcode, clusterTH_pca, fill = !!sym(cluster_col)), show.legend = FALSE) +
    geom_hline(yintercept = mean_clusterTH_pca, color = 'red', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score (PCA)') +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plotting for UMAP
  mean_clusterTH_umap <- mean(seurat_obj@meta.data$clusterTH_umap)
  df_umap <- seurat_obj@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(!!sym(cluster_col), -clusterTH_umap) %>%
    mutate(barcode = factor(barcode, levels = barcode))

  p2 <- ggplot(df_umap) +
    geom_col(aes(barcode, clusterTH_umap, fill = !!sym(cluster_col)), show.legend = FALSE) +
    geom_hline(yintercept = mean_clusterTH_umap, color = 'red', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score (UMAP)') +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plotting for t-SNE
  mean_clusterTH_tsne <- mean(seurat_obj@meta.data$clusterTH_tsne)
  df_tsne <- seurat_obj@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(!!sym(cluster_col), -clusterTH_tsne) %>%
    mutate(barcode = factor(barcode, levels = barcode))

  p3 <- ggplot(df_tsne) +
    geom_col(aes(barcode, clusterTH_tsne, fill = !!sym(cluster_col)), show.legend = FALSE) +
    geom_hline(yintercept = mean_clusterTH_tsne, color = 'red', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score (t-SNE)') +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Feature plots
  fp1 <- FeaturePlot(seurat_obj, reduction = "pca", features = "clusterTH_pca") +
    scale_color_viridis(discrete = FALSE)

  fp2 <- FeaturePlot(seurat_obj, reduction = "umap", features = "clusterTH_umap") +
    scale_color_viridis(discrete = FALSE)

  fp3 <- FeaturePlot(seurat_obj, reduction = "tsne", features = "clusterTH_tsne") +
    scale_color_viridis(discrete = FALSE)

  # Print plots
  print(p1)
  print(p2)
  print(p3)
  print(fp1)
  print(fp2)
  print(fp3)
}
