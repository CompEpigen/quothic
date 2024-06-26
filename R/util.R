#' @import Seurat ggplot2 dplyr
library(Seurat)
library(ggplot2)
library(dplyr)

#' Scatter Plot from Seurat Object Metadata
#'
#' @param seurat_obj A Seurat object.
#' @param x_col The column name to be used for the x-axis.
#' @param y_col The column name to be used for the y-axis.
#' @param cluster_col The column name to be used for coloring the dots.
#' @description This function generates a scatter plot using two columns from the metadata of a Seurat object and colors the dots based on a specified column of clusters.
#' @export
#' @examples
#' library(Seurat)
#' # Assuming seurat_obj is a Seurat object with appropriate metadata
#' scatter_plot_seurat(seurat_obj, "UMAP_1", "UMAP_2", "Cell_Type")
scatter_plot_seurat <- function(seurat_obj, x_col, y_col, cluster_col) {
  # Ensure the columns exist in metadata
  if (!(x_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", x_col, "not found in metadata"))
  }
  if (!(y_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", y_col, "not found in metadata"))
  }
  if (!(cluster_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", cluster_col, "not found in metadata"))
  }

  seurat_obj@meta.data[,cluster_col] = as.character(seurat_obj@meta.data[,cluster_col])
  # Calculate means of x and y columns
  mean_x <- mean(seurat_obj@meta.data[[x_col]], na.rm = TRUE)
  mean_y <- mean(seurat_obj@meta.data[[y_col]], na.rm = TRUE)

  # Create the scatter plot
  p <- ggplot(seurat_obj@meta.data, aes_string(x = x_col, y = y_col, color = cluster_col)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "black") +
    labs(color = cluster_col) +
    xlab(paste0("Dedifferentiation (", x_col, ")")) +
    ylab(paste0("Transcriptome heterogeneity (", y_col, ")")) +  # fixed to use y_col for the y-axis label
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  # Print the plot
  print(p)

}
