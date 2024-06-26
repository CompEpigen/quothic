#' @import Seurat dplyr tidyr ggplot2 mclust moments
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mclust)
library(moments)

#' Calculate the Default Reference Vector (Mean Expression) with optional log2 transformation
#'
#' @param data_matrix A matrix containing single-cell RNA-seq data
#' @return A numeric vector of mean expression values across all cells
get_default_ref <- function(data_matrix) {
  #print(data_matrix[1:5,1:5])
  ref_vec <- rowMeans(as.matrix(data_matrix))
  return(ref_vec)
}


#' Calculate Shannon Entropy for a Gene Expression Vector
#'
#' @param data_series Numeric vector of gene expression data
#' @param left Numeric, left boundary value for histogram bins
#' @param right Numeric, right boundary value for histogram bins
#' @param step Numeric, interval for histogram bins
#' @return Numeric value representing the Shannon entropy
calc_entropy <- function(data_series) {
  if (length(data_series) == 0) {
    return(NA)  # Return NA for empty data series
  }
  if (!is.numeric(data_series)) {
    stop("Data series must be numeric")
  }

  bins_cnt <- hist(data_series, breaks = "Sturges", plot = FALSE)$counts
  if (sum(bins_cnt > 0)==0){
    return(NA)
    warning("No sufficient values for entropy calculation")
  }else{
    bins_cnt <- bins_cnt[bins_cnt > 0]
    prob_vec <- bins_cnt / sum(bins_cnt)
    res_sum <- -sum(prob_vec * log2(prob_vec))
    return(res_sum)
  }
}

#' RCSA to Find the Intrinsic Reference Cell
#'
#' @param seurat_obj A Seurat object
#' @param components Integer, number of Gaussian components in the model
#' @return Numeric vector, refined reference vector after RCSA
RCSA <- function(data_matrix, components = 2) {
  idx = 0
  ref_vec <- get_default_ref(data_matrix)
  error <- sum((ref_vec^2)) / length(ref_vec)
  updated_error = error

  while (updated_error > 1e-6) {
    if (idx >= 30){
      cat('Unable to converge in 30 iterations, reverting to average as reference.\n')
      ref_vec <- get_default_ref(data_matrix)
      return(ref_vec)
    }
    cat('Calculating reference entropy...\n')
    series_entropy <- apply(data_matrix, 2, function(col) {
      calc_entropy(col - ref_vec)
    })

    gmm <- Mclust(series_entropy, G = components)
    comp <- which.max(gmm$parameters$mean)

    gau_wt <- gmm$parameters$pro[comp]
    gau_mean <- gmm$parameters$mean[comp]
    gau_sigma <- sqrt(gmm$parameters$variance$sigmasq)

    entropy_den <- gau_wt * dnorm(series_entropy, mean = gau_mean, sd = gau_sigma)
    cell_weight <- entropy_den / sum(entropy_den)

    pre_ref <- ref_vec
    ref_vec <- rowSums(as.matrix(data_matrix) * cell_weight) / sum(cell_weight)
    error <- sum((pre_ref - ref_vec)^2) / length(ref_vec)
    idx = idx + 1
    if (updated_error > error){ updated_error = error }
    cat(sprintf('\r%g: Stablility = %f\n', idx, updated_error))
  }

  return(ref_vec)
}

#' Generate Shannon Entropy for intrinsic transcription heterogeneity in a cell
#'
#' This function was amended from the scEntropy python codes (https://github.com/jzlei/scEntropy) provided by https://www.worldscientific.com/doi/abs/10.1142/S1793048020500010
#'
#' @param seurat_obj A Seurat object
#' @param ref_vec Optional numeric vector, reference cell expression vector
#' @param method Character string, method for generating reference vector. It must be either "predefined" or "RCSA".
#' @param check_log2 Logical, whether to check and perform log2 transformation if necessary
#' @param scale Normalize the entropy values to 0-1 range
#' @return Numeric vector, Shannon entropy values relative to the reference
#' @export
shannon_entropy <- function(seurat_obj, ref_vec = NULL, method = 'predefined', check_log2=TRUE, scale=FALSE) {
  data_matrix <- seurat_obj@assays$RNA@data
  # Optionally check if data needs log2 transformation
  if (check_log2) {
    # Example criterion: check if the skewness of the data is high
    skewness_values <- apply(data_matrix, 1, skewness)
    if (mean(skewness_values, na.rm = TRUE) > 2) {  # Arbitrary skewness threshold
      message("Performing log2 transformation due to high skewness")
      data_matrix <- log2(data_matrix + 1)
    }
  }

  if (method == 'predefined') {
    if (is.null(ref_vec)) {
      ref_vec <- get_default_ref(data_matrix)
    }
  } else if (method == 'RCSA') {
    ref_vec <- RCSA(data_matrix)
  } else {
    stop('Parameter error in scEntropy function')
  }

  if (length(ref_vec) != nrow(seurat_obj@assays$RNA@counts)) {
    stop('The dimension of ref_vec is not consistent with df_content!')
  }

  series_entropy <- apply(seurat_obj@assays$RNA@counts, 2, function(col) {
    calc_entropy(col - ref_vec)
  })

  # Normalize the values to c(-1,1) range
  if (scale){
    series_entropy = 2 * (series_entropy - min(series_entropy)) / (max(series_entropy) - min(series_entropy)) - 1
  }
  return(series_entropy)
}

