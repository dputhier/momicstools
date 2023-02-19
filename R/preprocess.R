#' Preprocess scRNA-seq or spatial transcriptomics data with Seurat
#'
#' This function takes a raw dataset in Seurat object format and preprocesses it using
#' NormalizeData or SCTransform, RunPCA, FindNeighbors, FindClusters, and RunUMAP.
#'
#' @param seurat_obj A raw dataset in Seurat object format.
#' @param normalization_method A character string specifying the normalization method to use,
#' either "NormalizeData" or "SCTransform". Defaults to "NormalizeData".
#' @param n.features The number of features to use in the PCA analysis. Defaults to 2000.
#' @param dims The number of dimensions to use in the PCA analysis. Defaults to 50.
#' @param resolution The resolution parameter for clustering. Defaults to 0.6.
#' @param n.neighbors The number of neighbors for the FindNeighbors step. Defaults to 30.
#' @param seed An optional seed for reproducibility.
#'
#' @return A preprocessed Seurat object.
#'
#' @examples
#' # Load the stxKidney dataset (spatial transcriptomics)
#' 
#' #InstallData('stxBrain')
#' LoadData("stxBrain", type = "anterior1")
#' 
#' # Preprocess the data
#' anterior1_preprocessed <- preprocess_seurat_data(anterior1, normalization_method = "SCTransform")
#'
#' # Load the pbmc3k dataset (scRNA-seq)
#' data(pbmc3k)
#' 
#' # Preprocess the data
#' pbmc3k_preprocessed <- preprocess_seurat_data(pbmc3k)
#'
#' @import Seurat
#' @importFrom Seurat NoLegend
#' @importFrom ggplot2 ggplot geom_point aes scale_color_gradientn 
preprocess_seurat_data <- function(seurat_obj, 
                                   normalization_method = c("NormalizeData", "SCTransform"),
                                   n.features = 2000, 
                                   dims = 50, 
                                   resolution = 0.6, 
                                   n.neighbors = 30, 
                                   seed = NULL) {
  
  normalization_method <- match.args(normalization_method)
  
  # Normalize the data
  if("Spatial" %in% names(seurat_obj)){
    assay <- "Spatial"
  }else if("RNA" %in% names(seurat_obj)){
    assay <- "RNA"
  }

  if (normalization_method == "SCTransform") {
    seurat_obj <- Seurat::SCTransform(seurat_obj, 
                                      assay=assay, 
                                      verbose = FALSE)
  } else {
    seurat_obj <- Seurat::NormalizeData(seurat_obj, 
                                        normalization.method = "LogNormalize", 
                                        scale.factor = 10000)
  }
  
  # Run PCA
  seurat_obj <- Seurat::RunPCA(seurat_obj, 
                               npcs = n.features, 
                               verbose = FALSE)
  
  # Find neighbors
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dims, 
                                      verbose = FALSE, 
                                      reduction.use = "pca", 
                                      k.param = n.neighbors, 
                                      force.recalc = TRUE, 
                                      seed.use = seed)
  
  # Find clusters
  seurat_obj <- Seurat::FindClusters(seurat_obj, 
                             resolution = resolution, 
                             verbose = FALSE)
  
  # Run UMAP
  seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                        dims = 1:dims, 
                        verbose = FALSE, 
                        seed.use = seed)
  
  return(seurat_obj)
} 