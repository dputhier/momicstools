#' @name infer_cell_dist_to_border
#' @title Given a list of vectors containing cell markers, infer the relative location of cell relative to border.
#' @param seurat_object a Seurat object
#' @param dist_to_border a dataframe produced by compute_dist_to_border().
#' @keywords distance, spots, spatial transcriptomics, visium
#' @return a dataframe with coordinates x, y, xend, yend (see geom_segments).
#' @examples
#' ## Install and process the brain dataset
#' library(SeuratData)
#' library(ggplot2)
#' #InstallData("stxBrain")
#' brain <- LoadData("stxBrain", type = "anterior1")
#' brain <- preprocess(brain, normalization="sct")
#' brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#' brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
#' brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
#' brain <- FindClusters(brain, verbose = FALSE)
#' brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
#' # Select some points (just to give an example)
#' coord_spot_brain <- getFlippedTissueCoordinates(brain)
#' coord_spot_brain$k <- 0
#' coord_spot_brain[SeuratObject::WhichCells(brain, idents=0), ]$k <- 1 # class 0 is the class of interest (labeled 1 against 0 for others)
#' border_segments <- compute_visium_ortho_hull(coord_spot_brain, size_x=3.6, size_y=3.4, delta=0.5)
#' dist_to_border <- compute_dist_to_border(coord_spot_brain, border_segments)
#' markers <- list(c("A", "B", "C"), c("D", "E", "F" ))
#' infer_cell_dist_to_border(brain, dist_to_border, markers)
#' @export infer_cell_dist_to_border
infer_cell_dist_to_border <- function(seurat_object,
                                      dist_to_border,
                                      markers,
                                      name="dist2border",
                                      dist_quantiles=seq(from=0, to=100, by=10)){


}

# ...
#dist_to_border[dist_to_border$k==0, 'dist2_inter'] <- -dist_to_border[dist_to_border$k==0, 'dist2_inter']
#cut(dist_to_border$dist2_inter, 6)
#brain[["dist2border"]] <- dist_to_border$dist2_inter
#brain[["dist2border_class"]] <-ggplot2::cut_interval(sort(dist_to_border$dist2_inter), 10)
#plot(sort(dist_to_border$dist2_inter))
