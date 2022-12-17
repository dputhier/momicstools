#' @name infer_cell_dist_to_border
#' @title Given a list of vectors containing cell markers, infer the relative location
#' @param data a data.frame.
#' @param x the column name storing the x coord
#' @param y the column name storing the y coord
#' @param k the column name storing the classes (0 not part of the class of interest, 1 part of the class of interest)
#' @param size_y the size of the square (y axis)
#' @param size_x the size of the square (x axis)
#' @param step_y the distance between two points on the y axis.
#' @param step_x the distance between two points on the x axis.
#' @param delta add more or less flexibility to search for neighbor points
#' @param verbose whether the function should print (debug) info during processing.
#' @keywords hull, spatial transcriptomics, visium
#' @return a dataframe with coordinates x, y, xend, yend (see geom_segments).
#' @examples
#' ## Install and process the brain dataset
#' library(Seurat)
#' library(SeuratData)
#' library(ggplot2)
#' InstallData("stxBrain")
#' library(ohmiki)
#' brain <- LoadData("stxBrain", type = "anterior1")
#' brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#' brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
#' brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
#' brain <- FindClusters(brain, verbose = FALSE)
#' brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
#' spatial_plot <- SpatialDimPlot(brain, label = TRUE, label.size = 3, pt.size.factor = 1.5)
#' ## Retrieve x/y coordinates and group from ggplot object
#' coord_st_data <- ggplot_build(spatial_plot)$data[[1]][,c("x", "y", "group")]
#' coord_st_data$group <- coord_st_data$group - 1 # group are 1-based in ggplot compare to seurat
#' ## Cluster 0 is, for instance, the cluster of interest.
#' cluster_to_show <- 0 # Could be also c(a, b)
#' coord_st_data$k <- ifelse(coord_st_data$group %in% cluster_to_show, 1, 0)
#' ## Compute the segments of the hull.
#' path <- compute_visium_ortho_hull(coord_st_data, size_x=3.2, size_y=3.6, delta=0.5)
#' ## Add the segments to the ggplot diagram
#' spatial_plot +
#'   theme_bw() +
#'   geom_segment(data=path,
#'                mapping=aes(x=x1,
#'                            y=y1,
#'                            xend=x2,
#'                            yend=y2),
#'                inherit.aes = F,
#'                color="white",
#'                size=0.7)
#' @export compute_visium_ortho_hull
compute_visium_ortho_hull <- function(data,
