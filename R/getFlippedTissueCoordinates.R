#' @name getFlippedTissueCoordinates
#' @title Seurat object internally store spot coordinates (see Seurat::GetTissueCoordinates()). 
#' However, at least in the case of Visium, data are flipped and rotated before SpatialDimPlot. 
#' This function  return the rotated/flipped tissue Coordinates from a Seurat object.
#' @param seurat_obj a seurat object with tissue coordinates.
#' @param as_data_frame return x/y coords as data.frame. Default to SeuratObject.
#' @return a seurat object with slots $x_coord and $y_coord.
#' @examples
#' #' ## Install and process the brain dataset
#' library(SeuratData)
#' library(Seurat)
#' #InstallData("stxBrain")
#' anterior1 <- LoadData("stxBrain", type = "anterior1")
#' anterior1_df <- getFlippedTissueCoordinates(anterior1, as_data_frame=TRUE)
#' plot(anterior1_df)
#' @export getFlippedTissueCoordinates
getFlippedTissueCoordinates <- function(seurat_obj, 
                                        feature=NULL,
                                        as_data_frame=FALSE){
  # Get coord and flip
  coord_spot <- GetTissueCoordinates(seurat_obj)[,2:1] # rotation
  colnames(coord_spot) <- c("x", "y")
  min_coord_y <- min(coord_spot$y)
  max_coord_y <- max(coord_spot$y)
  coord_spot$y <- -coord_spot$y + 2*min_coord_y + max_coord_y-min_coord_y
  
  # prepare return
  coord_spot_x <- coord_spot$x
  names(coord_spot_x) <- colnames(seurat_obj)
  coord_spot_y <- coord_spot$x
  names(coord_spot_y) <- colnames(seurat_obj)
  if(!as_data_frame){
    seurat_obj$x_coord <- coord_spot$x 
    seurat_obj$y_coord <-coord_spot$x   
    return(seurat_obj)
  }else{
    m <- data.frame(x=coord_spot$x, 
                    y=coord_spot$y,
                    row.names = colnames(seurat_obj))
  }

  
  
}