#' Given a Seurat Spatial object, this function creates a scatter plot of the 
#' spatial expression of a gene across spots, where the X and Y coordinates 
#' represent the spatial location of spots and the color represents the 
#' expression level of the gene.
#'
#' @param seurat_obj A Seurat object containing spatial expression data.
#' @param gene_name The name of the gene to plot.
#' @param metadata Provide the name of a metadata that will be used instead of genes (i.e. from 
#' meta.data) slot of a seurat object.
#' @param intensity_slot The assay slot to use for the gene expression values.
#'        Must be one of "sct", "counts", or "data". Default is "sct".
#' @param title The title of the plot. Default is an empty string.
#' @param legend Whether to display a legend for the color scale. Default is FALSE.
#' @param pt_size The size of the points in the plot. Default is 2.1.
#' @param pt_shape The shape of the points in the plot. Default is 16 (a circle).
#'
#' @return A ggplot2 object containing the scatter plot.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_gradientn theme_void
#'              ggtitle element_text margin
#' @importFrom Seurat NoLegend
#'
#' @examples
#' library(SeuratData)
#' library(Seurat)
#' #InstallData("stxBrain")
#' anterior1 <- LoadData("stxBrain", type = "anterior1")
#' plot_xy_intensity(seurat_obj = anterior1, gene_name = "Hpca")
#' anterior1 <- SCTransform(anterior1, assay = "Spatial")
#' plot_xy_intensity(seurat_obj = anterior1, gene_name = "Hpca", intensity_slot="sct")
#' plot_xy_intensity(seurat_obj = anterior1, metadata = "nCount_SCT")
plot_xy_intensity <- function(seurat_obj=NULL,
                              gene_name=NULL,
                              metadata=NULL,
                              intensity_slot=c("data", "counts", "sct"),
                              title="",
                              legend=TRUE,
                              pt_size=2.1,
                              pt_shape=16,
                              colours=c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A")){
  
  intensity_slot <- match.arg(intensity_slot)

  if(is.null(seurat_obj))
    print_msg("Please provide a seurat object...", msg_type = "STOP")
  
  if(is.null(gene_name) & is.null(metadata))
    print_msg("Please provide a value for gene_name or metadata...", msg_type = "STOP")
  
  print_msg("Getting x/y coordinates", msg_type = "DEBUG")
  
  xy_coord <- getFlippedTissueCoordinates(seurat_obj, 
                                          as_data_frame = TRUE)

  print_msg("Extracting expression values", msg_type = "DEBUG")
  
  if(!is.null(metadata)){

    if(!metadata %in% colnames(seurat_obj@meta.data))
      print_msg("The metadata was not found in the object", msg_type = "STOP")
    intensities <- as.vector(seurat_obj@meta.data[,metadata])
  }else{
    if(intensity_slot=="sct"){
      slot_intensity <- seurat_obj@assays$SC
    }else if(intensity_slot=="counts"){
      slot_intensity <- seurat_obj@assays$Spatial@counts
    }else if(intensity_slot=="data"){
      slot_intensity <- seurat_obj@assays$Spatial@data
    }
    
    if(is.null(slot_intensity))
      print_msg("Slot is empty.", msg_type = "STOP")
    
    if(!gene_name %in% rownames(slot_intensity))
      print_msg("The gene_name was not found in the object", msg_type = "STOP")
    intensities <- as.vector(slot_intensity[gene_name, ])
  }
  print_msg("Creating a ggplot diagram.", msg_type = "DEBUG")
  
  df <- cbind(xy_coord, intensities)
  colnames(df) <- c("x", "y", "intensity")
  p<- ggplot(df, aes(x=x, y=y, color=intensity)) +
    geom_point(size=pt_size, shape = pt_shape) +
    theme_void() +
    scale_color_gradientn(colours=colours) +
    ggtitle(title) +
    theme(legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.margin = margin(0,0,0,0),
          plot.title = element_text(face = "bold", size=20))
  if(!legend)
    p <- p + NoLegend()

  return(p)
}





  
