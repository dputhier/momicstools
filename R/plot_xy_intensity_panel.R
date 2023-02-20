#' Plot a panel of spatial expression/metadata scatter plots
#'
#' This function creates a panel of scatter plots for the spatial expression of a list of genes across spots, 
#' where the X and Y coordinates represent the spatial location of spots and the color represents the 
#' expression level of the gene.
#'
#' @param seurat_obj A Seurat object containing spatial expression data.
#' @param genes A vector of gene names to plot.
#' @param metadata Provide a vector of metadata that will be used instead of genes (i.e. from 
#' meta.data) slot of a seurat object.
#' @param panel_names A vector of panel names to use for each gene plot. Default is "A", 
#' "B", "C", etc.
#' @param ncol_layout Number of columns to use for the panel layout. Default is the ceiling 
#' of the number of genes divided by 2.
#' @param intensity_slot The assay slot to use for the gene expression values.
#'        Must be one of "sct", "counts", or "data". Default is "sct".
#' @param legend Whether to display a legend for the color scale. Default is FALSE.
#' @param pt_size The size of the points in the plot. Default is 2.1.
#' @param pt_shape The shape of the points in the plot. Default is 16 (a circle).
#' @param colours A vector of colors.
#' @importFrom ggplot2 ggplot theme_void
#' @importFrom patchwork plot_layout
#'
#' @return A ggplot2 object containing the panel of scatter plots.
#'
#' @examples
#' library(SeuratData)
#' library(Seurat)
#' #InstallData("stxBrain")
#' set_verbosity(3)
#' anterior1 <- LoadData("stxBrain", type = "anterior1")
#' anterior1 <- SCTransform(anterior1, assay = "Spatial")
#' genes_to_plot <- c("Hpca", "Olig1", "Klf2")
#' panel <- c("Panel A", "Panel B", "Panel C")
#' plot_xy_intensity_panel(seurat_obj = anterior1, genes = genes_to_plot, panel_names = panel)
#' metadata <- c("nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT")
#' plot_xy_intensity_panel(seurat_obj = anterior1, metadata = metadata)
#' @export
plot_xy_intensity_panel <- function(seurat_obj=NULL,
                                    genes=NULL,
                                    metadata=NULL,
                                    intensity_slot=c("data", "counts", "sct"),
                                    panel_names=NULL,
                                    ncol_layout=NULL,
                                    legend=TRUE,
                                    pt_size=2.1,
                                    pt_shape=16,
                                    colours=colors_for_gradient("J1")
                                    ){
  
  if(is.null(panel_names)){
    if(is.null(genes)){
      panel_names <- LETTERS[1:length(metadata)]
    }else{
      panel_names <- LETTERS[1:length(genes)]
    }
  }
   
  if(is.null(ncol_layout)){
    if(is.null(genes)){
      ncol_layout <- ceiling(length(metadata)/2)
    }else{
      ncol_layout <- ceiling(length(genes)/2)
    }
  }
  

  if(is.null(seurat_obj))
    print_msg("Please provide a seurat object...", msg_type = "STOP")

  if(is.null(genes) & is.null(metadata))
    print_msg("Please provide a list of genes or metadata...", msg_type = "STOP")

  
  if(is.null(genes) & is.null(metadata))
    print_msg("Please provide a list of genes or metadata...", msg_type = "STOP") 
  
  print_msg(paste0("Panel names : ",  panel_names), msg_type = "DEBUG")

  plot_panels <- NULL
  
  if(is.null(metadata)){
    if(length(panel_names) != length(genes))
      print_msg("panel_names and genes should have same length.", msg_type = "STOP") 
    
    for(i in 1:length(genes)){
      
      gene_curr <- genes[i]
      panel_curr <- panel_names[i]
      
      print_msg(paste0("Creating diagram for gene: ",  gene_curr), msg_type = "DEBUG")
      
      plot_cur <- plot_xy_intensity(seurat_obj=seurat_obj,
                                    gene_name=gene_curr,
                                    intensity_slot=intensity_slot,
                                    title=panel_curr,
                                    legend=TRUE)
      
      if(is.null(plot_panels)){
        plot_panels <- plot_cur
      }else{
        plot_panels <- plot_panels + plot_cur
      }
    }
  }else{
    for(i in 1:length(metadata)){
      
      
      metadata_curr <- metadata[i]
      panel_curr <- panel_names[i]
      
      print_msg(paste0("Creating diagram for metadata: ",  metadata_curr), msg_type = "DEBUG")
      
      plot_cur <- plot_xy_intensity(seurat_obj=seurat_obj,
                                    metadata=metadata_curr,
                                    title=panel_curr,
                                    legend=legend,
                                    pt_size=pt_size,
                                    pt_shape=pt_shape,
                                    colours=colours)
      
      if(is.null(plot_panels)){
        plot_panels <- plot_cur
      }else{
        plot_panels <- plot_panels + plot_cur
      }
    }
  }
  
  print_msg("Preparing diagram layout", msg_type = "DEBUG")
  
  plot_panels + patchwork::plot_layout(ncol=ncol_layout)
}