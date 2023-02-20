test_that("test for plot_xy_intensity_panel...", {
  
  set_verbosity(1)
  if (! "stxBrain" %in% InstalledData()$Dataset) InstallData('stxBrain')
  
  LoadData("stxBrain", type = "anterior1")
  
  # Preprocess the data
  anterior1_preprocessed <- preprocess_seurat_data(anterior1, 
                                                   normalization_method = "NormalizeData", 
                                                   approx=FALSE)
  
  ## Retrieve x/y coordinates and group from ggplot object
  coord_st_data <- getFlippedTissueCoordinates(anterior1, as_data_frame=TRUE)
  coord_st_data$k <- ifelse(Idents(anterior1_preprocessed)==0, 1, 0)
  
  ## Compute the segments of the hull.
  path <- visium_hull(coord_st_data, size_x=3.2, size_y=3.6, delta=0.5)
  
  ## Add the segments to the ggplot diagram
  spatial_plot <- SpatialDimPlot(anterior1_preprocessed, 
                                 label = TRUE, 
                                 label.size = 3, 
                                 colorhull="white",
                                 pt.size.factor = 1.5)
  spatial_plot +
   theme_bw() +
   geom_segment(data=path,
                mapping=aes(x=x1,
                            y=y1,
                            xend=x2,
                            yend=y2),
                inherit.aes = F,
                color=colorhull,
                size=0.7)
}
