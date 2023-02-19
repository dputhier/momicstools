test_that("test for plot_xy_intensity_panel...", {
  
  # Define a Seurat object to use for testing
  data("stxBrain", package = "SeuratData")
  anterior1 <- stxBrain[["anterior1"]]
  anterior1 <- SCTransform(anterior1, assay = "Spatial")
  
  # Test 1: function returns a ggplot object
  test_that("function returns ggplot object", {
    gg <- plot_xy_intensity_panel(anterior1, c("Hpca", "Gja1", "Dlx2"), panel_names = c("A", "B", "C"))
    expect_s3_class(gg, "gg")
  })
  
  # Test 2: function throws error if no Seurat object is provided
  test_that("function throws error if no Seurat object is provided", {
    expect_error(plot_xy_intensity_panel(genes = c("Hpca", "Gja1", "Dlx2"), panel_names = c("A", "B", "C")), "Please provide a seurat object...")
  })
  
  # Test 3: function throws error if gene list and panel name list are not the same length
  test_that("function throws error if gene list and panel name list are not the same length", {
    expect_error(plot_xy_intensity_panel(anterior1, c("Hpca", "Gja1", "Dlx2"), panel_names = c("A", "B")), "panel_names and genes should have same length.")
  })
  
}