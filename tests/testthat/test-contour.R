test_that("visium_hull works", {
  data(coord_spot_brain)
  cluster_to_show <- 7 # Could be also c(a, b)
  coord_spot_brain$k <- ifelse(coord_spot_brain$group %in% cluster_to_show, 1, 0)
  path <- visium_hull(coord_spot_brain, size_x=3.2, size_y=3.6, delta=0.5)
  #ggplot(coord_spot_brain, aes(x = x, y = y, col=group)) + geom_point() +
  #  geom_segment(data=path, mapping=aes(x=x1, y=y1, xend=x2, yend=y2),
  #               inherit.aes = F, color="red", size=0.7)
  expect_equal(sum(apply(path, 2, sum)), 591666.5, tolerance = 0.0000001)

  data(coord_spot_brain)
  cluster_to_show <- 4 # Could be also c(a, b)
  coord_spot_brain$k <- ifelse(coord_spot_brain$group %in% cluster_to_show, 1, 0)
  path <- visium_hull(coord_spot_brain, size_x=3.2, size_y=3.6, delta=0.5)
  expect_equal(sum(apply(path, 2, sum)), 701231.8, tolerance = 0.0000001)

  data(coord_spot_brain)
  cluster_to_show <- 0 # Could be also c(a, b)
  coord_spot_brain$k <- ifelse(coord_spot_brain$group %in% cluster_to_show, 1, 0)
  path <- visium_hull(coord_spot_brain, size_x=3.2, size_y=3.6, delta=0.5)
  expect_equal(sum(apply(path, 2, sum)), 420257.7, tolerance = 0.0000001)
})
