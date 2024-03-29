% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visium_hull.R
\name{visium_hull}
\alias{visium_hull}
\title{Internal function. Compute location of segments used to create an orthogonal
hull around points related to a particular class.}
\usage{
visium_hull(
  data,
  x = "x",
  y = "y",
  k = "k",
  size_x = 4.6,
  size_y = 5,
  step_x = 2.6,
  step_y = 2.4,
  delta = 0.3
)
}
\arguments{
\item{data}{a data.frame.}

\item{x}{The column name storing the x coord}

\item{y}{The column name storing the y coord}

\item{k}{The column name storing the classes (0 not part of the class of interest, 1 part of the class of interest)}

\item{size_x}{The size of the square (x axis)}

\item{size_y}{The size of the square (y axis)}

\item{step_x}{The distance between two points on the x axis.}

\item{step_y}{The distance between two points on the y axis.}

\item{delta}{Add more or less flexibility to search for neighbor points}
}
\value{
A dataframe with coordinates x, y, xend, yend (see geom_segments).
}
\description{
Internal function. Compute location of segments used to create an orthogonal
hull around points related to a particular class.
}
\examples{
if (! "stxBrain" \%in\% InstalledData()$Dataset) InstallData('stxBrain')
LoadData("stxBrain", type = "anterior1")
# Preprocess the data
anterior1_preprocessed <- preprocess_seurat_data(anterior1, normalization_method = "NormalizeData", approx=F)
## Retrieve x/y coordinates and group from ggplot object
coord_st_data <- getFlippedTissueCoordinates(anterior1_preprocessed, as_data_frame=TRUE)
# Create a hull around cluster 0
coord_st_data$k <- ifelse(Idents(anterior1_preprocessed)==0, 1, 0)
## Compute the segments of the hull.
path <- visium_hull(coord_st_data, size_x=3.2, size_y=3.6, delta=0.5)
## Add the segments to the ggplot diagram
spatial_plot <- SpatialDimPlot(anterior1_preprocessed, label = TRUE, label.size = 3, pt.size.factor = 1.5)
spatial_plot +
  theme_bw() +
  geom_segment(data=path,
               mapping=aes(x=x1,
                           y=y1,
                           xend=x2,
                           yend=y2),
               inherit.aes = F,
               color="white",
               size=0.7)
}
\keyword{hull,}
\keyword{spatial}
\keyword{transcriptomics,}
\keyword{visium}
