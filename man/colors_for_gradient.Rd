% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{colors_for_gradient}
\alias{colors_for_gradient}
\title{Generate a vector of colors for a gradient}
\usage{
colors_for_gradient(
  palette = c("Je1", "Seurat_like", "Ju1", "De1", "De2", "De3", "De4", "De5", "De6",
    "De7", "De8")
)
}
\arguments{
\item{palette}{A character vector specifying the palette to use. One of: "Je1",
"Seurat_like", "Ju1", "De1",  "De2", "De3", "De4", "De5", "De6", "De7", "De8".}
}
\value{
A character vector of color codes.
}
\description{
This function generates a vector of colors for a gradient, given
a specified palette name.
}
\examples{
colors_for_gradient()
colors_for_gradient(palette = "Seurat_like")

}
