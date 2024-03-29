% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_gallery.R
\name{displayGallery_imgObj}
\alias{displayGallery_imgObj}
\title{displayGallery_imgObj}
\usage{
displayGallery_imgObj(
  imgObj,
  cell_ids,
  imgType = c("labels", "morph"),
  morphToShow = c(1, 2, 3),
  method = c("raster", "browser"),
  allFrame = TRUE,
  colormode = c("color", "grayscale"),
  imagesPerRow = NULL
)
}
\arguments{
\item{imgObj}{multi-frame image object for 2D cell label matrix of 3D rgb array}

\item{cell_ids}{a vector of cell_id for each frame in input image obj (default = NULL if not to display cell_ids on plot)}

\item{imgType}{type of input image in \code{imgObj}, labels for 2D matrix, morph for multi-channel 3D array}

\item{morphToShow}{a vector of 3 elements for the channel index of each RGB channel in \code{imgObj}.}

\item{method}{display method, 'raster' as gallery in plot, 'browser' as one cell per frame in viewer}

\item{allFrame}{flag to show all frames, refer to \code{EBImage::display}}

\item{colormode}{color mode to display, default = color}

\item{imagesPerRow}{number of images per row, default = NULL to pick automatically}
}
\description{
display multi-frame image object as cell gallery in plot window
}
