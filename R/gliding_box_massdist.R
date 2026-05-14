# R/gliding_box_massdist.R

#' Gliding-box mass distribution
#'
#' Computes the frequency distribution of box masses for a 3D binary
#' voxel array using the gliding-box method with an integral image.
#'
#' @param C A 3D logical or integer array (voxels).
#' @param box_sizes Integer vector of box side lengths to evaluate.
#'
#' @return A list of data frames, one per box size, each with columns:
#'   \describe{
#'     \item{box_size}{The box side length.}
#'     \item{mass}{Number of filled voxels in the box.}
#'     \item{frequency}{How many boxes had that mass.}
#'   }
#'
#' @export
gliding_box_massdist <- function(C, box_sizes) {

  # Input checks in R before touching C++
  if (!is.array(C) || length(dim(C)) != 3L) {
    stop("`C` must be a 3D array.")
  }

  box_sizes <- as.integer(box_sizes)

  if (any(is.na(box_sizes)) || any(box_sizes <= 0L)) {
    stop("`box_sizes` must be a vector of positive integers.")
  }

  # Convert to unsigned integer cube for Rcpp
  storage.mode(C) <- "integer"

  .gliding_box_massdist(C, box_sizes)
}