// A script for extracting 3D gliding-box mass distributions
// from voxel arrays using the integral image algorithm
// Extension for lacunr package

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export(.gliding_box_massdist)]]
Rcpp::List gliding_box_massdist(
    arma::ucube C,
    IntegerVector box_sizes) {

  // ---------------------------------------
  // Validate inputs up front
  // ---------------------------------------
  if (C.n_elem == 0) {
    stop("Input voxel array is empty.");
  }

  if (box_sizes.size() == 0) {
    stop("box_sizes vector is empty.");
  }

  for (int i = 0; i < box_sizes.size(); i++) {
    if (box_sizes[i] <= 0) {
      stop("All box sizes must be positive integers.");
    }
  }

  // ---------------------------------------
  // Define dimensions
  // ---------------------------------------
  unsigned int Xdim = C.n_rows;
  unsigned int Ydim = C.n_cols;
  unsigned int Zdim = C.n_slices;

  // ---------------------------------------
  // Initialize integral image
  // (+1 padding for summed-volume table,
  //  zero-initialized so border is 0)
  // ---------------------------------------
  arma::dcube int_img(
    Xdim + 1,
    Ydim + 1,
    Zdim + 1,
    arma::fill::zeros
  );

  // ---------------------------------------
  // Build integral image
  //
  // Standard 3D inclusion-exclusion:
  //
  //   I(x,y,z) = C(x-1,y-1,z-1)      <- voxel value (0-indexed)
  //     + I(x-1, y,   z  )            <- left face
  //     + I(x,   y-1, z  )            <- front face
  //     + I(x,   y,   z-1)            <- bottom face
  //     - I(x-1, y-1, z  )            <- left-front edge (double counted)
  //     - I(x-1, y,   z-1)            <- left-bottom edge (double counted)
  //     - I(x,   y-1, z-1)            <- front-bottom edge (double counted)
  //     + I(x-1, y-1, z-1)            <- origin corner (triple subtracted)
  // ---------------------------------------
  for (unsigned int z = 1; z <= Zdim; ++z) {
    for (unsigned int y = 1; y <= Ydim; ++y) {
      for (unsigned int x = 1; x <= Xdim; ++x) {

        int_img(x, y, z) =
          (double) C(x - 1, y - 1, z - 1)
          + int_img(x - 1, y,     z    )
          + int_img(x,     y - 1, z    )
          + int_img(x,     y,     z - 1)
          - int_img(x - 1, y - 1, z    )
          - int_img(x - 1, y,     z - 1)
          - int_img(x,     y - 1, z - 1)
          + int_img(x - 1, y - 1, z - 1);
      }
    }
  }

  // ---------------------------------------
  // Create output list
  // (one dataframe per box size)
  // ---------------------------------------
  Rcpp::List out(box_sizes.size());

  // ---------------------------------------
  // Loop through requested box sizes
  // ---------------------------------------
  for (int i = 0; i < box_sizes.size(); i++) {

    // Cast here, after the <= 0 check above
    unsigned int box_size = (unsigned int) box_sizes[i];

    // -------------------------------------
    // Validate box size against dimensions
    // -------------------------------------
    if (box_size > Xdim ||
        box_size > Ydim ||
        box_size > Zdim) {

      stop("box_size %d exceeds one or more voxel array dimensions (%d x %d x %d).",
           box_sizes[i], Xdim, Ydim, Zdim);
    }

    // -------------------------------------
    // Gliding-box position counts and
    // corresponding integral-image spans.
    //
    // For a dimension of size D and box of
    // size B, there are (D - B + 1) positions.
    // The integral image spans from index B
    // to D (inclusive), giving exactly that
    // many entries. The "low" span covers
    // 0 to (D - B), i.e. the top-left corner
    // of each box in the integral image.
    //
    // Named clearly to avoid confusion:
    //   n_pos_* = number of gliding positions
    //   lo_end_* = last index of low span (0..lo_end)
    //   hi_start_* = first index of high span (hi_start..dim)
    // -------------------------------------
    unsigned int n_pos_x  = Xdim - box_size;   // = D - B (span length is n_pos + 1)
    unsigned int n_pos_y  = Ydim - box_size;
    unsigned int n_pos_z  = Zdim - box_size;

    unsigned int hi_start_x = box_size;         // = B
    unsigned int hi_start_y = box_size;
    unsigned int hi_start_z = box_size;

    // -------------------------------------
    // Compute box-mass cube using
    // inclusion-exclusion on integral image.
    //
    // Each entry M(i,j,k) is the sum of the
    // voxels in the box whose top-left corner
    // is at (i, j, k) in the original array.
    // -------------------------------------
    arma::dcube M =
      // +++ far corner of box (x1,y1,z1)
      int_img(arma::span(hi_start_x, Xdim),
              arma::span(hi_start_y, Ydim),
              arma::span(hi_start_z, Zdim))
      // --- subtract the 3 faces outside the box
      - int_img(arma::span(0, n_pos_x),
                arma::span(hi_start_y, Ydim),
                arma::span(hi_start_z, Zdim))
      - int_img(arma::span(hi_start_x, Xdim),
                arma::span(0, n_pos_y),
                arma::span(hi_start_z, Zdim))
      - int_img(arma::span(hi_start_x, Xdim),
                arma::span(hi_start_y, Ydim),
                arma::span(0, n_pos_z))
      // +++ add back the 3 edges (double-subtracted)
      + int_img(arma::span(0, n_pos_x),
                arma::span(0, n_pos_y),
                arma::span(hi_start_z, Zdim))
      + int_img(arma::span(0, n_pos_x),
                arma::span(hi_start_y, Ydim),
                arma::span(0, n_pos_z))
      + int_img(arma::span(hi_start_x, Xdim),
                arma::span(0, n_pos_y),
                arma::span(0, n_pos_z))
      // --- subtract the origin corner (triple-added)
      - int_img(arma::span(0, n_pos_x),
                arma::span(0, n_pos_y),
                arma::span(0, n_pos_z));

    // -------------------------------------
    // Flatten box masses to 1D vector
    // and convert to unsigned integer
    // -------------------------------------
    arma::vec masses   = arma::vectorise(M);
    arma::uvec mass_int = arma::conv_to<arma::uvec>::from(masses);

    if (mass_int.n_elem == 0) {
      stop("Empty mass vector for box_size %d.", box_sizes[i]);
    }

    // -------------------------------------
    // Build histogram: mass value -> count.
    //
    // arma::histc bins by edge values
    // (replaces deprecated arma::hist).
    // Bins are 0, 1, ..., max_mass.
    // -------------------------------------
    arma::uword max_mass = mass_int.max();

    arma::uvec bins = arma::regspace<arma::uvec>(0, max_mass);

    // histc counts elements in [bins(k), bins(k+1))
    // with the last bin catching == bins.back().
    // For integer masses this is exact.
    arma::uvec freq = arma::histc(mass_int, bins);

    // -------------------------------------
    // Assemble output dataframe
    // -------------------------------------
    Rcpp::IntegerVector box_col(bins.n_elem, box_sizes[i]);

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
      Rcpp::Named("box_size")   = box_col,
      Rcpp::Named("mass")       = Rcpp::wrap(bins),
      Rcpp::Named("frequency")  = Rcpp::wrap(freq)
    );

    out[i] = df;
  }

  return out;
}