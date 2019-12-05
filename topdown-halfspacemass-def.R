#' Half-space mass
#'
#' Data depth is a statistical method which models data distribution
#' in terms of center outward ranking rather than density or linear ranking.
#'
#' A random subsample is projected onto a random direction t times.
#' For each projection, a split point <score> is randomly selected
#' between a range adjusted by <scope>;
#' then the fraction of points that fall in either side of <score> are recorded.
#' @param data d-dimensional numeric matrix of data points
#' @param n_halfspace number of directions in the d-dimensional space
#' @param subsample fraction of data to use for computing
#' @param scope size of the convex region covering the density of data
#' @param seed
#'
#' @return named list of directions, scores, mass_left and mass_right
#' @export
#'
#' @examples
train_depth <- function(data, n_halfspace, subsample = 1, scope = 1, seed) {
  data <- validate_data(data)
  validate_args(n_halfspace, subsample, scope)

  set.seed(seed)

  # n_halfspace x ncol(data)
  directions <- generate_directions(n_halfspace, ncol(data))
  # ceiling(subsample * nrow(data)) x n_halfspace
  scalar_projections <- compute_scalar_projections(data,
    directions,
    subsample_frac = subsample
  )
  # compute min, max columnwise
  min_vals <- apply(scalar_projections, 2, FUN = min)
  max_vals <- apply(scalar_projections, 2, FUN = max)

  mid_vals <- (min_vals + max_vals) * 0.5
  offsets <- (max_vals - min_vals) * scope * 0.5

  # map uniformly generated random numbers from [0, 1]
  # to [mid_vals - offsets, mid_vals + offsets]
  # using mins + (maxs - mins) * runif(n, 0, 1)
  scores <- (mid_vals - offsets) + 2 * offsets * runif(n_halfspace, 0, 1)

  # compare scalar projections of each direction with each score
  # and compute fraction of datapoints to the left of this direction
  mass_left <- rowMeans(apply(scalar_projections, 1, `<`, scores))
  # mass to the right is simply 1 - left fraction
  mass_right <- 1 - mass_left

  list(
    "directions" = directions,
    "scores" = scores,
    "mass_left" = mass_left,
    "mass_right" = mass_right
  )
}

#' Half-space mass
#'
#' Data depth is a statistical method which models data distribution
#' in terms of center outward ranking rather than density or linear ranking.
#'
#' Given query points, project each data point onto each of the directions,
#' and the number of training points that fall on the same
#' side as the query points are accumulated and output
#' as estimated value of the half-space mass (or depth) for the query points.
#'
#' @param data d-dimensional numeric matrix of data points
#' @param halfspaces named list of directions, scores, mass_left and mass_right
#' @param metric one of ("mass", "depth") for either mass or tukey depth
#'
#' @return computed mass or depth for each query point
#' @export
#'
#' @examples
evaluate_depth <- function(data, halfspaces, metric = c("mass", "depth")) {
  directions <- halfspaces[["directions"]]
  scores <- halfspaces[["scores"]]
  mass_left <- halfspaces[["mass_left"]]
  mass_right <- halfspaces[["mass_right"]]

  data <- validate_data(data)

  # make sure dimensions are equal between training and test data
  if (ncol(data) != ncol(directions)) {
    stop(sprintf(
      "Inconsistent dims between training data (%sD) and test data (%sD)",
      ncol(directions), ncol(data)
    ))
  }

  scalar_projections <- compute_scalar_projections(data, directions)

  # n_halfspace x nrow(data)
  left_idx <- apply(scalar_projections, 1, `<`, scores)

  depths <- numeric(nrow(data))

  # define aggregation function of left and right mass
  # sum of both sides for metric == "mass" (or simply not "depth")
  # or the least number of points in either side for metric == "depth"
  # if more than two options upgrade to switch
  agg_fun <- sum
  rescale_factor <- 1 / ncol(scalar_projections)
  if (metric == "depth") {
    agg_fun <- min
    # TODO: rescale_factor <- nrow(original_data)
    # no idea how to get that unless data is also passed as result
    # or encode it inside left/right mass...probably not worth it
    rescale_factor <- 1
  }

  for (i in seq_len(nrow(data))) {
    # aggregate each side and then aggregate both aggregations
    depths[i] <- agg_fun(
      agg_fun(mass_left[left_idx[, i]]),
      agg_fun(mass_right[!left_idx[, i]])
    )
  }

  depths * rescale_factor
}

########## Validation ##########################################################
validate_data <- function(data) {
  if (is.data.frame(data)) {
    # keep only numerical data
    data <- data[, unlist(lapply(data, is.numeric)), drop = FALSE]

    tryCatch(
      expr = {
        data <- data.matrix(data)
      },
      error = function(e) {
        stop("Cannot convert dataframe to matrix")
      }
    )
  }
  # TODO: should I also allow (one dimensional) vectors?
  checkmate::assert_matrix(data, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
  data
}

validate_args <- function(n_halfspace, subsample, scope) {
  checkmate::assert_count(n_halfspace, positive = TRUE)
  # actually 0 should not be allowed?
  checkmate::assert_number(subsample, lower = 0, upper = 1)
  checkmate::assert_number(scope, lower = 1, finite = TRUE)
}

########## Helper functions ####################################################
# normalise 1D vector or normalise matrix row-wise
normalise <- function(vecs) {
  # assume matrix
  sum_ <- rowSums
  # if vector use sum instead
  if (is.null(dim(vecs))) sum_ <- sum
  vecs / sqrt(sum_(vecs^2))
}

# generate random vectors that fall uniformly on the unit circle
generate_directions <- function(nrow, ncol) {
  rand_directions <- matrix(rnorm(nrow * ncol, 0, 1), ncol = ncol)
  normalise(rand_directions)
}

# compute scalar projections of a matrix (data) onto multiple directions
# if necessary, subsample first
compute_scalar_projections <- function(data, directions, subsample_frac = 1.0) {
  nrows <- nrow(data)
  subsample_size <- ceiling(nrows * subsample_frac)

  # if subsampling needed
  if (subsample_size < nrows) {
    sample_rows <- sample(nrows, size = subsample_size, replace = FALSE)
    return(data[sample_rows, , drop = FALSE] %*% t(directions))
  }

  data %*% t(directions)
}

# train_data <- data.frame(z1 = c(1, 2, 3), z2 = c(2, 3, 1), z3 = c(3, 2, 1))
# test_data <- data.frame(z1 = c(0), z2 = c(0), z3 = c(0))
# n_halfspace <- 1e3
# subsample <- 1
# seed <- 20191125
#
#
# halfspaces <- train_depth(train_data,
#   n_halfspace = n_halfspace,
#   subsample = subsample,
#   seed = seed
# )
#
# evaluate_depth(test_data, halfspaces, metric = "depth") * nrow(train_data)
# evaluate_depth(test_data, halfspaces, metric = "mass")


# bench::mark(
#   slow = train_depth(train_data,
#     n_halfspace = n_halfspace,
#     subsample = subsample,
#     seed = seed
#   ),
#   fast = train_depth_fast(train_data,
#     n_halfspace = n_halfspace,
#     subsample = subsample,
#     seed = seed
#   ),
#   relative = TRUE
# )[, 1:8]