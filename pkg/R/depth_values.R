#' Depth functions wrapper
#'
#' Computes the depth values with respect to the distribution provided in `data`
#' of either all the coordinates given in `u`, or all observations from `data`
#' if `u` is not provided. Depth computations rely on the packages
#' \href{https://CRAN.R-project.org/package=depth}{`depth`} and
#' \href{https://CRAN.R-project.org/package=ddalpha}{`ddalpha`}.
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param u Matrix of numerical values containing the coordinates for which
#'   depth values are to be computed (in rows, with two columns for X and Y
#'   coordinates). When missing or NULL, depth values will be computed for all
#'   observations from `data`
#' @param method Depth function used. Valid options are "L2", "Lui",
#'   "Mahalanobis", "Oja", "Projection" (default), "Spatial", "Tukey"
#' @param warnings Logical value, to display the warnings and error raised by
#'   the underlying depth functions
#'
#' @return Array of numerical values containing the depth of each observation
#'   from `data`, or from `u` if provided (these values are all set to 0 in the
#'   occurrence of errors)
#'
#' @seealso [bril()], [filter_outliers()], [median_rec()], [median_mv()]
#'
#' @examples
#'
#' # Illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2)),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' # Compute depths
#' D <- depth_values(XY, method = "L2", warnings = TRUE)
#'
#' # Plot distribution with depth color scale
#' plotColors <- colorRampPalette(c("maroon4", "steelblue4", "green4", "gold"))(20)
#' plot(XY, pch = 20, asp = 1, col = plotColors[as.numeric(cut(D, breaks = 20))],
#'      xlab = "X", ylab = "Y")
#'
#' # Plot depth values
#' plot(1:nrow(XY), D, pch = 20, col = plotColors[as.numeric(cut(D, breaks = 20))],
#'      xlab = "Index", ylab = "Depth")
#'
#' # Compute depth for a single point
#' depth_values(XY, c(10, 3), method = "L2")
#'
#' # Compute depth for three sets of coordinates
#' depth_values(XY, rbind(c(10, 3), c(65, 8), c(0, 1)), method = "L2")
#'
#' @export
#'
depth_values <- function(data, u = NULL, method = "Projection", warnings = FALSE) {
  if (missing(data)) {
    stop("Parameter `data` containing the observations is mandatory")
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Parameter `data` must contain numerical values only")
  }
  if (is.null(ncol(data)) || ncol(data) != 2) {
    stop("Parameter`data` must contain exactly 2 columns, and more than one row")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  if (is.character(u)) {
    method <- u
    u <- NULL
  }

  opts_method <- c("L2", "Liu", "Mahalanobis", "Oja", "Projection", "Spatial", "Tukey", "Zonoid", "Qhpeeling", "Betaskeleton", "Potential")
  if (is.character(method) && (tolower(method) %in% tolower(opts_method))) {
    method <- opts_method[which(tolower(method) == tolower(opts_method))]
  } else {
    stop("Parameter `method` should be one of \"", paste(opts_method, collapse = "\", "), "\"")
  }

  if (is.data.frame(u)) {
    data <- as.matrix(u)
  }
  if (is.null(u)) {
    u <- data
  } else {
    if (!is.numeric(u)) {
      stop("`u` must be `NULL` or contain numerical values only")
    } else {
      if (is.null(dim(u))) {
        u <- t(as.matrix(u))
      }
      # when one row only is provided
      if (ncol(u) != 2) {
        stop("`u`, when provided, must contain exactly 2 columns")
      }
    }
  }

  # Exit if not enough observations to compute depth values
  minSampleSize <- ncol(data) + 1
  if (nrow(data) < minSampleSize) {
    if (warnings) {
      warning(sprintf(
        "depth_values(method=\"%s\"): sample size of `data``is %d but should be greater than %d, zeros returned",
        method, nrow(data), minSampleSize
      ))
    }
    return(rep(0, dim(u)[1]))
  }

  # Set the wrapper for each depth function
  switch(method,
         Oja = {
           depthFun <- function(x, y) {
             apply(x, 1, function(x, y) {
               # faster than ddalpha
               depth::depth(x, y, method = "Oja")
             }, y)
           }
         },
         Tukey = {
           depthFun <- function(x, y) {
             ddalpha::depth.halfspace(x, y, exact = TRUE)
           }
         },
         Liu = {
           depthFun <- function(x, y) {
             ddalpha::depth.simplicial(x, y)
           }
         },
         Projection = {
           depthFun <- function(x, y) {
             ddalpha::depth.projection(x, y)
           }
         },
         L2 = {
           depthFun <- function(x, y) {
             ddalpha::depth.L2(x, y)
           }
         },
         Spatial = {
           depthFun <- function(x, y) {
             ddalpha::depth.spatial(x, y)
           }
         },
         Mahalanobis = {
           depthFun <- function(x, y) {
             ddalpha::depth.Mahalanobis(x, y)
           }
         },
         Betaskeleton = {
           depthFun <- function(x, y) {
             ddalpha::depth.betaSkeleton(x, y)
           }
         },
         Zonoid = {
           depthFun <- function(x, y) {
             ddalpha::depth.zonoid(x, y)
           }
         },
         Qhpeeling = {
           depthFun <- function(x, y) {
             ddalpha::depth.qhpeeling(x, y)
           }
         },
         Potential = {
           depthFun <- function(x, y) {
             ddalpha::depth.potential(x, y)
           }
         },
         stop(sprintf("depth_values() - unsupported value for parameter: method = \"%s\" ", method))
  )

  # Apply the depth function, catching the exception in case of too few samples or singular matrix
  depthValues <- withCallingHandlers(
    withRestarts({
      depthFun(u, data)
    }, muffleError = function() NULL),
    warning = function(w) {
      if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
        warning("Warning caught in depth_values() with method=\"", method, "\" and ",
                nrow(data), " samples:\n* ", w, call. = FALSE)
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      if (warnings) {
        warning("Error caught in depth_values() with method=\"", method, "\" and ",
                nrow(data), " samples:\n* ", e, call. = FALSE)
      }
      invokeRestart("muffleError")
    }
  )

  # Check result obtained and return
  if (!exists("depthValues") || is.null(depthValues) || anyNA(depthValues)) {
    if (warnings) {
      warning("In depth_values() with method=\"", method, "\" and ", nrow(data),
              " samples: depth estimator failed, zeros returned", call. = FALSE)
    }
    return(rep(0, dim(u)[1]))
  } else {
    return(depthValues)
  }
}
