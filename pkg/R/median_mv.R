
#' Multivariate Median
#'
#' Computes the Multivariate Median of the distribution provided in `data`
#' (depth computations rely on the packages
#' \href{https://CRAN.R-project.org/package=depth}{`depth`},
#' \href{https://CRAN.R-project.org/package=OjaNP}{`OjaNP`} and
#' \href{https://CRAN.R-project.org/package=ddalpha}{`ddalpha`},
#' while convex body minizers "MCD" and "MVE" use the package
#' \href{https://CRAN.R-project.org/package=MASS}{`MASS`}).
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param method Method to use. Valid options are "CW", "MCD", "MVE", "L2",
#'   "Lui", "Mahalanobis", "Oja", "Projection" (default), "Spatial" and "Tukey"
#' @param sampleMedian Logical value. If `TRUE` (default), the function will
#'   return the Sample Median (observation from `data` with the highest depth).
#'   If `FALSE`, it will return the classic Multivariate Median (point in space
#'   with the highest depth), when applicable (methods "Oja", "Turkey" and
#'   "Spatial")
#' @param warnings Logical value, to display the warnings and error caught by
#'   the underlying functions
#'
#' @return Coordinates of the Multivariate Median
#'
#' @seealso [median_rec()], [depth_values()], [bril()]
#'
#' @examples
#'
#' # Illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(200, c(0, 0), diag(2)),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(100, c(5, 5), diag(2) * 200)
#' )
#'
#' # Compute median
#' m <- median_mv(XY, method = "L2", warnings = TRUE)
#'
#' # Plot results
#' plot(XY, asp = 1, xlab = "X", ylab = "Y")
#' points(m[1], m[2], col = "red", pch = 3, cex = 1.5, lwd = 3)
#'
#' ## Others examples of medians
#' \dontrun{
#' median_mv(XY, method = "Oja")
#' median_mv(XY, method = "Tukey", sampleMedian = TRUE)
#' median_mv(XY, method = "Tukey", sampleMedian = FALSE)
#' }
#'
#' @export
#'
median_mv <- function(data, method = "Projection", sampleMedian = TRUE, warnings = FALSE) {
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
    stop("Parameter `data` must contain exactly 2 columns, and more than one row")
  }
  if (!is.logical(sampleMedian)) {
    stop("Parameter `sampleMedian` must be a logical value")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  opts_method <- c("MCD", "MVE", "L2", "Liu", "Mahalanobis", "Oja", "Projection", "Spatial", "Tukey", "Zonoid", "Qhpeeling", "Betaskeleton", "Potential")
  if (is.character(method) && (tolower(method) %in% tolower(opts_method))) {
    method <- opts_method[which(tolower(method) == tolower(opts_method))]
  } else {
    stop("Parameter `method` should be one of \"", paste(opts_method, collapse = "\", "), "\"")
  }

  # Internal Parameters
  methodNsamp <- "best" # for MVE and MCD
  alpha <- 0.5

  # Exit if not enough observations to compute depth median, returning the coordinate-wise median instead
  minSampleSize <- ncol(data) + 1
  if (nrow(data) < minSampleSize) {
    if (warnings) {
      warning(sprintf(
        "median_mv(method=\"%s\"): sample size of `data``is %d but should be greater than %d, coordinate-wise median returned.",
        method, nrow(data), minSampleSize
      ))
    }
    return(c(stats::median(data[, 1]), stats::median(data[, 2])))
  }

  # Compute the median
  result <- withCallingHandlers(
    withRestarts({
      if (method %in% c("MCD", "MVE")) {
        robustScatter <- MASS::cov.rob(data,
                                       cor = FALSE,
                                       max(minSampleSize, floor(nrow(data) * (1 - alpha))),
                                       tolower(method), nsamp = methodNsamp
        )

        if (!exists("robustScatter") || !("center" %in% names(robustScatter)) || anyNA(robustScatter$center)) {
          result <- NULL
        } else {
          result <- robustScatter$center
        }
      } else if (method %in% c("CW")) {
        result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
      } else if (method == "Oja") {
        if (sampleMedian) {
          result <- OjaNP::ojaMedian(data, alg = "evolutionary")
        } else {
          depths <- depth_values(data, method = method, warnings = warnings)
          if (all(depths == 0)) {
            result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
          } else {
            result <- data[which.max(depths), ]
          }
        }
      } else if (method == "Tukey") {
        if (sampleMedian) {
          result <- depth::med(data, method = method, approx = TRUE)$median
        } else {
          depths <- depth_values(data, method = method, warnings = warnings)
          if (all(depths == 0)) {
            result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
          } else {
            result <- data[which.max(depths), ]
          }
        }
      } else if (method %in% c("Liu", "Spatial")) {
        if (sampleMedian) {
          result <- depth::med(data, method = method)$median
        } else {
          depths <- depth_values(data, method = method, warnings = warnings)
          if (all(depths == 0)) {
            result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
          } else {
            result <- data[which.max(depths), ]
          }
        }

        # Only sample medians for the methods below
      } else {
        depths <- depth_values(data, method = method, warnings = warnings)
        if (all(depths == 0)) {
          result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
        } else {
          result <- data[which.max(depths), ]
        }
      }
      result
    },
    muffleError = function() NULL
    ),
    warning = function(w) {
      if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
        warning("Warning caught in median_mv() with method=\"", method, "\" and ",
                nrow(data), " samples:\n* ", w,
                call. = FALSE
        )
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      if (warnings) {
        warning("Error caught in median_mv() with method=\"", method, "\" and ",
                nrow(data), " samples:\n* ", e,
                call. = FALSE
        )
      }
      invokeRestart("muffleError")
    }
  )

  # Check output
  if (is.null(result) || length(result) != 2 || anyNA(result) || !is.numeric(result)) {
    if (warnings) {
      warning("In median_mv() with method=\"", method, "\" and ", nrow(data),
              " samples: median estimator failed, coordinate-wise median returned instead.",
              call. = FALSE
      )
    }
    result <- c(stats::median(data[, 1]), stats::median(data[, 2]))
  }
  return(result)

}
