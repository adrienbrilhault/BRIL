
#' Recursive estimate of central location
#'
#' Recursive estimate of central location based on depth measures (from the
#' packages
#' \href{https://cran.r-project.org/web/packages/depth/index.html}{`depth`} and
#' \href{https://cran.r-project.org/web/packages/ddalpha/index.html}{`ddalpha`})
#' or convex body minimizers (package
#' \href{https://cran.r-project.org/web/packages/MASS/index.html}{`MASS`}).
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param method Method to use. Valid options are "MCD" and "MVE" for convex
#'   body minimizers, or "L2", "Lui", "Mahalanobis", "Oja", "Projection"
#'   (default), "Spatial" and "Tukey" for depth functions
#' @param alpha Proportion of samples trimmed at each iteration (numerical value
#'   between 0 and 1, default: 0.5)
#' @param maxIterations Set to a positive integer to limit the number of
#'   iterations, to NULL or 0 (default) for no limits.
#' @param warnings Logical value, to display the warnings and error raised by
#'   the underlying functions
#'
#' @return The function returns an S3 object of type `BRIL.MedianRec`,
#'   containing the following values:
#'   \item{`median`}{Coordinate of the recursive median}
#'   \item{`max`}{Coordinate of the sample with the highest depth (or the
#'   center of the first iteration in the case of convex body minimizers)}
#'   \item{`iterations`}{List containing the indices from the samples of `data`
#'   selected at each iterations}
#'
#' @seealso [plot.BRIL.MedianRec()], [print.BRIL.MedianRec()], [median_mv()], [depth_values()], [bril()]
#'
#' @examples
#'
#' # illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2)),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' # Compute the recursive median
#' res <- median_rec(XY)
#'
#' print(res)
#' plot(res)
#'
#' @export
#'
median_rec <- function(data, method = "Projection", alpha = 0.5, maxIterations = NULL, warnings = FALSE) {

  ######## Check parameters & initialize

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
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Parameter `alpha` must be a numerical value between 0 and 1")
  }
  if (!is.null(maxIterations) && maxIterations < 1) {
    stop("Parameter `maxIterations` must be either NULL or an integer >= 1")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  opts_method <- c(
    "MCD", "MVE", "L2", "Liu", "Mahalanobis", "Oja", "Projection", "Spatial",
    "Tukey", "Zonoid", "Qhpeeling", "Betaskeleton", "Potential"
  )
  if (is.character(method) && (tolower(method) %in% tolower(opts_method))) {
    method <- opts_method[which(tolower(method) == tolower(opts_method))]
  } else {
    stop("Parameter `method` should be one of \"", paste(opts_method, collapse = "\", \""), "\"")
  }


  # Internal Parameters
  methodNsamp <- "best" # for MVE and MCD

  # Output of the function
  output <- list()
  class(output) <- "BRIL.MedianRec"
  output$call$data <- data
  output$call$method <- method
  output$call$alpha <- alpha
  output$call$maxIterations <- maxIterations
  output$median <- c(stats::median(data[, 1]), stats::median(data[, 2]))
  output$max <- output$median
  output$iterations <- list(1:nrow(data))

  nbPoints <- nrow(data)

  # Return the coordinate-wise median if insufficient number of samples
  minSampleSize <- ncol(data) + 1

  if (nrow(data) < minSampleSize) {
    if (warnings) {
      warning(sprintf("median_rec(%s): Sample Size of `data` is %d but should be greater than %d, Coordinate-Wise median returned", method, nrow(data), minSampleSize))
    }
    return(output)
  }

  iterationNB <- 1

  if (method %in% c("MCD", "MVE")) {
    while (!is.null(nbPoints) && nbPoints > minSampleSize && (is.null(maxIterations) || iterationNB <= maxIterations)) {
      robustScatter <- withCallingHandlers(
        withRestarts(
          {
            robustScatter <- MASS::cov.rob(data,
              cor = FALSE,
              max(minSampleSize, floor(nbPoints * (1 - alpha))),
              tolower(method), nsamp = methodNsamp
            )

            # Check values are correct if it restarts from a warning
            if (!exists("robustScatter") || !("center" %in% names(robustScatter)) ||
              !("best" %in% names(robustScatter)) || anyNA(robustScatter$center) ||
              anyNA(robustScatter$best)) {
              robustScatter <- NULL
              if (warnings) {
                warning("In median_rec() with method=\"", method, "\" and ", nrow(data),
                  " samples: robust estimator failed",
                  call. = FALSE
                )
              }
            }
            robustScatter
          },
          muffleError = function() NULL
        ),
        warning = function(w) {
          if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
            warning("Warning caught in median_rec() with method=\"", method, "\" and ",
              nrow(data), " samples:\n* ", w,
              call. = FALSE
            )
          }
          invokeRestart("muffleWarning")
        },
        error = function(e) {
          if (warnings) {
            warning("Error caught in median_rec() with method=\"", method, "\" and ",
              nrow(data), " samples:\n* ", e,
              call. = FALSE
            )
          }

          invokeRestart("muffleError")
        }
      )

      if (is.null(robustScatter)) {
        if (warnings) {
          warning("Exiting median_rec() after ", iterationNB, "iterations (", nbPoints, " samples left)")
        }
        nbPoints <- 0
        break
      }

      # Use the center of the first iteration as `max`
      if (iterationNB == 1) {
        output$max <- robustScatter$center
      }

      indRec <- robustScatter$best
      output$iterations[[iterationNB + 1]] <- output$iterations[[iterationNB]][indRec]
      output$median <- robustScatter$center
      data <- data[indRec, ]

      if (is.null(nrow(data))) {
        break
      } else {
        nbPoints <- nrow(data)
        iterationNB <- iterationNB + 1
      }
    }
  } else {
    while (!is.null(nbPoints) && nbPoints > minSampleSize &&
      (is.null(maxIterations) || iterationNB <= maxIterations)) {
      depthValues <- depth_values(data, method = method, warnings = warnings)
      if (all(depthValues == 0)) {
        if (warnings) {
          warning(
            "Exiting median_rec() after", iterationNB,
            "iterations (", nbPoints, " samples left)"
          )
        }

        break
      }

      indRec <- sort.int(depthValues, decreasing = TRUE, index.return = TRUE)$ix

      # Return the deepest sample of the first iteration as `max`
      if (iterationNB == 1) {
        output$max <- data[indRec[1], ]
      }

      output$iterations[[iterationNB + 1]] <- output$iterations[[iterationNB]][
        indRec[1:max(minSampleSize, floor(nbPoints * (1 - alpha)))]
      ]
      data <- data[indRec[1:max(minSampleSize, floor(nbPoints * (1 - alpha)))], ]

      if (is.null(nrow(data))) {
        output$median <- data
        break
      } else {
        nbPoints <- nrow(data)
        output$median <- colMeans(data)
        iterationNB <- iterationNB + 1
      }
    }
  }

  # check our output is valid
  if (!is.numeric(output$median) || length(output$median) != 2 || anyNA(output$median)) {
    if (warnings) {
      warning(
        "Issue in median_rec(method=\"", method,
        "\"): `median` is not numeric, coordinate-wise median returned instead"
      )
    }
    output$median <- c(
      stats::median(output$call$data[, 1]),
      stats::median(output$call$data[, 2])
    )
  }
  if (!is.numeric(output$max) || length(output$max) != 2 || anyNA(output$max)) {
    if (warnings) {
      warning(
        "Issue in median_rec(method=\"", method,
        "\"): `median` is not numeric, coordinate-wise median returned instead"
      )
    }
    output$max <- c(
      stats::median(output$call$data[, 1]),
      stats::median(output$call$data[, 2])
    )
  }

  return(output)
}


#' Print method for `BRIL.MedianRec` objects
#'
#' @param x An object of class `BRIL.MedianRec` (see [median_rec()])
#' @param ... Other arguments passed to or from other methods
#'
#' @seealso [median_rec()], [plot.BRIL.MedianRec()]
#'
#' @export
#'
print.BRIL.MedianRec <- function(x, ...) {
  if (class(x) != "BRIL.MedianRec") {
    stop("the object provided is not of class \"BRILfiltering\"")
  }

  cat(sprintf(
    "\n=> Results for median_rec() using \"%s\" method (alpha=%.2f)\n",
    x$call$method, x$call$alpha
  ))
  cat(sprintf("  (%d samples, %d iterations)\n\n", nrow(x$call$data), length(x$iterations)))

  cat("Recursive median (median):\n")
  print(x$median, ...)
  cat("\nSample median (max):\n")
  print(x$max, ...)
  cat("\n")
  invisible(x)
}



#' Plot method for `BRIL.MedianRec` objects
#'
#' @param x An object of class `BRIL.MedianRec` (see [median_rec()])
#' @param nbIterations Number of iterations to display, or 0 to show all of them
#'   (default: 5)
#' @param showMedian Logical, to show the final recursive median (indicated by a
#'   "+")
#' @param showMax Logical, to show the overall deepest point, or the center of
#'   the first MCD/MVE iteration (indicated by a "x")
#' @param ... Other arguments passed to or from other methods
#'
#' @examples
#'
#' # illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2)),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' # Process the data
#' res <- median_rec(XY)
#'
#' # Default plot
#' plot(res)
#'
#' # Adjust axis
#' plot(res, asp = 1, xlim = c(-20, 20), ylim = c(-20, 20))
#'
#' # Change other graphical options
#' plot(res, showMedian = TRUE, pch = 16, main = "Recursive Median")
#' @seealso [median_rec()], [median_mv()], [bril()]
#'
#' @export
#'
plot.BRIL.MedianRec <- function(x, nbIterations = 5, showMedian = FALSE, showMax = FALSE, ...) {
  if (!is.logical(showMedian) || !is.logical(showMax)) {
    stop("showMedian and showMax must be logical value")
  }
  if (class(x) != "BRIL.MedianRec") {
    stop("the object provided is not of class \"BRIL.MedianRec\"")
  }
  if (!is.numeric(nbIterations) || nbIterations < 1) {
    stop("nbIteration must be a positive integer")
  }


  # Read extra parameters
  otherArgs <- list(...)

  do.call(plot, utils::modifyList(
    list(
      x = x$call$data, ylab = "X", xlab = "Y", pch = 1,
      col = "grey60", cex = 1, lwd = 1
    ),
    otherArgs[names(otherArgs) %in% c("col") == FALSE]
  ))

  myColours <- rev(grDevices::rainbow(min(nbIterations, length(x$iterations))))

  for (i in 2:min(nbIterations + 1, length(x$iterations))) {
    do.call(graphics::points, utils::modifyList(
      list(x$call$data[x$iterations[[i]], 1],
        x$call$data[x$iterations[[i]], 2],
        col = myColours[i - 1], pch = 1
      ),
      otherArgs[names(otherArgs) %in% c("col") == FALSE]
    ))
  }
  if (showMedian) {
    graphics::points(x$median[1], x$median[2], col = "black", pch = 3, cex = 1.5, lwd = 3)
  }
  if (showMax) {
    graphics::points(x$max[1], x$max[2], col = "black", pch = 4, cex = 1.5, lwd = 3)
  }
}
