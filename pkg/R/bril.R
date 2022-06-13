
suppressPackageStartupMessages(require("CompQuadForm", quietly = TRUE))

#' Bootstrap and Refine Iterative Locator
#'
#' Robust Estimate of Mode in Multivariate Distribution
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param maxIterations Maximum number of iterations performed by the algorithm
#'   (set to NULL or 0 for unlimited number)
#' @param minUnassigned Numerical value between 0 and 1 (default: 0.1),
#'   providing the proportion of unassigned samples from `data` below which the
#'   algorithm will terminate
#' @param method Method to use. Valid options are "MCD" and "MVE" for convex
#'   body minimizers, or "L2", "Lui", "Mahalanobis", "Oja", "Projection"
#'   (default), "Spatial" and "Tukey" for depth functions
#' @param alpha Proportion of samples trimmed at each iteration of the recursive
#'   median estimate (numerical value between 0 and 1, default: 0.5), see
#'   [median_rec()]
#' @param testUnimodal Statistical test used for unimodality. Valid options are
#'   "DIP" only, see [filter_outliers()]
#' @param threshUnimodal Threshold of significance for the unimodality test
#'   (numerical value between 0 and 1, default: 0.05)
#' @param distUnimodal Distance metric used for ordering the samples in the
#'   unimodal filtering. Valid options are "Euclidean" (default), or "MCD",
#'   "MVE", and "OGK" for robust distances. "Euclidean" is strongly advised for
#'   unimodality tests.
#' @param testNormal Statistical test used for normality. Valid options are
#'   "Mardia" (default), "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro",
#'   "Lillie" and "Chisq", see [filter_outliers()]
#' @param threshNormal Threshold of significance for the normality test
#'   (numerical value between 0 and 1, default: 0.05)
#' @param distNormal Distance metric used for ordering the samples in the normal
#'   filtering. Valid options are "Euclidean", or "MCD" (default), "MVE", and
#'   "OGK" for robust distances. Robust distances are strongly advised for
#'   normality tests.
#' @param trimmedPerFilteringIteration Number of samples trimmed at each
#'   iteration of the unimodality and normality filtering (default: 1), see
#'   [filter_outliers()]
#' @param exitWhenUnimodal Logical value. `TRUE` will terminate the execution of
#'   the algorithm as soon as an unimodal subset is encountered on the start of
#'   a global iteration. `FALSE` (default) will let that last iteration proceed
#'   before terminating
#' @param debug Logical value. `TRUE` will compute all p.values in the filtering
#'   steps (even after they exceed the selection threshold, see
#'   [plot.BRIL.Filtering()])
#' @param warnings Logical value, to display the warnings and errors caught
#'
#' @return The function returns an S3 object of type `BRIL` containing the
#'   following values:
#'   \item{`call`}{Parameters of the call (contains `data`,
#'   `maxIterations`, `minUnassigned`, `method`, `alpha`, `testUnimodal`,
#'   `threshUnimodal`, `distUnimodal`, `testNormal`, `threshNormal`,
#'   `distNormal`, `trimmedPerFilteringIteration`, and `exitWhenUnimodal`)}
#'   \item{`iterations`}{A list with every global iteration of the algorithm,
#'   each containing the two filtering procedures performed: `filteringUnimodal`
#'   and `filteringNormal` (both being S3 object of class `BRIL.Filtering`, see
#'   [filter_outliers()])}
#'   \item{`nbClusters`}{Number of groups encountered}
#'   \item{`labels`}{Labels of the groups encountered (corresponding to the
#'   number of the iteration they were identified in)}
#'   \item{`clustersCenters`}{Matrix containing the coordinates of the centers
#'   of each group (row-wise)}
#'   \item{`clustersSizes`}{Array with the number of samples in each group}
#'   \item{`mainCluster`}{Index of the group identified as main mode}
#'   \item{`mode`}{Coordinates of the main mode}
#'
#' @references Adrien Brilhault, Sergio Neuenschwander, and Ricardo Rios - A New
#'   Robust Multivariate Mode Estimator for Eye-tracking Calibration - Behavior
#'   Research Methods, 2022 - \href{https://rdcu.be/cI9Pf}{rdcu.be/cI9Pf}
#'
#' @seealso [plot.BRIL()], [print.BRIL()], [filter_outliers()], [median_rec()],
#'   [median_mv()], [depth_values()]
#'
#' @examples
#'
#' # Create a sample distribution and run bril() function
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' res <- bril(XY, debug = TRUE)
#' print(res)
#'
#' # Plot the mode and groups encountered
#' plot(res)
#'
#' # Plot the different iterations (interactive)
#' \dontrun{
#' plot(res, contents = "iterations", asp = 1)
#' }
#'
#' # See ?plot.BRIL() for other plotting examples
#'
#' @export
#'
bril <- function(data, maxIterations = NULL, minUnassigned = 0.1, method = "Projection", alpha = 0.5,
                 testUnimodal = "DIP", threshUnimodal = 0.05, distUnimodal = "Euclidean",
                 testNormal = "Mardia", threshNormal = 0.05, distNormal = "MCD",
                 trimmedPerFilteringIteration = 1, exitWhenUnimodal = FALSE,
                 debug = FALSE, warnings = FALSE) {

  DEV_DEBUG <- FALSE

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
    stop("Parameter`data` must contain exactly 2 columns, and more than one row")
  }
  if (!is.null(maxIterations) && (!is.numeric(maxIterations) || maxIterations < 0)) {
    stop("Parameter `maxIterations` must be either NULL or a positive entire number")
  }
  if (!is.numeric(minUnassigned) || minUnassigned < 0 || minUnassigned > 1) {
    stop("Parameter `minUnassigned` must be a numerical value between 0 and 1 (included)")
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Parameter `alpha` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(threshUnimodal) || threshUnimodal < 0 || threshUnimodal > 1) {
    stop("Parameter `threshUnimodal` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(threshNormal) || threshNormal < 0 || threshNormal > 1) {
    stop("Parameter `threshNormal` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(trimmedPerFilteringIteration) || trimmedPerFilteringIteration < 1) {
    stop("Parameter `trimmedPerFilteringIteration` must be an entire number superior or equal to 1")
  }
  trimmedPerFilteringIteration <- round(trimmedPerFilteringIteration)
  if (!is.logical(debug)) {
    stop("Parameter `debug` must be a logical value")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  opts_test <- c("DIP", "Mardia", "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro", "Lillie", "Chisq")
  if (is.character(testUnimodal) && (tolower(testUnimodal) %in% tolower(opts_test))) {
    testUnimodal <- opts_test[which(tolower(testUnimodal) == tolower(opts_test))]
  } else {
    stop("Parameter `testUnimodal` should be one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  if (is.character(testNormal) && (tolower(testNormal) %in% tolower(opts_test))) {
    testNormal <- opts_test[which(tolower(testNormal) == tolower(opts_test))]
  } else {
    stop("Parameter `testNormal` should be one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  opts_dist <- c("Euclidean", "MCD", "MVE", "OGK")
  if (missing(distUnimodal) || is.null(distUnimodal)) {
    distUnimodal <- "Euclidean"
  } else {
    if (is.character(distUnimodal) && tolower(distUnimodal) %in% tolower(opts_dist)) {
      distUnimodal <- opts_dist[which(tolower(distUnimodal) == tolower(opts_dist))]
    } else {
      stop("Parameter `distUnimodal` should be NULL or one of \"", paste(opts_dist, collapse = "\", \""), "\"")
    }
  }

  if (missing(distNormal) || is.null(distNormal)) {
    distNormal <- "MCD"
  } else {
    if (is.character(distNormal) && tolower(distNormal) %in% tolower(opts_dist)) {
      distNormal <- opts_dist[which(tolower(distNormal) == tolower(opts_dist))]
    } else {
      stop("Parameter `distNormal` should be NULL or one of \"", paste(opts_dist, collapse = "\", \""), "\"")
    }
  }

  if (missing(maxIterations) || is.null(maxIterations)) {
    maxIterations <- 0
  } else {
    maxIterations <- round(maxIterations)
  }

  ## Internal Parameter
  reorderClustersBySize <- FALSE
  recenterAfterUnimodal <- TRUE
  exitWhenUnimodal <- FALSE

  # Output of the function
  output <- list()
  class(output) <- "BRIL"
  output$call$data <- data
  output$call$maxIterations <- maxIterations
  output$call$minUnassigned <- minUnassigned
  output$call$method <- method
  output$call$alpha <- alpha
  output$call$testUnimodal <- testUnimodal
  output$call$threshUnimodal <- threshUnimodal
  output$call$distUnimodal <- distUnimodal
  output$call$testNormal <- testNormal
  output$call$threshNormal <- threshNormal
  output$call$distNormal <- distNormal
  output$call$trimmedPerFilteringIteration <- trimmedPerFilteringIteration
  output$call$exitWhenUnimodal <- exitWhenUnimodal

  # Exit if not enough observations to compute depth values
  minSampleSize <- ncol(data) + 1
  if (nrow(data) < minSampleSize) {
    if (warnings) {
      warning(sprintf(
        "depth_values(method=\"%s\"): sample size of `data``is %d but should be greater than %d, coordinate-wise median returned",
        method, nrow(data), minSampleSize
      ))
    }
    output$mode <- c(stats::median(data[, 1]), stats::median(data[, 2]))
    return(output)
  }

  # Initialization
  iterationNB <- 1
  output$iterations <- list()
  output$nbClusters <- 0
  output$clustersCenters <- c()
  output$clustersSizes <- c()
  output$labels <- rep(0, dim(data)[1])

  ######## Main loop, selecting a new group at each iteration

  while ((maxIterations == 0 || iterationNB <= maxIterations) &&
         (sum(output$labels == 0) > max(minSampleSize, minUnassigned * nrow(data)))) {
    unassignedIndices <- which(output$labels == 0)

    if (DEV_DEBUG)
      cat("BRIL - Iteration",iterationNB, "(nb samples:",nrow(output$call$data),"- unprocessed:",length(unassignedIndices),"\n")

    ## Bootstrap : compute first estimate

    bootstrapEstimate <- median_rec(data[unassignedIndices, ], method = method, alpha = alpha, warnings = warnings)$median

    if (DEV_DEBUG)
      cat("* Bootstrap estimate (median_rec): ",bootstrapEstimate,"\n")

    ## Refine estimate

    # Filter furthest outliers until reaching an unimodal subset
    filteringUnimodal <- filter_outliers(
      data[unassignedIndices, ],
      center = bootstrapEstimate, test = testUnimodal,
      threshold = threshUnimodal, distType = distUnimodal, debug = debug, warnings = warnings
    )
    output$iterations[[iterationNB]] <- list("filteringUnimodal" = filteringUnimodal)
    unimodalIndices <- unassignedIndices[filteringUnimodal$selected]

    if (DEV_DEBUG)
      cat("* Filtering unimodal:", length(filteringUnimodal$selected), "selected","\n")

    # If the current distribution, before any filtering, was already unimodal,
    # terminate the overall recursive search now or at the end of this iteration
    if (utils::head(filteringUnimodal$p.values, 1) > threshUnimodal) {
      if (exitWhenUnimodal == TRUE && iterationNB != 1) {
        maxIterations <- iterationNB
        break
      } else {
        maxIterations <- iterationNB
      }
    }

    # Re-estimate central location within this subset
    if (recenterAfterUnimodal == TRUE) {
      bootstrapEstimate <- median_rec(
        data[unimodalIndices, ], method = method, alpha = alpha, warnings = warnings)$median
    }

    if (DEV_DEBUG)
      cat("* Adjusted center:", bootstrapEstimate,"\n")

    # Filter furthest outliers until reaching a normal subset
    filteringNormal <- filter_outliers(
      data[unimodalIndices, ],
      center = bootstrapEstimate, test = testNormal,
      threshold = threshNormal, distType = distNormal, debug = debug, warnings = warnings
    )
    output$iterations[[iterationNB]]$filteringNormal <- filteringNormal
    normalIndices <- unimodalIndices[filteringNormal$selected]

    # Estimate the center of this group (refined estimate)
    center <- colMeans(data[normalIndices, ])

    if (DEV_DEBUG) {
      cat("* Filtering normal:", length(filteringNormal$selected), "selected","\n")
      cat("* Cluster final center:", center,"\n")
    }

    # Update the groups structures
    output$labels[normalIndices] <- iterationNB
    output$nbClusters <- output$nbClusters + 1
    output$clustersCenters <- rbind(output$clustersCenters, center)
    output$clustersSizes <- c(output$clustersSizes, length(normalIndices))
    rownames(output$clustersCenters) <- NULL

    iterationNB <- iterationNB + 1

    # # If the size of the remaining samples is inferior to the largest
    # # group encountered so far, it will not be the main mode, we can exit
    # if (sum(output$labels == 0) < max(output$clustersSizes))
    #   break
  }

  ######## Select main mode and return

  if (output$nbClusters == 1) {
    output$mainCluster <- 1
    output$mode <- output$clustersCenters
  } else {
    output$mainCluster <- which.max(utils::head(output$clustersSizes, -1))
    output$mode <- output$clustersCenters[output$mainCluster, ]
  }

  return(output)
}


#' Print method for `BRIL` objects
#'
#' @param x An object of class `BRIL` (see [bril()])
#' @param maxDisplayed Number of elements to display in the output. Set to NULL
#'   (or 0) to show all values.
#' @param ... Other arguments passed to or from other methods
#'
#' @seealso [bril()], [plot.BRIL()]
#'
#' @export
#'
print.BRIL <- function(x, maxDisplayed = NULL, ...) {
  if (class(x) != "BRIL") {
    stop("the object provided is not of class \"BRIL\"")
  }

  if (is.null(maxDisplayed) || maxDisplayed == 0) {
    maxDisplayed <- length(x$labels)
  } else if (!is.numeric(maxDisplayed) || maxDisplayed < 0) {
    stop("Parameter `maxDisplayed` must be a numerical value superior or equal to 1")
  }
  maxDisplayed <- round(maxDisplayed)

  cat("\n=> Results for bril() using method \"", x$call$method, "\" (alpha=", x$call$alpha,
      "), ", x$call$testUnimodal, " Unimodality Test (> ", x$call$threshUnimodal,
      "), and ", x$call$testNormal, " Normality test (> ", x$call$threshNormal, ")\n",
      sep = ""
  )

  cat("   ", nrow(x$call$data), " samples: ", x$nbClusters, " clusters identified (sizes ",
      toString(x$clustersSize), ")",
      sep = ""
  )
  if (sum(x$labels == 0) > 0) {
    cat(", and ", sum(x$labels == 0), " samples unassigned", sep = "")
  }

  cat("\n\nMode:\n")
  print(x$mode, ...)

  cat("\n\nClusters Sizes:\n")
  print(x$clustersSizes, ...)

  cat("\n\nClusters Centers:\n")
  print(x$clustersCenters, ...)

  cat("\n\nLabels:\n")
  print(utils::head(x$labels, maxDisplayed), ...)
  if (length(x$labels) > maxDisplayed) {
    cat(paste0(" ... (", length(x$labels) - maxDisplayed, " hidden)"))
  }

  cat("\n\n")

  invisible(x)
}


#' Plot method for `BRIL` objects
#'
#' @param x An object of class `BRIL` (see [bril()])
#' @param contents Contents to be displayed, options are "scatterplot",
#'   or "iterations" (only one option possible)
#' @param showClusters Logical value used when `contents = "scatterplot"`, to
#'   show or not the different clusters
#' @param showMode Logical value used when `contents = "scatterplot"`, to show
#'   or not the main mode
#' @param col Default color of samples
#' @param colMode Color of the mode when `contents = "scatterplot"`
#' @param colClusters List (or array) of colors for each of the
#'   clusters/iterations (length must be at least equal to the number of
#'   groups identified by the function [bril()], i.e. `x$nbClusters`)
#' @param iterationsIndices Numerical value or array of numerical values, used
#'   when `contents = "iterations"`, which provides the indices of the
#'   iterations to be plotted. If more than one iteration is requested, an
#'   interactive menu in the console will be used for the selection. 0 or NULL
#'   (default) will include all the iterations. Values that are negative or
#'   superior to the number of iterations performed by the execution [bril()]
#'   will be ignored
#' @param iterationsOptions List of additional parameters to be passed to the
#'   [plot.BRIL.Filtering()] function when `contents = "iterations"` is selected
#'   (see [plot.BRIL.Filtering()] for details). Example: `iterationsOptions =
#'   list(xlab = NA, ylab = NA, contents = c("p.values", "scatterplot"), asp =
#'   1)`
#' @param ... Other arguments passed to or from other methods (such as pch for
#'   the symbols, main and sub for title and subtitle, xlab, xmin, ...)
#'
#' @seealso [bril()], [print.BRIL()], [filter_outliers()],
#'   [print.BRIL.Filtering()]
#'
#' @examples
#'
#' # Create a sample distribution and run bril() function
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#' res <- bril(XY, debug = TRUE)
#'
#' # Plot the mode and groups encountered (default)
#' plot(res)
#'
#' # Plot the mode only
#' plot(res, showClusters = FALSE)
#'
#' # Plot the mode only (with extra graphic options)
#' plot(res, showClusters = FALSE, main = "Multivariate Mode Estimate",
#'      col = "blue", colMode = "black", asp = 1, pch = 3 )
#'
#' # Plot the clusters without the mode
#' plot(res, showMode = FALSE, col = "gray",
#'      colClusters = c("yellow","cyan","purple","red"))
#'
#' # Plot the second iteration
#' plot(res, contents = "iteration", iterationsIndices = 2)
#'
#' # Plot the second iteration (with arguments to plot.BRIL.filtering())
#' plot(res, contents = "iteration", iterationsIndices = 2,
#'      iterationsOptions = list(
#'        contents = c("scatterplot", "p.values"),
#'        colSelection = "blue", mfrow = c(2,1), asp = 1))
#'
#' \dontrun{
#' # Plot all iterations (interactive mode)
#' plot(res, contents = "iterations")
#'
#' # Plot the 3 first iterations with options (interactive mode)
#' plot(res, contents = "iterations", iterationsIndices = c(1:3),
#'      iterationsOptions = list(
#'        contents = c("scatterplot"),
#'        xlim  = c(-50,50),  ylim  = c(-30,30), asp = 1))
#' }
#'
#' @importFrom graphics plot points
#' @export
#'
plot.BRIL <- function(x, contents = "plot", showClusters = TRUE, showMode = TRUE,
                      col, colMode, colClusters = NULL, iterationsIndices = NULL,
                      iterationsOptions = NULL, ...) {
  if (class(x) != "BRIL") {
    stop("the object provided is not of class \"BRIL\"")
  }

  # Read extra parameters
  otherArgs <- list(...)

  if (!is.null(iterationsOptions) && !is.list(iterationsOptions)) {
    warning("Parameter `iterationsOptions`, if provided, must be a list")
    iterationsOptions <- NULL
  }
  if (length(contents) > 1) {
    warning("Parameter `contents` can only have one value (among \"scatterplot\" and \"iterations\"")
    contents <- contents[1]
  }
  if (!is.character(contents)) {
    stop("Parameter `contents` must be a string (either \"scatterplot\" or \"iterations\"")
  }
  contents <- tolower(contents)

  if (!is.logical(showClusters)) {
    stop("Parameter `showClusters` must be a logical value")
  }
  if (!is.logical(showMode)) {
    stop("Parameter `showMode` must be a logical value")
  }
  if (is.null(iterationsIndices) || anyNA(iterationsIndices) || iterationsIndices == 0) {
    iterationsIndices <- 1:x$nbClusters
  }
  if (!is.numeric(iterationsIndices)) {
    stop("Parameter `iterationsIndices`, if provided, must be numerical value(s)")
  }
  iterationsIndices <- sort(unique(iterationsIndices))
  if (max(iterationsIndices) > x$nbClusters) {
    warning("Only ", x$nbClusters, " iterations, values above discarded")
    iterationsIndices <- iterationsIndices[iterationsIndices <= x$nbClusters]
  }
  if (min(iterationsIndices) < 1) {
    warning("Parameter `iterations` can not take negative values")
    iterationsIndices <- iterationsIndices[iterationsIndices > 0]
  }
  if (missing(col) || is.null(col)) {
    col <- "black"
  }
  if (missing(colMode) || is.null(colMode)) {
    if (showClusters == TRUE) {
      colMode <- "black"
    } else {
      colMode <- "orange"
    }
  }
  myColours <- grDevices::rainbow(x$nbClusters)
  if (!missing(colClusters) && !is.null(colClusters)) {
    if (length(colClusters) >= x$nbClusters) {
      myColours <- colClusters
    } else if (length(colClusters) >= length(iterationsIndices)) {
      myColours[iterationsIndices] <- colClusters
    } else {
      warning(
        "Parameter `colClusters` require at least ", length(iterationsIndices),
        " elements, reverting to default colors"
      )
    }
  }

  if (any(unlist(lapply(c("scater", "scatter", "plot"), grepl, contents)))) {
    do.call(plot, utils::modifyList(
      list(x = x$call$data, xlab = "X", ylab = "Y", col = col),
      otherArgs
    ))

    if (showClusters == TRUE) {
      for (i in 1:x$nbClusters) {
        selection <- x$call$data[x$labels == i, ]
        do.call(points, utils::modifyList(
          list(selection[, 1], selection[, 2], col = myColours[i]),
          otherArgs
        ))
      }
    }
    if (showMode == TRUE) {
      points(x$mode[1], x$mode[2], col = colMode, pch = 4, lwd = 2
      )
    }
  } else if (grepl("iter", contents)) {
    xlim <- c(min(x$call$data[, 1]), max(x$call$data[, 1]))
    ylim <- c(min(x$call$data[, 2]), max(x$call$data[, 2]))

    defaultOptions <- list(col = col)

    if (is.null(iterationsOptions) ||
        (is.list(iterationsOptions) && "contents" %in% names(iterationsOptions) &&
         (length(iterationsOptions$contents) > 1 || iterationsOptions$contents == "all"))) {
      defaultOptions <- append(defaultOptions, list(
        xlim = c(min(x$call$data[, 1]), max(x$call$data[, 1])),
        ylim = c(min(x$call$data[, 2]), max(x$call$data[, 2]))
      ))
    }

    if (!is.null(iterationsOptions)) {
      extraOptions <- utils::modifyList(iterationsOptions, otherArgs)
    } else {
      extraOptions <- otherArgs
    }

    if (length(iterationsIndices) == 1) {

      # Plot unimodal filtering
      do.call(plot, utils::modifyList(
        append(defaultOptions, list(
          x = x$iterations[[iterationsIndices]]$filteringUnimodal,
          colSelection = myColours[iterationsIndices]
        )),
        extraOptions
      ))

      # Plot normal filtering
      do.call(plot, utils::modifyList(
        append(defaultOptions, list(
          x = x$iterations[[iterationsIndices]]$filteringNormal,
          colSelection = myColours[iterationsIndices]
        )),
        extraOptions
      ))
    } else {

      # build options list
      options <- c()
      for (i in iterationsIndices) {
        options <- c(options, paste0("Iteration ", i, " - Unimodal Filtering"))
        options <- c(options, paste0("iteration ", i, " - Normal Filtering"))
      }

      # Plot the first figure and present menu menu
      choice <- 1
      repeat{
        if (choice %% 2 == 1) {
          do.call(plot, utils::modifyList(
            append(defaultOptions, list(
              x = x$iterations[[(choice + 1) %/% 2]]$filteringUnimodal,
              colSelection = myColours[((choice + 1) %/% 2)]
            )),
            extraOptions
          ))
        } else {
          do.call(plot, utils::modifyList(
            append(defaultOptions, list(
              x = x$iterations[[(choice + 1) %/% 2]]$filteringNormal,
              colSelection = myColours[((choice + 1) %/% 2)]
            )),
            extraOptions
          ))
        }

        # re-present menu waiting user choice
        cat("\nCURRENTLY DISPLAYED: ", options[choice], "\n\n")
        choice <- utils::menu(options, graphics = FALSE, title = "NEW PLOT SELECTION:")

        if(choice == 0){
          break
        }
      }
    }
  }
  else {
    stop("Value for parameter `contents` is not supported. Valid options are
         \"scatterplot\" or \"iterations\"")
  }
}
