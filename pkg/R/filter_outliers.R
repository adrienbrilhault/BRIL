#' Recursive outlier filtering based on robust distances and multinormality
#' tests
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param center Coordinates used to computes the distances of the samples and
#'   order them (array of numerical two values, for X and Y)
#' @param test Statistical test to use. Valid options are "DIP" for unimodality
#'   test, or "Mardia", "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro",
#'   "Lillie", and "Chisq" for multivariate normality test
#' @param threshold Threshold of significance for the statistical test (between
#'   0 and 1, default: 0.05)
#' @param distType Distance metric used to order the samples. Valid options are
#'   "Euclidean", "MCD", "MVE", and "OGK". If empty or NULL, "Euclidean" will be
#'   automatically selected for unimodality tests, and "MCD" for normality
#'   tests
#' @param trimmedPerIteration Number of samples trimmed at each iteration
#'   (positive integer, default: 1)
#' @param debug Logical value. `TRUE` will compute all p.values, even after
#'   exceeding the threshold, for plotting purpose (see [plot.BRIL.Filtering()])
#' @param warnings Logical value, to display the warnings and errors caught
#'
#' @details For unimodality tests parameter `distType` should be set to
#'   "Euclidean" (as the distribution might contains a large amount of
#'   outliers). For normality tests robust distances are preferable, using a
#'   robust estimate estimates of location and scatter ("MCD","MVE", or "OGK").
#'
#' @return The function returns an S3 object of type `BRIL.Filtering`
#'   containing the following values:
#'   \item{`call`}{Parameters of the call (contains `data`, `test`,
#'   `testType`, `center`, `threshold`, `trimmedPerIteration` and `distType`)}
#'   \item{`distances`}{Distances of each sample from `data` to the `center`
#'   provided}
#'   \item{`p.values`}{P.Values of the test at each iteration}
#'   \item{`index.p.values`}{Subset size corresponding to each P.Value, for
#'   plotting purpose}
#'   \item{`selected`}{Indices of the samples from `data` selected at the end of
#'   the filtering}
#'   \item{`cutoffDistance`}{Distance of the furthest inlier selected}
#'
#' @seealso [plot.BRIL.Filtering()], [print.BRIL.Filtering()], [bril()],
#'   [median_rec()], [median_mv()], [depth_values()]
#'
#' @examples
#'
#' ## Example 1
#'
#' # Illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' # Compute an estimate for the center
#' center <- median_rec(XY)$median
#'
#' # Remove non unimodal outliers from this location
#' filtering <- filter_outliers(XY, center, test = "DIP", debug = TRUE)
#' print(filtering, maxDisplayed = 200)
#' plot(filtering)
#'
#' ## Example 2
#'
#' # Illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 4 - 1.5),
#'   mvtnorm::rmvnorm(150, c(5, 5), diag(2) * 400)
#' )
#'
#' # Compute an estimate for the center
#' center <- median_rec(XY)$median
#'
#' # Remove non normal outliers from this location
#' filtering <- filter_outliers(XY, center, test = "Chisq", distType = "MVE", debug = TRUE)
#' print(filtering)
#' plot(filtering, asp = 1)
#'
#' @export
#'
filter_outliers <- function(data, center, test = "Mardia", threshold = 0.05, distType, trimmedPerIteration = 1, debug = FALSE, warnings = FALSE) {

  ######## Check parameters & initialize

  if (missing(data)) {
    stop("Parameter `data` containing the observations is mandatory")
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  # convert from dataframe to matrix
  if (!is.numeric(data)) {
    stop("Parameter `data` must contain numerical values only")
  }
  if (is.null(ncol(data)) || ncol(data) != 2) {
    stop("Parameter`data` must contain exactly 2 columns, and more than one row")
  }
  if (missing(center)) {
    stop("Parameter `center` is mandatory")
  }
  if (!is.numeric(center) || length(center) != 2) {
    stop("Parameter `center` must be a numerical array of two values corresponding to X et Y coordinates")
  }
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Parameter `alpha` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.logical(debug)) {
    stop("Parameter `debug` must be a logical value")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  if (!is.numeric(trimmedPerIteration) || trimmedPerIteration < 1) {
    stop("Parameter `trimmedPerIteration` must be an entire number superior or equal to 1")
  }

  opts_test <- c("DIP", "Mardia", "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro", "Lillie", "Chisq")
  if (is.character(test) && (tolower(test) %in% tolower(opts_test))) {
    test <- opts_test[which(tolower(test) == tolower(opts_test))]
  } else {
    stop("Parameter `test` should be one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  if (test %in% c("DIP")) {
    testType <- "Unimodality"
  } else {
    testType <- "Normality"
  }

  opts_dist <- c("Euclidean", "MCD", "MVE", "OGK")
  if (missing(distType) || is.null(distType)) {
    if (testType == "Unimodality") {
      distType <- "Euclidean"
    } else {
      distType <- "MCD"
    }
  } else {
    if (is.character(distType) && tolower(distType) %in% tolower(opts_dist)) {
      distType <- opts_dist[which(tolower(distType) == tolower(opts_dist))]
    } else {
      stop("Parameter `distType` should be NULL or one of \"", paste(opts_dist, collapse = "\", \""), "\"")
    }

    if (warnings && testType == "Unimodality" && distType != "Euclidean") {
      warning(
        "For ", testType, " Tests, Euclidean Distances are preferable (\"",
        distType, "\" selected), consider leaving `distType` empty or setting it to \"Euclidean\""
      )
    }
    if (warnings && testType == "Normality" && distType == "Euclidean") {
      warning(
        "For ", testType, " Tests, Robust Distances are preferable (\"",
        distType, "\" selected), consider leaving `distType` empty or setting it to one of \"",
        paste(opts_dist[grep("Euclidean", opts_dist, invert = TRUE)], collapse = "\", "), "\""
      )
    }
  }



  trimmedPerIteration <- round(trimmedPerIteration)
  #   if (debug)
  #     trimmedPerIteration <- 1

  # Output of the function
  output <- list()
  class(output) <- "BRIL.Filtering"
  output$call$data <- data
  output$call$test <- test
  output$call$testType <- testType
  output$call$center <- center
  output$call$threshold <- threshold
  output$call$trimmedPerIteration <- trimmedPerIteration
  output$call$distType <- distType

  ######## Order samples

  # Compute robust distances
  if (distType != "Euclidean") {
    distances <- withCallingHandlers(
      withRestarts({
        if (toupper(distType) == "MCD") {
          robustScatter <- robustbase::covMcd(data, tolSolve = 1e-20)
        } else if (toupper(distType) == "OGK") {
          # other sigmamu : s_mad, s_IQR, s_Sn, s_Qn
          robustScatter <- robustbase::covOGK(data, sigmamu = robustbase::scaleTau2)
        } else if (toupper(distType) == "MVE") {
          robustScatter <- MASS::cov.mve(data)
        }

        # Check values are correct if it restarts from a warning
        if (!exists("robustScatter") || !("cov" %in% names(robustScatter)) ||
            is.null(robustScatter$cov) || all(robustScatter$cov == 0)) {
          distances <- NULL
          if (warnings) {
            warning("In filter_outliers() with test=\"", test, "\" and ", nrow(data),
                    " samples: robust distances could not be computed by \"", distType,
                    "\", reverting to euclidean distances",
                    call. = FALSE
            )
          }
        } else {
          distances <- sqrt(stats::mahalanobis(data, center, robustScatter$cov))
        }
        distances
      },
      muffleError = function() NULL
      ),
      warning = function(w) {
        if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
          warning("Warning caught in filter_outliers() with distType=\"", distType, "\" and ",
                  nrow(data), " samples:\n*  ", w,
                  call. = FALSE
          )
        }
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        if (warnings) {
          warning("Error caught in filter_outliers() with distType=\"", distType, "\" and ",
                  nrow(data), " samples:\n*  ", e,
                  call. = FALSE
          )
        }
        invokeRestart("muffleError")
      }
    )
  }

  # If Euclidean distances selected, or if robust estimates failed
  if (distType == "Euclidean" || is.null(distances)) {

    # Compute the euclidean distances of each sample to the center
    distances <- apply(data, 1, function(x, y) {
      stats::dist(rbind(x, y), method = "euclidean")
    }, center)
  }

  # order these distances
  selected <- sort(distances, index.return = TRUE)$ix
  output$distances <- distances


  ######## Filter outliers recursively

  # compute first test for all values
  if (testType == "Unimodality") {
    output$p.values <- test_unimodality(distances, test = test, warnings = warnings)
  } else {
    output$p.values <- test_multinormality(distances, test = test, data = data, center = center, warnings = warnings)
  }
  output$index.p.values <- length(selected)

  # Remove the furthest sample(s) until reaching multinormality
  while ((length(selected) - trimmedPerIteration > 2) && (is.na(utils::tail(output$p.values, 1)) || utils::tail(output$p.values, 1) < threshold)) {

    # remove furthest sample(s)
    selected <- utils::head(selected, - trimmedPerIteration)

    # apply test and store p.value with sample index
    if (testType == "Unimodality") {
      output$p.values <- c(output$p.values, test_unimodality(distances[selected], test = test, warnings = warnings))
    } else {
      output$p.values <- c(output$p.values, test_multinormality(distances[selected],
                                                                test = test, data = data[selected, ],
                                                                center = center, warnings = warnings
      ))
    }
    output$index.p.values <- c(output$index.p.values, utils::tail(output$index.p.values, 1) - trimmedPerIteration)
  }

  output$selected <- selected
  output$cutoffDistance <- distances[utils::tail(selected, 1)]

  # when debugging, compute the p.values for all the samples, even after having reached normality
  while (debug && ((length(selected) - trimmedPerIteration)) > 2) {
    selected <- utils::head(selected, - trimmedPerIteration)
    if (testType == "Unimodality") {
      output$p.values <- c(output$p.values, test_unimodality(distances[selected], test = test, warnings = warnings))
    } else {
      output$p.values <- c(output$p.values, test_multinormality(distances[selected],
                                                                test = test, data = data[selected, ],
                                                                center = center, warnings = warnings
      ))
    }
    output$index.p.values <- c(output$index.p.values, utils::tail(output$index.p.values, 1) - trimmedPerIteration)
  }

  return(output)
}




#' Test the unimodality of a distribution
#
#' @param values Unidimensional array of numerical values (distances)
#' @param test Statistical test used (for now the only option is "DIP")
#' @param warnings Logical value, to display the warnings and errors caught
#
#' @return p-value of the test (lower values suggest multimodality)
#'
#' @seealso [filter_outliers()]
#'
test_unimodality <- function(values, test = "DIP", warnings = FALSE) {
  opts_test <- c("DIP")
  if (is.character(test) && tolower(test) %in% tolower(opts_test)) {
    test <- opts_test[which(tolower(test) == tolower(opts_test))]
  } else {
    stop("Parameter `test` should be one of \"", paste(opts_test, collapse = "\", "), "\"")
  }

  # Compute the test statistics
  p.value <- withCallingHandlers(
    withRestarts({
      if (test == "DIP") {
        p.value <- diptest::dip.test(values)$p.value
      }

      # Check values are correct if it restarts from a warning
      if (!exists("p.value") || is.na(p.value) || is.null(p.value)) {
        if (warnings) {
          warning("In unimodality_test() with test=\"", test, "\" and ", length(values),
                  " samples: p.value is NA\n",
                  call. = FALSE
          )
        }
        p.value <- NA
      }
      p.value
    },
    muffleError = function() NA
    ),
    warning = function(w) {
      if (warnings) {
        warning("Warning caught in unimodality_test() with test=\"", test, "\" and ",
                length(values), " samples:\n*  ", w,
                call. = FALSE
        )
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      if (warnings) {
        if (warnings) {
          warning("Error caught in unimodality_test() with test=\"", test, "\" and ",
                  length(values), " samples:\n*  ", e,
                  call. = FALSE
          )
        }
      }
      invokeRestart("muffleError")
    }
  )

  return(p.value)
}



#' Test the multivariate normality of a distribution
#
#' @param values Unidimensional array of numerical values (distances)
#' @param test Statistical test used, valid options are "Mardia", "Kurtosis",
#'   "Skewness", "KS", "KS-adj", "Shapiro", "Lillie", and "Chisq"
#' @param threshold Threshold of significance for the statistical test
#'   (default: 0.05)
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param center Center of the observations from `data`
#' @param warnings Logical value, to display the warnings and errors caught
#
#' @details Parameter `data` is only required for the tests "Mardia",
#' "Skewness", "Kurtosis" and "Chisq", while parameter `center` is only
#' required for "Chisq"
#
#' @return p-value of the test (lower values suggest non normality)
#'
#' @seealso [filter_outliers()]
#'
test_multinormality <- function(values, test = "Mardia", threshold = 0.05, data, center, warnings = FALSE) {

  ######## Check parameters

  opts_test <- c("Mardia", "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro", "Lillie", "Chisq")
  if (is.character(test) && (tolower(test) %in% tolower(opts_test))) {
    test <- opts_test[which(tolower(test) == tolower(opts_test))]
  } else {
    stop("Parameter `test` should be one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  if ((test %in% c("Mardia", "Skewness", "Kurtosis", "Chisq")) && (missing(data) || !is.numeric(data) || dim(data)[2] != 2)) {
    stop("For tests \"Mardia\", \"Skewness\", \"Kurtosis\" and \"Chisq\", parameter `data` with the observations (as a matrix of numerical values) is required.")
  }

  if ((test %in% c("Chisq")) && (missing(center) || !is.numeric(center) || length(center) != 2)) {
    stop("For test \"Chisq\", parameter `center` (a pair of numerical values) is required.")
  }

  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Parameter `threshold` must be a numerical value between 0 and 1")
  }

  ######## Compute p.values for the test

  p.value <- withCallingHandlers(
    withRestarts({
      if (test == "Shapiro") {
        p.value <- stats::shapiro.test(values)$p.value
      } else if (test == "Lillie") {
        p.value <- nortest::lillie.test(values)$p.value # adapted KS for unknown variance and expected value
      } else if (test == "KS-adj") {
        # adjust the confidence value in function of size n : confKS=1.36 / sqrt(n)
        p.value <- threshold / (1.36 / sqrt(length(values))) * stats::ks.test(values, "pnorm", mean(values), stats::sd(values))$p.value
      } else if (test == "KS") {
        p.value <- stats::ks.test(values, "pnorm", mean(values), stats::sd(values))$p.value
      } else if (test == "Kurtosis") {
        p.value <- ICS::mvnorm.kur.test(data, method = "integration")$p.value
      } else if (test == "Skewness") {
        p.value <- ICS::mvnorm.skew.test(data)$p.value
      } else if (test == "Chisq") {
        p.value <- stats::ks.test(stats::mahalanobis(data, center, stats::cov(data), tol = 1e-20), "pchisq", df = 2)$p.value
      } else if (test == "Mardia") {
        p.value <- ICS::mvnorm.kur.test(data, method = "integration")$p.value
        if (p.value > threshold) {
          # Mardia test rejects the null hypothesis only if both kurtosis and skewness reject them
          p.value <- ICS::mvnorm.skew.test(data)$p.value
        }
      }
      # Check values are correct if it restarts from a warning
      if (!exists("p.value") || is.na(p.value) || is.null(p.value)) {
        if (warnings) {
          warning("In multinormality_test() with test=\"", test, "\" and ", length(values),
                  " samples: p.value is NA\n",
                  call. = FALSE
          )
        }
        p.value <- NA
      }
      p.value
    },
    muffleError = function() NA
    ),
    warning = function(w) {
      if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
        warning("Warning caught in multinormality_test() with test=\"", test, "\" and ",
                length(values), " samples:\n*  ", w,
                call. = FALSE
        )
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      if (warnings) {
        warning("Error caught in multinormality_test() with test=\"", test, "\" and ",
                length(values), " samples:\n*  ", e,
                call. = FALSE
        )
      }
      invokeRestart("muffleError")
    }
  )

  return(p.value)
}




#' Print method for `BRIL.Filtering` objects
#'
#' @param x An object of class `BRIL.Filtering` (see [filter_outliers()])
#' @param maxDisplayed Number of elements to display in the output (default:
#'   200). Set to NULL (or 0) to show all values.
#' @param ... Other arguments passed to or from other methods
#'
#' @seealso [filter_outliers()], [plot.BRIL.Filtering()], [bril()]
#'
#' @export
#'
print.BRIL.Filtering <- function(x, maxDisplayed = 200, ...) {
  if (class(x) != "BRIL.Filtering") {
    stop("the object provided is not of class \"BRIL.Filtering\"")
  }

  if (is.null(maxDisplayed) || maxDisplayed == 0) {
    maxDisplayed <- length(x$distances)
  } else if (!is.numeric(maxDisplayed) || maxDisplayed < 0) {
    stop("Parameter `maxDisplayed` must be a numerical value superior or equal to 1")
  }
  maxDisplayed <- round(maxDisplayed)

  cat(sprintf(
    "\n=> Results for filter_outlier() using %s %s Test (p.value > %.2f) and %s Distances\n",
    x$call$test, x$call$testType, x$call$threshold,
    ifelse(x$call$distType == "Euclidean", "Euclidean", paste(x$call$distType, "-based Robust", sep = ""))
  ))
  cat(sprintf(
    "   %d samples: %d selected, %d filtered (%d trimmed per iteration)",
    nrow(x$call$data), length(x$selected), nrow(x$call$data) - length(x$selected), x$call$trimmedPerIteration
  ))


  cat("\n\nSelected indices:\n")
  print(utils::head(x$selected, maxDisplayed), ...)
  if (length(x$selected) > maxDisplayed) {
    cat(paste0(" ... (", length(x$selected) - maxDisplayed, " hidden)"))
  }


  cat("\n\n", toupper(x$call$test), " Test p.values:\n", sep = "")
  print(utils::head(x$p.values, maxDisplayed), ...)
  if (length(x$p.values) > maxDisplayed) {
    cat(paste0(" ... (", length(x$p.values) - maxDisplayed, " hidden)"))
  }

  cat("\n\nOutliers cutoff distance from center (", toString(x$call$center), "):\n", sep = "")
  print(x$cutoffDistance, ...)

  cat("\n\n", ifelse(x$call$distType == "Euclidean", "Euclidean ", paste(x$call$distType, "-based Robust ", sep = "")),
      "Distances:\n",
      sep = ""
  )
  print(utils::head(x$distances, maxDisplayed), ...)
  if (length(x$distances) > maxDisplayed) {
    cat(paste0(" ... (", length(x$distances) - maxDisplayed, " hidden)"))
  }

  cat("\n\n")

  invisible(x)
}


#' Plot method for `BRIL.Filtering` objects
#'
#' @param x An object of class `BRIL.Filtering` (see [filter_outliers()])
#' @param contents Contents to be displayed, options are "p.values",
#'   "scatterplot", "dist", "hist" and "all"
#' @param showCenter Logical value, to show the center used in the filtering
#' @param showSelection Logical value, to highlight the samples selected by the
#'   filtering process
#' @param col Default color for non-selected samples (default: "black")
#' @param colSelection Color of the selected samples (default: "red")
#' @param colCenter Color for the center in "scatterplot" (default: "orange")
#' @param mtextTitles Logical value, `TRUE` to set smaller titles/subtitles on
#' top, `FALSE` to use the default plot title options.
#' @param mfrow Number of rows and columns of the figure (example:
#'   c(4,1))
#' @param ... Other arguments passed to or from other methods (such as pch for
#'   the symbols, main and sub for title and subtitle, xlab, xmin, ...)
#'
#' @details
#'   \emph{Red intercept lines correspond to the
#'   selection based on the p.values exceeding the given threshold.} \cr
#'   \emph{To display all the p.values, rerun the function [filter_outliers()]
#'   with the parameter `debug = TRUE`} \cr \cr
#'
#' `contents` options:
#'   - \bold{"p.values"} provides a plot of the test p.values (in function of
#'   the subset size)
#'   - \bold{"scatterplot"} displays the data in cartesian coordinates. Selected
#'   samples are displayed in red, and the center used to compute distances as
#'   an orange cross
#'   - \bold{"dist"} shows the distances of each sample to the center provided
#'   in [filter_outliers()] (in function of sample index)
#'   - \bold{"hist"} draws an histogram of the distances of the samples to
#'   the center provided in [filter_outliers()]
#'   - \bold{"all"} displays a figure with all of the options above
#'
#
#' @seealso [filter_outliers()], [print.BRIL.Filtering()], [bril()],
#'   [median_rec()], [median_mv()], [depth_values()]
#'
#' @examples
#'
#' # Illustrative data
#' XY <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' # Process the data
#' filtering <- filter_outliers(XY, median_rec(XY)$median, test = "DIP", debug = TRUE)
#'
#' # Plot all default figures
#' plot(filtering)
#'
#' # Plot P.Values and Scatterplot only
#' plot(filtering, contents = c("pvalues", "scatterplot"))
#'
#' # Change the layout to vertical
#' plot(filtering, contents = c("pvalues", "scatterplot"), mfrow = c(2, 1))
#'
#' # Remove title, subtitle, and axis labels
#' plot(filtering, contents = "scatterplot", main = "", sub = "",
#'   ylab = "", xlab = "")
#'
#' # Other graphical options
#' plot(filtering, contents = "scatterplot", asp = 1,
#'   xlim = c(-30, 30), ylim = c(-30, 30))
#'
#' plot(filtering,
#'   contents = "scatterplot", asp = 1, pch = 4, lwd = 2, col = "blue",
#'   colSelection = "green", showCenter = FALSE
#' )
#'
#' plot(filtering, contents = "hist", main = "My Histogram",
#'   showSelection = FALSE, breaks = 50)
#'
#' @export
#'
plot.BRIL.Filtering <- function(x, contents = c("p.values", "scatterplot", "dist", "hist"), showCenter = TRUE,
                                showSelection = TRUE, col = "black", colSelection = "red", colCenter = "orange",
                                mtextTitles = TRUE, mfrow, ...) {
  if (class(x) != "BRIL.Filtering") {
    stop("the object provided is not of class \"BRIL.Filtering\"")
  }

  if (!is.logical(mtextTitles)) {
    stop("Parameter `mtextTitles` must be a logical value")
  }

  contents <- unique(tolower(contents))
  if ("all" %in% contents) {
    contents <- c("p.values", "scatterplot", "dist", "hist")
  }
  nbPlots <- length(contents)


  # Number of row and lines of the figure
  if (!missing(mfrow) && is.numeric(mfrow) && length(mfrow) == 2) {
    graphics::par(mfrow = mfrow)
  } else {
    if (nbPlots != 1) {
      if (nbPlots <= 3) {
        graphics::par(mfrow = c(1, nbPlots))
      } else {
        graphics::par(mfrow = c(2, ceiling(nbPlots / 2)))
      }
    }
  }

  # Read extra parameters
  otherArgs <- list(...)

  # Plot
  while (length(contents) >= 1) {
    if (any(unlist(lapply(c("scater", "scatter", "plot"), grepl, contents[1])))) {
      title <- ifelse("main" %in% names(otherArgs), list(otherArgs$main), "Samples selected")
      subtitle <- ifelse("sub" %in% names(otherArgs), list(otherArgs$sub),
                         paste0("(", nrow(x$call$data), " samples, ", length(x$selected), " selected)")
      )

      do.call(plot, utils::modifyList(
        list(
          x = x$call$data, xlab = "X", ylab = "Y", col = col,
          main = ifelse(mtextTitles, list(NULL), title), sub = ifelse(mtextTitles, list(NULL), subtitle)
        ),
        otherArgs[names(otherArgs) %in% c("main", "sub") == FALSE]
      ))

      if (showSelection == TRUE) {
        selection <- x$call$data[x$selected, ]
        do.call(points, utils::modifyList(
          list(selection[, 1], selection[, 2], col = colSelection),
          otherArgs[names(otherArgs) %in% c("col") == FALSE]
        ))
      }
      if (showCenter == TRUE) {
        points(x$call$center[1], x$call$center[2],
                         col = colCenter, pch = 3, lwd = 2
        )
      }

      if (mtextTitles) {
        graphics::mtext(side = 3, line = 2, adj = 0, cex = 1, title)
        graphics::mtext(side = 3, line = 1, adj = 0, cex = 0.7, subtitle)
      }
    }
    else if (any(unlist(lapply(c("p.val", "pval", "value"), grepl, contents[1])))) {
      title <- ifelse("main" %in% names(otherArgs), list(otherArgs$main), paste(x$call$test, "Test"))
      subtitle <- ifelse("sub" %in% names(otherArgs), list(otherArgs$sub),
                         paste0("P.Values for the ", x$call$test, " Test of ", x$call$testType, " (p < ", x$call$threshold, ")")
      )

      do.call(plot, utils::modifyList(
        list(
          x = x$index.p.values, x$p.values, col = col,
          xlim = c(0, max(x$index.p.values, na.rm = TRUE)),
          ylim = c(0, max(0.3, max(x$p.values, na.rm = TRUE))),
          xlab = "Subset Size", ylab = "P.Value",
          main = ifelse(mtextTitles, list(NULL), title), sub = ifelse(mtextTitles, list(NULL), subtitle)
        ),
        otherArgs[names(otherArgs) %in% c("main", "sub", "asp", if (nbPlots > 1) c("xlim", "ylim")) == FALSE]
      ))

      if (showSelection == TRUE) {
        selection <- which(x$index.p.values <= length(x$selected))
        do.call(points, utils::modifyList(
          list(x$index.p.values[selection], x$p.values[selection], col = colSelection),
          otherArgs[names(otherArgs) %in% c("col") == FALSE]
        ))
        graphics::abline(h = x$call$threshold, col = colSelection)
      }
      if (mtextTitles) {
        graphics::mtext(side = 3, line = 2, adj = 0, cex = 1, title)
        graphics::mtext(side = 3, line = 1, adj = 0, cex = 0.7, subtitle)
      }
    }
    else if (any(unlist(lapply(c("dist"), grepl, contents[1])))) {
      title <- ifelse("main" %in% names(otherArgs), list(otherArgs$main),
                      ifelse(x$call$distType == "Euclidean", "Distances", "Robust Distances")
      )
      subtitle <- ifelse("sub" %in% names(otherArgs), list(otherArgs$sub),
                         paste0(
                           ifelse(x$call$distType == "Euclidean", "Euclidean Distances",
                                  paste0(x$call$distType, "-based Robust Distances")
                           ),
                           " to (", toString(lapply(x$call$center, function(x) sprintf("%.2f", x))), ")"
                         )
      )

      do.call(plot, utils::modifyList(
        list(
          x = x$distances, xlab = "Index", ylab = "Distance", col = col,
          main = ifelse(mtextTitles, list(NULL), title), sub = ifelse(mtextTitles, list(NULL), subtitle)
        ),
        otherArgs[names(otherArgs) %in% c("main", "sub", "asp", if (nbPlots > 1) c("xlim", "ylim")) == FALSE]
      ))

      if (showSelection == TRUE) {
        selection <- which(x$distances <= x$cutoffDistance)
        do.call(graphics::points, utils::modifyList(
          list(selection, x$distances[selection], col = colSelection),
          otherArgs[names(otherArgs) %in% c("col") == FALSE]
        ))
        graphics::abline(h = x$cutoffDistance, col = colSelection)
      }
      if (mtextTitles) {
        graphics::mtext(side = 3, line = 2, adj = 0, cex = 1, title)
        graphics::mtext(side = 3, line = 1, adj = 0, cex = 0.7, subtitle)
      }
    }

    else if (any(unlist(lapply(c("hist"), grepl, contents[1])))) {
      title <- ifelse("main" %in% names(otherArgs), list(otherArgs$main),
                      paste0("Histogram of ", ifelse(x$call$distType == "Euclidean", "Distances", "Robust Distances"))
      )
      subtitle <- ifelse("sub" %in% names(otherArgs), list(otherArgs$sub),
                         paste0(
                           "Histogram of ", ifelse(x$call$distType == "Euclidean", "Distances",
                                                   paste0(x$call$distType, "-based Robust Distances")
                           ),
                           " to (", toString(lapply(x$call$center, function(x) sprintf("%.2f", x))), ")"
                         )
      )

      do.call(graphics::hist, utils::modifyList(
        list(
          x = x$distances, breaks = 30, col = "grey90",
          xlab = "Distances", ylab = "Count",
          main = ifelse(mtextTitles, list(NULL), title), sub = ifelse(mtextTitles, list(NULL), subtitle)
        ),
        otherArgs[names(otherArgs) %in% c("main", "sub", "asp", if (nbPlots > 1) c("xlim", "ylim", "col")) == FALSE]
      ))
      if (showSelection == TRUE) {
        graphics::abline(v = x$cutoffDistance, col = colSelection)
      }
      if (mtextTitles) {
        graphics::mtext(side = 3, line = 2, adj = 0, cex = 1, title)
        graphics::mtext(side = 3, line = 1, adj = 0, cex = 0.7, subtitle)
      }
    }
    else {
      warning(sprintf(
        "Option \"%s\" unrecognized.\nValid options are \"p.values\", \"scatterplot\", \"dist\", \"hist\", and \"all\".",
        contents[1]
      ))
    }
    # parse next option
    contents <- utils::tail(contents, -1)
  }

  # Restore standard grids
  if (nbPlots != 1) {
    graphics::par(mfrow = c(1, 1))
  }
}
