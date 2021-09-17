
########################################################################### \
## Project: BRIL
## Script purpose: Wrapper for clustering methods
## Author: Adrien Brilhault
########################################################################### \

if (F) {
  install.packages("cluster")
  install.packages("mclust")
  install.packages("tclust")
  install.packages("dbscan")
}

suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("tclust"))
suppressPackageStartupMessages(library("dbscan"))


#' modeClustering Cluster the data, and returns the center of the main cluster
#'
#' @param data Matrix of numerical values containing the observations (one per
#'   row, with two columns for X and Y coordinates)
#' @param method Method to use. Valid options are "PAM", "TClust", "MClust",
#'   "HClust", "DBSCAN" body minimizers, or "L2", "Lui", "Mahalanobis", "Oja",
#'   "Projection" (default), "Spatial" and "Tukey" for depth functions
#' @param k.max Maximum number of clusters (default: 10)
#' @param alpha Fraction of outliers (only used with method "TClust")
#'
#' @return The function returns an S3 object of type `BRIL` containing the
#'   following values:
#'   \item{`nbClusters`}{Number of groups encountered}
#'   \item{`labels`}{Labels of the sample, corresponding to their cluster}
#'   \item{`clustersSizes`}{Array with the number of samples in each group}
#'   \item{`mainCluster`}{Index of the cluster identified as main mode}
#'   \item{`mode`}{Coordinates of the main mode}
#'
modeClustering <- function(data, method = "PAM", k.max = 10, alpha) {

  method <- match.arg(method, c("PAM", "TClust", "MClust", "HClust", "DBSCAN"))

  minSampleSize <- 4

  # Return the mean position in case other methods fail or if sample size is too small
  result <- list()
  result$mode <- colMeans(data)
  result$labels <- rep(1, times = nrow(data))
  result$nbClusters <- 1
  result$mainCluster <- 1

  if (is.null(nrow(data)) || nrow(data) < minSampleSize) {
    warning(sprintf(
      "modeClustering (%s): Sample Size is %d but should be greater than %d, mean returned",
      method, nrow(data), minSampleSize
    ))
    return(result)
  }

  silWidth <- -1
  bestk <- 1

  if (method == "PAM") {

    ## ========================================================================\
    ## -------- PAM  ---------
    ## ========================================================================\

    for (k in 2:k.max) {
      res <- pam(daisy(data, metric = "euclidean"), k, diss = TRUE)

      if (res$silinfo$avg.width > silWidth) { # silhouette test
        silWidth <- res$silinfo$avg.width
        bestk <- k
        # find label of the largest cluster
        mainCluster <- as.numeric(names(sort(table(res$clustering), decreasing = T, index.return = T))[1])
        # find index of the medoid of this cluster
        mode <- data[res$medoids[mainCluster], ]
        # store all labels
        labels <- res$clustering
      }
    }
    # cat("##### BEST K (Kmeans) = ", bestk,"\n")

    result$mode <- mode
    result$labels <- labels
    result$nbClusters <- bestk
    result$mainCluster <- mainCluster
  } else if (method == "HClust") {

    ## ========================================================================\
    ## -------- Hierarchical Clustering  ---------
    ## ========================================================================\


    dataDistances <- dist(data, method = "euclidean")
    res <- hclust(dataDistance, method = "centroid")

    for (k in 2:k.max) {
      clusters <- cutree(res, k = k)
      sil.avg.width <- mean(silhouette(clusters, dataDistances)[, 3])
      # cat("Silouette for ",k,"clusters :",sil.avg.width,"\n")
      if (sil.avg.width > silWidth) {
        silWidth <- sil.avg.width
        bestk <- k
        mainCluster <- names(sort(table(clusters), decreasing = T, index.return = T))[1]
        mode <- colMeans(data[which(clusters == mainCluster), ])
        labels <- clusters
      }
    }
    # cat("##### BEST K (HClust) = ",bestk,"\n")

    result$mode <- mode
    result$labels <- labels
    result$nbClusters <- bestk
    result$mainCluster <- mainCluster
  } else if (method == "MClust") {

    ## ========================================================================\
    ## -------- MClust  ---------
    ## ========================================================================\

    res <- Mclust(data, G = 1:k.max, verbose = FALSE)
    resSummary <- summary(res, parameters = TRUE)

    result$mode <- res$parameters$mean[, 1]
    result$labels <- res$classification
    result$nbClusters <- res$G
    result$mainCluster <- 1
  } else if (method == "DBSCAN") {

    ## ========================================================================\
    ## -------- DBScan  ---------
    ## ========================================================================\

    dataScaled <- scale(data)

    # eps.list <- seq(0.1, 2, 0.1)
    # eps.list <- seq(0.1, 5, 2)
    eps.list <- c(seq(0.1, 1, 0.1), seq(2, 10, 1))

    # minPts.list <- seq(5, 20, 1)
    # minPts.list<-c(5, 20, 50, 100)
    # minPts.list<-c(5,10)
    minPts.list <- c(5, 10, 15, 20)

    distances <- daisy(dataScaled, metric = "euclidean")

    for (eps in eps.list) {
      for (minPts in minPts.list) {
        res <- dbscan::dbscan(dataScaled, eps = eps, minPts = minPts)
        # res <- fcp::dbscan(data, eps = eps, MinPts = minPts)
        # res <- fpc::dbscan(distances, eps = eps, MinPts = minPts, method = "dist")

        if (length(unique(res$cluster)) >= 2) {
          silWidth.temp <- summary(silhouette(res$cluster, distances))$avg.width
        } else {
          silWidth.temp <- -1
        }

        if (silWidth.temp > silWidth) {
          bestEPS <- eps
          silWidth <- silWidth.temp
          bestk <- length(unique(res$cluster[res$cluster != 0]))
          clusters <- res$cluster
          # Select the cluster with the highest number of element
          mainCluster <- names(sort(table(clusters[clusters != 0]), decreasing = T, index.return = T))[1]
          mode <- colMeans(data[which(clusters == mainCluster), ])
          labels <- clusters
        }
        # print(res)
      }
    }
    # cat("##### BEST K (HClust) = ",cutClustering," BEST EPS = ", bestEPS, "\n")

    result$mode <- mode
    result$labels <- labels
    result$nbClusters <- bestk
    result$mainCluster <- mainCluster
  } else if (method == "TClust") {

    ## ========================================================================\
    ## -------- Tclust  ---------
    ## ========================================================================\

    distances <- daisy(data, metric = "euclidean")

    if (!missing(alpha)) {
      alpha.list <- alpha
    } else {
      alpha.list <- seq(0, 0.5, 0.1)
      # alpha.list <- c(0)
    }

    for (k in 2:k.max) {
      for (alpha in alpha.list) {
        res <- tclust(data, k = k, alpha = alpha)
        sil <- silhouette(res$cluster, as.matrix(distances))

        if (is.na(sil)) {
          silWidth.temp <- -1
        } else {
          silWidth.temp <- summary(sil)$avg.width
        }

        # cat("TClust : k =",k," - alpha = ",alpha," - sil= ",silWidth.temp,"\n")

        if (silWidth.temp >= silWidth) {
          silWidth <- silWidth.temp
          bestk <- length(unique(res$cluster[res$cluster != 0]))
          clusters <- res$cluster
          bestAlpha <- alpha
          mainCluster <- names(sort(table(clusters[clusters != 0]), decreasing = T, index.return = T))[1]
          mode <- colMeans(data[which(clusters == mainCluster), ])
          labels <- clusters
        }
      }
    }

    # print(bestAlpha)
    result$mode <- mode
    result$labels <- labels
    result$nbClusters <- bestk
    result$mainCluster <- mainCluster
  }
  else {
    stop(sprintf("modeClustering - unsuported value for parameter: method = \"%s\" ", method))
  }

  if (length(result$mode) != ncol(data) || any(is.nan(result$mode) || !is.numeric(result$mode))) {
    print(result$mode)
    stop(sprintf("Error in modeClustering (%s): mode is not numeric", method))
  }

  return(result)
}
