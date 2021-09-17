
########################################################################### \
## Project: BRIL
## Script purpose: Create, test and save artificial benchmarks for BRIL
## Author: Adrien Brilhault
########################################################################### \


## ========================================================================\
## -------- Initialization  --------
## ========================================================================\


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (F) {
  install.packages("doParallel")
  install.packages("doSNOW")
  install.packages("ggplot2")
  install.packages("mvtnorm")
  install.packages("CompQuadForm")
  install.packages("remotes")
  install.packages("tcltk")
  remotes::install_github("adrienbrilhault/BRIL", subdir = "pkg")
}


library(BRIL)
library(mvtnorm)
library(ggplot2)
library(doParallel)
library(doSNOW)
library(tcltk) # For progress bar widget


## ========================================================================\
## -------- Parameters  --------
## ========================================================================\

params = list()

######### Path to save the benchmark resultss

params$folderData <- "results/Test"

######### Misc

params$writeRes <- T # Output resultss to file
params$plotRes <- T # Display plots
params$savePlots <- F # Save plots de jpg files in result folder
params$consoleMode <- F # disables plots, progress widget, and all graphic outputs
params$progressBar <- T # display progress bar
params$verbose <- T # Verbose
params$debug <- T # Run some additional debugging code

######### Number of repetitions

params$oneRun <- T # Only one simulation (for debug)
params$montecarloRepetitions <- 200 # Number of simulations

######### Parallel simulations

params$multithread <- T # Set to true to parallelize the simulations
params$nbCores <- parallel::detectCores()
# params$nbCores <- 2


######### Number of clusters and samples

params$nbClusterValues <- 2:5
# params$nbClusterValues <- 3:5
# params$nbClusterValues <- 3

# params$sampleSizeValues <- c(250,500,1000,2000)
# params$sampleSizeValues <- seq(250,2000,250)
params$sampleSizeValues <- c(500)

######### Inliers/Outliers ratios

params$randomClustersSize <- F # for the outliers clusters
params$inliersRatioValues <- list(
  "1" = c(100),
  "2" = c(50:55, 60, 70, 80, 100),
  "3" = c(33.3, 34, 35, 36, 38, 40, 45, 50, 55, 60, 70, 80, 100),
  "4" = c(25, 26, 28, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100),
  "5" = c(20, 22, 24, 28, 30, 35, 40, 50, 55, 60, 70, 80, 100)
)
# params$inliersRatioValues <- c(25, 35)
# params$inliersRatioValues <- 40

######### Multivariate Gaussians parameters

params$randomVariance <- T
# params$randomVarianceInterval <- c(20,300)
# params$randomVarianceInterval <- c(20,150)
params$randomVarianceInterval <- c(50, 150) # in percents of the main cluster variance
params$randomCovariance <- T
params$plotArea <- matrix(c(0, 50, 0, 50), 2, 2)
# params$plotArea <- matrix(c(0, 30, 0, 50), 2, 2)

######### Uniform Noise

params$uniformNoisePercent <- T # If false values are expressed as absolute number of samples, if true as a percentage of sample size
params$uniformNoiseValues <- c(0, 25)
# params$uniformNoiseValues <- seq(0, 50, 5)
# params$uniformNoiseValues <- 25

######### Clustering algorithms

params$doClustering <- T # Run clustering algorithms
params$clusteringKMax <- 10
# params$plotClustering <- params$plotRes
params$plotClustering <- T
# params$clusteringMethods <- c("PAM", "MClust", "DBSCAN", "TClust")
params$clusteringMethods <- c("PAM", "MClust", "DBSCAN")

######### BRIL tuning parameters

params$bril.unimodalThresh <- 0.05
params$bril.testNormal <- "Mardia"
params$bril.normalThresh <- 0.10
params$bril.distNormal <- "MCD"
params$bril.alpha <- 0.5
params$bril.trimmedPerFilteringIteration <- 1
params$bril.warnings <- TRUE
# params$bril.methods =  c("MVE","MCD","Projection","L2")
# params$bril.methods =  c("Projection")
params$bril.methods <- c("MVE", "MCD", "Tukey", "Liu", "Oja", "Projection", "Spatial", "L2", "Mahalanobis")

######### END PARAMETERS

if (params$consoleMode) {
  params$plotRes <- FALSE
  params$plotClustering <- FALSE
}
if (params$oneRun) {
  params$writeRes <- FALSE
  params$multithread <- FALSE
  params$progressBar <- FALSE
}
if (params$writeRes && !file.exists(params$folderData)) {
  dir.create(file.path(params$folderData))
}
if (params$writeRes) {
  save(params, file = file.path(params$folderData,paste0("params_",format(Sys.time(), "%Y%m%d_%H%M%S"),".Rdata")))
}

## ========================================================================\
## -------- Simulations  --------
## ========================================================================\


if (is.list(params$inliersRatioValues)) {
  totalIterations = sum(sapply(params$nbClusterValues, function(x) length(params$inliersRatioValues[[as.character(x)]]))) *
    length(params$sampleSizeValues) * length(params$uniformNoiseValues) * params$montecarloRepetitions
} else {
  totalIterations = length(params$nbClusterValues) * length(params$sampleSizeValues) * length(params$uniformNoiseValues) * params$montecarloRepetitions
}
iterationsCompleted = 0

if (params$progressBar && !params$consoleMode ) {
  globalProgressBar <- tkProgressBar(max = totalIterations, width = 550, initial = 1,
                                     title = sprintf("%d/%d - Unit %d/%d (Repetition %d/%d)",
                                                     1+iterationsCompleted, totalIterations,
                                                     1+iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                                                     1+iterationsCompleted, totalIterations))
}

start.time <- Sys.time()

for (nbClusters in params$nbClusterValues) {

  if (is.list(params$inliersRatioValues)) {
    if (toString(nbClusters) %in% names(params$inliersRatioValues)) {
      inliersRatios <- params$inliersRatioValues[[toString(nbClusters)]]
    } else {
      stop("Inliers ratio values not found for ",nbClusters," clusters")
    }
  } else {
    inliersRatios <- params$inliersRatioValues
  }

  for (inliersRatio in inliersRatios) {
    for (sampleSize in params$sampleSizeValues) {
      for (noiseRatio in params$uniformNoiseValues) {

        cat(sprintf(
          "\n#### Start of unit %d/%d: NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d (%d repetitions)\n##################################################################################################\n\n",
          1+iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions, nbClusters, inliersRatio, sampleSize, noiseRatio, params$montecarloRepetitions
        ))

        ## Structure to store the results
        results.all <- data.frame()

        ## Progress bar
        opts <- list()
        if (params$progressBar) {
          if (params$consoleMode) {

            progressBar <- utils::txtProgressBar(min = 1, max = params$montecarloRepetitions, style = 3)
            progress <- function(n) utils::setTxtProgressBar(progressBar, n)

          } else {

            setTkProgressBar(globalProgressBar,iterationsCompleted+1,
                             title = sprintf("%d/%d - Unit %d/%d (Repetition %d/%d)",
                                             1+iterationsCompleted, totalIterations,
                                             1+iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                                             1,params$montecarloRepetitions),
                             label = sprintf("Current Unit: NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d",
                                             nbClusters, inliersRatio, sampleSize, noiseRatio))

            progress <- function(n) tcltk::setTkProgressBar(globalProgressBar,iterationsCompleted+n,
                                                            title = sprintf("%d/%d - Unit %d/%d (Repetition %d/%d)",
                                                                            1+n+iterationsCompleted, totalIterations,
                                                                            1+iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                                                                            1+n,params$montecarloRepetitions))
          }
          # Function called to update progress in the foreach loop
          opts <- list(progress = progress)
        }

        ## -------- Multithreading  --------

        if (params$multithread) {

          # cl <- makeCluster(params$nbCores, type = "SOCK", outfile = "")
          cl <- makeCluster(params$nbCores)

          # registerDoParallel(cl)
          registerDoSNOW(cl)

          `%mydopar%` <- `%dopar%`
        } else {
          `%mydopar%` <- `%do%`
        }

        listPackages <- c("ggplot2")

        FOREACHLOOP <- T # set this to false if you switch below the foreach loop for a classic one

        # for (montecarlo.i in 1:params$montecarloRepetitions) {
        results.all <- foreach(montecarlo.i = (1:ifelse(params$oneRun, 1, params$montecarloRepetitions)), .combine = rbind, .errorhandling = "remove", .options.snow = opts, .packages = listPackages, .inorder = F) %mydopar% {

          ## ========================================================================\
          ## --- GENERATE DATA  ----
          ## ========================================================================\

          ## ========= Population components sizes  ============

          ## Fist cluster are the inliers

          nbSamplesPerCluster <- c(inliersRatio)
          if (params$randomClustersSize == FALSE) {
            nbSamplesPerCluster <- c(nbSamplesPerCluster, rep((100 - inliersRatio) / (nbClusters - 1), times = nbClusters - 1))
          } else {
            if (nbClusters > 2) {
              nbSamplesPerCluster  <- c(nbSamplesPerCluster, 100 - inliersRatio)
            } else {
              for (i in 1:(nbClusters - 2)) {
                minSize <- max(0, (100 - sum(nbSamplesPerCluster)) - (inliersRatio - 1) * (nbClusters - 1 - i))
                size <- runif(1, minSize, inliersRatio - 1)
                nbSamplesPerCluster <- c(nbSamplesPerCluster, size)
              }
              nbSamplesPerCluster <- c(nbSamplesPerCluster, max(0, 100 - sum(nbSamplesPerCluster)))
            }
          }


          if (params$uniformNoisePercent == FALSE) {

            noise.nb <- noiseRatio
            samples.nb <- sampleSize + noise.nb

          } else {

            noise.nb <- round(sampleSize * noiseRatio / 100)
            samples.nb <- sampleSize
          }


          nbSamplesPerCluster <- round((samples.nb - noise.nb) * nbSamplesPerCluster / sum(nbSamplesPerCluster))

          if (params$verbose) {
            print(sprintf("Number of samples=%d, noise=%d, total=%d", samples.nb - noise.nb, noise.nb, samples.nb))
            print(sprintf("Number of samples per cluster: %s", toString(nbSamplesPerCluster)))
          }

          ## ========= Population positions  ============

          if (params$verbose) {
            print("Drawing positions")
          }

          ## Pick random positions for clusters

          clusters.position <- matrix(c(
            sample(params$plotArea[1, 1]:params$plotArea[2, 1], nbClusters, replace = T),
            sample(params$plotArea[1, 2]:params$plotArea[2, 2], nbClusters, replace = T)
          ), ncol = 2)

          ## Pick new random position if clusters are too close to each other or colinears
          if (nbClusters > 1) {
            while (min(dist(clusters.position)) < 10) {
              # while (min(dist(clusters.position)) < 10 || min(testAllAngles(clusters.position)) < 30) {
              clusters.position <- matrix(c(
                sample(params$plotArea[1, 1]:params$plotArea[2, 1], nbClusters, replace = T),
                sample(params$plotArea[1, 2]:params$plotArea[2, 2], nbClusters, replace = T)
              ), ncol = 2)
            }
          }

          sigma <- matrix(c(1, 0, 0, 1), nc = 2)


          ## ========= Gaussian variances  ============

          clusters.variance <- rep(1, nbClusters)
          if (params$randomVariance && nbClusters > 1) {
            clusters.variance[2:nbClusters] <- runif(nbClusters - 1, params$randomVarianceInterval[1] / 100, params$randomVarianceInterval[2] / 100)
          }


          ## ========= Draw Population  ===========================================
          ## ========================================================================\

          if (params$verbose) {
            print("Drawing samples")
          }

          ## ========= Mixture of bivariate distributions

          # mixbivnorm <- mvtnorm::rmvnorm(nbSamplesPerCluster[1], clusters.position[1, ], sigmaL[[1]])
          mixbivnorm <- mvtnorm::rmvnorm(nbSamplesPerCluster[1], clusters.position[1, ], sigma)

          if (nbClusters > 1) {
            for (i in 2:nbClusters) {
              if (nbSamplesPerCluster[i] != 0) {

                if (params$randomCovariance == T) {
                  sigmaT <- 1.5 * sigma * clusters.variance[i]
                  sigmaT[1, 1] <- sigmaT[1, 1] / runif(1, min = 1, max = 8)
                  theta <- runif(1, min = 0, max = 2 * pi)
                  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
                  sigmaT <- rot %*% sigmaT %*% t(rot)
                } else {
                  sigmaT <- sigma * clusters.variance[i]
                  # sigmaT <- sigmaL[[i]] * clusters.variance[i]
                }

                mixbivnorm <- rbind(mixbivnorm, mvtnorm::rmvnorm(nbSamplesPerCluster[i], clusters.position[i, ], sigmaT))
              }
            }
          }

          samples.XY <- mixbivnorm


          ## ========= Uniform noise populations

          if (noise.nb > 0) {

            ### Set the intervals of the area with uniform noise
            if (T) {

              areaNoise <- matrix(c(
                min(mixbivnorm[, 1] - 10), max(mixbivnorm[, 1] + 10),
                min(mixbivnorm[, 2] - 10), max(mixbivnorm[, 2]) + 10
              ), 2)
            } else {
              # areaNoise <- rbind(params$plotArea[1, ] - 30, params$plotArea[2, ] + 30)
              areaNoise <- params$plotArea
            }

            ### Generate the noise

            # mixuniform <- cbind(runif(noise.nb, areaNoise[1, 1], areaNoise[2, 1]),runif(noise.nb, areaNoise[1, 2], areaNoise[2, 2]))

            i <- 0
            mixuniform <- matrix(0, nrow = 0, ncol = 2)
            while (i < noise.nb) {
              ### Create a random point in the interval
              noisePoint <- cbind(
                runif(1, areaNoise[1, 1], areaNoise[2, 1]),
                runif(1, areaNoise[1, 2], areaNoise[2, 2])
              )

              ### Compute distances of this point to the clusters
              distancesToPoint <- apply(clusters.position, 1, function(x, y) {
                dist(rbind(x, y))
              }, noisePoint)
              # split(x, rep(1:nrow(clusters.position), each = ncol(clusters.position)))

              ### Keep if it's above 3 standard deviations away from each cluster
              outsideClusters <- mapply(function(x, y) {
                x < 3 * sqrt(y)
              }, distancesToPoint, clusters.variance)
              if (!any(outsideClusters)) {
                mixuniform <- rbind(mixuniform, noisePoint)
                i <- i + 1
              }
            }

            samples.XY <- rbind(samples.XY, mixuniform)
          }


          ## ========================================================================\
          ## -------- COMPUTE MEDIANS  --------
          ## ========================================================================\

          if (params$verbose) {
            print("Computing estimates")
          }

          ## ========= Simple medians  ============

          results <- data.frame(
            Method = character(), X = double(), Y = double(), nbClust = integer(),
            Group1 = character(), Group2 = character(), Group3 = character(), Time = integer()
          )

          ## Mean of all inliers
          sTime <- as.numeric(Sys.time()) * 1000
          MeanInliers <- colMeans(samples.XY[1:nbSamplesPerCluster[1], ])
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = "MeanInliers", Group1 = "Reference", Group2 = "Reference", Group3 = "Reference",
            X = MeanInliers[1], Y = MeanInliers[2], nbClust = 0, Time = execTime
          ))

          ## Mean of all samples
          sTime <- as.numeric(Sys.time()) * 1000
          Mean <- c(mean(samples.XY[, 1]), mean(samples.XY[, 2]))
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = "Mean", Group1 = "Classic", Group2 = "Classic", Group3 = "Mean",
            X = Mean[1], Y = Mean[2], nbClust = 0, Time = execTime
          ))

          ## Coordinate-wise median
          sTime <- as.numeric(Sys.time()) * 1000
          MedCW <- c(median(samples.XY[ ,1]), median(samples.XY[ ,2]))
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = "Med-CW", Group1 = "Classic", Group2 = "Classic", Group3 = "CW",
            X = MedCW[1], Y = MedCW[2], nbClust = 0, Time = execTime
          ))


          ## ========================================================================\
          ## -------- Depth Medians  ---------
          ## ========================================================================\


          for (method in params$bril.methods) {

            if (method %in% c("MVE","MCD")) {
              Group1 <- "Cov"
            } else {
              Group1 <- "Depth"
            }

            if (! method %in% c("MVE","MCD")) {

              ## Traditional Multivariate Median

              sTime <- as.numeric(Sys.time()) * 1000
              medMM <- BRIL::median_mv(samples.XY, method = method, sampleMedian = F)
              execTime <- as.numeric(Sys.time()) * 1000 - sTime

              results <- rbind(results, data.frame(
                Method = paste0("Med-", method),
                Group1 = Group1, Group2 = "Med", Group3 = method, nbClust = 0,
                X = medMM[1], Y = medMM[2], Time = execTime
              ))

              ## Sample Median

              sTime <- as.numeric(Sys.time()) * 1000
              medMax <- BRIL::median_mv(samples.XY, method = method, sampleMedian = T)
              execTime <- as.numeric(Sys.time()) * 1000 - sTime

              results <- rbind(results, data.frame(
                Method = paste0("Max-", method),
                Group1 = Group1, Group2 = "Max", Group3 = method, nbClust = 0,
                X = medMax[1], Y = medMax[2], Time = execTime
              ))

              ## Sup Median

              sTime <- as.numeric(Sys.time()) * 1000
              depths <- BRIL::depth_values(samples.XY, method = method)
              deepestSamples <- samples.XY[sort.int(depths, decreasing = TRUE, index.return = TRUE)$ix[1:floor(nrow(samples.XY) * 0.5)], ]
              MedSup <- BRIL::median_mv(deepestSamples, method = method, sampleMedian = T)
              execTime <- as.numeric(Sys.time()) * 1000 - sTime

              results <- rbind(results, data.frame(
                Method = paste0("Sup-", method),
                Group1 = Group1, Group2 = "Sup", Group3 = method, nbClust = 0,
                X = MedSup[1], Y = MedSup[2], Time = execTime
              ))
            }

            ## Rec Median

            sTime <- as.numeric(Sys.time()) * 1000
            medRec <- BRIL::median_rec(samples.XY, method = method, alpha = 0.5)
            execTime <- as.numeric(Sys.time()) * 1000 - sTime

            results <- rbind(results, data.frame(
              Method = paste0("Rec-", method),
              Group1 = Group1, Group2 = "Rec", Group3 = method, nbClust = 0,
              X = medRec$median[1], Y = medRec$median[2], Time = execTime
            ))

            ## BRL (refined)

            sTime <- as.numeric(Sys.time()) * 1000
            medBRL <- BRIL::bril(samples.XY, method = method, maxIterations = 1, alpha = params$bril.alpha, threshUnimodal = params$bril.unimodalThresh,
                                 testNormal = params$bril.testNormal, threshNormal = params$bril.normalThresh, distNormal = params$bril.distNormal,
                                 trimmedPerFilteringIteration = params$bril.trimmedPerFilteringIteration, warnings = params$bril.warnings)
            execTime <- as.numeric(Sys.time()) * 1000 - sTime

            results <- rbind(results, data.frame(
              Method = paste0("BRL-", method),
              Group1 = Group1, Group2 = "BRL", Group3 = method, nbClust = 0,
              X = medBRL$mode[1], Y = medBRL$mode[2], Time = execTime
            ))

            ## BRIL

            sTime <- as.numeric(Sys.time()) * 1000
            medBRIL <- BRIL::bril(samples.XY, method = method, alpha = params$bril.alpha, threshUnimodal = params$bril.unimodalThresh,
                                  testNormal = params$bril.testNormal, threshNormal = params$bril.normalThresh, distNormal = params$bril.distNormal,
                                  trimmedPerFilteringIteration = params$bril.trimmedPerFilteringIteration, warnings = params$bril.warnings)
            execTime <- as.numeric(Sys.time()) * 1000 - sTime
            eval(parse(text = paste0("ind", "BRIL", method, "<- which(medBRIL$labels == medBRIL$mainCluster)")))

            results <- rbind(results, data.frame(
              Method = paste0("BRIL-", method),
              Group1 = Group1, Group2 = "BRIL", Group3 = method, nbClust = medBRIL$nbClusters,
              X = medBRIL$mode[1], Y = medBRIL$mode[2], Time = execTime
            ))

            ## ========= Plot

            if (params$plotRes || params$savePlots ) {

              res <- medBRIL
              samples.df <- as.data.frame(samples.XY)
              colnames(samples.df) <- c("X", "Y")
              samples.df$cluster <- 2
              samples.df$category <- c(rep("Inlier",nbSamplesPerCluster[1]),rep("Outlier",sum(nbSamplesPerCluster)-nbSamplesPerCluster[1]), rep("Noise",noise.nb))
              samples.df$category <- factor(samples.df$category,levels = c("Inlier","Outlier","Noise"))
              samples.df$labels <- factor(res$labels)

              labelsCluster <- c()
              for (i in sort.int(unique(res$labels))) {
                labelsCluster[[as.character(i)]] <- sprintf("%d   (%d pts)", i, length(which(res$labels == i)))
              }
              labelsCluster[[as.character(res$mainCluster)]] <- paste(labelsCluster[[as.character(res$mainCluster)]], " *")


              g <- ggplot(samples.df, aes(X, Y, color = labels, group = labels, shape = category)) +
                geom_point(size=1.5, alpha=0.7) +
                geom_point(aes(x = res$mode[1], y = res$mode[2]), # Black X for BRIL
                           colour = "black", shape = 4, size = 2.5, stroke = 2)
              if (!method %in% c("MVE","MCD")) {
                g <- g + geom_point(aes(x = medMax[1], y = medMax[2]),  # Gray x for Max
                                    colour = "gray40", shape = 4, size = 1.5, stroke = 1.5, alpha=1)
              }
              g <- g + geom_point(aes(x = medRec$median[1], y = medRec$median[2]), # Gray + for Rec
                                  colour = "gray60", shape = 3, size = 1.5, stroke = 1.5, alpha=0.8) +
                geom_point(aes(x = medBRL$mode[1], y = medBRL$mode[2]),  # Gray o for BRL
                           colour = "gray80", shape = 1, size = 2, stroke = 1, alpha=0.6) +
                scale_color_discrete(labels = labelsCluster, name = "clusters") +
                coord_fixed(ratio = 1, expand = TRUE) +
                labs(title = paste0("BRIL-", method),
                     subtitle = sprintf("NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d",
                                        nbClusters, inliersRatio, sampleSize, noiseRatio))
              if (params$plotRes) {
                plot(g)
              }
              if (params$savePlots) {
                filename <- sprintf("Sim %d - BRIL-%s - Clusters %i Inliers %.1f SampleSize %i Noise %i.jpg",
                                    iterationsCompleted + montecarlo.i, method,  nbClusters, inliersRatio, sampleSize, noiseRatio)
                ggsave(g, filename = file.path(params$folderData, filename), width = 7, height = 6)
              }
            }
          }

          ## ========= Clustering  ===========================================
          ## ========================================================================\


          if (params$doClustering) {

            if (params$verbose) {
              print("Clustering")
            }

            source("utils_clustering.R", local = TRUE)

            for (method in params$clusteringMethods) {

              sTime <- as.numeric(Sys.time()) * 1000
              resClust <- modeClustering(samples.XY, method = method, k.max = params$clusteringKMax)
              execTime <- as.numeric(Sys.time()) * 1000 - sTime
              eval(parse(text = paste0("ind", method, "<- which(resClust$labels == resClust$mainCluster)")))

              results <- rbind(results, data.frame(
                Method = method,
                Group1 = "Clustering", Group2 = "Classic", Group3 = "Classic", nbClust = resClust$nbClusters,
                X = resClust$mode[1], Y = resClust$mode[2], Time = execTime
              ))


              ## ========= Plot Clustering resultss

              if (params$savePlots || (params$plotRes && params$plotClustering)) {

                res <- resClust
                samples.df <- as.data.frame(samples.XY)
                colnames(samples.df) <- c("X", "Y")
                samples.df$category <- c(rep("Inlier",nbSamplesPerCluster[1]),rep("Outlier",sum(nbSamplesPerCluster)-nbSamplesPerCluster[1]), rep("Noise",noise.nb))
                samples.df$category <- factor(samples.df$category,levels = c("Inlier","Outlier","Noise"))
                samples.df$labels <- factor(res$labels)

                labelsCluster <- c()
                for (i in sort.int(unique(res$labels))) {
                  labelsCluster[[as.character(i)]] <- sprintf("%d   (%d pts)", i, length(which(res$labels == i)))
                }
                labelsCluster[[as.character(res$mainCluster)]] <- paste(labelsCluster[[as.character(res$mainCluster)]], " *")


                g <- ggplot(samples.df, aes(X, Y, color = labels, group = labels, shape = category)) +
                  geom_point(size=1.5, alpha=0.7) +
                  geom_point(aes(x = res$mode[1], y = res$mode[2]), # Black X for mode
                             colour = "black", shape = 4, size = 2.5, stroke = 2) +
                  scale_color_discrete(labels = labelsCluster, name = "clusters") +
                  coord_fixed(ratio = 1, expand = TRUE) +
                  labs(title = method,
                       subtitle = sprintf("NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d",
                                          nbClusters, inliersRatio, sampleSize, noiseRatio))

                if (params$plotClustering) {
                  plot(g)
                }
                if (params$savePlots) {
                  filename <- sprintf("Sim %d - %s - Clusters %i Inliers %.1f SampleSize %i Noise %i.jpg",
                                      iterationsCompleted + montecarlo.i, method,  nbClusters, inliersRatio, sampleSize, noiseRatio)
                  ggsave(g, filename = file.path(params$folderData, filename), width = 7, height = 6)
                }
              }
            }
          }

          ## ========================================================================\
          ## -------- PROCESS RESULTS   ---------
          ## ========================================================================\

          if (params$verbose) {
            print("Storing resultss")
          }

          ## Add extra fields to the results

          results$error <- 0
          results$squared_error <- 0
          results$close <- TRUE
          results$inside <- FALSE
          results$fdr <- 0
          results$fnr <- 0

          results$inliersRatio <- inliersRatio
          results$nbClusters <- nbClusters
          results$samples.nb <- sampleSize
          results$noise.nb <- noise.nb
          results$noise <- noiseRatio
          results$montecarloRepetitions <- params$montecarloRepetitions
          results$montecarlo.i <- montecarlo.i
          results$Time <- sapply(results$Time, round)

          inliers.exist <- 1
          if (as.integer(inliersRatio) < as.integer(100 / nbClusters)) {
            inliers.exist <- -1
          } else if (inliersRatio != 100 && as.integer(inliersRatio) == as.integer(100 / nbClusters)) {
            inliers.exist <- 0
          }
          results$inliers.exist <- inliers.exist

          referencePosition <- clusters.position[1, ]

          for (i in 1:nrow(results)) {

            ## Add distances
            results$error[i] <- dist(rbind(referencePosition, results[i, c("X", "Y")]))[1]
            results$squared_error[i] <- results$error[i] * results$error[i]
            if (nbClusters > 1) {
              results$close[i] <- results$error[i] <= min(dist(rbind(as.numeric(results[i, c("X", "Y")]), clusters.position))[1:nbClusters])
            }

            results$inside[i] <- (results$error[i] < 3 * sqrt(clusters.variance[1]))

            if (results$Group1[i] == "Clustering" || results$Group2[i] == "BRIL") {
              listIndices <- eval(parse(text = paste0("ind", gsub("-", "", results$Method[i], fixed = TRUE))))
              results$fdr[i] <- 100.0 * length(which(listIndices > nbSamplesPerCluster[1])) / length(listIndices)
              results$fnr[i] <- 100.0 * (nbSamplesPerCluster[1] - length(which(listIndices <= nbSamplesPerCluster[1]))) / nbSamplesPerCluster[1]
            }
          }

          ## Convert Logicals Columns to Numeric
          cols <- sapply(results, is.logical)
          results[, cols] <- as.numeric(as.matrix(results[, cols]))

          # Remove values not required anymore
          results <- subset(results, select = -c(X, Y))

          ######################## Follow the Progress

          timeElapsed <- Sys.time() - start.time
          nbIterDone <- iterationsCompleted + montecarlo.i
          avgTimePerIteration <- timeElapsed / nbIterDone
          timeRemaining <- (totalIterations - nbIterDone) * avgTimePerIteration
          if (timeRemaining < as.difftime(60, format = "%S", units = "secs")) {
            units(timeRemaining) <- "secs"
          } else if (timeRemaining < as.difftime(60, format = "%M", units = "mins")) {
            units(timeRemaining) <- "mins"
          } else {
            units(timeRemaining) <- "hours"
          }

          cat(sprintf(
            "\n## End of repetiton %d/%d  : NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d (Simulation %d/%d - %s elapsed - %s remaining)\n\n",
            montecarlo.i, params$montecarloRepetitions, nbClusters, inliersRatio, sampleSize, noiseRatio, nbIterDone, totalIterations, format(timeElapsed, digits = 0), format(timeRemaining, digits = 0)
          ))

          if (params$progressBar && !params$multithread) {
            progress(montecarlo.i)
          }


          if (params$oneRun && !params$multithread) {
            # readline(prompt='Press [enter] to exit')
            # stop("Aborting due to 'oneRun' parameter")
            cat("\n---- Aborting due to 'oneRun' parameter ----\n\n")
            return()
          }

          # Sys.sleep(0.3)

          ## Aggregate results in traditional loop
          if (!FOREACHLOOP) {
            results.all <- rbind(results.all, results)
          }

          results
        }

        ## ========================================================================\
        ## -------- FINALIZE UNIT AND WRITE FILE   ---------
        ## ========================================================================\

        iterationsCompleted <- iterationsCompleted + params$montecarloRepetitions

        timeElapsed <- Sys.time() - start.time
        nbIterDone <- iterationsCompleted
        avgTimePerIteration <- timeElapsed / nbIterDone
        timeRemaining <- (totalIterations - nbIterDone) * avgTimePerIteration
        if (timeRemaining < as.difftime(60, format = "%S", units = "secs")) {
          units(timeRemaining) <- "secs"
        } else if (timeRemaining < as.difftime(60, format = "%M", units = "mins")) {
          units(timeRemaining) <- "mins"
        } else {
          units(timeRemaining) <- "hours"
        }

        cat(sprintf(
          "\n#### End of unit %d/%d: NbCluster=%i, Inliers=%.1f, Size=%d, Noise=%d (%d repetitions)",
          iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions, nbClusters, inliersRatio, sampleSize, noiseRatio, params$montecarloRepetitions),
          "\n#### Time elapsed:", format(timeElapsed, digits = 2),", remaining: ", format(timeRemaining, digits = 2),"\n")


        ## Stop the cluster
        if (params$multithread) {
          stopCluster(cl)
        }

        ## Output to file
        if (params$writeRes) {
          filename = sprintf("results Clusters %i Inliers %.1f SampleSize %i Noise %i.csv",
                             nbClusters, inliersRatio, sampleSize, noiseRatio)
          write.table(results.all, file = file.path(params$folderData, filename),
                      row.names = F, col.names = T, sep = ";"
          )
          cat("\nSaving results to file \"",filename,"\"...\n\n", sep = "")
        }

        if (params$oneRun && !params$multithread) {
          cat("\n---- Aborting due to 'oneRun' parameter ----\n\n")
          return()
        }
      }
    }
  }
}

if (params$progressBar && !params$consoleMode ) {
  setTkProgressBar(globalProgressBar,iterationsCompleted,
                   title = sprintf("%d/%d - Unit %d/%d (Repetition %d/%d)",
                                   iterationsCompleted, totalIterations,
                                   iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                                   params$montecarloRepetitions,params$montecarloRepetitions),
                   label = "Benchmark completed.")
}
