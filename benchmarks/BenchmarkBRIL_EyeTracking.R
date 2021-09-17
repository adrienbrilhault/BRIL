
########################################################################### \
## Project: BRIL
## Script purpose: Create, test and save BRIL benchmarks with eye-tracking data
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
  install.packages("CompQuadForm")
  install.packages("remotes")
  remotes::install_github("adrienbrilhault/BRIL", subdir = "pkg")
}

library(BRIL)
library(ggplot2)
# library(doParallel)
library(doSNOW)
library(tcltk) # For progress bar widget

set.seed(123)
closeAllConnections()

## ========================================================================\
## -------- Parameters  --------
## ========================================================================\

params = list()

######### Path to save the benchmark resultss

params$folderOutput <- "Results/RunEyetrackingData"
params$datasetFilepath <- "../data/eyeTrackingCalibrations.Rdata"

######### Misc

params$writeRes <- T # Output results to file
params$plotRes <- F # Display plots
params$savePlots <- F # Save plots de jpg files in result folder
params$consoleMode <- F # disables plots, progress widget, and all graphic outputs
params$progressBar <- T # display progress bar
params$verbose <- T # Verbose
params$debug <- T # Run some additional debugging code

######### Number of repetitions

params$oneRun <- F # Only one simulation (for debug)
params$montecarloRepetitions <- 20 # Number of simulations
params$sampleSizeValues <- 1000

######### Parallel simulations

params$multithread <- T # Set to true to parallelize the simulations
params$outputToFiles <- T # Redirect standard output to file in the worker
params$nbCores <- parallel::detectCores()
# params$nbCores <- 4

######### Selected session

params$sessionsBlacklist <- c()
params$sessionsWhitelist <- c()
params$sessionsWhitelist <- c("juj011a00","ded00800","ded005a01")

######### Clustering algorithms

params$doClustering <- T # Run clustering algorithms
params$clusteringKMax <- 10
params$plotClustering <- T
params$clusteringMethods <- c("PAM", "MClust", "DBSCAN")

######### BRIL tuning parameters

params$bril.unimodalThresh <- 0.05
params$bril.testNormal <- "Chisq"
params$bril.normalThresh <- 0.05
params$bril.distNormal <- "MCD"
params$bril.alpha <- 0.5
params$bril.trimmedPerFilteringIteration <- 1
params$bril.warnings <- FALSE
# params$bril.methods =  c("MVE","MCD","Tukey", "Liu", "Oja", "Projection", "Spatial")
# params$bril.methods =  c("MVE","MCD","Oja","Projection")
# params$bril.methods =  c("MVE","MCD","Projection","L2")
# params$bril.methods =  c("MVE","L2")
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

## ========================================================================\
## -------- BENCHMARK  --------
## ========================================================================\

## -------- Load data  --------

load(file = params$datasetFilepath)

sessions <- unique(metadata$session)
if (!is.null(params[["sessionsWhitelist"]])){
  sessions <- intersect(sessions, params$sessionsWhitelist)
}
if (!is.null(params[["sessionsBlacklist"]])){
  sessions <- setdiff(sessions, params$sessionsBlacklist)
}
params$sessions <- sessions

if (params$writeRes && !file.exists(params$folderOutput)) {
  dir.create(file.path(params$folderOutput))
}
if (params$writeRes) {
  save(params, file = file.path(params$folderOutput,paste0("params_",format(Sys.time(), "%Y%m%d_%H%M%S"),".Rdata")))
}


totalIterations <- sum(sapply(sessions, function(x) length(unique(metadata$target[metadata$session == x])))) *
  length(params$sampleSizeValues) * params$montecarloRepetitions
iterationsCompleted <- 0

if (params$progressBar && !params$consoleMode ) {
  globalProgressBar <- tkProgressBar(max = totalIterations, width = 550, initial = 1,
                                     title = sprintf("%d/%d - Session %d/%d - Target %d/%d (Repetition %d/%d)",
                                                     1, totalIterations, 1, length(sessions),
                                                     1, metadata$nTargets[metadata$session == sessions[1]][1],
                                                     1, params$montecarloRepetitions))
}

start.time <- Sys.time()


for (sampleSize in params$sampleSizeValues) {

  for (session in sessions) {

    targets <- sort(unique(metadata$target[metadata$session == session]))

    result.session <- data.frame()

    cat(sprintf("\n@@@@@@@@@@@@@@ Start of Session %d/%d: \"%s\" (sample size %d) @@@@@@@@@@@@@@\n",
                which(session == sessions), length(sessions), session, sampleSize ))

    for (target in targets) {

      cat("\n##################################################################################################\n",
          sprintf("#### Start of unit %d/%d : Session %s (%d/%d) - Sample Size %d (%d/%d) - Target %d (%d/%d)",
                  1+iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                  session, which(session==sessions), length(sessions),
                  sampleSize, which(params$sampleSizeValues==sampleSize), length(params$sampleSizeValues),
                  target, which(targets==target), length(targets)),
          "\n##################################################################################################\n\n", sep = "")

      # Filter the corresponding trials
      selectedIndices <- metadata$session == session & metadata$target == target
      samples <- cbind(as.vector(t(x[selectedIndices, ])), as.vector(t(y[selectedIndices, ])))

      # Remove trailing NA fill
      samplesFiltered <- samples[!is.na(samples[ ,1]) | !is.na(samples[ ,2]), ]

      # Remove overflowed values and error codes from eyetracker
      samplesFiltered <- samplesFiltered[samplesFiltered[ ,1] != -32768 & samplesFiltered[ ,1] != -32768, ]
      samplesFiltered <- samplesFiltered[samplesFiltered[ ,1] != 32767 & samplesFiltered[ ,1] != 32767, ]

      ## Structure to store the results
      results.block<- data.frame()

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
                           label = sprintf("Current Unit: SampleSize=%i, Session=%s, Target=%d",
                                           sampleSize, session, target))

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

        cl <- makeCluster(params$nbCores, outfile = "")

        # registerDoParallel(cl)
        registerDoSNOW(cl)

      } else {
        registerDoSEQ()
      }

      listPackages <- c("ggplot2")

      FOREACHLOOP <- T # set this to false if you switch below the foreach loop for a classic one

      # for (montecarlo.i in 1:params$montecarloRepetitions) {
      results.block <- foreach(montecarlo.i = (1:ifelse(params$oneRun, 1, params$montecarloRepetitions)), .combine = rbind, .errorhandling = "remove", .options.snow = opts, .packages = listPackages, .inorder = F) %dopar% {

        ## No output on parrallel worker, redirect output to file to monitor progress
        if (params$multithread & params$outputToFiles) {
          sinkFilename <- sprintf("Simulation %03d of %03d - Session %s - Target %d - SampleSize %d - Repetiton %d of %d.txt",
                                  iterationsCompleted+montecarlo.i, totalIterations, session, target, sampleSize, montecarlo.i, params$montecarloRepetitions)
          sink(file.path(params$folderOutput, file=sinkFilename))

        }

        cat(sprintf(
          "\n#### Session %s (%d/%d) -  Sample Size %d (%d/%d) - Target %d (%d/%d) - Repetition %d/%d\n\n",
          session, which(sessions==session), length(sessions),
          sampleSize, which(params$sampleSizeValues==sampleSize), length(params$sampleSizeValues),
          target, which(targets==target), length(targets),
          montecarlo.i, ifelse(params$oneRun, 1, params$montecarloRepetitions)
        ))

        ## ========= Draw random samples  ============

        # Draw set of random samples of size sampleSize
        samples.XY <- samplesFiltered[sample(seq_len(nrow(samplesFiltered)), sampleSize), ]

        ## Save data for debug
        if (params$multithread & params$outputToFiles) {
          save(samples.XY,session,target,montecarlo.i,iterationsCompleted,
               file = file.path(params$folderOutput, sprintf("%s - DATA.Rdata", sub(".txt","",sinkFilename))))
        }

        ## ========================================================================\
        ## -------- COMPUTE MEDIANS  --------
        ## ========================================================================\

        if (params$verbose) {
          cat("Computing estimates\n")
        }

        ## ========= Simple medians  ============

        results <- data.frame(
          Method = character(), X = double(), Y = double(),
          Group1 = character(), Group2 = character(), Group3 = character(), Time = integer()
        )

        ## Mean of all samples
        sTime <- as.numeric(Sys.time()) * 1000
        Mean <- c(mean(samples.XY[, 1]), mean(samples.XY[, 2]))
        execTime <- as.numeric(Sys.time()) * 1000 - sTime

        results <- rbind(results, data.frame(
          Method = "Mean", Group1 = "Classic", Group2 = "Classic", Group3 = "Mean",
          X = Mean[1], Y = Mean[2], Time = execTime
        ))

        ## Coordinate-wise median
        sTime <- as.numeric(Sys.time()) * 1000
        MedCW <- c(median(samples.XY[ ,1]), median(samples.XY[ ,2]))
        execTime <- as.numeric(Sys.time()) * 1000 - sTime

        results <- rbind(results, data.frame(
          Method = "Med-CW", Group1 = "Classic", Group2 = "Classic", Group3 = "CW",
          X = MedCW[1], Y = MedCW[2], Time = execTime
        ))


        ## ========================================================================\
        ## -------- Depth Medians  ---------
        ## ========================================================================\

        for (method in params$bril.methods) {

          cat("-- BRIL", method, "\n")

          if (method %in% c("MVE","MCD")) {
            Group1 <- "Cov"
          } else {
            Group1 <- "Depth"
          }

          if (! method %in% c("MVE","MCD")) {

            ## Sample Median
            cat("---- Sample Median\n")

            sTime <- as.numeric(Sys.time()) * 1000
            medMax <- BRIL::median_mv(samples.XY, method = method, sampleMedian = T)
            execTime <- as.numeric(Sys.time()) * 1000 - sTime

            results <- rbind(results, data.frame(
              Method = paste0("Max-", method),
              Group1 = Group1, Group2 = "Max", Group3 = method,
              X = medMax[1], Y = medMax[2], Time = execTime
            ))

          }

          ## Rec Median
          cat("---- Rec Median\n")

          sTime <- as.numeric(Sys.time()) * 1000
          medRec <- BRIL::median_rec(samples.XY, method = method, alpha = 0.5)
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = paste0("Rec-", method),
            Group1 = Group1, Group2 = "Rec", Group3 = method,
            X = medRec$median[1], Y = medRec$median[2], Time = execTime
          ))

          ## BRL Median (refined)
          cat("---- BRL Median\n")

          sTime <- as.numeric(Sys.time()) * 1000
          medBRL <- BRIL::bril(samples.XY, method = method, maxIterations = 1, alpha = params$bril.alpha, threshUnimodal = params$bril.unimodalThresh,
                               testNormal = params$bril.testNormal, threshNormal = params$bril.normalThresh, distNormal = params$bril.distNormal,
                               trimmedPerFilteringIteration = params$bril.trimmedPerFilteringIteration, warnings = params$bril.warnings)
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = paste0("BRL-", method),
            Group1 = Group1, Group2 = "BRL", Group3 = method,
            X = medBRL$mode[1], Y = medBRL$mode[2], Time = execTime
          ))

          ## BRIL Median
          cat("---- BRIL Median\n")

          sTime <- as.numeric(Sys.time()) * 1000
          medBRIL <- BRIL::bril(samples.XY, method = method, maxIterations = 5, alpha = params$bril.alpha, threshUnimodal = params$bril.unimodalThresh,
                                testNormal = params$bril.testNormal, threshNormal = params$bril.normalThresh, distNormal = params$bril.distNormal,
                                trimmedPerFilteringIteration = params$bril.trimmedPerFilteringIteration, warnings = params$bril.warnings)
          execTime <- as.numeric(Sys.time()) * 1000 - sTime

          results <- rbind(results, data.frame(
            Method = paste0("BRIL-", method),
            Group1 = Group1, Group2 = "BRIL", Group3 = method,
            X = medBRIL$mode[1], Y = medBRIL$mode[2], Time = execTime
          ))

        } # END LOOP BRIL METHOD

        ## ========= Clustering  ===========================================
        ## ========================================================================\

        if (params$doClustering) {

          if (params$verbose) {
            cat("Clustering\n")
          }

          source("utils_clustering.R", local = TRUE)

          for (method in params$clusteringMethods) {

            cat("-- Clustering", method, "\n")

            sTime <- as.numeric(Sys.time()) * 1000
            resClust <- modeClustering(samples.XY, method = method, k.max = params$clusteringKMax)
            execTime <- as.numeric(Sys.time()) * 1000 - sTime

            results <- rbind(results, data.frame(
              Method = method,
              Group1 = "Clustering", Group2 = "Classic", Group3 = "Classic",
              X = resClust$mode[1], Y = resClust$mode[2], Time = execTime
            ))
          }
        }

        ## ========================================================================\
        ## -------- PROCESS RESULTS   ---------
        ## ========================================================================\

        if (params$verbose) {
          cat("Storing results\n")
        }

        ## Add extra fields to the results
        results$errorPx <- 0
        results$errorDeg <- 0
        results$errorRaw <- 0
        results$squared_error <- 0

        results$session <- session
        results$target <- target
        results$samples.nb <- sampleSize
        results$montecarloRepetitions <- params$montecarloRepetitions
        results$montecarlo.i <- montecarlo.i
        results$Time <- sapply(results$Time, round)

        ## Compute the errors
        referencePosition <- groundTruth[groundTruth$session == session & groundTruth$target == target,3:4]
        sessionGroundTruth <- groundTruth[groundTruth$session == session,]
        distances <- as.matrix(stats::dist(sessionGroundTruth[order(sessionGroundTruth$target), 3:4]))
        correction <- sqrt(100^2 + 100^2) / mean(distances[2:nrow(distances),1])
        pixeltodegrees <- 1 / 25.470

        for (i in 1:nrow(results)) {

          ## Add distances
          results$errorRaw[i] <- as.numeric(dist(rbind(referencePosition, results[i, c("X", "Y")]))[1])
          results$errorPx[i] <- results$errorRaw[i] * correction
          results$errorDeg[i] <- results$errorPx[i]*pixeltodegrees
          results$squared_error[i] <- results$errorPx[i] * results$errorPx[i]
        }

        ######################## Display the Progress

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
          "\n -> End of repetiton %d/%d  : SampleSize %d, Session \"%s\", Target %d (Simulation %d/%d - %s elapsed - %s remaining)\n\n",
          montecarlo.i, params$montecarloRepetitions, sampleSize, session, target, nbIterDone, totalIterations, format(timeElapsed, digits = 0), format(timeRemaining, digits = 0)
        ))

        if (params$progressBar && !params$multithread) {
          progress(montecarlo.i)
        }

        ## Abort if only one run
        if (params$oneRun && !params$multithread) {
          # readline(prompt='Press [enter] to exit')
          # stop("Aborting due to 'oneRun' parameter")
          cat("\n---- Aborting due to 'oneRun' parameter ----\n\n")
          return()
        }

        ## Monitor and call garbage collector
        if (params$outputToFiles == TRUE && params$multithread == TRUE) {
          sink()
          mem_usage <- toString(gc())
          cat(mem_usage, file = file.path(params$folderOutput, file=sprintf("%s FINISHED ( %s elapsed - %s remaining).txt",
                                                                            sub(".txt","",sinkFilename),format(timeElapsed, digits = 0), format(timeRemaining, digits = 0))))
        }
        Sys.sleep(0.3)

        ## Aggregate results if standard "for" loop
        if (!FOREACHLOOP) {
          results.block <- rbind(results.block, results)
        }

        ## Return results if "foreach" loop
        results

      } # END FOREACH LOOP ON MONTE CARLO REPETITIONS


      ## ========================================================================\
      ## -------- FINALIZE UNIT AND WRITE FILE   ---------
      ## ========================================================================\

      iterationsCompleted = iterationsCompleted + params$montecarloRepetitions

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

      cat("\n##################################################################################################\n",
          sprintf("#### End of unit %d/%d : Session %s (%d/%d) -  Sample Size %d (%d/%d) - Target %d (%d/%d)",
                  iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                  session, which(sessions==session), length(sessions),
                  sampleSize, which(params$sampleSizeValues==sampleSize), length(params$sampleSizeValues),
                  target, which(targets==target), length(targets)),
          "\n####     Time elapsed: ", format(timeElapsed, digits = 2),", remaining: ", format(timeRemaining, digits = 2),
          "\n##################################################################################################\n", sep = "")


      ## Stop the cluster
      if (params$multithread) {
        stopCluster(cl)
      }

      ## Output to file
      if (params$writeRes) {
        filename = sprintf("results SampleSize %i Session %s Target %d.csv",
                           sampleSize, session, target)
        write.table(results.block, file = file.path(params$folderOutput, filename),
                    row.names = F, col.names = T, sep = ";"
        )
        cat("\nSaving results to file \"",filename,"\"...\n\n", sep = "")
      }

      if (params$oneRun && !params$multithread) {
        cat("\n---- Aborting due to 'oneRun' parameter ----\n\n")
        return()
      }

    } # END LOOP TARGETS
  } # END LOOP SESSIONS
} # END LOOP SAMPLE SIZES



## End the progress bar when the benchmark is completed

if (params$progressBar && !params$consoleMode ) {
  setTkProgressBar(globalProgressBar,iterationsCompleted,
                   title = sprintf("%d/%d - Unit %d/%d (Repetition %d/%d)",
                                   iterationsCompleted, totalIterations,
                                   iterationsCompleted/params$montecarloRepetitions, totalIterations/params$montecarloRepetitions,
                                   params$montecarloRepetitions,params$montecarloRepetitions),
                   label = "Benchmark completed.")
}


