
########################################################################### \
## Project: BRIL
## Script purpose: Analyze results from BRIL benchmark on eye-tracking data
## Author: Adrien Brilhault
########################################################################### \


if (F) {
  install.packages("tidyverse")
  install.packages("xlsx")

}

library(tidyverse)
library(xlsx)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(ggnewscale)
library(xlsx)

theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## ========================================================================\
## -------- Read data  --------
## ========================================================================\

######### Load dataset

load("../Data/eyeTrackingCalibrations.Rdata")

######### Load data from csv files or Rdata

folderResult <- "Results/RunEyetrackingData"

forceReload <- FALSE

if (!file.exists(file.path(folderResult, "resultsBenchmarkEyeTracking.Rdata")) || forceReload) {

  ######### Read each csv file in folder

  filelist <- list.files(path = folderResult, pattern = "csv$")

  df <- data.frame()

  for (filename in filelist) {

    filename.full <- file.path(folderResult, filename)
    cat("Reading",filename.full,"\n")

    ## Read data from file
    dfTemp <- read.csv(filename.full, header = T, sep = ";",stringsAsFactors = TRUE)

    ## Merge
    df <- rbind(df, dfTemp)
  }

  ## Save to Rdata
  results <- df
  save(results, file = file.path(folderResult, "resultsBenchmarkEyeTracking.Rdata"))


} else {

  ######### Load results benchmark

  load(file.path(folderResult, "resultsBenchmarkEyeTracking.Rdata"))
}




## ========================================================================\
## -------- Analyze Results  --------
## ========================================================================\


if (FALSE) {

  ## Check contents
  sort(unique(results$session))
  sort(unique(results$target))
  sort(unique(results$samples.nb))
  sort(unique(results$montecarlo.i))
  sort(unique(results$Method))
}

######### Compute summary statistics



## Results for each method and session, and across sessions
resultAnalysis <- bind_rows(
  results %>%
    filter(Group2 %in% c("BRIL","Classic")) %>%
    group_by(Method, session) %>%
    summarise(error = mean(errorPx), sd_error = sd(errorPx), sum_squarred_error = as.integer(sum(squared_error))),
  results %>%
    filter(Group2 %in% c("BRIL","Classic")) %>%
    group_by(Method) %>%
    summarise(session= "all", error = mean(errorPx), sd_error = sd(errorPx), sum_squarred_error = as.integer(sum(squared_error)))
)


## Table for the 3 main datasets expanded in columns
resultAnalysis <- resultAnalysis %>% arrange(Method)
finalSummary <- data.frame(Method = unique(resultAnalysis$Method))
for (i in 1:nrow(finalSummary)) {

  finalSummary$All_error[i] <- resultAnalysis$error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "all"]
  finalSummary$All_sd[i] <- resultAnalysis$sd_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "all"]
  finalSummary$All_sse[i] <- resultAnalysis$sum_squarred_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "all"]

  finalSummary$Set1_error[i] <- resultAnalysis$error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "juj011a00"]
  finalSummary$Set1_sd[i] <- resultAnalysis$sd_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "juj011a00"]
  finalSummary$Set1_sse[i] <- resultAnalysis$sum_squarred_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "juj011a00"]

  finalSummary$Set2_error[i] <- resultAnalysis$error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded00800"]
  finalSummary$Set2_sd[i] <- resultAnalysis$sd_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded00800"]
  finalSummary$Set2_sse[i] <- resultAnalysis$sum_squarred_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded00800"]

  finalSummary$Set3_error[i] <- resultAnalysis$error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded005a01"]
  finalSummary$Set3_sd[i] <- resultAnalysis$sd_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded005a01"]
  finalSummary$Set3_sse[i] <- resultAnalysis$sum_squarred_error[resultAnalysis$Method == finalSummary$Method[i] & resultAnalysis$session == "ded005a01"]
}

finalSummary <- finalSummary %>% arrange(All_error)
.Last.value %>% View()
finalSummary

## Save to file
if (F) {
  write.xlsx(as.data.frame(finalSummary), file.path(folderResult, "results_eyetracking_fixations.xlsx"), sheetName = "Summary",
             col.names = TRUE, row.names = FALSE, append = FALSE)
}

## -------- Extra analysis

if (F) {

  ## Results per Session
  resultAnalysis <- results %>%
    filter(Group2 %in% c("BRIL","Classic")) %>%
    group_by(Method, session) %>%
    summarise(montecarlo.i = "all", group1 = first(Group1),group2 = first(Group2),group3 = first(Group3),
              error = mean(errorPx), sd_error = sd(errorPx), sum_squarred_error = as.integer(sum(squared_error))) %>%
    arrange(session, sum_squarred_error)
  .Last.value %>% View()


  ## Results per session, target and Monte-Carlo repetition
  resultAnalysis <- results %>%
    filter(Method %in% c("PAM","BRIL-Projection")) %>%
    group_by(Method, session, target, montecarlo.i) %>%
    summarise(group1 = first(Group1),group2 = first(Group2),group3 = first(Group3),
              error = mean(errorPx))
  .Last.value %>% View()

  ## Results per session and target
  resultAnalysis <- results %>%
    filter(Method %in% c("PAM","BRIL-Projection")) %>%
    # group_by(Method, session, montecarlo.i) %>%
    group_by(Method, session, target) %>%
    summarise(x=mean(X),y=mean(Y),
              error = mean(errorPx), sd_error = sd(errorPx), sum_squarred_error = as.integer(sum(squared_error))) %>%
    arrange(session,target)
  .Last.value %>% View()
}


## ========================================================================\
## -------- Plot Final Calibrations  --------
## ========================================================================\

if (F) {
  showFigures <- T
  saveFigures <- T

  plotIndividualCalibs <- TRUE
  plotGrouplCalibs <- TRUE

  resultsPlot <- results

  sessions <- sort(unique(resultsPlot$session))
  sessions <- c("ded005a01")
  sessions <- c("juj011a00", "ded00800", "ded005a01")

  listmethods <- c("BRIL-Projection")
  iteration <- 1 # montecarlo repetition used

  ## Load if needed
  if (F) {
    load("Results/RunEyetrackingData/resultsBenchmarkEyeTracking.Rdata")
    load("../Data/eyeTrackingCalibrations.Rdata")
  }

  dfCalibration <- list()
  dfCalibrationPlot <- list()
  groundTruthPlot <- list()
  xLimGlobal <- list()
  yLimGlobal <- list()

  for (session in sessions) {

    sessionGroundTruth <- groundTruth[groundTruth$session == session, ]
    distances <- as.matrix(stats::dist(sessionGroundTruth[order(sessionGroundTruth$target), 3:4]))
    correction <- sqrt(100^2 + 100^2) / mean(distances[2:nrow(distances), 1])
    groundTruthPlot[[session]] <- sessionGroundTruth

    ## -------- retrieve samples from the calibration session

    dfCalibration[[session]] <- data.frame()
    targets <- sort(unique(metadata$target[metadata$session == session]))

    ## Pool all samples of the session, per target, filtering error values
    for (target in targets) {
      selectedIndices <- which(metadata$session == session & metadata$target == target)
      df <- data.frame(
        X = as.vector(t(x[selectedIndices, ])),
        Y = as.vector(t(y[selectedIndices, ]))
      )
      df <- df[!is.na(df$X) | !is.na(df$Y), ]
      df <- df[df$X != -32768 & df$Y != -32768, ]
      df <- df[df$X != 32767 & df$Y != 32767, ]

      df$Target <- target
      dfCalibration[[session]] <- rbind(dfCalibration[[session]], df)
    }

    dfCalibration[[session]]$Target <- factor(dfCalibration[[session]]$Target, targets)

    ## correct coordinates to screen referencial
    dfCalibration[[session]]['X'] <- (dfCalibration[[session]]['X'] - sessionGroundTruth[sessionGroundTruth$target==1,'X'])*correction
    dfCalibration[[session]]['Y'] <- (dfCalibration[[session]]['Y'] - sessionGroundTruth[sessionGroundTruth$target==1,'Y'])*correction
    groundTruthPlot[[session]]['X'] <- (groundTruthPlot[[session]]['X'] - sessionGroundTruth[sessionGroundTruth$target==1,'X'])*correction
    groundTruthPlot[[session]]['Y'] <- (groundTruthPlot[[session]]['Y'] - sessionGroundTruth[sessionGroundTruth$target==1,'Y'])*correction
    resultsPlot[resultsPlot$session == session, 'Y'] <- (resultsPlot[resultsPlot$session == session, 'Y'] - sessionGroundTruth[sessionGroundTruth$target==1,'Y'])*correction
    resultsPlot[resultsPlot$session == session, 'X'] <- (resultsPlot[resultsPlot$session == session, 'X'] - sessionGroundTruth[sessionGroundTruth$target==1,'X'])*correction


    set.seed(123)
    dfCalibrationPlot[[session]] <- sample_n(dfCalibration[[session]], min(nrow(dfCalibration[[session]]), 6000))

    ## -------- Set the plotting area for the session (based on samples and on groundtruth)

    xLim <- c(
      min(groundTruthPlot[[session]]$X) - (max(groundTruthPlot[[session]]$X) - min(groundTruthPlot[[session]]$X)) * 0.25,
      max(groundTruthPlot[[session]]$X) + (max(groundTruthPlot[[session]]$X) - min(groundTruthPlot[[session]]$X)) * 0.25
    )
    yLim <- c(
      min(groundTruthPlot[[session]]$Y) - (max(groundTruthPlot[[session]]$Y) - min(groundTruthPlot[[session]]$Y)) * 0.25,
      max(groundTruthPlot[[session]]$Y) + (max(groundTruthPlot[[session]]$Y) - min(groundTruthPlot[[session]]$Y)) * 0.25
    )
    xLimGlobalSession <- as.numeric(quantile(dfCalibration[[session]]$X, c(0.001, 0.999)))
    yLimGlobalSession <- as.numeric(quantile(dfCalibration[[session]]$Y, c(0.001, 0.999)))

    widthPlot <- 600
    percentile <- 0
    while (diff(as.numeric(quantile(dfCalibration[[session]]$X, c(percentile, 1-percentile)))) > widthPlot-200)
      percentile <- percentile + 0.001
    xLimGlobal[[session]] <- c(round(as.numeric(quantile(dfCalibration[[session]]$X,percentile)-100)/50)*50,
                               widthPlot + round(as.numeric(quantile(dfCalibration[[session]]$X,percentile))/50)*50)

    heightPlot <- 800
    percentile <- 0
    while (diff(as.numeric(quantile(dfCalibration[[session]]$Y, c(percentile, 1-percentile)))) > heightPlot-200)
      percentile <- percentile + 0.001
    yLimGlobal[[session]] <- c(round(as.numeric(quantile(dfCalibration[[session]]$Y,percentile)-100)/50)*50,
                               heightPlot + round(as.numeric(quantile(dfCalibration[[session]]$Y,percentile))/50)*50)


    ## -------- Individual Plot

    if (plotIndividualCalibs && (showFigures || saveFigures)) {

      for (method in listmethods) {

        index <- which(resultsPlot$session == session & resultsPlot$Method == method & resultsPlot$montecarlo.i == iteration )
        resultModes <- resultsPlot[index, c("X", "Y", "target")]
        names(resultModes)[names(resultModes) == "target"] <- "Target"
        resultModes$Target <- factor(resultModes$Target, levels(dfCalibration[[session]]$Target))
        posX <- xLimGlobal[[session]][2] - (xLimGlobal[[session]][2] %% 100)
        posY <- yLimGlobal[[session]][1] + 25
        # posY <- yLimGlobal[[session]][1] + (yLimGlobal[[session]][1] %% 100)

        g <- ggplot(dfCalibrationPlot[[session]], aes(X, Y, color = Target, group = Target)) +
          geom_point(shape = 1, alpha = 0.9) +
          scale_colour_hue(c = 40, l = 80) +
          geom_point(data = groundTruthPlot[[session]], aes(X, Y), colour = "gray30", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          new_scale_color() +
          scale_colour_hue(c = 60, l = 40) +
          geom_point(data = resultModes, aes(X, Y, colour = Target), size = 3, stroke = 2, shape = 1, inherit.aes = FALSE) +
          # coord_fixed(ratio = 1, expand = TRUE, xlim = xLim, ylim = yLim) +
          # coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobalSession, ylim = yLimGlobalSession) +
          coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobal[[session]], ylim = yLimGlobal[[session]]) +
          # coord_fixed(ratio=1, expand = TRUE) +
          theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
          theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
          scale_x_continuous(breaks = seq(-2000,2000, by = 200))+
          scale_y_continuous(breaks = seq(-2000,2000, by = 200))+
          # ggtitle(paste0("Session ", session, " - ", method))
          annotate("errorbarh", xmin = posX-100, xmax = posX, y = posY, colour = "gray40",
                   size = 1, height = 20) +
          annotate("text", label = "100 px",x = posX-50, y = posY+10,
                   colour = "black", vjust = "inward", hjust = "center", fontface = "italic", size = 3.5)+
          ggtitle(paste0("Session ", session))

        if (showFigures)
          print(g)

        ## Save figure
        if (saveFigures) {
          sizeFig = c(6,6)

          ggsave(paste0("Figures/pdf/Calibration ", session, " - ", method, ".pdf"), g, width = sizeFig[1], height =  sizeFig[2])
          ggsave(paste0("Figures/Calibration ", session, " - ", method, ".jpg"), g, width = sizeFig[1], height =  sizeFig[2])
          ggsave(paste0("Figures/pdf/Zoom - Calibration ", session, " - ", method, ".pdf"), g + coord_fixed(xlim = xLim, ylim = yLim), width = sizeFig[1], height =  sizeFig[2])
          ggsave(paste0("Figures/Zoom - Calibration ", session, " - ", method, ".jpg"), g + coord_fixed(xlim = xLim, ylim = yLim), width = sizeFig[1], height =  sizeFig[2])
        }
      }
    }
  }

  ## -------- Global Plot

  if (plotGrouplCalibs && (showFigures || saveFigures)) {

    for (method in listmethods) {

      groupedPlots <- list()

      for (session in sessions) {

        posX <- xLimGlobal[[session]][2] - (xLimGlobal[[session]][2] %% 100)
        posY <- yLimGlobal[[session]][1] + 25
        # posY <- yLimGlobal[[session]][1] + (yLimGlobal[[session]][1] %% 100)

        groupedPlots[[session]] <-
          ggplot(dfCalibrationPlot[[session]], aes(X, Y, color = Target, group = Target)) +
          geom_point(shape = 1, alpha = 0.9) +
          scale_colour_hue(c = 40, l = 80) +
          geom_point(data = groundTruthPlot[[session]], aes(X, Y), colour = "gray30", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          new_scale_color() +
          scale_colour_hue(c = 60, l = 40) +
          geom_point(data = resultModes, aes(X, Y, colour = Target), size = 3, stroke = 2, shape = 1, inherit.aes = FALSE) +
          # coord_fixed(ratio = 1, expand = TRUE, xlim = xLim, ylim = yLim) +
          # coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobalSession, ylim = yLimGlobalSession) +
          coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobal[[session]], ylim = yLimGlobal[[session]]) +
          # coord_fixed(ratio=1, expand = TRUE) +
          theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
          theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
          scale_x_continuous(breaks = seq(-2000,2000, by = 200))+
          scale_y_continuous(breaks = seq(-2000,2000, by = 200))+
          # ggtitle(paste0("Session ", session, " - ", method)) +
          ggtitle(paste0("Session ", session)) +
          annotate("errorbarh", xmin = posX-100, xmax = posX, y = posY, colour = "gray40",
                   size = 1, height = 20) +
          annotate("text", label = "100 px",x = posX-50, y = posY+10,
                   colour = "black", vjust = "inward", hjust = "center", fontface = "italic", size = 3.5)+
          theme(legend.position="none")

        if (showFigures) {
          print(groupedPlots[[session]])
        }

        ## Save figure
        if (saveFigures) {
          sizeFig = c(6,6)

          sessionPlot <- groupedPlots[[session]] + ggtitle("")
          if (session == sessions[length(sessions)]) {
            sessionPlot <- sessionPlot + theme(legend.position="right")
            sizeFig = c(7,6)
          }
          ggsave(paste0("Figures/pdf/Final Calibration Session ", session," - ",method, ".pdf"),
                 sessionPlot, width = sizeFig[1], height = sizeFig[2])
          ggsave(paste0("Figures/Final Calibration Session ", session," - ",method, ".jpg"),
                 sessionPlot, width = sizeFig[1], height = sizeFig[2])
        }
      }

      sizeFig = c(18,6)

      fullFigure <- grid.arrange(
        groupedPlots[[1]],groupedPlots[[2]],groupedPlots[[3]],ncol = 3)

      ## Save figure
      if (saveFigures) {
        ggsave(paste0("Figures/pdf/Final Calibration All Sessions - ", method,".pdf"),
               fullFigure, width = sizeFig[1], height = sizeFig[2])
        ggsave(paste0("Figures/Final Calibration All Sessions - ", method,".jpg"),
               fullFigure, width = sizeFig[1], height = sizeFig[2])
      }
    }
  }
}

