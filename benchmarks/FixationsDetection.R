########################################################################### \
## Project: BRIL
## Script purpose: Process eyetracking calibrations based on fixations detections
## Author: Adrien Brilhault
########################################################################### \


## ========================================================================\
## -------- Initialization  --------
## ========================================================================\


if (F) {
  install.packages("remotes")
  remotes::install_github("adrienbrilhault/BRIL", subdir = "pkg")
  install.packages("tidyverse")
  install.packages("plotly")
  install.packages("grid")
  install.packages("gridExtra")
  install.packages("viridis")
  install.packages("ggnewscale")
  install.packages("xlsx")
}

library(BRIL)
library(tidyverse)
library(plotly)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(ggnewscale)
library(xlsx)

theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils_fixations.R") # Saccade/Fixation detection and other helpers

## Load data
load("../Data/eyeTrackingCalibrations.Rdata")

## Select sessions to process
# sessions <- c("ded005a01")
# sessions <- c("ded00800")
# sessions <- unique(metadata$session)
sessions <- c("juj011a00", "ded00800", "ded005a01")

## Signal recorded frequency
signalFreq <- 1000
## Converting factor (pixels -> degrees) based on monitor size, distance and resolution
pixeltodegrees <- 1 / 25.470
## Converting factor (raw coordinates -> pixels), computed automatically for each session (based on groundTruth annotations)
convertionPixels <- list()
## Discard fixations in the first 200ms (set to 0 to keep all coordinates)
skipFirstMs <- round(200 * 1000 / signalFreq)
## Discard fixations starting after 1000ms (set to 0 to keep all coordinates)
skipAfterMs <- 0

## Other params
showFigures <- TRUE
saveFigures <- FALSE

plotTrialSamples <- FALSE
plotTrialFixations <- TRUE
plotUnfilteredData <- TRUE


## ========================================================================\
## --------  FIXATIONS  DETECTION--------
## ========================================================================\

## ========================================================================\
## -------- o Parse Sessions  --------
## ========================================================================\

resultsFixations <- data.frame()

for (session in sessions) {
  cat("\n\nSession", session, "\n\n")

  ## -------- compute conversion raw coordinates -> pixels in the session

  sessionGroundTruth <- groundTruth[groundTruth$session == session, ]
  distances <- as.matrix(stats::dist(sessionGroundTruth[order(sessionGroundTruth$target), 3:4]))
  correction <- sqrt(100^2 + 100^2) / mean(distances[2:nrow(distances), 1])
  convertionPixels[session] <- correction


  ## -------- retrieve samples from the calibration session

  dfCalibration <- data.frame()
  duration <- unique(metadata$duration[metadata$session == session])
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
    dfCalibration <- rbind(dfCalibration, df)
  }
  dfCalibration$Target <- factor(dfCalibration$Target, targets)

  ## -------- Set the plotting area for the session (based on samples and on groundtruth)

  xLim <- c(
    min(sessionGroundTruth$X) - (max(sessionGroundTruth$X) - min(sessionGroundTruth$X)) * 0.25,
    max(sessionGroundTruth$X) + (max(sessionGroundTruth$X) - min(sessionGroundTruth$X)) * 0.25
  )
  yLim <- c(
    min(sessionGroundTruth$Y) - (max(sessionGroundTruth$Y) - min(sessionGroundTruth$Y)) * 0.25,
    max(sessionGroundTruth$Y) + (max(sessionGroundTruth$Y) - min(sessionGroundTruth$Y)) * 0.25
  )
  xLimGlobal <- as.numeric(quantile(dfCalibration$X, c(0.001, 0.999)))
  yLimGlobal <- as.numeric(quantile(dfCalibration$Y, c(0.001, 0.999)))

  ## -------- Compute velocities for the whole session to fix the saccade velocity threshold

  velSession <- velocities(dfCalibration$X, dfCalibration$Y,
    errors = c(-32768, 32767), freq = 1000, lowpassfilter = TRUE, filterFreq = 30)

  minSaccadeVelocitySession <- sqrt(mean(velSession^2))

  velSessionRaw <- velocities(dfCalibration$X, dfCalibration$Y,
    errors = c(), freq = 1000, lowpassfilter = FALSE)

  minSaccadeVelocitySessionRaw <- sqrt(mean(velSessionRaw^2))


  ## ========================================================================\
  ## -------- o Parse Targets   --------
  ## ========================================================================\

  dfplot <- sample_n(dfCalibration, 5000)

  for (target in targets) {

    selectedIndices <- which(metadata$session == session & metadata$target == target)

    ## ========================================================================\
    ## -------- o Parse Trials    --------
    ## ========================================================================\

    for (i in selectedIndices) {
      print(paste("Session", session, "- Target", target, "- Trial", metadata$trial[i], ":", duration, "samples"))

      ## ========================================================================\
      ## -------- Plot positions per axis and index  --------
      ## ========================================================================\

      ## -------- Raw data

      if (showFigures | saveFigures) {
        g1 <- ggplot(data.frame(Index = 1:duration, X = x[i, 1:duration], Y = y[i, 1:duration]) %>% gather(Axis, Position, X:Y),
          aes(Index, Position, color = Axis)) +
          geom_point(shape = 1) +
          geom_line() +
          scale_color_manual(values = viridisLite::plasma(3)) +
          labs(x = "") +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
          theme(legend.title = element_blank(), legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
                legend.background = element_rect(fill = "white", color = "black")) +
          scale_y_continuous(position = "right", sec.axis = dup_axis()) +
          theme(axis.text.y.right = element_blank(),
                axis.ticks.y.right = element_blank(),
                axis.title.y.left = element_blank())
        g1
        if (plotUnfilteredData) {
          g1 <- g1 + ggtitle("Raw") +
            theme(plot.title = element_text(size = 11)) +
            theme(plot.title = element_text(hjust = 0.5))
        }
      }

      ## -------- Filtered data

      lowpassfreq <- 30
      xFilt <- filter_eyedata(x[i, 1:duration], freq = 1000, lowpassfilter = T, filterFreq = lowpassfreq)
      yFilt <- filter_eyedata(y[i, 1:duration], freq = 1000, lowpassfilter = T, filterFreq = lowpassfreq)

      if (showFigures | saveFigures) {
        g1b <- ggplot(data.frame(Index = 1:duration, X = xFilt, Y = yFilt) %>% gather(Axis, Position, X:Y),
          aes(Index, Position, color = Axis)) +
          geom_point(shape = 1) +
          geom_line() +
          scale_color_manual(values = viridisLite::plasma(3)) +
          labs(x = "") +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
          theme(legend.title = element_blank()) +
          theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
            legend.background = element_rect(fill = "white", color = "black")) +
          scale_y_continuous(position = "right", sec.axis = dup_axis()) +
          theme(axis.text.y.right = element_blank(),
            axis.ticks.y.right = element_blank(),
            axis.title.y.left = element_blank())
        g1b

        if (plotUnfilteredData) {
          g1b <- g1b + ggtitle("Filtered") +
            theme(plot.title = element_text(size = 11)) +
            theme(plot.title = element_text(hjust = 0.5))

          ## Use the same vertical scale for both plots
          yLimPos <- c(
            min(layer_scales(g1)$y$range$range, layer_scales(g1b)$y$range$range),
            max(layer_scales(g1)$y$range$range, layer_scales(g1b)$y$range$range)
          )
          g1 <- g1 + ylim(yLimPos[1], yLimPos[2])
          g1b <- g1b + ylim(yLimPos[1], yLimPos[2])
        }
      }

      ## ========================================================================\
      ## -------- Detect Fixations and Saccades   --------
      ## ========================================================================\

      ## -------- Compute velocities with filters XY data

      velRaw <- velocities(x[i, 1:duration], y[i, 1:duration], errors = c(), lowpassfilter = FALSE)
      vel <- velocities(x[i, 1:duration], y[i, 1:duration], lowpassfilter = T, freq = signalFreq, filterFreq = 30)

      ## -------- Detect Saccades and Fixations

      missingValues <- which((x[i, 1:duration] %in% c(-32768, 32767)) |
        (x[i, 1:duration] %in% c(-32768, 32767)))

      events <- find_fixations(
        vel = vel,
        minFixationDuration = 100,
        minSaccadeVelocity = minSaccadeVelocitySession,
        missingValues = missingValues
      )
      if (events$nbfixations == 0) {
        warning("No fixation encountered in Session ", session, " - Target ", target, " - Trial ", metadata$trial[i])
      }
      dfVel <- data.frame(Index = 1:duration, Velocity = vel)
      dfVel$Type <- "Fixation"
      dfVel$Type[events$saccade] <- "Saccade"
      dfVel$Type[missingValues] <- "Skipped"

      ## -------- Remove initial data if required

      if (skipFirstMs > 0) {

        dfVel$Type[1:skipFirstMs] <- "Skipped"

        if (events$nbfixations > 0) {

          ## Remove fixations ending within the first "skipFirstMs"
          events$fixations <- matrix(events$fixations[events$fixations[, 2] >= skipFirstMs, ], ncol = 2)

          ## Shorten the fixations starting within the first "skipFirstMs"
          events$fixations[events$fixations[, 1] < skipFirstMs, 1] <- skipFirstMs

          events$durations <- apply(events$fixations, 1, diff) + 1
          events$nbfixations <- length(events$durations)
        }
      }

      ## -------- Remove late fixations if required

      if (skipAfterMs > 0) {

        events$fixations <- matrix(events$fixations[events$fixations[, 1] < skipAfterMs, ], ncol = 2)
        events$durations <- apply(events$fixations, 1, diff) + 1
        events$nbfixations <- length(events$durations)
        dfVel$Type[skipAfterMs:duration] <- "Skipped"
      }


      ## -------- Mark the fixation of longest duration

      if (events$nbfixations > 0) {
        mainFixation <- events$fixations[which.max(events$durations), ]
        dfVel$Type[mainFixation[1]:mainFixation[2]] <- "MainFixation"
      }
      dfVel$Type <- factor(dfVel$Type, levels = c("MainFixation", "Fixation", "Saccade", "Skipped"))

      ## -------- Estimate reference coordinates from the fixation

      if (events$nbfixations > 0) {

        ## Save Longest Fixation

        averagePosition <- c(
          mean(x[i, which(dfVel$Type == "MainFixation")]),
          mean(y[i, which(dfVel$Type == "MainFixation")])
        )
        error <- as.numeric(dist(rbind(
          as.matrix(sessionGroundTruth[sessionGroundTruth$target == target, c("X", "Y")]),
          averagePosition
        )))
        longestDurations <- sort(events$durations, decreasing = TRUE)
        longestFixation <- events$fixations[which.max(events$durations), ]
        resultsFixations <- rbind(resultsFixations, data.frame(
          Method = "LongestFixation", Session = session, Target = target, Trial = metadata$trial[i],
          Duration = duration, X = averagePosition[1], Y = averagePosition[2],
          ErrorRaw = error, ErrorPx = error * convertionPixels[[session]], ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
          NbFixations = events$nbfixations, Fix1Start = longestFixation[1], Fix1End = longestFixation[2],
          Fix1Dur = longestDurations[1], Fix2Dur = longestDurations[2], Fix3Dur = longestDurations[3]
        ))

        ## Save First Fixation

        averagePositionFirst <- c(
          mean(x[i, seq(events$fixations[1, 1], events$fixations[1, 2])]),
          mean(y[i, seq(events$fixations[1, 1], events$fixations[1, 2])])
        )

        error <- as.numeric(dist(rbind(
          as.matrix(sessionGroundTruth[sessionGroundTruth$target == target, c("X", "Y")]),
          averagePosition
        )))
        firstDurations <- events$durations
        firstFixation <- events$fixations[1, ]
        resultsFixations <- rbind(resultsFixations, data.frame(
          Method = "FirstFixation", Session = session, Target = target, Trial = metadata$trial[i],
          Duration = duration, X = averagePositionFirst[1], Y = averagePositionFirst[2],
          ErrorRaw = error, ErrorPx = error * convertionPixels[[session]], ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
          NbFixations = events$nbfixations, Fix1Start = firstFixation[1], Fix1End = firstFixation[2],
          Fix1Dur = firstDurations[1], Fix2Dur = firstDurations[2], Fix3Dur = firstDurations[3]
        ))
      } 

      ## ========================================================================\
      ## -------- Plot Velocities   --------
      ## ========================================================================\

      if (showFigures | saveFigures) {

        ## -------- Colors Saccades/Fixations

        colorsFixations <- c("#21908CFF", "#440154FF", "#a39300ff", "gray40")

        ## -------- Velocities UnFiltered

        dfVelRaw <- data.frame(Index = 1:duration, Velocity = velRaw)
        minSaccadeVelocityTrial <- sqrt(mean(dfVel$Velocity^2))
        dfVelRaw$Type <- "Fixation"
        dfVelRaw$Type[dfVelRaw$Velocity > minSaccadeVelocitySessionRaw] <- "Saccade"
        dfVelRaw$Type <- factor(dfVelRaw$Type, levels = c("MainFixation", "Fixation", "Saccade", "Skipped"))

        g2 <- ggplot(dfVelRaw, aes(Index, Velocity, color = Type, group = NA)) +
          geom_line(size = 1) +
          scale_color_manual(values = colorsFixations[sort(unique(as.numeric(dfVelRaw$Type)))]) +
          geom_hline(yintercept = minSaccadeVelocitySessionRaw, linetype = "dashed", colour = "black") +
          annotate("text", label = paste0("Thresh. (", round(minSaccadeVelocitySessionRaw), ")"),
            x = nrow(dfVel) - 10, y = minSaccadeVelocitySessionRaw - 2,
            colour = "black", vjust = 1, hjust = "inward", fontface = "italic", size = 3.5) +
          theme(legend.title = element_blank()) +
          theme(legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
            legend.background = element_rect(fill = "white", color = "black")) +
          scale_y_continuous(position = "right", sec.axis = dup_axis()) +
          theme(axis.text.y.right = element_blank(),
            axis.ticks.y.right = element_blank(),
            axis.title.y.left = element_blank())
        g2
        # ggplotly(g2)

        ## -------- Velocities Filtered

        g2b <- ggplot(dfVel, aes(Index, Velocity, color = Type, group = NA)) +
          geom_line(size = 1) +
          scale_color_manual(values = colorsFixations[sort(unique(as.numeric(dfVel$Type)))]) +
          geom_hline(yintercept = minSaccadeVelocitySession, linetype = "dashed", colour = "black") +
          annotate("text", label = paste0("Thresh. (", round(minSaccadeVelocitySession), ")"),
            x = nrow(dfVel) - 10, y = minSaccadeVelocitySession - 2,
            colour = "black", vjust = 1, hjust = "inward", fontface = "italic", size = 3.5) +
          theme(legend.title = element_blank()) +
          theme(legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
            legend.background = element_rect(fill = "white", color = "black")) +
          scale_y_continuous(position = "right", sec.axis = dup_axis()) +
          theme(axis.text.y.right = element_blank(),
            axis.ticks.y.right = element_blank(),
            axis.title.y.left = element_blank()
          )
        g2b
        # ggplotly(g2)
      }

      ## ========================================================================\
      ## -------- Plot samples in cartesian coordinates   --------
      ## ========================================================================\

      ## -------- Simple Global Plot (no fixation/saccades)

      if (plotTrialSamples & (showFigures | saveFigures)) {
        dfTrial <- data.frame(X = x[i, 1:duration], Y = y[i, 1:duration])
        g <- ggplot(dfplot, aes(X, Y, color = Target, group = Target)) +
          geom_point(shape = 1) +
          scale_colour_brewer(palette = "Pastel1") +
          geom_point(data = dfTrial, aes(X, Y), colour = "black", shape = 1, inherit.aes = FALSE) +
          geom_point(data = groundTruth[groundTruth$session == session, ], aes(X, Y), colour = "red", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          geom_point(data = groundTruth[groundTruth$session == session & groundTruth$target == target, ], aes(X, Y), colour = "green1", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobal, ylim = yLimGlobal) +
          theme(legend.position = "right") +
          ggtitle(paste0("Session ", session, " - Target ", target, " - Trial ", sprintf("%02d", metadata$trial[i])))

        if (showFigures) {
          print(g)
        }

        ## Save figure
        if (saveFigures) {
          ggsave(paste0("Figures/pdf/Full - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".pdf"), g)
          ggsave(paste0("Figures/Full - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".jpg"), g)
          ggsave(paste0("Figures/pdf/Zoom - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".pdf"), g + coord_fixed(xlim = xLim, ylim = yLim))
          ggsave(paste0("Figures/Zoom - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".jpg"), g + coord_fixed(xlim = xLim, ylim = yLim))
        }
      }

      ## -------- Cartesian Plot with fixations and saccades

      if (showFigures | saveFigures) {

        dfTrial <- data.frame(X = x[i, 1:duration], Y = y[i, 1:duration], Type = dfVel$Type)
        gf <- ggplot(dfTrial, aes(X, Y, color = Type, group = factor(Type))) +
          geom_point(data = dfplot, aes(X, Y), colour = "gray80", shape = 1, inherit.aes = FALSE) +
          geom_point(shape = 1) +
          scale_color_manual(values = colorsFixations[sort(unique(as.numeric(dfVel$Type)))])

        if (events$nbfixations > 0) {
          gf <- gf + geom_point(data = data.frame(X = averagePosition[1], Y = averagePosition[2]), aes(X, Y), colour = "cyan", size = 1.5, stroke = 1, shape = 3, inherit.aes = FALSE)
        }

        gf <- gf + geom_point(data = groundTruth[groundTruth$session == session, ], aes(X, Y), colour = "red", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          geom_point(data = groundTruth[groundTruth$session == session & groundTruth$target == target, ], aes(X, Y), colour = "black", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
          coord_fixed(ratio = 1, expand = TRUE, xlim = xLimGlobal, ylim = yLimGlobal) +
          guides(color = FALSE) +
          labs(x = "", y = "") +
          ggtitle("Gaze Coordinates") +
          theme(plot.title = element_text(size = 11))
        gf
      }

      ## ========================================================================\
      ## -------- Final Figure   --------
      ## ========================================================================\

      if (plotTrialFixations & (showFigures | saveFigures)) {

        ## Figure combining velocities, XY position per index, and cartesian plot
        mainTittle <- paste0("Session ", session, " - Target ", target, " - Trial ", sprintf("%02d", metadata$trial[i]))

        if (plotUnfilteredData) {
          w <- 14
          h <- 7
          fullFigure <- grid.arrange(gf, rbind(ggplotGrob(g1b), ggplotGrob(g2b), size = "last"),
            rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"),
            ncol = 3, widths = c(2, 1, 1),
            top = textGrob(mainTittle, gp = gpar(fontsize = 12))
          )
        } else {
          w <- 10
          h <- 7
          fullFigure <- grid.arrange(gf, rbind(ggplotGrob(g1b), ggplotGrob(g2b), size = "last"),
            ncol = 2, widths = c(2, 1.5),
            top = textGrob(mainTittle, gp = gpar(fontsize = 12))
          )
        }

        ## Save figure
        if (saveFigures) {
          ggsave(paste0("Figures/pdf/Fixations - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".pdf"),
            fullFigure,
            width = w, height = h
          )
          ggsave(paste0("Figures/Fixations - Session ", session, " - Target", target, " - Trial", sprintf("%02d", metadata$trial[i]), ".jpg"),
            fullFigure,
            width = w, height = h
          )
        }
      }
    } # LOOP TRIAL
  } # LOOP TARGET
} # LOOP SESSION

## Save results
if (F) {
  save(resultsFixations, file = "./Results/Fixations/resultsFixations.Rdata")
}

## ========================================================================\
## -------- RUN BRIL FOR COMPARAISON   --------
## ========================================================================\

if (F) {
  sessions <- c("juj011a00", "ded00800", "ded005a01")
  # sessions <- c("ded005a01")

  ## Load if needed
  if (F) {
    load("../data/eyeTrackingCalibrations.Rdata")
  }

  convertionPixels <- list()
  pixeltodegrees <- 1 / 25.470
  resultsBril <- data.frame()

  for (session in sessions) {

    ## -------- compute conversion raw coordinates -> pixels for the session

    sessionGroundTruth <- groundTruth[groundTruth$session == session, ]
    distances <- as.matrix(stats::dist(sessionGroundTruth[order(sessionGroundTruth$target), 3:4]))
    correction <- sqrt(100^2 + 100^2) / mean(distances[2:nrow(distances), 1])
    convertionPixels[session] <- correction

    ## -------- retrieve samples from the calibration session

    dfCalibration <- data.frame()
    duration <- unique(metadata$duration[metadata$session == session])
    targets <- sort(unique(metadata$target[metadata$session == session]))

    ## Pool all samples of the session, per target, filtering error values
    for (target in targets) {
      cat("Session", session, "- Target", target, "\n")

      selectedIndices <- which(metadata$session == session & metadata$target == target)
      df <- data.frame(
        X = as.vector(t(x[selectedIndices, ])),
        Y = as.vector(t(y[selectedIndices, ]))
      )
      df <- df[!is.na(df$X) | !is.na(df$Y), ]
      df <- df[df$X != -32768 & df$Y != -32768, ]
      df <- df[df$X != 32767 & df$Y != 32767, ]
      df$Target <- target
      dfCalibration <- rbind(dfCalibration, df)

      ## -------- Compute BRIL results

      data <- sample_n(df[, c("X", "Y")], 1000)

      res <- bril(data, maxIterations = 5, method = "Projection", testNormal = "Chisq")

      error <- as.numeric(dist(rbind(
        as.matrix(groundTruth[groundTruth$session == session & groundTruth$target == target, c("X", "Y")]),
        res$mode
      )))

      resultsBril <- rbind(resultsBril, data.frame(
        Method = "BRIL", Session = session, Target = target, X = res$mode[1], Y = res$mode[2],
        ErrorRaw = error, ErrorPx = error * convertionPixels[[session]], ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees
      ))
    }
  } # LOOP SESSION

  resultsBril$Target <- factor(resultsBril$Target)

  ## Save results
  if (F) {
    save(resultsBril, file = "./Results/Fixations/resultsBril.Rdata")
  }
}

## ========================================================================\
## -------- RESULTS ANALYSIS  --------
## ========================================================================\

if (F) {

  ## Load if needed
  if (F) {
    load("../Data/eyeTrackingCalibrations.Rdata")
    load("./Results/Fixations/resultsBril.Rdata")
    load("./Results/Fixations/resultsFixations.Rdata")
  }

  resultsAnalysis <- data.frame()
  convertionPixels <- list()
  pixeltodegrees <- 1 / 25.470

  for (session in unique(resultsFixations$Session)) {

    ## -------- compute conversion raw coordinates -> pixels for the session

    duration <- unique(metadata$duration[metadata$session == session])
    sessionGroundTruth <- groundTruth[groundTruth$session == session, ]
    distances <- as.matrix(stats::dist(sessionGroundTruth[order(sessionGroundTruth$target), 3:4]))
    correction <- sqrt(100^2 + 100^2) / mean(distances[2:nrow(distances), 1])
    convertionPixels[session] <- correction

    for (target in unique(resultsFixations$Target[resultsFixations$Session == session])) {

      ## -------- For each target, average the position of the first fixation in each trial

      selectedIndices <- which(resultsFixations$Session == session & resultsFixations$Target == target & resultsFixations$Method == "FirstFixation")
      averagePosition <- c(
        mean(resultsFixations$X[selectedIndices], na.rm = TRUE),
        mean(resultsFixations$Y[selectedIndices], na.rm = TRUE)
      )
      error <- as.numeric(dist(rbind(
        as.matrix(groundTruth[groundTruth$session == session & groundTruth$target == target, c("X", "Y")]),
        averagePosition
      )))
      resultsAnalysis <- rbind(resultsAnalysis, data.frame(
        Method = "AvgFirstFixations", Session = session, Target = target,
        Duration = resultsFixations$Duration[selectedIndices[1]], X = averagePosition[1], Y = averagePosition[2],
        ErrorRaw = error, ErrorPx = error * convertionPixels[[session]], ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
        NbFixations = mean(resultsFixations$NbFixations[selectedIndices], na.rm = TRUE),
        Fix1Dur = mean(resultsFixations$Fix1Dur[selectedIndices], na.rm = TRUE),
        Fix2Dur = mean(resultsFixations$Fix2Dur[selectedIndices], na.rm = TRUE),
        Fix3Dur = mean(resultsFixations$Fix3Dur[selectedIndices], na.rm = TRUE)
      ))

      ## -------- For each target, average the position of the longest fixation in each trial

      selectedIndices <- which(resultsFixations$Session == session & resultsFixations$Target == target & resultsFixations$Method == "LongestFixation")
      averagePosition <- c(
        mean(resultsFixations$X[selectedIndices], na.rm = TRUE),
        mean(resultsFixations$Y[selectedIndices], na.rm = TRUE)
      )
      error <- as.numeric(dist(rbind(
        as.matrix(groundTruth[groundTruth$session == session & groundTruth$target == target, c("X", "Y")]),
        averagePosition
      )))
      resultsAnalysis <- rbind(resultsAnalysis, data.frame(
        Method = "AvgLongestFixations", Session = session, Target = target,
        Duration = resultsFixations$Duration[selectedIndices[1]], X = averagePosition[1], Y = averagePosition[2],
        ErrorRaw = error, ErrorPx = error * convertionPixels[[session]], ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
        NbFixations = mean(resultsFixations$NbFixations[selectedIndices], na.rm = TRUE),
        Fix1Dur = mean(resultsFixations$Fix1Dur[selectedIndices], na.rm = TRUE),
        Fix2Dur = mean(resultsFixations$Fix2Dur[selectedIndices], na.rm = TRUE),
        Fix3Dur = mean(resultsFixations$Fix3Dur[selectedIndices], na.rm = TRUE)
      ))


      ## -------- For each target, select the trial with the longest fixation

      selectedIndices <- which(resultsFixations$Session == session & resultsFixations$Target == target)
      longestFixations <- selectedIndices[order(resultsFixations$Fix1Dur[selectedIndices], decreasing = TRUE)]
      error <- as.numeric(dist(rbind(
        as.matrix(groundTruth[groundTruth$session == session & groundTruth$target == target, c("X", "Y")]),
        c(resultsFixations$X[longestFixations[1]], resultsFixations$Y[longestFixations[1]])
      )))

      resultsAnalysis <- rbind(resultsAnalysis, data.frame(
        Method = "LongestFixationAcrossTrials", Session = session, Target = target,
        Duration = resultsFixations$Duration[selectedIndices[1]], X = resultsFixations$X[longestFixations[1]], Y = resultsFixations$Y[longestFixations[1]],
        ErrorRaw = error,
        ErrorPx = error * convertionPixels[[session]],
        ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
        NbFixations = resultsFixations$NbFixations[longestFixations[1]],
        Fix1Dur = resultsFixations$Fix1Dur[longestFixations[1]],
        Fix2Dur = resultsFixations$Fix2Dur[longestFixations[2]],
        Fix3Dur = resultsFixations$Fix3Dur[longestFixations[3]]
      ))

      ## -------- Include BRIL results if needed, as comparaison

      if (exists("resultsBril") & T) {

        brilRes <- resultsBril[resultsBril$Session == session & resultsBril$Target == target, ]
        # Recompute Error (optional)
        error <- as.numeric(dist(rbind(
          as.matrix(groundTruth[groundTruth$session == session & groundTruth$target == target, c("X", "Y")]),
          brilRes[c("X", "Y")]
        )))

        resultsAnalysis <- rbind(resultsAnalysis, data.frame(
          Method = "BRIL", Session = session, Target = target,
          Duration = duration, X = brilRes$X, Y = brilRes$Y,
          ErrorRaw = error,
          ErrorPx = error * convertionPixels[[session]],
          ErrorDeg = error * convertionPixels[[session]] * pixeltodegrees,
          NbFixations = NA,
          Fix1Dur = NA,
          Fix2Dur = NA,
          Fix3Dur = NA
        ))
      }
    } # Loop TARGET

    ## -------- Summary per calibration session

    for (method in unique(resultsAnalysis$Method)) {

      selectedIndices <- which(resultsAnalysis$Session == session & resultsAnalysis$Method == method)

      resultsAnalysis <- rbind(resultsAnalysis, data.frame(
        Method = method, Session = session, Target = "All",
        Duration = resultsAnalysis$Duration[selectedIndices[1]],
        X = NA, Y = NA,
        ErrorRaw = mean(resultsAnalysis$ErrorRaw[selectedIndices], na.rm = TRUE),
        ErrorPx = mean(resultsAnalysis$ErrorPx[selectedIndices], na.rm = TRUE),
        ErrorDeg = mean(resultsAnalysis$ErrorDeg[selectedIndices], na.rm = TRUE),
        NbFixations = mean(resultsAnalysis$NbFixations[selectedIndices], na.rm = TRUE),
        Fix1Dur = mean(resultsAnalysis$Fix1Dur[selectedIndices], na.rm = TRUE),
        Fix2Dur = mean(resultsAnalysis$Fix2Dur[selectedIndices], na.rm = TRUE),
        Fix3Dur = mean(resultsAnalysis$Fix3Dur[selectedIndices], na.rm = TRUE)
      ))
    }
  } # LOOP SESSION

  ## -------- Summary across all sessions

  for (method in unique(resultsAnalysis$Method)) {

    selectedIndices <- which(resultsAnalysis$Method == method & resultsAnalysis$Target == "All")

    resultsAnalysis <- rbind(resultsAnalysis, data.frame(
      Method = method, Session = "All", Target = "All",
      Duration = NA, X = NA, Y = NA,
      ErrorRaw = mean(resultsAnalysis$ErrorRaw[selectedIndices], na.rm = TRUE),
      ErrorPx = mean(resultsAnalysis$ErrorPx[selectedIndices], na.rm = TRUE),
      ErrorDeg = mean(resultsAnalysis$ErrorDeg[selectedIndices], na.rm = TRUE),
      NbFixations = mean(resultsAnalysis$NbFixations[selectedIndices], na.rm = TRUE),
      Fix1Dur = mean(resultsAnalysis$Fix1Dur[selectedIndices], na.rm = TRUE),
      Fix2Dur = mean(resultsAnalysis$Fix2Dur[selectedIndices], na.rm = TRUE),
      Fix3Dur = mean(resultsAnalysis$Fix3Dur[selectedIndices], na.rm = TRUE)
    ))
  }

  resultsAnalysis[is.na(resultsAnalysis)] <- NA # Uniforming NA and NaN


  ## -------- Global Rankings Summary

  globalSummary <- bind_rows(
    resultsAnalysis[resultsAnalysis$Target != "All", ] %>%
      group_by(Method, Session) %>%
      summarise(
        errorPx = mean(ErrorPx, na.rm = TRUE), errorPx_sd = sd(ErrorPx, na.rm = TRUE),
        errorDeg = mean(ErrorDeg, na.rm = TRUE), errorDeg_sd = sd(ErrorDeg, na.rm = TRUE),
        sum_squarred_errorPx = sum(ErrorPx)
      ) %>% arrange(desc(Session), errorPx),
    resultsAnalysis[resultsAnalysis$Target != "All", ] %>%
      group_by(Method) %>%
      summarise(
        Session = "All", errorPx = mean(ErrorPx, na.rm = TRUE), errorPx_sd = sd(ErrorPx, na.rm = TRUE),
        errorDeg = mean(ErrorDeg, na.rm = TRUE), errorDeg_sd = sd(ErrorDeg, na.rm = TRUE),
        sum_squarred_errorPx = sum(ErrorPx)
      ) %>% arrange(errorPx)
  )
  print(globalSummary)
  view(globalSummary)


  ## Save results
  if (F) {
    save(resultsAnalysis, file = "./Results/Fixations/resultsAnalysis.Rdata")
  }

  ## Write results
  if (F) {
    write.xlsx(as.data.frame(globalSummary), "./Results/Fixations/results_summary.xlsx",
      sheetName = "Summary",
      col.names = TRUE, row.names = FALSE, append = FALSE
    )
    write.xlsx(as.data.frame(resultsAnalysis), "./Results/Fixations/results_detailled.xlsx",
      sheetName = "Details",
      col.names = TRUE, row.names = FALSE, append = FALSE
    )
    write.xlsx(as.data.frame(resultsFixations), "./Results/Fixations/results_fixations.xlsx",
      sheetName = "LongestFixations",
      col.names = TRUE, row.names = FALSE, append = FALSE
    )
    write.xlsx(as.data.frame(resultsBril), "./Results/Fixations/results_bril.xlsx",
      sheetName = "Bril",
      col.names = TRUE, row.names = FALSE, append = FALSE
    )
  }
}


## ========================================================================\
## -------- PLOT FINAL CALIBRATIONS   --------
## ========================================================================\

if (F) {

  sessions <- c("juj011a00", "ded00800", "ded005a01")
  # sessions <- c("ded005a01")

  ## Load if needed
  if (F) {
    load("../Data/eyeTrackingCalibrations.Rdata")
    load("./Results/Fixations/resultsFixations.Rdata")
    load("./Results/Fixations/resultsBril.Rdata")
    load("./resultsAnalysis.Rdata")
  }

  for (session in unique(resultsFixations$Session)) {
    sessionGroundTruth <- groundTruth[groundTruth$session == session, ]

    ## -------- retrieve samples from the calibration session

    dfCalibration <- data.frame()
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
      dfCalibration <- rbind(dfCalibration, df)
    }
    dfCalibration$Target <- factor(dfCalibration$Target, targets)
    set.seed(123)
    dfCalibrationPlot <- sample_n(dfCalibration, min(nrow(dfCalibration), 6000))

    ## -------- Set the plotting area for the session (based on samples and on groundtruth)

    xLim <- c(
      min(sessionGroundTruth$X) - (max(sessionGroundTruth$X) - min(sessionGroundTruth$X)) * 0.25,
      max(sessionGroundTruth$X) + (max(sessionGroundTruth$X) - min(sessionGroundTruth$X)) * 0.25
    )
    yLim <- c(
      min(sessionGroundTruth$Y) - (max(sessionGroundTruth$Y) - min(sessionGroundTruth$Y)) * 0.25,
      max(sessionGroundTruth$Y) + (max(sessionGroundTruth$Y) - min(sessionGroundTruth$Y)) * 0.25
    )
    xLimGlobal <- as.numeric(quantile(dfCalibration$X, c(0.001, 0.999)))
    yLimGlobal <- as.numeric(quantile(dfCalibration$Y, c(0.001, 0.999)))

    ## -------- Plot each method

    # listmethods <- c("AvgLongestFixations", "LongestFixationAcrossTrials", "BRIL")
    listmethods <- c("AvgFirstFixations", "AvgLongestFixations", "LongestFixationAcrossTrials", "BRIL")
    # listmethods <- c("AvgLongestFixations","LongestFixationAcrossTrials")
    # listmethods <- c("BRIL")

    for (method in listmethods) {
      if (method == "BRIL") {
        resultModes <- resultsBril[resultsBril$Session == session, c("X", "Y", "Target")]
      } else {
        resultModes <- resultsAnalysis[resultsAnalysis$Session == session & resultsAnalysis$Target != "All" & resultsAnalysis$Method == method, c("X", "Y", "Target")]
      }

      g <- ggplot(dfCalibrationPlot, aes(X, Y, color = Target, group = Target)) +
        geom_point(shape = 1, alpha = 0.9) +
        scale_colour_hue(c = 40, l = 80) +
        geom_point(data = groundTruth[groundTruth$session == session, ], aes(X, Y), colour = "gray30", size = 2, stroke = 1.5, shape = 4, inherit.aes = FALSE) +
        new_scale_color() +
        scale_colour_hue(c = 60, l = 40) +
        geom_point(data = resultModes, aes(X, Y, colour = Target), size = 3, stroke = 2, shape = 1, inherit.aes = FALSE) +
        coord_fixed(ratio = 1, expand = TRUE, xlim = xLim, ylim = yLim) +
        theme(legend.position = "right") +
        theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
        ggtitle(paste0("Session ", session))

      print(g)

      ## Save figure
      if (saveFigures) {
        sizeFig = c(6, 4)

        ggsave(paste0("Figures/pdf/Calibration ", session, " - ", method, ".pdf"), g, width = sizeFig[1], height =  sizeFig[2])
        ggsave(paste0("Figures/Calibration ", session, " - ", method, ".jpg"), g, width = sizeFig[1], height =  sizeFig[2])
        ggsave(paste0("Figures/pdf/Full - Calibration ", session, " - ", method, ".pdf"), g + coord_fixed(xlim = xLimGlobal, ylim = yLimGlobal), width = sizeFig[1], height =  sizeFig[2])
        ggsave(paste0("Figures/Full - Calibration ", session, " - ", method, ".jpg"), g + coord_fixed(xlim = xLimGlobal, ylim = yLimGlobal), width = sizeFig[1], height =  sizeFig[2])
      }
    }
  } # LOOP SESSION
}

## ========================================================================\
## -------- MOVIE  --------
## ========================================================================\

## Code to create a movie clip from the plot images

if (F) {

  ## Rename files
  old_files <- list.files("Figures/tmp", pattern = "*.jpg", full.names = TRUE)
  new_files <- paste0("Figures/tmp/img", sprintf("%04d", 1:length(old_files)), ".jpg")
  file.copy(from = old_files, to = new_files)

  ## Create clip

  system("ffmpeg -r 1/2 -i Figures/tmp/img%04d.jpg -vcodec libx264 -vf fps=25 -crf 15 -pix_fmt yuv420p movie.mp4")
  # system("ffmpeg -r 25 -qscale 2 -i Figures/Zoom/tmp/img%04d.jpg output.mp4")
  # system("ffmpeg -r 1/2 -pattern_type glob -i Figures/Zoom/tmp/*.jpg -c:v libx264 -vf fps=25 -pix_fmt yuv420p movie.mp4")
}
