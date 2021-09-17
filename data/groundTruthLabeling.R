########################################################################### \
## Project: BRIL
## Script purpose: Annotate the reference coordinates in eyetracking calibrations.
## A video of the procedure is available at: https://www.youtube.com/watch?v=10ZapuMvK1s7
## Another video, showing the raw data trial by trial, is available at
## https://www.youtube.com/watch?v=ZQzThht0VMw
##
## Author: Adrien Brilhault
########################################################################### \


## ========================================================================\
## -------- Initialization  --------
## ========================================================================\


if (F) {
  install.packages("remotes")
  remotes::install_github("adrienbrilhault/BRIL", subdir = "pkg")
  install.packages("ggplot2")
}

library(BRIL)
library(ggplot2)

## Load eyetracking data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("eyeTrackingCalibrations.Rdata")

## Select one or several sessions to process
sessions <- unique(metadata$session)
# sessions <- c("ded005a01")

## Uncomment next line to create a new window for the plots
# X11()

newGroundTruth <- data.frame()

# Analyze each session
for (session in sessions) {

  ## ========================================================================\
  ## -------- Read data  --------
  ## ========================================================================\

  dfCalibration <- data.frame()
  targets <- sort(unique(metadata$target[metadata$session == session]))

  # Iterate through each target within the session
  for (target in targets) {

    # Filter the corresponding trials
    selectedIndices <- metadata$session == session & metadata$target == target
    df <- data.frame(
      X = as.vector(t(x[selectedIndices, ])),
      Y = as.vector(t(y[selectedIndices, ]))
    )

    # Remove trailing NA fill
    df <- df[!is.na(df$X) | !is.na(df$Y), ]

    # Remove overflowed values and error codes from eye-tracker
    df <- df[df$X != -32768 & df$Y != -32768, ]
    df <- df[df$X != 32767 & df$Y != 32767, ]

    # Update structures
    df$Target <- target
    dfCalibration <- rbind(dfCalibration, df)
  }

  ## ========================================================================\
  ## -------- Interactive annotation  --------
  ## ========================================================================\

  dfCalibration$Target <- factor(dfCalibration$Target, targets)
  plotColors <- grDevices::rainbow(length(unique(dfCalibration$Target)))


  ############## Select the overall bounding box

  dfplot <- sample_n(dfCalibration, 5000)

  plot(dfplot$X, dfplot$Y, asp = 1, col = plotColors[dfplot$Target], xlab = "", ylab = "")
  mtext(side = 3, line = 2, adj = 0, cex = 1, paste("Calibration", session))
  mtext(side = 3, line = 1, adj = 0, cex = 0.7, "Click the top left and bottom righ coordinates encompassing the correct data")

  clickCoords <- locator(n = 2)

  filteredSamples <- dfCalibration[dfCalibration$X > min(clickCoords$x) &
    dfCalibration$X < max(clickCoords$x) &
    dfCalibration$Y > min(clickCoords$y) &
    dfCalibration$Y < max(clickCoords$y), ]

  ############## Select the bounding box for each target

  dfplot <- sample_n(filteredSamples, min(2000, nrow(filteredSamples)))

  selectedSamples <- list()
  groundTruthCoordinates <- matrix(data = NA, nrow = length(targets), ncol = 2)
  i <- 1

  while (i <= length(targets)) {
    target <- targets[i]
    targetSamples <- filteredSamples[filteredSamples$Target == target, ]

    ############## Select the bounding box for current target

    cat("\n\n## Calibration", session, "- Selecting samples of Target", target, "\n")

    plot(dfplot$X, dfplot$Y, asp = 1, col = plotColors[dfplot$Target], xlab = "", ylab = "")
    points(sample_n(targetSamples[, 1:2], min(2000, nrow(targetSamples))), col = "black")
    points(targetSamples[, 1:2], col = "black")

    mtext(side = 3, line = 2, adj = 0, cex = 1, paste("Calibration", session, "- Target", target))
    mtext(side = 3, line = 1, adj = 0, cex = 0.7,
      paste("Click the top left and bottom righ coordinates encompassing the valid samples for Target", target, "(in black)"))
    clickCoords <- locator(n = 2)

    targetSamples <- targetSamples[targetSamples$X > min(clickCoords$x) &
      targetSamples$X < max(clickCoords$x) &
      targetSamples$Y > min(clickCoords$y) &
      targetSamples$Y < max(clickCoords$y), ]

    ############## Select the min density for current target

    cat("## Calibration", session, "- Selecting minimal density for Target", target, "\n")

    targetSamples <- sample_n(targetSamples, min(3000, nrow(targetSamples)))
    D <- depth_values(targetSamples[, 1:2], method = "Potential", warnings = TRUE)
    plotDensityColors <- rev(grDevices::rainbow(23)[0:20])
    plot(targetSamples[, 1:2], pch = 20, asp = 1, col = plotDensityColors[as.numeric(cut(D, breaks = 20))], xlab = "", ylab = "")
    mtext(side = 3, line = 2, adj = 0, cex = 1, paste("Calibration", session, "- Target", target))
    mtext(side = 3, line = 1, adj = 0, cex = 0.7,
      paste("Select the valid sample of minimal density (all warmer points will be averaged as ground-truth)"))
    clickCoords <- locator(n = 1)
    minDepth <- depth_values(targetSamples[, 1:2], c(clickCoords$x, clickCoords$y), method = "Potential")

    ############## Show result for confirming target

    cat("## Calibration", session, "- Confirming results for Target", target, "\n")

    selectedSamples[[i]] <- targetSamples[D > minDepth, 1:2]
    groundTruthCoordinates[i, ] <- colMeans(selectedSamples[[i]])

    plot(dfplot$X, dfplot$Y, asp = 1, col = plotColors[dfplot$Target], xlab = "", ylab = "")
    points(selectedSamples[[i]], col = "black")
    points(groundTruthCoordinates[i, 1], groundTruthCoordinates[i, 2], col = "white", pch = 3, cex = 1, lwd = 3)
    mtext(side = 3, line = 2, adj = 0, cex = 1, paste("Calibration", session, "- Target", target))
    mtext(side = 3, line = 1, adj = 0, cex = 0.7,
      paste("Samples selected for Target", target, "in black, with their center as a white cross. Confirm in the command-line."))

    cat("\n", nrow(selectedSamples[[i]]), "samples for Target", target, "of Calibration", session,
        "\n Mean coordinates: ", toString(groundTruthCoordinates[i, ]), "\n\n")

    choice <- ""
    while (!tolower(choice) %in% c("y", "n", "a")) {
      choice <- readline(prompt = "To confirm press \"y\", to redo the current target press \"n\", to abord press \"a\": ")
    }

    if (tolower(choice) == "y") {
      i <- i + 1
    }
    if (tolower(choice) == "a") {
      return()
    }

    ############## Global confirmation after all targets have been done

    if (i == length(targets) + 1) {
      cat("\n\n#### Calibration", session, "- Final Results\n\nReference coordinates:\n")
      print(groundTruthCoordinates)

      plot(dfplot$X, dfplot$Y, asp = 1, col = plotColors[dfplot$Target], xlab = "", ylab = "")
      for (j in seq_along(targets)) {
        points(selectedSamples[[j]], col = "black")
        points(groundTruthCoordinates[j, 1], groundTruthCoordinates[j, 2], col = "white", pch = 3, cex = 1, lwd = 3)
      }
      mtext(side = 3, line = 2, adj = 0, cex = 1, paste("Calibration", session, "- Final Result"))
      mtext(
        side = 3, line = 1, adj = 0, cex = 0.7,
        paste("Samples selected for each targeti n black, with their center as a white cross. Confirm in the command-line.")
      )

      choice <- ""
      while (!tolower(choice) %in% c("y", "n", "a")) {
        choice <- readline(prompt = "To confirm the calibration press \"y\", to restart press \"n\", to abord press \"a\": ")
      }

      if (tolower(choice) == "n") {
        i <- 1
      }
      if (tolower(choice) == "a") {
        return()
      }
    }
  }

  ## ========================================================================\
  ## -------- Update structures                                       --------
  ## ========================================================================\

  for (i in seq_along(targets)) {
    groundTruth <- rbind(groundTruth, data.frame(
      session = session, target = targets[i], X = groundTruthCoordinates[i, 1], Y = groundTruthCoordinates[i, 2]
    ))
  }
}

## Save results
if (F) {
  save(groundTruth, file = "groundTruth.Rdata")
  save(groundTruth, metadata, x, y, file = "eyeTrackingCalibrations.Rdata")
}
