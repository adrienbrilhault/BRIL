########################################################################### \
## Project: BRIL
## Script purpose: Helper functions to process eye-tracking data (filtering,
## computing velocities, detecting saccades and fixations)
## Author: Adrien Brilhault
########################################################################### \

library(signal)

#' Filter raw eye-tracking data
#'
#' Process raw eye-tracking data by replacing invalid values (provided in
#' `errors`) with the last valid coordinate received, and applying an optional
#' low-pass filter
#'
#' @param data raw eye-tracking data (numerical away of one dimensions)
#' @param errors List of invalid values to be discaded
#' @param lowpassfilter Logical value (to apply or not the low-pass filter)
#' @param freq Sampling frequency of the signal provided in `data`
#' @param filterFreq Cutoff frequency of the lowpass filter (in Hz)
#'
#' @return Filtered data (same dimensions as `data`)
#'
filter_eyedata <- function(data, errors = c(-32768, 32767), lowpassfilter = TRUE, freq = 1000, filterFreq = 30) {

  ## Fill incorrect values with the last valid position
  for (i in 2:length(data)) {
    if (data[i] %in% errors) {
      data[i] <- data[i - 1]
    }
  }

  ## Apply a low pass filter
  if (lowpassfilter) {

    ## Low pass filter specifications
    # flt <- fir2(60, c(0, filterFreq / (freq / 2), filterFreq / (freq / 2), 1), c(1, 1, 0, 0))
    # flt <- fir1(60, filterFreq / (freq / 2), "low")
    # flt <- ellip(5, 0.1, 60, filterFreq / (freq / 2), "low")
    flt <- butter(3, filterFreq / (freq / 2), "low")

    ## Add buffer on both ends
    data <- c(data[seq(50, 1, -1)], data, data[seq(length(data), length(data) - 49, -1)])
    ## Apply the filter
    data <- filtfilt(flt, data)
    ## Remove the buffer
    data <- data[51:(length(data) - 50)]
  }

  return(data)
}


#' Computes a signal's velocities
#'
#' @param x Input signal (horizontal coordinates)
#' @param y Input signal (vetical coordinates, of same size as `x`)
#' @param errors List of invalid values to be discaded
#' @param lowpassfilter Logical value (to apply or not the low-pass filter)
#' @param freq Sampling frequency of the signal provided in `x` and `y`
#' @param filterFreq Cutoff frequency of the lowpass filter (in Hz)
#'
#' @return an array of same size than `x` and `y`, containing the velocities
#'   computed
#'
velocities <- function(x, y, errors = c(-32768, 32767), lowpassfilter = TRUE, freq = 1000, filterFreq = 30) {
  xFilt <- filter_eyedata(x, errors, lowpassfilter, freq, filterFreq)
  yFilt <- filter_eyedata(x, errors, lowpassfilter, freq, filterFreq)

  return(c(0, sqrt(diff(xFilt)^2 + diff(yFilt)^2)))
}

#' Detect Saccades and Fixations
#'
#' @param vel Numerical array containing the signal velocity
#' @param minFixationDuration Minimum durations of fixations
#' @param minSaccadeVelocity Minimum saccade velocity
#' @param missingValues Indices of the samples in `vel`
#'
#' @return
#'
find_fixations <- function(vel, minFixationDuration = NULL, minSaccadeVelocity = NULL, missingValues = NULL) {

  ## If not provided, compute the RMS value of the signal to be used as detection threshold
  if (missing(minSaccadeVelocity) || is.null(minSaccadeVelocity) || is.na(minSaccadeVelocity)) {
    minSaccadeVelocity <- sqrt(mean(vel^2))
  }

  ## Classify samples above the threshold as saccades
  saccade <- vel > minSaccadeVelocity


  ## ========================================================================\
  ## -------- Propagate saccades on both ends  --------
  ## ========================================================================\

  i <- 1
  insideSaccade <- F
  while (i <= length(saccade)) {

    ## Propagate the saccade at its end while velocities are strictly decreasing
    if (insideSaccade) {
      if (saccade[i] == T) {
        i <- i + 1
        next
      } else {
        j <- i
        while (j <= length(saccade) && saccade[j] == F && vel[j] < vel[j - 1] && !(j %in% missingValues)) {
          saccade[j] <- T
          j <- j + 1
        }
        i <- j
        insideSaccade <- F
        next
      }
    }

    ## Propagate the saccade before its start while velocities are strictly decreasing
    if (!insideSaccade) {
      if (saccade[i] == F) {
        i <- i + 1
        next
      } else {
        j <- i - 1
        while (j >= 1 && saccade[j] == F && vel[j] < vel[j + 1] && !(j %in% missingValues)) {
          saccade[j] <- T
          j <- j - 1
        }
        insideSaccade <- T
        i <- i + 1
        next
      }
    }
  }

  ## Discard first and last incongruent samples
  saccade[1] <- saccade[2]
  saccade[length(saccade)] <- saccade[length(saccade)]

  ## Missing values are considered as saccades
  if (!missing(missingValues) && !is.null(missingValues)) {
    saccade[missingValues] <- T
  }

  nbfixations <- 0
  fixations <- matrix(NA, 1, 2)

  ## No saccade in the whole signal
  if (!any(saccade)) {
    fixations[1, ] <- c(1, length(vel))
    return(list(saccade = saccade, nbfixations = 1, fixations = fixations, durations = diff(fixations[1, ]) + 1))
  }

  ## No fixation in the whole signal
  if (all(saccade)) {
    return(list(saccade = saccade, nbfixations = 0, fixations = fixations, durations = 0))
  }

  ## ========================================================================\
  ## -------- Find fixations and discard short ones  --------
  ## ========================================================================\

  fixationStart <- min(which(saccade == FALSE))
  i <- fixationStart
  insideFixation <- T
  while (i <= length(saccade)) {

    ## Mark new fixation start
    if (!insideFixation) {
      if (saccade[i] == T) {
        i <- i + 1
        next
      } else {
        fixationStart <- i
        insideFixation <- T
        i <- i + 1
        next
      }
    }

    ## Mark fixation end, and discard if it was shorter than specified
    if (insideFixation) {
      if (saccade[i] == F) {
        i <- i + 1
        next
      } else {
        fixationEnd <- i - 1

        ## Discard the fixation is it was shorter than the min fixation duration
        if (is.numeric(minFixationDuration) && (fixationEnd - fixationStart < minFixationDuration)) {
          saccade[fixationStart:fixationEnd] <- T
        } else {

          ## Add the fixation to the list
          nbfixations <- nbfixations + 1
          if (nbfixations == 1) {
            fixations[nbfixations, ] <- c(fixationStart, fixationEnd)
          } else {
            fixations <- rbind(fixations, c(fixationStart, fixationEnd))
          }
        }

        insideFixation <- F
        i <- i + 1
        next
      }
    }
  }

  ## -------- If we where still inside a fixation at the end of the signal

  if (insideFixation) {
    fixationEnd <- length(saccade)

    ## Discard the last fixation is it was shorter than the min fixation duration
    if (is.numeric(minFixationDuration) && (fixationEnd - fixationStart < minFixationDuration)) {
      saccade[fixationStart:fixationEnd] <- T
    } else {

      ## Add the fixation to the list
      nbfixations <- nbfixations + 1
      if (nbfixations == 1) {
        fixations[nbfixations, ] <- c(fixationStart, fixationEnd)
      } else {
        fixations <- rbind(fixations, c(fixationStart, fixationEnd))
      }
    }
  }

  ## -------- Return results

  return(list(
    saccade = saccade, nbfixations = nbfixations,
    fixations = fixations, durations = apply(fixations, 1, diff) + 1
  ))
}
