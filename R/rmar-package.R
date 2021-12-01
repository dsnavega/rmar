# This file is part of rmar
#
# Copyright (C) 2021, David Senhora Navega
#
# rmar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rmar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with rmar. If not, see <http://www.gnu.org/licenses/>.
#
# David Senhora Navega
# Laboratory of Forensic Anthropology
# Department of Life Sciences
# University of Coimbra
# Calçada Martim de Freitas, 3000-456, Coimbra
# Portugal

#' Median Absolute Error (MedAE)
#'
#' Evaluates MedAE on predictions of a regression model
#'
#' @author David Senhora Navega
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return MedAE
#'
medae <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))
  medae <- stats::median(abs(data$x - data$y))

  return(medae)

}

#' Mean Absolute Error (MAE)
#'
#' Evaluates MAE on predictions of a regression model
#'
#' @author David Senhora Navega
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return MAE
#'
mae <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))
  mae <- mean(abs(data$x - data$y))

  return(mae)

}

#' Root Mean Squared Error (RMSE)
#'
#' Evaluates RMSE on predictions of a regression model
#'
#' @author David Senhora Navega
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return RMSE
#'
rmse <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))
  rmse <- sqrt(mean((data$x - data$y) ^ 2))

  return(rmse)

}

#' R Squared
#'
#' Evaluates R Squared (Explained Variance)
#'
#' @author David Senhora Navega
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return R Squared value
#'
#' @references
#' Gelman A, Goodrich B, Gabry J, Vehtari A. R-squared for Bayesian regression
#' models. Am Stat 2019;73(3):307–9.
#'
rsquared <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))

  # Predictions Variance
  pss <- sum((data$y - mean(data$y)) ^ 2)
  # Residual Variance
  rss <- sum((data$x - data$y) ^ 2)
  # R Squared
  rsquared <- pss / (pss + rss)

  return(rsquared)

}

#' Regression Prediction Bias
#'
#' Evaluate Regression Prediction Bias by the slope of the regression model of
#' residuals on known values.
#'
#' @author David Senhora Navega
#' @noRd
#'
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return Regression prediction bias value
#'
#' @details
#' A positive bias means that lower known values are systematically
#' overestimated and upper known values are underestimated. A value near 0
#' indicates no bias
#'
prediction_bias <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = known - predicted))

  Sxy <- sum((data$x - mean(data$x)) * (data$y - mean(data$y)))
  Sxx <- sum((data$x - mean(data$x)) ^ 2)
  slope <- Sxy / Sxx
  return(slope)

}

#' Prediction Coverage
#'
#' Evaluates the coverage of prediction intervals.
#'
#' @author David Senhora Navega
#' @noRd
#' @md
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric matrix. See Details
#'
#'
#' @return Coverage
#'
#' @details
#'
#' predicted must be a matrix organized as one the follows options:
#' * 1) 3 columns where the first column stores the point estimate of the model
#' and the second and third columns store the lower and upper bounds of the
#' predictive interval if the algorithm is able to produce one.
#' * 2) 2 columns storing the lower and upper bounds of the predictive interval.
#'
prediction_coverage <- function(known, predicted) {

  n <- length(known) - sum(rowMeans(is.na(predicted)))
  predicted <- t(apply(predicted, MARGIN = 1, FUN = range))

  hits <- sapply(seq_len(n), function(i) {
    known[i] >= predicted[i, 1] & known[i] <= predicted[i, 2]
  })
  coverage <- sum(hits, na.rm = T) / n
  return(coverage)

}

#' Prediction Interval Width
#'
#' Evaluates the width of prediction intervals as a measure of efficiency
#'
#' @author David Senhora Navega
#' @noRd
#' @md
#'
#' @import stats
#' @param predicted a numeric matrix. See Details
#' @param tau quantiles to evaluate the distribution the  predictive intervals.
#' Default is 0.5, 0.025, and 0.975 which return the median, and 95% CI of PIW
#' @return Prediction Width
#'
#' @details
#' predicted must be a matrix organized as one the follows options:
#' * 1) 3 columns where the first column stores the point estimate of the model
#' and the second and third columns store the lower and upper bounds of the
#' predictive interval if the algorithm is able to produce one.
#'
#' * 2) 2 columns storing the lower and upper bounds of the predictive interval.
#'
prediction_interval_width <- function(predicted, tau = c(0.5, 0.025, 0.975)) {

  piw.range <- function(x) {
    if(any(is.infinite(x) | is.na(x))) {
      return(NA)
    } else {
      return(diff(range(x)))
    }
  }

  data <- na.omit(data.frame(y = predicted))
  piw <- quantile(x = apply(data, MARGIN = 1, FUN = piw.range), probs = tau)
  names(piw) <- paste0("PIW(", as.character(tau), ")")

  return(piw)

}
