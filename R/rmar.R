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

#' rmar: Regression Model Assessment in R
#'
#' @author David Senhora Navega
#'
#' @export
#'
#' @param known a numeric vector
#' @param predicted a numeric vector or matrix. See Details.
#' @param interval minimum and maximum value of prediction domain. Should be
#' supplied as vector of with two elements. If not supplied is assess from the
#' known argument.
#' @param digits an integer defining the precision of the round function.
#'
#' @details
#' If predicted is a matrix, it should be organized in way that the first column
#' represents the point estimate and the second and third columns represent the
#' lower and upper bound of the predictive interval.
rmar <- function(known, predicted, interval = NULL, digits = 3) {

  if (!(is.vector(known) & is.numeric(known)))
    stop("\n(-) y must be a numeric vector.")

  if ((!(is.vector(predicted) & is.numeric(predicted))) & !is.matrix(predicted))
    stop("\n(-) predicted must be a numeric vector or matrix.")

  if(NROW(known) !=  NROW(predicted))
    stop("\n(-) Number of known and predicted observations do not match.")

  # Define Prediction Domain
  if (is.null(interval)) {
    n <- sum(!is.na(known))
    bw <- sd(x = known, na.rm = T)  * n ^ -0.25
    interval <- range(x = known, na.rm = T) + c(-bw, bw) * 0.5
  } else {
    interval <- range(interval, na.rm = T)
  }

  lower <- interval[1]
  upper <- interval[2]

  # Define Valid Known and Predicted Data
  predicted <- cbind(predicted)
  complete <- apply(predicted, 1, function(x) !all(is.na(x)))
  known <- known[complete]
  predicted <- predicted[complete,, drop = F]

  # Compute Measures
  medae <- medae(known = known, predicted = predicted[, 1])
  mae <- mae(known = known, predicted = predicted[, 1])
  rmse <- rmse(known = known, predicted = predicted[, 1])
  rsquared <- rsquared(known = known, predicted = predicted[, 1])
  bias <- prediction_bias(known = known, predicted = predicted[, 1])

  vector.measures <- c(medae,  mae, rmse, rsquared, bias)
  names(vector.measures) <- c("MedAE", "MAE", "RMSE", "R Squared", "Bias")

  if (ncol(predicted) == 3) {
    coverage <- prediction_coverage(known = known, predicted = predicted)
    piw <- prediction_interval_width(predicted = predicted)
    matrix.measures <- c(coverage, piw)
    names(matrix.measures) <- c("P(alpha)",names(piw))
  }

  if (exists("matrix.measures"))
    return(round(x = c(vector.measures, matrix.measures), digits = digits))

  return(round(x = vector.measures, digits = digits))

}
