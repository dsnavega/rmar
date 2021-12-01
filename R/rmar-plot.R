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

#' Regression Accuracy Plot
#'
#' @author David Senhora Navega
#' @export
#'
#' @import ggplot2 scales
#'
#' @param known a numeric vector
#' @param predicted a numeric vector or matrix. See Details.
#' @param interval minimum and maximum value of prediction domain. Should be
#' supplied as vector of with two elements. If not supplied is assess from the
#' known argument.
#' @param label a character defining an alternative label for the x axis.
#' @param breaks breaks on plot axis. Default is 10.
#' @param base_size ggplot parameter to control text size in plots.
#' Default is 12
#' @param digits number of digits controlling floats precision. Default is 3
#'
#' @details
#' If predicted is a matrix organized it should be organized as follows:
#' 3 columns where the first column stores the point estimate of the model
#' and the second and third columns store the lower and upper bounds of the
#' predictive interval if the algorithm is able to produce one.
#'
#' @return a ggplot2 graphical object
#'
rma_accuracy_plot <- function(
  known, predicted, interval = NULL, label = NULL,
  breaks = 10, base_size = 12, digits = 3
) {

  if (!(is.vector(known) & is.numeric(known)))
    stop("\n(-) y must be a numeric vector.")

  if ((!(is.vector(predicted) & is.numeric(predicted))) & !is.matrix(predicted))
    stop("\n(-) predicted must be a numeric vector or matrix.")

  if(NROW(known) !=  NROW(predicted))
    stop("\n(-) Number of known and predicted observations do not match.")


  if (is.null(interval)) {
    n <- sum(!is.na(known))
    bw <- sd(x = known, na.rm = T)  * n ^ -0.25
    interval <- range(x = known, na.rm = T) + c(-bw, bw) * 0.5
  } else {
    interval <- range(interval, na.rm = T)
  }

  # Plot Data
  lower <- interval[1]
  upper <- interval[2]

  predicted <- cbind(predicted)[,1]

  identity_data <- data.frame(Known = lower:upper, Predicted = lower:upper)
  scatter_data <- data.frame(Known = known, Predicted = predicted)

  medae <- round(
    x = medae(known = known, predicted = predicted),
    digits = digits
  )

  mae <- round(
    x = mae(known = known, predicted = predicted),
    digits = digits
  )

  rmse <- round(
    x = rmse(known = known, predicted = predicted),
    digits = digits
  )

  rsquared <- round(
    x = rsquared(known = known, predicted = predicted),
    digits = digits
  )

  # Build Graphical Object
  title <- "Prediction Accuracy Analysis"
  subtitle <- paste0(
    "MedAE: ", medae,
    "; MAE: ", mae,
    "; RMSE: ", rmse,
    "\nR Squared: ", rsquared
  )


  if (is.null(label)) {
    label <- "Known"
  } else {
    if (!is.character(label))
      stop("\n label must be a character vector")
  }

  grob <- with(identity_data,
    ggplot2::ggplot(data = identity_data) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(x = Known, y = Predicted),
        linetype = "dashed", colour = "red", show.legend = F
      ) +
      ggplot2::geom_point(
        data = scatter_data, mapping = ggplot2::aes(x = Known, y = Predicted)
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = breaks)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = breaks)) +
      ggplot2::labs(x = label) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::ggtitle(label = title, subtitle = subtitle)
  )

  class(grob) <- c("ggsilent", class(grob))
  return(grob)

}

#' Regression Residual and Bias Plot
#'
#' @author David Senhora Navega
#' @export
#'
#' @import ggplot2 scales
#'
#' @param known a numeric vector
#' @param predicted a numeric vector or matrix. See Details.
#' @param interval minimum and maximum value of prediction domain. Should be
#' supplied as vector of with two elements. If not supplied is assess from the
#' known argument.
#' @param label a character defining an alternative label for the x axis.
#' @param scaling a numeric between 0 and 1 defining a scaling of the y axis.
#' Default is 0.5
#' @param breaks breaks on plot axis. Default is 10.
#' @param base_size ggplot parameter to control text size in plots.
#' Default is 12
#' @param digits number of digits controlling floats precision. Default is 3
#'
#' @details
#' 3 columns where the first column stores the point estimate of the model
#' and the second and third columns store the lower and upper bounds of the
#' predictive interval if the algorithm is able to produce one.
#'
#' @return a ggplot2 graphical object
#'
rma_bias_plot <- function(
  known, predicted, interval = NULL, label = NULL,
  scaling = 0.5, breaks = 10, base_size = 12, digits = 3
) {

  if (!(is.vector(known) & is.numeric(known)))
    stop("\n(-) y must be a numeric vector.")

  if ((!(is.vector(predicted) & is.numeric(predicted))) & !is.matrix(predicted))
    stop("\n(-) predicted must be a numeric vector or matrix.")

  if(NROW(known) !=  NROW(predicted))
    stop("\n(-) Number of known and predicted observations do not match.")


  if (is.null(interval)) {
    n <- sum(!is.na(known))
    bw <- sd(x = known, na.rm = T)  * n ^ -0.25
    interval <- range(x = known, na.rm = T) + c(-bw, bw) * 0.5
  } else {
    interval <- range(interval, na.rm = T)
  }

  predicted <- cbind(predicted)[,1]
  bias <- prediction_bias(known = known, predicted = predicted)
  title <- "Prediction Residuals & Bias Analysis"
  subtitle <- paste0(
    "Prediction Bias: ", round(bias, digits = digits),
    "\nResidual = Known - Predicted"
  )

  if (is.null(label)) {
    label <- "Known"
  } else {
    if (!is.character(label))
      stop("\n label must be a character vector")
  }

  residual_data <- data.frame(Known = known, Residual = known - predicted)
  grob <- with(residual_data,
    ggplot2::ggplot(data = residual_data) +
      ggplot2::geom_point(mapping = aes(x = Known, y = Residual)) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dotted",
        colour = "red", size = 1.25, show.legend = F) +
      ggplot2::stat_smooth(
        mapping = aes(x = Known, y = Residual), method = "lm", se = F
      ) +
      ggplot2::scale_x_continuous(breaks = pretty_breaks(n = breaks)) +
      ggplot2::scale_y_continuous(
        breaks = pretty_breaks(n = breaks),
        limits = scaling * c(-diff(interval), diff(interval))
      ) +
      ggplot2::labs(x = label) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::ggtitle(label = title, subtitle = subtitle)
  )

  class(grob) <- c("ggsilent", class(grob))
  return(grob)

}

#' Regression Efficiency Plot
#'
#' In regression, efficiency can be evaluated through the analysis of the width
#' of the predictive interval. The function evaluates and plot the width of the
#' predictive interval and its coverage, the probability of an estimated
#' predictive interval contain the known value.
#'
#' @author David Senhora Navega
#' @export
#' @import ggplot2 scales
#'
#' @param known a numeric vector
#' @param predicted a numeric matrix. See Details.
#' @param interval minimum and maximum value of prediction domain. Should be
#' supplied as vector of with two elements. If not supplied is assess from the
#' known argument.
#' @param breaks breaks on plot axis. Default is 10.
#' @param base_size ggplot parameter to control text size in plots.
#' Default is 12
#' @param digits number of digits controlling floats precision. Default is 3
#'
#' @details
#' predicted matrix must be organized as follows:
#' 3 columns where the first column stores the point estimate of the model
#' and the second and third columns store the lower and upper bounds of the
#' predictive interval if the algorithm is able to produce one.
#'
#' @return a ggplot2 graphical object
#'
rma_efficiency_plot <- function(
  known, predicted, interval = NULL,
  breaks = 10, base_size = 12, digits = 3
) {

  if (!(is.vector(known) & is.numeric(known)))
    stop("\n(-) y must be a numeric vector.")

  if (!is.matrix(predicted))
    stop("\n(-) predicted must be a matrix. See Details of help file.")

  if (NROW(known) !=  NROW(predicted))
    stop("\n(-) Number of known and predicted observations do not match.")


  if (is.null(interval)) {
    n <- sum(!is.na(known))
    bw <- sd(x = known, na.rm = T)  * n ^ -0.25
    interval <- range(x = known, na.rm = T) + c(-bw, bw) * 0.5
  } else {
    interval <- range(interval, na.rm = T)
  }

  # Sort Data According to Predicted
  index <- seq_len(nrow(predicted))
  new.order <- order(predicted[,1])
  known <- known[new.order]
  predicted <- predicted[new.order, ]

  piw_data <- na.omit(object = data.frame(index, known, predicted))
  colnames(piw_data) <- c("Index", "Known", "Center", "Lower", "Upper")

  coverage <- prediction_coverage(known = known, predicted = predicted)
  piw <- prediction_interval_width(predicted = predicted)
  piw <- round(x = piw, digits = digits)
  piw <- paste0(piw[1], ", 95% CI: ",piw[2]," - ", piw[3])

  title <- "Predictive Interval Width Analysis (Efficiency)"
  subtitle <- paste0(
    "Predictive Interval Width: ", piw,
    "\nCoverage Probability: ", round(x = coverage, digits = digits)
  )

  grob <- with(piw_data,
    ggplot2::ggplot(data = piw_data) +
      ggplot2::geom_point(mapping = aes(x = Index, y = Known)) +
      ggplot2::geom_smooth(
        mapping = aes(
          x = Index, y = Center, ymin = Lower, ymax = Upper, alpha = 0.25
        ),
        show.legend = FALSE, stat = "identity"
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks(n = breaks),limits = interval
      ) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::ggtitle(label = title, subtitle = subtitle)
  )

  class(grob) <- c("ggsilent", class(grob))
  return(grob)

}

#' @author David Senhora Navega
#' @noRd
#' @details
#' Suppress ggplot messages on print
#' See https://tinyurl.com/yxw2bov4
print.ggsilent <- function(object) {
  suppressMessages(print(object))
}
