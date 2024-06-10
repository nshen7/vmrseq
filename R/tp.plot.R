#' Plot transition probability distribution
#'
#' @param tp 'transitProbs' object storing information about trained transition
#' probabilities. Can be obtained from function `estimTransitProbs`
#' @param line_size size of fitted loess smooth line. Default value is 0.2.
#' @param plot_train logical value indiating whether to plot training data.
#' Default is TRUE.
#' @param point_size size of training data points. Only applicable when
#' `plot_train = T`. Default value is 0.2.
#'
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom stats na.omit var
#' @importFrom tidyr pivot_longer
#'
#' @return A plot of the transition probability distribution
#' @export
#'
#' @examples
#' tp.plot(tp0)
#' 
tp.plot <- function(tp,
                    line_size = 0.2,
                    plot_train = T,
                    point_size = 0.2) {

  # compatibility checks within transProbs object
  if(nrow(tp@transit_probs) != tp@max_dist_bp)
    stop("Number of rows in `transit_probs` not equal to `maxdist_bp`.")
  if(nrow(tp@buffer_probs) != tp@buffer_bp)
    stop("Number of rows in `buffer_probs` not equal to `buffer_bp`.")
  if(nrow(tp@transit_probs) + nrow(tp@buffer_probs) != nrow(tp@train))
    stop("Number of rows in `transit_probs` and `buffer_probs` not adding up to number of rows in `train`.")

  plot_df <-
    data.frame(
      dist_bp = 1:nrow(tp@train),
      rbind(tp@transit_probs, tp@buffer_probs),
      tp@train
    ) %>%
    tidyr::pivot_longer(
      cols = -1,
      names_to = c(".value", "type"),
      names_pattern = "(.*)_(.*)"
    )

  type_labs <- c("P(0|0)","P(0|1)","P(1|0)","P(1|1)")
  names(type_labs) <- c("00","01","10","11")

  if (plot_train) {
    plot_df %>%
      ggplot() +
      geom_ribbon(aes(dist_bp, ymin = pbar-sqrt(var), ymax = pbar+sqrt(var)),
                  fill = "lightblue", size = point_size, alpha = 0.4) +
      geom_point(aes(dist_bp, pbar), color = "grey", size = point_size) +
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      geom_path(aes(dist_bp, phat), color = "red", size = line_size) +
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs))
  } else {
    plot_df %>%
      ggplot() +
      geom_path(aes(dist_bp, phat), color = "red", size = line_size) +
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs))
  }
}

