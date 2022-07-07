#' Plot fitted transition probabilities
#'
#' @param tp
#' @param lineSize size of fitted loess smooth line. Default value is 0.2.
#' @param plotTrain logical value indiating whether to plot training data.
#' Default is TRUE.
#' @param pointSize size of training data points. Only applicable when
#' `plotTrain = T`. Default value is 0.2.
#'
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return
#' @export
#'
#' @examples
#'
plotTransitProbs <- function(tp,
                             lineSize = 0.2,
                             plotTrain = T,
                             pointSize = 0.2) {

  # compatibility checks within transProbs object
  if(nrow(tp@transit_probs) != tp@max_dist_bp)
    stop("Number of rows in `transit_probs` not equal to `maxdist_bp`.")
  if(nrow(tp@buffer_probs) != tp@buffer_bp)
    stop("Number of rows in `buffer_probs` not equal to `buffer_bp`.")
  if(nrow(tp@transit_probs) + nrow(tp@buffer_probs) != nrow(tp@train))
    stop("Number of rows in `transit_probs` and `buffer_probs` not adding up to number of rows in `train`.")

  plot_df <- data.frame(dist_bp = 1:nrow(tp@train), rbind(tp@transit_probs, tp@buffer_probs), tp@train) %>%
    dplyr::select(-starts_with("var")) %>%
    pivot_longer(cols = -1, names_to = c(".value", "type"),
                 names_pattern = "(.*)_(.*)")

  type_labs <- c("P(0|0)","P(0|1)","P(1|0)","P(1|1)")
  names(type_labs) <- c("00","01","10","11")

  if (plotTrain) {
    plot_df %>%
      ggplot() +
      geom_point(aes(dist_bp, pbar), color = "grey", size = pointSize) +
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      geom_path(aes(dist_bp, phat), color = "red", size = lineSize) +
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs))
  } else {
    plot_df %>%
      ggplot() +
      geom_path(aes(dist_bp, phat), color = "red", size = lineSize) +
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs))
  }
}

