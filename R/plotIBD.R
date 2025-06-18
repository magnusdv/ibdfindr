
#' Plot IBD segments and posteriors
#'
#' @param data A data frame, typically produced by [ibdPosteriors()].
#' @param segments A data frame with IBD segments, typically produced by [findSegments()].
#' @param chrom A vector of chromosomes to plot (default: all).
#' @param ncol Number of columns in the plot. By default a suitable layout
#'   is chosen automatically.
#'
#' @returns A `ggplot2` plot.
#'
#' @seealso [findSegments()], [ibdPosteriors()]
#'
#' @examples
#' x = subset(cousinsDemo, CHROM %in% 11:13)
#'
#' post = ibdPosteriors(x, k1 = 0.2, a = 10)
#' segs = findSegments(x, k1 = 0.2, a = 10)
#' plotIBD(post, segs)
#'
#' @importFrom ggplot2 aes
#' @export
plotIBD = function(data, segments = NULL, chrom = NULL, ncol = NULL) {
  if(!is.null(chrom)) {
    data = data[data$chrom %in% chrom, , drop = FALSE]
    segments = segments[segments$chrom %in% chrom, , drop = FALSE] # NULL ok
  }
  chrs = unique(data$chrom)
  ncol = ncol %||% ceiling(sqrt(length(chrs)))

  p = ggplot2::ggplot(data, aes(x = cm)) +
    ggplot2::facet_wrap("chrom", ncol = ncol, scales = "free_x",
                        labeller = ggplot2::labeller(chrom = \(i) paste0("Chr", i))) +
    ggplot2::geom_jitter(aes(y = ibs/2), show.legend = FALSE,  color = "gray",
                         width = 0, height = .075, size = 1, alpha = 1) +
    ggplot2::geom_line(aes(y = post), col = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(
      name = "Posterior",
      limits = c(-0.21, 1.1), expand = c(0,0),
      breaks = c(0, 1),
      sec.axis = ggplot2::sec_axis(~., name = NULL, breaks = c(0, 0.5, 1), labels = paste0("IBS", 0:2))
    ) +
    ggplot2::labs(
      x = "Position (cM)",
      y = "Posterior IBD probability",
      title = "IBD segments and posteriors"
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(margin = ggplot2::margin(b=2), hjust = 0),
      axis.text.y.right  = ggplot2::element_text(color = 8, size = 6),
      axis.ticks.y.right = ggplot2::element_line(color = 8)
    )

  if(!is.null(segments)) {
    p = p + ggplot2::geom_segment(
      data = as.data.frame(segments), size = 1.5, color = "red",
      aes(x = start, xend = end, y = -0.15, yend = -0.15))
  }
  p
}


utils::globalVariables(
  c("cm","ibs","post","start","end")
)
