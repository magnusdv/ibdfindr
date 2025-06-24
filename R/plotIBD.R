
#' Plot IBD segments and posteriors
#'
#' @param x A data frame, typically produced by [ibdPosteriors()], or a list
#'   produced by [findIBD()].
#' @param segments A data frame with IBD segments, typically produced by
#'   [findSegments()].
#' @param chrom A vector of chromosomes to plot (default: all).
#' @param ncol Number of columns in the plot. By default a suitable layout is
#'   chosen automatically.
#' @param title Plot title. Generated automatically if `NA` (default).
#' @param base_size Base font size.
#'
#' @returns A `ggplot2` plot.
#'
#' @seealso [findSegments()], [ibdPosteriors()]
#'
#' @examples
#' x = subset(cousinsDemo, CHROM %in% 3:4)
#' ibd = findIBD(x, k1 = 0.2, a = 5)
#' plotIBD(ibd)
#'
#' @importFrom ggplot2 aes
#' @export
plotIBD = function(x, segments = NULL, chrom = NULL, ncol = NULL,
                   title = NA, base_size = 12) {
  if(is.data.frame(x)) {
    data = x
    ids = NULL
  }
  else { # x is a list from findIBD()
    if(!"posteriors" %in% names(x))
      stop2("If `x` is a list, it must contain an element named 'posteriors'")
    data = x$posteriors
    segments = x$segments
    ids = x$ids
  }

  if(!"ibs" %in% names(data))
    stop2("Expected the input to contain a column named `ibs`")

  if(!is.null(chrom)) {
    data = data[data$chrom %in% chrom, , drop = FALSE]
    segments = segments[segments$chrom %in% chrom, , drop = FALSE] # NULL ok
  }
  chrs = unique(data$chrom)
  ncol = ncol %||% min(4, ceiling(sqrt(length(chrs))))

  # Title
  if(is.na(title) && is.null(ids))
    title = NULL

  if(is.na(title) && !is.null(ids)) {
    idsString = paste(ids, collapse = " vs. ")
    if(!is.null(segments))
      title = paste("IBD segments and posteriors for", idsString)
    else
      title = paste("IBD posteriors for", idsString)
  }

  p = ggplot2::ggplot(data, aes(x = cm)) +
    ggplot2::facet_wrap("chrom", ncol = ncol, scales = "free_x",
                        labeller = ggplot2::labeller(chrom = \(i) paste0("Chr", i))) +
    ggplot2::geom_jitter(aes(y = ibs/2), show.legend = FALSE,  color = "gray",
                         width = 0, height = .075, size = 1, alpha = 1) +
    ggplot2::geom_line(aes(y = post), col = 1) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::scale_y_continuous(
      name = "IBD posterior",
      limits = c(-0.21, 1.1), expand = c(0,0),
      breaks = c(0, 1),
      sec.axis = ggplot2::sec_axis(~., name = NULL, breaks = c(0, 0.5, 1), labels = paste0("IBS", 0:2))
    ) +
    ggplot2::labs(
      x = "Position (cM)",
      y = "Posterior IBD probability",
      title = title
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(margin = ggplot2::margin(b=2), hjust = 0),
      axis.text.y.right  = ggplot2::element_text(color = 8, size = 7),
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
