# ============================================================
# TEXT LABELS THAT SHRINK TO AVOID OVERLAPS
# ============================================================
#
# The hierarchical group labels of MoPlot_multpred() are centred on their
# blocks and have a fixed size in points, whereas the distance between the
# blocks depends on the size of the plot. Their size is therefore fitted
# when the plot is drawn, and again whenever the device is resized.


#' Common scale that keeps neighbouring labels apart
#'
#' Labels in the same row (same `y`) are checked pairwise in order of `x`,
#' and each row is checked against the row above it, since the labels sit
#' on their baseline (`vjust = 0`). The scale is the largest value not
#' above one at which every pair of neighbours is separated by at least
#' `padding`, bounded below by `min_scale`. A single scale is returned for
#' all rows, so that the hierarchy keeps one text size.
#'
#' @param x Horizontal centres, in inches.
#' @param y Baselines, in inches; labels with the same value form a row.
#' @param widths,heights Label widths and heights at full size, in inches.
#' @param padding Minimum gap between neighbouring labels, in inches.
#' @param min_scale Smallest scale returned.
#'
#' @return A single number between `min_scale` and one.
#'
#' @keywords internal
#' @noRd
.fit_label_scale <- function(x,
                             y,
                             widths,
                             heights = rep(0, length(x)),
                             padding = 0.06,
                             min_scale = 0.6) {

  scale <- 1

  # Vertical: the tallest label of a row must fit below the next row.
  rows <- sort(unique(y))

  if (length(rows) > 1L) {

    row_heights <- vapply(rows, function(row) max(heights[y == row]), numeric(1))

    row_scales <- (diff(rows) - padding / 2) / utils::head(row_heights, -1L)

    scale <- min(scale, row_scales[is.finite(row_scales)])
  }

  # Horizontal: neighbours within a row must not touch.
  for (row in rows) {

    in_row <- which(y == row)

    if (length(in_row) < 2L) {
      next
    }

    in_row <- in_row[order(x[in_row])]

    centre_gaps <- diff(x[in_row])

    half_width_sums <- (utils::head(widths[in_row], -1L) + utils::tail(widths[in_row], -1L)) / 2

    pair_scales <- (centre_gaps - padding) / half_width_sums

    scale <- min(scale, pair_scales[is.finite(pair_scales)])
  }

  max(min_scale, min(1, scale))
}


#' Geom for group labels fitted at drawing time
#'
#' Behaves as `ggplot2::GeomText`, including its default aesthetics, but
#' draws a grob whose text size is fitted by `.fit_label_scale()` when it
#' is rendered.
#'
#' @keywords internal
#' @noRd
GeomMoplotFitText <- ggplot2::ggproto(
  "GeomMoplotFitText",
  ggplot2::GeomText,

  draw_panel = function(data,
                        panel_params,
                        coord,
                        padding = 0.06,
                        min_scale = 0.6) {

    data <- coord$transform(data, panel_params)

    grid::gTree(
      label_data = data,
      padding = padding,
      min_scale = min_scale,
      cl = "moplot_fit_text"
    )
  }
)


#' Fit the group labels to the space available on the device
#'
#' @param x A `"moplot_fit_text"` grob.
#'
#' @return The grob with its text children.
#'
#' @exportS3Method grid::makeContent
#' @noRd
makeContent.moplot_fit_text <- function(x) {

  data <- x$label_data

  font_sizes <- data$size * ggplot2::.pt

  font_faces <- if (is.null(data$fontface)) rep(1, nrow(data)) else data$fontface

  full_size_labels <- lapply(
    seq_len(nrow(data)),
    function(row) {
      grid::textGrob(
        data$label[row],
        gp = grid::gpar(fontsize = font_sizes[row], fontface = font_faces[row])
      )
    }
  )

  label_widths <- vapply(
    full_size_labels,
    function(label) grid::convertWidth(grid::grobWidth(label), "inches", valueOnly = TRUE),
    numeric(1)
  )

  label_heights <- vapply(
    full_size_labels,
    function(label) grid::convertHeight(grid::grobHeight(label), "inches", valueOnly = TRUE),
    numeric(1)
  )

  scale <- .fit_label_scale(
    x = grid::convertX(grid::unit(data$x, "npc"), "inches", valueOnly = TRUE),
    y = grid::convertY(grid::unit(data$y, "npc"), "inches", valueOnly = TRUE),
    widths = label_widths,
    heights = label_heights,
    padding = x$padding,
    min_scale = x$min_scale
  )

  text <- grid::textGrob(
    label = data$label,
    x = grid::unit(data$x, "npc"),
    y = grid::unit(data$y, "npc"),
    vjust = data$vjust,
    hjust = if (is.null(data$hjust)) 0.5 else data$hjust,
    gp = grid::gpar(
      col = ggplot2::alpha(data$colour, data$alpha),
      fontsize = font_sizes * scale,
      fontface = font_faces,
      fontfamily = if (is.null(data$family)) "" else data$family
    )
  )

  grid::setChildren(x, grid::gList(text))
}
