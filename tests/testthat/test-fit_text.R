test_that(".fit_label_scale() shrinks labels only when neighbours would touch", {

  # Two labels 1.5 inches wide whose centres are 2 inches apart fit.
  expect_identical(.fit_label_scale(x = c(0, 2), y = c(1, 1), widths = c(1.5, 1.5)), 1)

  # One inch apart they need (1 - 0.06) / 1.5 of their size.
  expect_equal(.fit_label_scale(x = c(0, 1), y = c(1, 1), widths = c(1.5, 1.5)), 0.94 / 1.5)

  # Labels in different rows do not constrain each other horizontally.
  expect_identical(.fit_label_scale(x = c(0, 0.1), y = c(1, 2), widths = c(1.5, 1.5)), 1)

  # A row must fit below the next one: 0.25 inches tall, 0.2 inches apart.
  expect_equal(
    .fit_label_scale(x = c(0, 0), y = c(1, 1.2), widths = c(1, 1), heights = c(0.25, 0.25)),
    (0.2 - 0.03) / 0.25
  )

  # The scale never drops below min_scale.
  expect_identical(.fit_label_scale(x = c(0, 0.1), y = c(1, 1), widths = c(2, 2)), 0.6)
})


test_that("group labels are fitted to the size of the device", {

  labels <- grid::gTree(
    label_data = data.frame(
      x = c(0.30, 0.55),
      y = c(0.9, 0.9),
      label = c("Delayed", "Mean(Delayed,Immediate)"),
      size = 3.87,
      fontface = "bold",
      colour = "black",
      alpha = NA,
      vjust = 0
    ),
    padding = 0.06,
    min_scale = 0.6,
    cl = "moplot_fit_text"
  )

  drawn_font_size <- function(width) {
    grDevices::pdf(NULL, width = width, height = 4)
    on.exit(grDevices::dev.off())
    grid::pushViewport(grid::viewport())
    fitted <- grid::makeContent(labels)
    unique(fitted$children[[1L]]$gp$fontsize)
  }

  full_size <- 3.87 * ggplot2::.pt

  expect_equal(drawn_font_size(12), full_size)
  expect_lt(drawn_font_size(5), full_size)
  expect_gte(drawn_font_size(5), 0.6 * full_size)
})


test_that("MoPlot_multpred() draws its group labels with the fitted geom", {

  set.seed(6)
  data <- expand.grid(
    encoding = factor(c("Shallow", "Deep", "Self-ref")),
    delay = factor(c("Immediate", "Delayed")),
    id = 1:10
  )
  data$recall <- stats::rnorm(nrow(data))
  stats::contrasts(data$encoding) <- stats::contr.helmert(3)
  stats::contrasts(data$delay) <- stats::contr.sum(2)

  model <- stats::lm(recall ~ encoding * delay, data = data)

  for (style in c("minimal", "tree")) {

    plot <- suppressMessages(suppressWarnings(
      MoPlot_multpred(model, data = data, style = style, caption = FALSE, interaction_info = FALSE)
    ))

    geoms <- vapply(plot[[1L]]$layers, function(layer) class(layer$geom)[1L], character(1))

    expect_true("GeomMoplotFitText" %in% geoms, info = style)

    # Drawing on a small device does not fail.
    grDevices::pdf(NULL, width = 5, height = 4)
    expect_no_error(print(plot))
    grDevices::dev.off()
  }
})
