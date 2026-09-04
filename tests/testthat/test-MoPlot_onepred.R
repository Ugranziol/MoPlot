test_that("MoPlot_onepred() handles Reverse Helmert contrasts", {

  set.seed(123)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 20))
  )

  dat$y <- c(
    rnorm(20, 10, 2),
    rnorm(20, 12, 2),
    rnorm(20, 15, 2),
    rnorm(20, 18, 2)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  plot <- MoPlot_onepred(
    model = mod,
    predictor = "group",
    caption = FALSE
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("MoPlot_onepred() displays raw data with show_data = TRUE", {

  set.seed(123)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 20))
  )

  dat$y <- c(
    rnorm(20, 10, 2),
    rnorm(20, 12, 2),
    rnorm(20, 15, 2),
    rnorm(20, 18, 2)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  plot <- MoPlot_onepred(
    model = mod,
    predictor = "group",
    show_data = TRUE,
    caption = FALSE
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("MoPlot_onepred() publication style assigns distinct colours to contrasts", {

  set.seed(123)

  dat <- data.frame(
    group = factor(
      rep(c("A", "B", "C", "D"), each = 20),
      levels = c("A", "B", "C", "D")
    )
  )

  dat$y <- c(
    rnorm(20, 10, 1),
    rnorm(20, 12, 1),
    rnorm(20, 14, 1),
    rnorm(20, 16, 1)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  prep <- .prepare_moplot_onepred(
    model = mod,
    predictor = "group",
    estimate_type = "observed",
    prediction_scale = "response",
    coefficient_scale = "model"
  )

  upper <- .plot_moplot_onepred_upper(
    prep = prep,
    style = "publication"
  )

  build <- ggplot2::ggplot_build(
    upper
  )

  segment_layer_index <- which(
    vapply(
      build$data,
      function(layer) {
        all(c("x", "y", "xend", "yend") %in% names(layer))
      },
      logical(1)
    )
  )

  expect_length(
    segment_layer_index,
    1L
  )

  segment_layer <- build$data[[segment_layer_index]]

  expect_gt(
    length(unique(segment_layer$colour)),
    1L
  )
})

test_that("MoPlot_onepred() minimal style keeps contrast colours in grayscale", {

  set.seed(123)

  dat <- data.frame(
    group = factor(
      rep(c("A", "B", "C"), each = 20)
    )
  )

  dat$y <- rnorm(
    nrow(dat)
  )

  contrasts(dat$group) <- reverse_helmert(3)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  prep <- .prepare_moplot_onepred(
    model = mod,
    predictor = "group",
    estimate_type = "observed",
    prediction_scale = "response",
    coefficient_scale = "model"
  )

  upper <- .plot_moplot_onepred_upper(
    prep = prep,
    style = "minimal"
  )

  build <- ggplot2::ggplot_build(upper)

  segment_layer_index <- which(
    vapply(
      build$data,
      function(layer) {
        all(c("x", "y", "xend", "yend") %in% names(layer))
      },
      logical(1)
    )
  )

  expect_length(
    segment_layer_index,
    1L
  )

  segment_layer <- build$data[[segment_layer_index]]

  expect_equal(
    length(unique(segment_layer$colour)),
    1L
  )
})
