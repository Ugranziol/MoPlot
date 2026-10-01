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


test_that("MoPlot_onepred() handles sliding-difference contrasts", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C", "D"), each = 20)))
  dat$y <- stats::rnorm(80, c(10, 12, 15, 18)[dat$group], 2)

  contrasts(dat$group) <- MASS::contr.sdif(4)

  plot <- suppressMessages(
    MoPlot_onepred(
      model = stats::lm(y ~ group, data = dat),
      predictor = "group",
      caption = FALSE
    )
  )

  expect_s3_class(plot, "patchwork")
})


test_that("per_difference expresses coefficients as the displayed differences", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C", "D"), each = 20)))
  dat$y <- stats::rnorm(80, c(10, 12, 15, 18)[dat$group], 2)

  contrasts(dat$group) <- stats::contr.helmert(4)

  model <- stats::lm(y ~ group, data = dat)

  prepare <- function(per_difference) {
    suppressMessages(
      .prepare_moplot_onepred(
        model = model,
        predictor = "group",
        estimate_type = "observed",
        per_difference = per_difference
      )
    )
  }

  coefficients <- prepare(FALSE)
  differences <- prepare(TRUE)

  segments <- differences$segment_data
  displayed <- segments$yend - segments$y

  expect_equal(differences$parameter_data$Display_estimate, displayed, tolerance = 1e-8)
  expect_equal(coefficients$parameter_data$Display_estimate, displayed / (2:4), tolerance = 1e-8)

  # Rescaling does not change the inferential quantities.
  expect_equal(differences$parameter_data$d, coefficients$parameter_data$d)

  expect_equal(
    differences$parameter_data$Display_CI_high - differences$parameter_data$Display_CI_low,
    (coefficients$parameter_data$Display_CI_high - coefficients$parameter_data$Display_CI_low) * (2:4)
  )
})


test_that("per_difference rescales ratios on the log scale", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C"), each = 60)))
  dat$success <- stats::rbinom(180, 1, c(0.2, 0.4, 0.7)[dat$group])

  contrasts(dat$group) <- stats::contr.helmert(3)

  model <- stats::glm(success ~ group, family = stats::binomial, data = dat)

  prep <- suppressMessages(
    .prepare_moplot_onepred(
      model = model,
      predictor = "group",
      estimate_type = "predicted",
      coefficient_scale = "interpreted",
      per_difference = TRUE
    )
  )

  expect_equal(
    prep$parameter_data$Display_estimate,
    exp(unname(stats::coef(model)[-1]) / c(1 / 2, 1 / 3)),
    tolerance = 1e-8
  )
})
