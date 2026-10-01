simulate_one_predictor_data <- function(k = 4, n_per_level = 30) {
  set.seed(11)
  data <- data.frame(
    group = factor(rep(LETTERS[seq_len(k)], each = n_per_level)),
    group_size = stats::rnorm(k * n_per_level)
  )
  data$y <- stats::rnorm(nrow(data), mean = c(0, 0.3, 1.5, 2)[seq_len(k)][data$group])
  data$success <- stats::rbinom(nrow(data), 1, c(0.05, 0.03, 0.5, 0.6)[seq_len(k)][data$group])
  data
}

quiet_contrasts <- function(...) {
  suppressMessages(ContrastsOnePred(...))
}

# Each coefficient must equal Scale_factor times the displayed difference.
displayed_coefficients <- function(result) {
  positive <- result[seq(1, nrow(result), by = 2), ]
  negative <- result[seq(2, nrow(result), by = 2), ]
  positive$Scale_factor * (positive$Estimate - negative$Estimate)
}


test_that("displayed differences reconstruct lm coefficients for every coding", {

  data <- simulate_one_predictor_data()

  codings <- list(
    treatment = stats::contr.treatment(4),
    treatment_base_3 = stats::contr.treatment(4, base = 3),
    sum = stats::contr.sum(4),
    sum_first_omitted = rbind(c(-1, -1, -1), diag(3)),
    helmert = stats::contr.helmert(4),
    reverse_helmert = reverse_helmert(4),
    sliding = MASS::contr.sdif(4),
    custom = cbind(c(-3, 1, 2, 0) / 6, c(0, -1, 1, 0) / 2, c(1, 1, 1, -3))
  )

  for (coding_name in names(codings)) {

    stats::contrasts(data$group) <- codings[[coding_name]]

    model <- stats::lm(y ~ group + group_size, data = data)

    result <- quiet_contrasts(
      model = model,
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "predicted"
    )

    expect_equal(
      displayed_coefficients(result),
      # The three group coefficients follow the intercept.
      unname(stats::coef(model)[2:4]),
      tolerance = 1e-8,
      info = coding_name
    )
  }
})


test_that("coefficients are matched by term, not by name prefix", {

  data <- simulate_one_predictor_data(k = 3)

  model <- stats::lm(y ~ group_size + group, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  p_values <- summary(model)$coefficients[c("groupB", "groupC"), "Pr(>|t|)"]

  expect_identical(
    result$Significance[c(1, 3)],
    ifelse(p_values < 0.05, "sig", "notsig"),
    ignore_attr = TRUE
  )
})


test_that("treatment contrasts with a non-first reference are returned", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- stats::contr.treatment(3, base = 3)

  result <- quiet_contrasts(
    model = stats::lm(y ~ group, data = data),
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_identical(result$Levels, c("A", "C", "B", "C"))
  expect_identical(attr(result, "contrast_type"), "Treatment")
})


test_that("sum contrasts with any omitted level are compared with the Grand Mean", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- rbind(c(-1, -1), c(1, 0), c(0, 1))

  model <- stats::lm(y ~ group, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_identical(result$Levels, c("B", "Grand Mean", "C", "Grand Mean"))
  expect_equal(displayed_coefficients(result), unname(stats::coef(model)[-1]))
  expect_identical(attr(result, "contrast_type"), "Sum")
})


test_that("non-orthogonal custom contrasts show weights and a divergence note", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- cbind(c(-3, 1, 2) / 6, c(0, -1, 1) / 2)

  expect_message(
    result <- ContrastsOnePred(
      model = stats::lm(y ~ group, data = data),
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "observed"
    ),
    "suggests C vs B"
  )

  expect_identical(result$Levels[3:4], c("W(A .20, C .80)", "B"))
  expect_equal(unique(result$Scale_factor), c(4 / 3, 10 / 9))
  expect_match(attr(result, "disclaimer"), "not orthogonal")
})


test_that("single-level estimates keep their link-scale confidence intervals", {

  data <- simulate_one_predictor_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)

  result <- quiet_contrasts(
    model = stats::glm(success ~ group + group_size, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "response"
  )

  single_levels <- result$Levels %in% LETTERS[1:4]

  expect_true(all(result$CI_low[single_levels] >= 0))
  expect_true(all(result$CI_high[single_levels] <= 1))
})


test_that("glm coefficients are reconstructed on the link scale", {

  data <- simulate_one_predictor_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)

  model <- stats::glm(success ~ group + group_size, family = stats::binomial, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "link"
  )

  expect_equal(displayed_coefficients(result), unname(stats::coef(model)[2:4]), tolerance = 1e-8)
})


test_that("an interaction with a numeric predictor triggers a warning", {

  data <- simulate_one_predictor_data(k = 3)

  expect_warning(
    quiet_contrasts(
      model = stats::lm(y ~ group * group_size, data = data),
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "predicted"
    ),
    "do not reconstruct"
  )
})


test_that("contrasts in data must match the fitted model", {

  data <- simulate_one_predictor_data(k = 3)
  model <- stats::lm(y ~ group, data = data)

  data$group <- factor(data$group, levels = c("A", "B", "C", "D"))

  expect_error(
    quiet_contrasts(model = model, predictor = "group", outcome = "y", data = data),
    "must be those used to fit the model"
  )
})


test_that("polynomial output carries the same attributes as other codings", {

  data <- simulate_one_predictor_data()
  data$ordered_group <- ordered(data$group)

  result <- quiet_contrasts(
    model = stats::lm(y ~ ordered_group, data = data),
    predictor = "ordered_group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  expect_identical(attr(result, "contrast_type"), "Polynomial")
  expect_identical(attr(result, "estimate_type"), "predicted")
  expect_true(all(c("Weight", "Polynomial_fit", "Display_Levels") %in% names(result)))
})


test_that("ModelContrastTable() uses the shared classification", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- rbind(c(-1, -1), c(1, 0), c(0, 1))

  variables <- ModelContrastTable(stats::lm(y ~ group, data = data), type = "x")

  expect_identical(
    variables$Contrast_Matrix_Type[variables$Variable == "group"],
    "Sum"
  )
})


test_that("other factors are averaged or held where the coefficients are defined", {

  set.seed(31)

  data <- expand.grid(
    group = factor(c("A", "B", "C")),
    other = factor(c("u", "v")),
    replicate = 1:25
  )
  data$x <- stats::rnorm(nrow(data), mean = data$other == "v")
  data$y <- stats::rnorm(
    nrow(data),
    as.integer(data$group) + 2 * (data$other == "v") + (data$group == "C") * (data$other == "v") + data$x
  )

  for (coding in list(stats::contr.treatment(2), stats::contr.sum(2))) {

    stats::contrasts(data$other) <- coding

    model <- stats::lm(y ~ group * other + x, data = data)

    expect_silent(
      result <- quiet_contrasts(
        model = model,
        predictor = "group",
        outcome = "y",
        data = data,
        estimate_type = "predicted"
      )
    )

    expect_equal(displayed_coefficients(result), unname(stats::coef(model)[2:3]), tolerance = 1e-8)
  }

  # Without an interaction, the other factor is averaged with equal weights.
  model <- stats::lm(y ~ group + other, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  balanced <- stats::predict(
    model,
    newdata = expand.grid(group = factor("A", levels = c("A", "B", "C")), other = factor(c("u", "v")))
  )

  expect_equal(result$Estimate[result$Levels == "A"][1], mean(balanced))
})
