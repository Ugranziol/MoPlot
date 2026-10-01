# ============================================================
# CONFIDENCE INTERVALS FOR OBSERVED MEANS WITHIN THE NATURAL RANGE
# ============================================================
#
# Symmetric intervals for observed means can leave the natural range
# of the response: a proportion close to 0 may receive a negative
# lower bound. Observed (raw) estimates therefore use intervals that
# respect the range of the response:
#
#   - proportions (binary outcomes or binomial models): Wilson score
#     interval; for averages of several cells, the Wilson interval is
#     computed with the effective sample size p (1 - p) / SE^2
#     (Korn and Graubard), which reduces to the classical Wilson
#     interval for a single cell;
#   - positive responses (count and positive-continuous families):
#     interval on the log scale, exp(log(m) +/- z SE / m); a mean of
#     zero receives the exact Poisson upper bound -log(alpha / 2) / n;
#   - other responses: symmetric interval, mean +/- critical * SE.
# ============================================================


#' Natural range of the response for observed intervals
#'
#' @param family_object Model family, or `NULL`.
#' @param values Observed response values.
#'
#' @return One of `"proportion"`, `"positive"`, `"unbounded"`.
#'
#' @keywords internal
#' @noRd
.response_range_type <- function(family_object,
                                  values) {

  values <- values[!is.na(values)]

  if (length(values) == 0L || !is.numeric(values)) {
    return("unbounded")
  }

  family_name <- if (is.null(family_object)) "" else tolower(family_object$family)

  is_binomial_family <- family_name %in% c("binomial", "quasibinomial")

  if ((is_binomial_family || all(values %in% c(0, 1))) &&
      all(values >= 0 & values <= 1)) {

    return("proportion")
  }

  is_positive_family <- grepl(
    "poisson|negative binomial|nbinom|gamma|inverse.gaussian|tweedie|genpois|compois",
    family_name
  )

  if (is_positive_family && all(values >= 0)) {
    return("positive")
  }

  "unbounded"
}


#' Confidence intervals for observed means within the natural range
#'
#' @param estimate Observed means (or weighted averages of means).
#' @param se Their standard errors.
#' @param n Number of observations; for a single cell, the cell size. For an
#'   average of cells, the total number of observations, used only when the
#'   effective sample size cannot be computed (estimate at a bound).
#' @param range_type Output of `.response_range_type()`.
#' @param conf_level Confidence level.
#' @param critical Critical value for the symmetric interval of unbounded
#'   responses (for example a t quantile). Defaults to the normal quantile.
#' @param single_cell Logical. For proportions, should `n` be used directly
#'   (single cells) rather than the effective sample size (averages)?
#'
#' @return A list with `ci_low` and `ci_high`.
#'
#' @keywords internal
#' @noRd
.observed_interval <- function(estimate,
                               se,
                               n,
                               range_type,
                               conf_level = 0.95,
                               critical = NULL,
                               single_cell = TRUE) {

  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (is.null(critical)) {
    critical <- z_value
  }

  n <- rep_len(n, length(estimate))
  critical <- rep_len(critical, length(estimate))

  ci_low <- estimate - critical * se
  ci_high <- estimate + critical * se

  if (identical(range_type, "proportion")) {

    # Effective sample size of an average of proportions.
    sample_size <- n

    if (!single_cell) {

      effective <- estimate * (1 - estimate) / se^2

      usable <- is.finite(effective) & estimate > 0 & estimate < 1 & se > 0

      sample_size[usable] <- effective[usable]
    }

    valid <- is.finite(estimate) & is.finite(sample_size) & sample_size > 0

    denominator <- 1 + z_value^2 / sample_size
    centre <- (estimate + z_value^2 / (2 * sample_size)) / denominator
    half_width <- z_value / denominator *
      sqrt(estimate * (1 - estimate) / sample_size + z_value^2 / (4 * sample_size^2))

    ci_low[valid] <- pmax(0, centre[valid] - half_width[valid])
    ci_high[valid] <- pmin(1, centre[valid] + half_width[valid])

    # Exact bounds at the edges, avoiding floating-point residues.
    ci_low[valid & estimate == 0] <- 0
    ci_high[valid & estimate == 1] <- 1
  }

  if (identical(range_type, "positive")) {

    on_log_scale <- is.finite(estimate) & estimate > 0 & is.finite(se)

    ci_low[on_log_scale] <- exp(log(estimate[on_log_scale]) - z_value * se[on_log_scale] / estimate[on_log_scale])
    ci_high[on_log_scale] <- exp(log(estimate[on_log_scale]) + z_value * se[on_log_scale] / estimate[on_log_scale])

    # All observations equal to zero: exact Poisson upper bound.
    at_zero <- is.finite(estimate) & estimate == 0 & is.finite(n) & n > 0

    ci_low[at_zero] <- 0
    ci_high[at_zero] <- -log((1 - conf_level) / 2) / n[at_zero]
  }

  list(
    ci_low = ci_low,
    ci_high = ci_high
  )
}
