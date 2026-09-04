#' Extract and interpret model parameters
#'
#' @description
#' `extract_param()` extracts model coefficients and inferential statistics
#' using `parameters::model_parameters()`. It returns a compact,
#' non-redundant parameter table and, when appropriate, transforms coefficients
#' and confidence intervals into a more directly interpretable metric.
#'
#' Depending on the model family and link function, the interpreted estimate
#' may represent an odds ratio, incidence rate ratio, hazard ratio, or another
#' multiplicative ratio. For models for which no meaningful automatic
#' transformation is available, the original coefficient is retained.
#'
#' The function also attempts to compute an approximate Cohen's d from each
#' coefficient-level t or z statistic.
#'
#' @param model A fitted statistical model supported by
#'   `parameters::model_parameters()`.
#'
#' @param ci A single numeric value between 0 and 1 specifying the confidence
#'   level. The default is `0.95`.
#'
#' @param include_intercept Logical. Should the intercept be included in the
#'   returned table? The default is `FALSE`.
#'
#' @param effects Character string specifying the model effects to extract for
#'   models that distinguish between fixed and random effects. The default is
#'   `"fixed"`.
#'
#' @param component Character string specifying the model component to extract
#'   for models with multiple components, such as zero-inflated or hurdle
#'   models. The default is `"conditional"`.
#'
#' @param transform Character string specifying how model coefficients should
#'   be transformed for interpretation. Possible values are:
#'
#'   \describe{
#'     \item{`"auto"`}{Automatically exponentiate coefficients from models
#'       with log or logit links and from Cox proportional-hazards models.}
#'     \item{`"none"`}{Retain the original coefficient and confidence interval
#'       as the interpreted estimate.}
#'     \item{`"exponentiate"`}{Exponentiate the coefficient and confidence
#'       interval regardless of the detected model family or link.}
#'   }
#'
#'   The default is `"auto"`.
#'
#' @details
#' The `Estimate` column always contains the coefficient on the scale on which
#' the model was fitted. For generalized models, this is generally the scale of
#' the linear predictor.
#'
#' When `transform = "auto"`, the following transformations are applied:
#'
#' \describe{
#'   \item{Identity link}{No transformation. `Interpreted_Estimate` is equal
#'     to `Estimate`.}
#'   \item{Binomial model with logit link}{The coefficient is exponentiated
#'     and interpreted as an odds ratio.}
#'   \item{Poisson or negative-binomial model with log link}{The coefficient
#'     is exponentiated and interpreted as an incidence rate ratio.}
#'   \item{Other models with a log link}{The coefficient is exponentiated and
#'     interpreted as a multiplicative ratio.}
#'   \item{Cox proportional-hazards model}{The coefficient is exponentiated
#'     and interpreted as a hazard ratio.}
#' }
#'
#' No automatic exponentiation is performed for links such as `probit`,
#' `cauchit`, or `cloglog`, because exponentiating an individual coefficient
#' does not generally produce a standard effect measure with a direct
#' interpretation.
#'
#' The interpreted estimates are not necessarily predictions on the response
#' scale. For example, an odds ratio is a ratio of odds rather than a predicted
#' probability. Predicted probabilities or expected responses require the
#' intercept and values of the other model predictors.
#'
#' Approximate Cohen's d values are computed as follows:
#'
#' \itemize{
#'   \item t statistics are converted using `effectsize::t_to_d()`;
#'   \item z statistics are converted using `effectsize::z_to_d()`.
#' }
#'
#' For t statistics, coefficient-specific degrees of freedom are used when
#' available. For z statistics, the number of model observations is obtained
#' using `stats::nobs()`.
#'
#' These conversions are approximations and assume conditions implicit in the
#' conversion formulas. In particular, the conversion to Cohen's d assumes two
#' groups of approximately equal size. Directly estimated effect sizes should
#' be preferred when the necessary data and an appropriate estimator are
#' available.
#'
#' @return
#' A data frame with one row for each extracted model parameter and the
#' following columns:
#'
#' \describe{
#'   \item{Parameter}{Name of the model parameter.}
#'   \item{Estimate}{Untransformed model coefficient on the model scale.}
#'   \item{SE}{Standard error of the untransformed coefficient.}
#'   \item{CI}{Confidence level used for the interval.}
#'   \item{CI_low}{Lower confidence limit of `Estimate`.}
#'   \item{CI_high}{Upper confidence limit of `Estimate`.}
#'   \item{Statistic}{Coefficient-level test statistic.}
#'   \item{Statistic_type}{Type of test statistic, usually `"t"` or `"z"`.}
#'   \item{df}{Coefficient-level degrees of freedom, when available.}
#'   \item{p}{Coefficient-level p-value, when available.}
#'   \item{Interpreted_Estimate}{Coefficient transformed into the
#'     interpretation indicated by `Estimate_type`.}
#'   \item{Interpreted_CI_low}{Lower confidence limit of the interpreted
#'     estimate.}
#'   \item{Interpreted_CI_high}{Upper confidence limit of the interpreted
#'     estimate.}
#'   \item{Estimate_type}{Interpretation of `Interpreted_Estimate`, such as
#'     `"Coefficient"`, `"Odds ratio"`, `"Incidence rate ratio"`,
#'     `"Hazard ratio"`, or `"Multiplicative ratio"`.}
#'   \item{d}{Approximate Cohen's d derived from the test statistic.}
#'   \item{d_CI_low}{Lower confidence limit of approximate Cohen's d.}
#'   \item{d_CI_high}{Upper confidence limit of approximate Cohen's d.}
#' }
#'
#' @note
#' The coefficient transformation is based on the model family and link
#' detected through `stats::family()` whenever available. Some model classes
#' may use non-standard structures or parameterizations and may therefore
#' require `transform = "none"` or `transform = "exponentiate"` to override
#' automatic behaviour.
#'
#' In models with multiple components, the interpretation may depend on the
#' component selected through `component`.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' # Linear regression
#' dat_lm <- data.frame(
#'   y = rnorm(120),
#'   x = rnorm(120),
#'   group = factor(rep(c("A", "B", "C"), each = 40))
#' )
#'
#' contrasts(dat_lm$group) <- contr.sum(3)
#'
#' mod_lm <- lm(y ~ x + group, data = dat_lm)
#'
#' extract_param(mod_lm)
#'
#' # Binomial logistic regression
#' dat_glm <- data.frame(
#'   y = rbinom(200, size = 1, prob = 0.40),
#'   x = rnorm(200),
#'   group = factor(sample(c("A", "B", "C"), 200, replace = TRUE))
#' )
#'
#' contrasts(dat_glm$group) <- contr.sum(3)
#'
#' mod_glm <- glm(
#'   y ~ x + group,
#'   data = dat_glm,
#'   family = binomial(link = "logit")
#' )
#'
#' extract_param(mod_glm)
#'
#' # Disable automatic coefficient transformation
#' extract_param(mod_glm, transform = "none")
#' }
#'
#' @importFrom parameters model_parameters
#' @importFrom effectsize t_to_d z_to_d
#' @importFrom stats family nobs
#' @export

extract_param <- function(
    model,
    ci = 0.95,
    include_intercept = FALSE,
    effects = "fixed",
    component = "conditional",
    transform = c("auto", "none", "exponentiate")
) {

  # ==========================================================================
  # Argument checks
  # ==========================================================================

  transform <- match.arg(transform)

  if (
    !is.numeric(ci) ||
    length(ci) != 1L ||
    is.na(ci) ||
    !is.finite(ci) ||
    ci <= 0 ||
    ci >= 1
  ) {
    stop(
      "`ci` must be a single numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.logical(include_intercept) ||
    length(include_intercept) != 1L ||
    is.na(include_intercept)
  ) {
    stop(
      "`include_intercept` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.character(effects) ||
    length(effects) != 1L ||
    is.na(effects) ||
    !nzchar(effects)
  ) {
    stop(
      "`effects` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  if (
    !is.character(component) ||
    length(component) != 1L ||
    is.na(component) ||
    !nzchar(component)
  ) {
    stop(
      "`component` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Internal helpers
  # ==========================================================================

  first_existing_column <- function(data, candidates) {

    available <- candidates[candidates %in% names(data)]

    if (length(available) == 0L) {
      return(NULL)
    }

    available[1L]
  }

  extract_numeric_column <- function(data, column, n) {

    if (is.null(column)) {
      return(rep(NA_real_, n))
    }

    suppressWarnings(
      as.numeric(data[[column]])
    )
  }

  remove_intercept <- function(data) {

    if (
      !include_intercept &&
      "Parameter" %in% names(data)
    ) {
      data <- data[
        data$Parameter != "(Intercept)",
        ,
        drop = FALSE
      ]
    }

    data
  }

  # ==========================================================================
  # Detect model family and link
  # ==========================================================================

  family_object <- tryCatch(
    stats::family(model),
    error = function(e) NULL
  )

  family_name <- if (
    !is.null(family_object) &&
    !is.null(family_object$family)
  ) {
    tolower(as.character(family_object$family)[1L])
  } else {
    NA_character_
  }

  link_name <- if (
    !is.null(family_object) &&
    !is.null(family_object$link)
  ) {
    tolower(as.character(family_object$link)[1L])
  } else {
    NA_character_
  }

  is_cox_model <- inherits(model, "coxph")

  # ==========================================================================
  # Determine transformation and interpretation
  # ==========================================================================

  exponentiate_estimate <- switch(
    transform,
    none = FALSE,
    exponentiate = TRUE,
    auto = {
      is_cox_model ||
        isTRUE(link_name %in% c("log", "logit"))
    }
  )

  if (!exponentiate_estimate) {

    estimate_type <- "Coefficient"

  } else if (is_cox_model) {

    estimate_type <- "Hazard ratio"

  } else if (identical(link_name, "logit")) {

    estimate_type <- "Odds ratio"

  } else if (
    identical(link_name, "log") &&
    grepl(
      "poisson|negative binomial|negative-binomial|nbinom",
      family_name
    )
  ) {

    estimate_type <- "Incidence rate ratio"

  } else if (identical(link_name, "log")) {

    estimate_type <- "Multiplicative ratio"

  } else {

    estimate_type <- "Exponentiated coefficient"
  }

  # ==========================================================================
  # Extract untransformed parameters
  # ==========================================================================

  raw_parameters <- tryCatch(
    parameters::model_parameters(
      model,
      ci = ci,
      effects = effects,
      component = component,
      exponentiate = FALSE,
      pretty_names = FALSE,
      verbose = FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "Model parameters could not be extracted. ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  raw_parameters <- as.data.frame(raw_parameters)
  raw_parameters <- remove_intercept(raw_parameters)
  rownames(raw_parameters) <- NULL

  if (!"Parameter" %in% names(raw_parameters)) {
    stop(
      paste0(
        "The table returned by `parameters::model_parameters()` ",
        "does not contain a `Parameter` column."
      ),
      call. = FALSE
    )
  }

  n_parameters <- nrow(raw_parameters)

  # ==========================================================================
  # Identify columns returned by model_parameters()
  # ==========================================================================

  estimate_column <- first_existing_column(
    raw_parameters,
    c("Coefficient", "Estimate")
  )

  se_column <- first_existing_column(
    raw_parameters,
    c("SE", "Std_Error", "Std. Error")
  )

  ci_column <- first_existing_column(
    raw_parameters,
    "CI"
  )

  ci_low_column <- first_existing_column(
    raw_parameters,
    c("CI_low", "CI_Low", "Conf.Low")
  )

  ci_high_column <- first_existing_column(
    raw_parameters,
    c("CI_high", "CI_High", "Conf.High")
  )

  p_column <- first_existing_column(
    raw_parameters,
    c("p", "p_value", "P")
  )

  df_column <- first_existing_column(
    raw_parameters,
    c("df_error", "df", "df_residual")
  )

  # ==========================================================================
  # Identify test statistic and statistic type
  # ==========================================================================

  if ("t" %in% names(raw_parameters)) {

    statistic_column <- "t"
    statistic_type <- rep("t", n_parameters)

  } else if ("z" %in% names(raw_parameters)) {

    statistic_column <- "z"
    statistic_type <- rep("z", n_parameters)

  } else if ("Statistic" %in% names(raw_parameters)) {

    statistic_column <- "Statistic"
    statistic_type <- rep(NA_character_, n_parameters)

    statistic_attribute <- attr(
      raw_parameters$Statistic,
      "statistic"
    )

    if (
      !is.null(statistic_attribute) &&
      grepl("^t$", statistic_attribute, ignore.case = TRUE)
    ) {
      statistic_type[] <- "t"
    }

    if (
      !is.null(statistic_attribute) &&
      grepl("^z$", statistic_attribute, ignore.case = TRUE)
    ) {
      statistic_type[] <- "z"
    }

    if (all(is.na(statistic_type))) {

      if (!is.null(df_column)) {

        df_values <- suppressWarnings(
          as.numeric(raw_parameters[[df_column]])
        )

        finite_df <- is.finite(df_values)

        statistic_type[finite_df] <- "t"
        statistic_type[!finite_df] <- "z"

      } else {

        statistic_type[] <- "z"
      }
    }

  } else {

    statistic_column <- NULL
    statistic_type <- rep(NA_character_, n_parameters)
  }

  # ==========================================================================
  # Build the non-redundant output table
  # ==========================================================================

  output <- data.frame(
    Parameter = as.character(raw_parameters$Parameter),

    Estimate = extract_numeric_column(
      raw_parameters,
      estimate_column,
      n_parameters
    ),

    SE = extract_numeric_column(
      raw_parameters,
      se_column,
      n_parameters
    ),

    CI = if (!is.null(ci_column)) {
      extract_numeric_column(
        raw_parameters,
        ci_column,
        n_parameters
      )
    } else {
      rep(ci, n_parameters)
    },

    CI_low = extract_numeric_column(
      raw_parameters,
      ci_low_column,
      n_parameters
    ),

    CI_high = extract_numeric_column(
      raw_parameters,
      ci_high_column,
      n_parameters
    ),

    Statistic = extract_numeric_column(
      raw_parameters,
      statistic_column,
      n_parameters
    ),

    Statistic_type = statistic_type,

    df = extract_numeric_column(
      raw_parameters,
      df_column,
      n_parameters
    ),

    p = extract_numeric_column(
      raw_parameters,
      p_column,
      n_parameters
    ),

    stringsAsFactors = FALSE
  )

  # ==========================================================================
  # Extract interpreted estimates
  # ==========================================================================

  if (exponentiate_estimate) {

    interpreted_parameters <- tryCatch(
      parameters::model_parameters(
        model,
        ci = ci,
        effects = effects,
        component = component,
        exponentiate = TRUE,
        pretty_names = FALSE,
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(interpreted_parameters)) {

      interpreted_parameters <- as.data.frame(
        interpreted_parameters
      )

      interpreted_parameters <- remove_intercept(
        interpreted_parameters
      )

      interpreted_estimate_column <- first_existing_column(
        interpreted_parameters,
        c("Coefficient", "Estimate")
      )

      interpreted_ci_low_column <- first_existing_column(
        interpreted_parameters,
        c("CI_low", "CI_Low", "Conf.Low")
      )

      interpreted_ci_high_column <- first_existing_column(
        interpreted_parameters,
        c("CI_high", "CI_High", "Conf.High")
      )

      if ("Parameter" %in% names(interpreted_parameters)) {

        parameter_match <- match(
          output$Parameter,
          interpreted_parameters$Parameter
        )

        output$Interpreted_Estimate <-
          if (!is.null(interpreted_estimate_column)) {
            suppressWarnings(
              as.numeric(
                interpreted_parameters[[interpreted_estimate_column]][parameter_match]
              )
            )
          } else {
            rep(NA_real_, nrow(output))
          }

        output$Interpreted_CI_low <-
          if (!is.null(interpreted_ci_low_column)) {
            suppressWarnings(
              as.numeric(
                interpreted_parameters[[interpreted_ci_low_column]][parameter_match]
              )
            )
          } else {
            rep(NA_real_, nrow(output))
          }

        output$Interpreted_CI_high <-
          if (!is.null(interpreted_ci_high_column)) {
            suppressWarnings(
              as.numeric(
                interpreted_parameters[[interpreted_ci_high_column]][parameter_match]
              )
            )
          } else {
            rep(NA_real_, nrow(output))
          }

      } else {

        output$Interpreted_Estimate <- NA_real_
        output$Interpreted_CI_low <- NA_real_
        output$Interpreted_CI_high <- NA_real_
      }

    } else {

      output$Interpreted_Estimate <- NA_real_
      output$Interpreted_CI_low <- NA_real_
      output$Interpreted_CI_high <- NA_real_
    }

  } else {

    output$Interpreted_Estimate <- output$Estimate
    output$Interpreted_CI_low <- output$CI_low
    output$Interpreted_CI_high <- output$CI_high
  }

  output$Estimate_type <- rep(
    estimate_type,
    nrow(output)
  )

  # ==========================================================================
  # Extract sample size for z-to-d conversions
  # ==========================================================================

  n_observations <- tryCatch(
    as.numeric(stats::nobs(model)),
    error = function(e) NA_real_
  )

  if (
    length(n_observations) != 1L ||
    is.na(n_observations) ||
    !is.finite(n_observations) ||
    n_observations <= 0
  ) {
    n_observations <- NA_real_
  }

  # ==========================================================================
  # Compute approximate Cohen's d
  # ==========================================================================

  output$d <- NA_real_
  output$d_CI_low <- NA_real_
  output$d_CI_high <- NA_real_

  for (i in seq_len(nrow(output))) {

    statistic_value <- output$Statistic[i]
    current_statistic_type <- output$Statistic_type[i]
    df_value <- output$df[i]

    if (
      is.na(statistic_value) ||
      !is.finite(statistic_value) ||
      is.na(current_statistic_type)
    ) {
      next
    }

    d_result <- tryCatch(
      {
        if (
          identical(current_statistic_type, "t") &&
          !is.na(df_value) &&
          is.finite(df_value) &&
          df_value > 0
        ) {

          effectsize::t_to_d(
            t = statistic_value,
            df_error = df_value,
            ci = ci
          )

        } else if (
          identical(current_statistic_type, "z") &&
          !is.na(n_observations) &&
          is.finite(n_observations)
        ) {

          effectsize::z_to_d(
            z = statistic_value,
            n = n_observations,
            ci = ci
          )

        } else {

          NULL
        }
      },
      error = function(e) NULL
    )

    if (is.null(d_result)) {
      next
    }

    d_result <- as.data.frame(d_result)

    d_column <- first_existing_column(
      d_result,
      c("d", "Cohens_d")
    )

    d_ci_low_column <- first_existing_column(
      d_result,
      c("CI_low", "CI_Low", "Conf.Low")
    )

    d_ci_high_column <- first_existing_column(
      d_result,
      c("CI_high", "CI_High", "Conf.High")
    )

    if (!is.null(d_column)) {
      output$d[i] <- suppressWarnings(
        as.numeric(d_result[[d_column]][1L])
      )
    }

    if (!is.null(d_ci_low_column)) {
      output$d_CI_low[i] <- suppressWarnings(
        as.numeric(d_result[[d_ci_low_column]][1L])
      )
    }

    if (!is.null(d_ci_high_column)) {
      output$d_CI_high[i] <- suppressWarnings(
        as.numeric(d_result[[d_ci_high_column]][1L])
      )
    }
  }

  # ==========================================================================
  # Final column order
  # ==========================================================================

  output <- output[
    ,
    c(
      "Parameter",
      "Estimate",
      "SE",
      "CI",
      "CI_low",
      "CI_high",
      "Statistic",
      "Statistic_type",
      "df",
      "p",
      "Interpreted_Estimate",
      "Interpreted_CI_low",
      "Interpreted_CI_high",
      "Estimate_type",
      "d",
      "d_CI_low",
      "d_CI_high"
    ),
    drop = FALSE
  ]

  rownames(output) <- NULL

  output
}
