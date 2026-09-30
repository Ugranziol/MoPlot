#' Summarize and Interpret Contrasts for One Categorical Predictor
#'
#' @description
#' ContrastsOnePred() returns observed or model-based estimates associated
#' with the contrasts used to code a single categorical predictor in a fitted
#' regression model.
#'
#' The function identifies the predictor's contrast coding and returns the
#' estimates for the levels or groups of levels involved in each contrast.
#'
#' The function supports standard linear and generalized linear models,
#' mixed-effects models fitted with lme4, mixed-effects models fitted with
#' nlme, and conditional models fitted with glmmTMB.
#'
#' @param model A fitted model object. Supported classes include objects
#' returned by [stats::lm()], [stats::glm()], [lme4::lmer()],
#' [lme4::glmer()], [nlme::lme()], and [glmmTMB::glmmTMB()].
#'
#' @param predictor Character string indicating the name of the categorical
#' predictor whose contrasts should be summarized. The corresponding
#' variable in data must be a factor.
#'
#' @param outcome Character string indicating the name of the dependent
#' variable. The corresponding variable in data must be numeric.
#'
#' @param data A data frame containing the variables used to fit `model`.
#'   The contrast matrix assigned to `data[[predictor]]` must correspond to
#'   the contrast matrix used when fitting the model.
#'
#' @param sig.level Numeric value indicating the significance threshold used
#' to classify contrasts as "sig" or "notsig". The default is .05.
#'
#' @param level Numeric value indicating the confidence level used for
#' confidence intervals. The default is .95.
#'
#' @param estimate_type Character string specifying how the level-specific
#'   estimates should be computed. Available options are:
#'
#'   \describe{
#'     \item{`"auto"`}{
#'       Uses observed means for models containing only the focal predictor
#'       and no interactions. Model-based predictions are used when the model
#'       contains additional predictors or interaction terms.
#'     }
#'     \item{`"observed"`}{
#'       Computes observed means directly from the supplied data.
#'     }
#'     \item{`"predicted"`}{
#'       Computes model-based predictions from the fixed effects and their
#'       variance-covariance matrix.
#'     }
#'   }
#'
#'   The default is `"auto"`.
#'
#' @param prediction_scale Character string indicating the scale used for
#' model-based predictions. Available options are:
#' \describe{
#' \item{"response"}{
#' Predictions are reported on the response scale. For generalized
#' models, the inverse link function is applied.
#' }
#' \item{"link"}{
#' Predictions are reported on the linear-predictor scale.
#' }
#' }
#' This argument is relevant only when estimate_type = "predicted".
#'
#' @details
#' The function retrieves the contrast matrix used for predictor and
#' identifies the levels contributing positively and negatively to each
#' contrast.
#'
#' Model-based predictions are computed manually from the fixed-effects model:
#'
#' \deqn{\hat{\eta} = X\hat{\beta}}
#'
#' Their variance on the linear-predictor scale is computed as:
#'
#' \deqn{\mathrm{Var}(\hat{\eta}) = X V_{\hat{\beta}} X^\mathsf{T}}
#'
#' where \eqn{X} is the fixed-effects model matrix and
#' \eqn{V_{\hat{\beta}}} is the variance-covariance matrix of the fixed-effect
#' coefficients.
#'
#' For estimate_type = "observed", estimates are observed means computed
#' directly from the data.
#'
#' For estimate_type = "predicted", estimates are computed from the
#' fixed-effect coefficients and their variance-covariance matrix.
#'
#' In mixed-effects models, random effects are excluded from model-based
#' predictions. Predicted estimates are therefore population-level or
#' fixed-effects predictions.
#'
#' Other numeric predictors are fixed at their observed means. Other
#' categorical predictors are fixed at their first factor level. Consequently,
#' model-based estimates are conditional predictions rather than averages over
#' all values of the remaining predictors.
#'
#' For generalized models with prediction_scale = "response", estimates are
#' transformed using the inverse link function. Standard errors on the response
#' scale are approximated using the delta method. Confidence intervals for
#' individual factor levels are first computed on the linear-predictor scale
#' and then transformed to the response scale.
#'
#' For non-polynomial contrasts involving more than one level on the positive
#' or negative side, the reported group estimate is the unweighted arithmetic
#' mean of the corresponding level-specific estimates.
#'
#' The Significance label refers to the p-value of the model coefficient
#' associated with the corresponding contrast. It does not represent a direct
#' test of the individual level-specific estimates or their confidence
#' intervals.
#'
#' For glmmTMB objects, only the conditional model component, cond, is
#' analyzed. Zero-inflation and dispersion components are not currently
#' summarized.
#'
#' @section Sum contrasts and the Grand Mean:
#' For sum-to-zero contrasts, each contrast is displayed as a comparison
#' between one factor level and a row labelled "Grand Mean".
#'
#' The Grand Mean is defined as the unweighted arithmetic mean of the
#' level-specific estimates. Each factor level therefore contributes equally,
#' regardless of its sample size.
#'
#' When estimate_type = "observed", the Grand Mean is:
#'
#' \deqn{
#' \frac{\bar{Y}_1 + \bar{Y}_2 + \cdots + \bar{Y}_k}{k}
#' }
#'
#' where \eqn{\bar{Y}_j} is the observed mean for level \eqn{j}.
#'
#' When estimate_type = "predicted" and
#' prediction_scale = "link", the Grand Mean is the arithmetic mean of the
#' level-specific linear predictors:
#'
#' \deqn{
#' \frac{\hat{\eta}_1 + \hat{\eta}_2 + \cdots + \hat{\eta}_k}{k}
#' }
#'
#' With sum contrasts, this quantity generally corresponds to the model
#' intercept when the remaining predictors are fixed at the values used to
#' construct the prediction data.
#'
#' When estimate_type = "predicted" and
#' prediction_scale = "response", the Grand Mean is the arithmetic mean of
#' the response-scale predictions:
#'
#' \deqn{
#' \frac{
#' g^{-1}(\hat{\eta}_1) +
#' g^{-1}(\hat{\eta}_2) +
#' \cdots +
#' g^{-1}(\hat{\eta}_k)
#' }{k}
#' }
#'
#' This is not generally equal to applying the inverse link function to the
#' average linear predictor:
#'
#' \deqn{
#' g^{-1}\left(\frac{1}{k}\sum_j \hat{\eta}_j\right)
#' \neq
#' \frac{1}{k}\sum_j g^{-1}(\hat{\eta}_j)
#' }
#'
#' For a binomial model with a logit link, the response-scale Grand Mean is the
#' arithmetic mean of the predicted probabilities.
#'
#' For a Poisson or negative-binomial model with a log link, the response-scale
#' Grand Mean is the arithmetic mean of the predicted expected counts.
#'
#' The standard error of an observed Grand Mean is computed by treating the
#' level-specific observed means as independent estimates.
#'
#' The standard error of a predicted Grand Mean is computed using the complete
#' variance-covariance matrix of the fixed-effect coefficients. On the response
#' scale, the required gradient is obtained using the multivariate delta method.
#'
#' @section Contrast types:
#' The function provides specific handling for:
#' \itemize{
#' \item Treatment contrasts generated by [stats::contr.treatment()].
#' \item Sum contrasts generated by [stats::contr.sum()].
#' \item Polynomial contrasts generated by [stats::contr.poly()].
#' \item Helmert contrasts generated by [stats::contr.helmert()].
#' \item Reverse Helmert contrasts generated by [reverse_helmert()].
#' \item Sliding-difference contrasts generated by [MASS::contr.sdif()].
#' \item Other user-defined contrast matrices.
#' }
#'
#' For treatment and sliding-difference contrasts, the coding matrix is
#' converted into the corresponding hypothesis-matrix representation using
#' the hypr package. This conversion makes the compared levels or groups of
#' levels directly interpretable.
#'
#' Polynomial contrasts are handled separately because all factor levels
#' contribute to each polynomial component through their polynomial weights.
#' For each component, the function returns one row per factor level and
#' reports both the polynomial weight and the fitted curve implied by that
#' component.
#'
#' @section P-values:
#' P-values are extracted from the coefficient table associated with the
#' fixed-effects or conditional component of the model.
#'
#' The function recognizes coefficient-table columns named:
#' \itemize{
#' \item "Pr(>|t|)";
#' \item "Pr(>|z|)";
#' \item "p-value".
#' }
#'
#' Only the main-effect coefficients associated with predictor are used.
#' Interaction coefficients involving predictor are excluded.
#'
#' For model classes that do not provide frequentist p-values, the
#' Significance column is set to NA.
#'
#' @return
#' A data frame with one or more rows per contrast.
#'
#' For non-polynomial contrasts, the output includes:
#'
#' \describe{
#' \item{Contrast}{
#' Name or sequential identifier of the contrast.
#' }
#' \item{Levels}{
#' Factor level or group of levels represented by the contrast.
#' }
#' \item{Estimate}{
#' Observed mean or model-based estimate.
#' }
#' \item{SE}{
#' Standard error of the reported estimate.
#' }
#' \item{CI_low}{
#' Lower confidence-interval boundary.
#' }
#' \item{CI_high}{
#' Upper confidence-interval boundary.
#' }
#' \item{Significance}{
#' Significance label associated with the corresponding model coefficient:
#' "sig", "notsig", or NA.
#' }
#' }
#'
#' For polynomial contrasts, the output additionally includes:
#'
#' \describe{
#' \item{Weight}{
#' Polynomial contrast weight associated with each factor level.
#' }
#' \item{Polynomial_fit}{
#' Fitted value implied by the intercept, the polynomial coefficient, and
#' the polynomial weight.
#' }
#' }
#'
#' The returned object also contains a "disclaimer" attribute describing how
#' estimates, standard errors, and confidence intervals were computed. The
#' same information is printed as a message when the function is executed.
#'
#' @examples
#' ###############################################################
#' ## Linear model with treatment contrasts
#' ###############################################################
#'
#' set.seed(123)
#'
#' dat <- data.frame(
#'   group = factor(rep(c("A", "B", "C"), each = 40))
#' )
#'
#' dat$y <- c(
#'   rnorm(40, mean = 10, sd = 2),
#'   rnorm(40, mean = 13, sd = 2),
#'   rnorm(40, mean = 17, sd = 2)
#' )
#'
#' contrasts(dat$group) <- contr.treatment(3)
#'
#' mod_lm <- lm(
#'   y ~ group,
#'   data = dat
#' )
#'
#' ContrastsOnePred(
#'   model = mod_lm,
#'   predictor = "group",
#'   outcome = "y",
#'   data = dat,
#'   estimate_type = "observed"
#' )
#'
#' ContrastsOnePred(
#'   model = mod_lm,
#'   predictor = "group",
#'   outcome = "y",
#'   data = dat,
#'   estimate_type = "predicted"
#' )
#'
#'
#' ###############################################################
#' ## Generalized linear model with sum contrasts
#' ###############################################################
#'
#' contrasts(dat$group) <- contr.sum(3)
#'
#' dat$binary_y <- rbinom(
#'   nrow(dat),
#'   size = 1,
#'   prob = plogis(
#'     -1 +
#'       0.5 * (dat$group == "B") +
#'       1.2 * (dat$group == "C")
#'   )
#' )
#'
#' mod_glm <- glm(
#'   binary_y ~ group,
#'   family = binomial,
#'   data = dat
#' )
#'
#' ContrastsOnePred(
#'   model = mod_glm,
#'   predictor = "group",
#'   outcome = "binary_y",
#'   data = dat,
#'   estimate_type = "predicted",
#'   prediction_scale = "response"
#' )
#'
#'
#' ###############################################################
#' ## Linear model with polynomial contrasts
#' ###############################################################
#'
#' dat$ordered_group <- ordered(
#'   dat$group,
#'   levels = c("A", "B", "C")
#' )
#'
#' contrasts(dat$ordered_group) <- contr.poly(3)
#'
#' mod_poly <- lm(
#'   y ~ ordered_group,
#'   data = dat
#' )
#'
#' ContrastsOnePred(
#'   model = mod_poly,
#'   predictor = "ordered_group",
#'   outcome = "y",
#'   data = dat,
#'   estimate_type = "predicted"
#' )
#'
#'
#' ###############################################################
#' ## Mixed-effects models
#' ###############################################################
#' \dontrun{
#'
#' library(lme4)
#' library(nlme)
#' library(glmmTMB)
#' library(MASS)
#'
#' set.seed(123)
#'
#' dat_mixed <- expand.grid(
#'   id = factor(1:60),
#'   time = 1:5
#' )
#'
#' dat_mixed$group <- factor(
#'   sample(
#'     c("A", "B", "C"),
#'     nrow(dat_mixed),
#'     replace = TRUE
#'   )
#' )
#'
#' random_intercept <- rnorm(
#'   nlevels(dat_mixed$id),
#'   mean = 0,
#'   sd = 1
#' )
#'
#' eta <- 10 +
#'   0.3 * dat_mixed$time +
#'   2 * (dat_mixed$group == "B") +
#'   5 * (dat_mixed$group == "C") +
#'   random_intercept[dat_mixed$id]
#'
#' dat_mixed$y <- eta +
#'   rnorm(
#'     nrow(dat_mixed),
#'     mean = 0,
#'     sd = 2
#'   )
#'
#'
#' ## Linear mixed model with Helmert contrasts
#'
#' contrasts(dat_mixed$group) <- contr.helmert(3)
#'
#' mod_lmer <- lmer(
#'   y ~ group + time + (1 | id),
#'   data = dat_mixed
#' )
#'
#' ContrastsOnePred(
#'   model = mod_lmer,
#'   predictor = "group",
#'   outcome = "y",
#'   data = dat_mixed,
#'   estimate_type = "predicted"
#' )
#'
#'
#' ## Linear mixed model fitted with nlme and sliding differences
#'
#' contrasts(dat_mixed$group) <- MASS::contr.sdif(3)
#'
#' mod_lme <- nlme::lme(
#'   y ~ group + time,
#'   random = ~1 | id,
#'   data = dat_mixed
#' )
#'
#' ContrastsOnePred(
#'   model = mod_lme,
#'   predictor = "group",
#'   outcome = "y",
#'   data = dat_mixed,
#'   estimate_type = "predicted"
#' )
#'
#'
#' ## Generalized mixed model with sum contrasts
#'
#' contrasts(dat_mixed$group) <- contr.sum(3)
#'
#' eta_count <- 1.5 +
#'   0.2 * dat_mixed$time +
#'   0.3 * (dat_mixed$group == "B") +
#'   0.7 * (dat_mixed$group == "C") +
#'   random_intercept[dat_mixed$id]
#'
#' dat_mixed$count <- rpois(
#'   nrow(dat_mixed),
#'   lambda = exp(eta_count)
#' )
#'
#' mod_tmb <- glmmTMB(
#'   count ~ group + time + (1 | id),
#'   family = poisson,
#'   data = dat_mixed
#' )
#'
#' ContrastsOnePred(
#'   model = mod_tmb,
#'   predictor = "group",
#'   outcome = "count",
#'   data = dat_mixed,
#'   estimate_type = "predicted",
#'   prediction_scale = "response"
#' )
#'
#' }
#' @export
ContrastsOnePred <- function(model = NULL,
                             predictor,
                             outcome,
                             data,
                             sig.level = .05,
                             level = .95,
                             estimate_type = c("auto", "observed", "predicted"),
                             prediction_scale = c("response", "link")) {

  # ============================================================
  # 1. ARGUMENT MATCHING AND BASIC VALIDATION
  # ============================================================

  estimate_type <- match.arg(estimate_type)
  prediction_scale <- match.arg(prediction_scale)

  if (is.null(model)) {
    stop("You must provide a fitted model.",call. = FALSE)
  }

  if (!.is_supported_model(model)) {
    stop(
      paste0(
        "Unsupported model class. Supported model classes include ",
        "lm, glm, lme4::lmer/glmer, lmerTest::lmer, nlme::lme, ",
        "and glmmTMB."
      ),
      call. = FALSE
    )
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.",call. = FALSE)
  }

  if (!is.character(predictor) ||
      length(predictor) != 1L ||
      is.na(predictor)) {

    stop("'predictor' must be a single character string.",call. = FALSE)
  }

  if (!is.character(outcome) ||
      length(outcome) != 1L ||
      is.na(outcome)) {

    stop("'outcome' must be a single character string.",call. = FALSE)
  }

  if (!predictor %in% names(data)) {
    stop("The predictor was not found in 'data'.",call. = FALSE)
  }

  if (!outcome %in% names(data)) {
    stop("The outcome was not found in 'data'.",call. = FALSE)
  }

  if (!is.factor(data[[predictor]])) {
    stop("The predictor must be represented by a factor variable.",call. = FALSE)
  }

  if (!is.numeric(data[[outcome]])) {
    stop("The dependent variable must be numeric.",call. = FALSE)
  }

  if (!is.numeric(sig.level) ||
      length(sig.level) != 1L ||
      is.na(sig.level) ||
      sig.level <= 0 ||
      sig.level >= 1) {

    stop("'sig.level' must be a single numeric value between 0 and 1.",call. = FALSE)
  }

  if (!is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1) {

    stop("'level' must be a single numeric value between 0 and 1.",call. = FALSE)
  }

  if (nlevels(data[[predictor]]) < 2L) {
    stop("The predictor must contain at least two factor levels.",call. = FALSE)
  }

  # ============================================================
  # 2. AUTOMATIC ESTIMATE-TYPE SELECTION
  # ============================================================

  if (identical(estimate_type, "auto")) {

    model_information <- ModelContrastTable(model,type = "summary")

    number_predictors <- as.integer(model_information$Value[model_information$Metric == "Predictors"])

    number_interactions <- as.integer(model_information$Value[model_information$Metric == "Interactions"])

    if (
      length(number_predictors) != 1L ||
      is.na(number_predictors)
    ) {
      stop(
        paste0(
          "The number of model predictors could not be determined ",
          "when `estimate_type = \"auto\"`."
        ),
        call. = FALSE
      )
    }

    if (
      number_predictors == 1L &&
      !is.na(number_interactions) &&
      number_interactions == 0L
    ) {
      estimate_type <- "observed"
    } else {
      estimate_type <- "predicted"
    }
  }


  z_value <- stats::qnorm(1 - (1 - level) / 2)

  # ============================================================
  # 2. INTERNAL FUNCTION: P-VALUES FOR THE FOCAL PREDICTOR
  # ============================================================
  #
  # The shared helper .extract_p_values() retrieves p-values by
  # coefficient name. Here, only the main-effect coefficients of
  # the focal predictor are selected. Interaction coefficients are
  # deliberately excluded.
  # ============================================================

  get_predictor_p_values <- function(model,
                                     predictor,
                                     n_contrasts) {

    fixed_effects <- .extract_fixef(model)
    coefficient_names <- names(fixed_effects)

    focal_indices <- startsWith(
      coefficient_names,
      predictor
    ) &
      !grepl(
        ":",
        coefficient_names,
        fixed = TRUE
      )

    focal_coefficient_names <- coefficient_names[focal_indices]

    # Ensure that the number of requested coefficient names matches
    # the number of columns in the predictor's contrast matrix.
    if (length(focal_coefficient_names) > n_contrasts) {

      focal_coefficient_names <- focal_coefficient_names[seq_len(n_contrasts)]
    }

    p_information <- .extract_p_values(
      model = model,
      coefficient_names = focal_coefficient_names,
      alpha = sig.level
    )

    p_values <- p_information$table$p_value

    if (length(p_values) < n_contrasts) {

      p_values <- c(
        p_values,
        rep(
          NA_real_,
          n_contrasts - length(p_values)
        )
      )
    }

    if (length(p_values) > n_contrasts) {

      p_values <- p_values[seq_len(n_contrasts)]
    }

    list(
      p = as.numeric(p_values),
      source = p_information$source
    )
  }

  # ============================================================
  # 3. INTERNAL FUNCTION: LEVEL-SPECIFIC ESTIMATES
  # ============================================================

  get_level_estimates <- function(model,
                                  predictor,
                                  outcome,
                                  data,
                                  estimate_type,
                                  prediction_scale,
                                  level) {

    level_names <- levels(data[[predictor]])

    n_levels <- length(level_names)

    z_value <- stats::qnorm(1 - (1 - level) / 2)

    # ----------------------------------------------------------
    # OBSERVED ESTIMATES
    # ----------------------------------------------------------

    if (estimate_type == "observed") {

      level_means <- tapply(
        data[[outcome]],
        data[[predictor]],
        mean,
        na.rm = TRUE
      )

      level_sds <- tapply(
        data[[outcome]],
        data[[predictor]],
        stats::sd,
        na.rm = TRUE
      )

      level_ns <- tapply(
        data[[outcome]],
        data[[predictor]],
        function(x) {
          sum(!is.na(x))
        }
      )

      level_ses <- as.numeric(level_sds / sqrt(level_ns))

      return(
        data.frame(
          Levels = level_names,
          Estimate = as.numeric(level_means),
          SE = level_ses,
          CI_low = as.numeric(level_means) -
            z_value * level_ses,
          CI_high = as.numeric(level_means) +
            z_value * level_ses,
          SD = as.numeric(level_sds),
          N = as.numeric(level_ns),
          stringsAsFactors = FALSE
        )
      )
    }

    # ----------------------------------------------------------
    # MODEL-BASED ESTIMATES
    # ----------------------------------------------------------
    #
    # Create one prediction row for each level of the focal
    # predictor. Other numeric predictors are fixed at their
    # observed means, whereas other factors are fixed at their
    # first level.
    #
    # Contrast matrices are not assigned manually here. They are
    # recovered from the fitted model and applied centrally by
    # .build_fixed_model_matrix().
    # ----------------------------------------------------------

    newdata <- data[rep(1L, n_levels), ,drop = FALSE]

    for (variable_name in names(newdata)) {

      if (is.factor(data[[variable_name]])) {

        newdata[[variable_name]] <- factor(
          rep(
            levels(data[[variable_name]])[1L],
            n_levels
          ),
          levels = levels(data[[variable_name]]),
          ordered = is.ordered(data[[variable_name]])
        )

      } else if (is.numeric(data[[variable_name]])) {

        newdata[[variable_name]] <- rep(
          mean(
            data[[variable_name]],
            na.rm = TRUE
          ),
          n_levels
        )
      }
    }

    # Replace the focal predictor with one row for each of its levels.
    newdata[[predictor]] <- factor(
      level_names,
      levels = level_names,
      ordered = is.ordered(data[[predictor]])
    )

    # The shared prediction helper guarantees the same fixed-effect
    # extraction, model-matrix construction, variance propagation,
    # inverse-link transformation, and delta-method logic used by
    # the other package functions.
    predictions <- .compute_fixed_predictions(
      model = model,
      newdata = newdata,
      prediction_scale = prediction_scale,
      conf_level = level
    )

    estimates <- data.frame(
      Levels = level_names,
      Estimate = predictions$estimate,
      SE = predictions$se,
      CI_low = predictions$ci_low,
      CI_high = predictions$ci_high,
      SD = NA_real_,
      N = NA_real_,
      stringsAsFactors = FALSE
    )

    # Preserve intermediate prediction quantities required for the
    # predicted Grand Mean in sum contrasts.
    attr(estimates, "model_matrix") <- attr(predictions, "model_matrix")

    attr(estimates, "vcov_beta") <- attr(predictions, "vcov_beta")

    attr(estimates, "eta") <- attr(predictions, "eta")

    attr(estimates, "family") <- attr(predictions, "family")

    estimates
  }

  # ============================================================
  # 4. INTERNAL FUNCTION: GRAND MEAN FOR SUM CONTRASTS
  # ============================================================

  get_sum_grand_mean <- function(level_estimates,
                                 estimate_type,
                                 prediction_scale,
                                 level) {

    estimate_values <- level_estimates$Estimate
    se_values <- level_estimates$SE

    z_value <- stats::qnorm(1 - (1 - level) / 2)

    # ----------------------------------------------------------
    # OBSERVED GRAND MEAN
    # ----------------------------------------------------------

    if (estimate_type == "observed") {

      valid_levels <- !is.na(estimate_values) & !is.na(se_values)

      number_valid_levels <- sum(valid_levels)

      if (number_valid_levels == 0L) {

        return(
          list(
            Estimate = NA_real_,
            SE = NA_real_,
            CI_low = NA_real_,
            CI_high = NA_real_
          )
        )
      }

      # Each factor level receives equal weight, irrespective of
      # the number of observations in that level.
      grand_mean <- mean(estimate_values[valid_levels])

      # Observed level means are based on disjoint groups and are
      # treated as independent for variance propagation.
      grand_mean_se <- sqrt(sum(se_values[valid_levels]^2)) / number_valid_levels

      return(
        list(
          Estimate = grand_mean,
          SE = grand_mean_se,
          CI_low = grand_mean -
            z_value * grand_mean_se,
          CI_high = grand_mean +
            z_value * grand_mean_se
        )
      )
    }

    # ----------------------------------------------------------
    # MODEL-BASED GRAND MEAN
    # ----------------------------------------------------------

    model_matrix <- attr(level_estimates,"model_matrix")

    variance_matrix <- attr(level_estimates,"vcov_beta")

    linear_predictors <- attr(level_estimates,"eta")

    family_object <- attr(level_estimates,"family")

    if (is.null(model_matrix) ||
        is.null(variance_matrix) ||
        is.null(linear_predictors)) {

      stop(
        paste0(
          "The information required to calculate the model-based ",
          "Grand Mean was not available."
        ),
        call. = FALSE
      )
    }

    # ----------------------------------------------------------
    # LINK-SCALE GRAND MEAN
    # ----------------------------------------------------------

    if (prediction_scale == "link" || is.null(family_object)) {

      mean_model_row <- matrix(colMeans(model_matrix),nrow = 1L)

      grand_mean <- mean(linear_predictors)

      grand_mean_variance <- as.numeric(
        mean_model_row %*%
          variance_matrix %*%
          t(mean_model_row)
      )

      grand_mean_se <- sqrt(max(grand_mean_variance,0))

      return(
        list(
          Estimate = grand_mean,
          SE = grand_mean_se,
          CI_low = grand_mean -
            z_value * grand_mean_se,
          CI_high = grand_mean +
            z_value * grand_mean_se
        )
      )
    }

    # ----------------------------------------------------------
    # RESPONSE-SCALE GRAND MEAN
    # ----------------------------------------------------------

    inverse_link <- family_object$linkinv
    derivative_inverse_link <- family_object$mu.eta

    response_values <- inverse_link(linear_predictors)

    grand_mean <- mean(response_values)

    response_derivatives <- as.numeric(derivative_inverse_link(linear_predictors))

    gradient <- colMeans(model_matrix * response_derivatives)

    gradient <- matrix(gradient,nrow = 1L)

    grand_mean_variance <- as.numeric( gradient %*% variance_matrix %*% t(gradient))

    grand_mean_se <- sqrt(max(grand_mean_variance,0))

    list(
      Estimate = grand_mean,
      SE = grand_mean_se,
      CI_low = grand_mean -
        z_value * grand_mean_se,
      CI_high = grand_mean +
        z_value * grand_mean_se
    )
  }

  # ============================================================
  # 5. INTERNAL FUNCTION: DISCLAIMER
  # ============================================================

  build_disclaimer <- function(estimate_type,
                               prediction_scale,
                               p_source,
                               is_mixed,
                               is_generalized,
                               is_polynomial,
                               is_sum) {

    disclaimer_text <- character()

    if (estimate_type == "observed") {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "Estimates are observed/raw means computed directly ",
          "from the data."
        )
      )

    } else {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "Estimates are model-based predictions computed from ",
          "the fixed effects and their variance-covariance matrix."
        )
      )

      if (is_mixed) {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "For mixed-effects models, random effects are not included; ",
            "predictions are population-level estimates."
          )
        )
      }

      if (is_generalized &&
          prediction_scale == "response") {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "For generalized models, predictions are reported on the ",
            "response scale. Individual confidence intervals are computed ",
            "on the linear-predictor scale and transformed using the ",
            "inverse link."
          )
        )
      }

      if (is_generalized &&
          prediction_scale == "link") {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "For generalized models, predictions are reported on the ",
            "link scale."
          )
        )
      }
    }

    if (is_sum) {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "For sum contrasts, Grand Mean is the unweighted arithmetic ",
          "mean of the displayed level-specific estimates on the ",
          "reported scale."
        )
      )

      if (estimate_type == "predicted" &&
          is_generalized &&
          prediction_scale == "response") {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "Its standard error is approximated using the multivariate ",
            "delta method."
          )
        )
      }
    }

    if (is_polynomial) {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "For polynomial contrasts, Weight reports the polynomial ",
          "weight and Polynomial_fit reports the curve implied by each ",
          "single polynomial contrast, computed as intercept plus the ",
          "contrast coefficient multiplied by the polynomial weight."
        )
      )
    }

    if (p_source == "lmerTest") {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "P-values for linear mixed-effects models were obtained ",
          "using lmerTest."
        )
      )
    }

    if (p_source == "not_available") {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "P-values were not available; Significance was set to NA."
        )
      )
    }

    paste(
      disclaimer_text,
      collapse = " "
    )
  }

  # ============================================================
  # 6. CONTRAST-MATRIX IDENTIFICATION
  # ============================================================

  original_contrast_matrix <- stats::contrasts(data[[predictor]])

  if (is.null(original_contrast_matrix)) {
    stop(
      "The predictor has no contrast matrix.",
      call. = FALSE
    )
  }

  contrast_matrix <- as.matrix(original_contrast_matrix)

  dimnames(contrast_matrix) <- NULL

  level_names <- levels(data[[predictor]])

  n_levels <- length(level_names)

  n_contrasts <- ncol(contrast_matrix)

  same_contrast_matrix <- function(first,
                                   second) {

    isTRUE(
      all.equal(
        unname(
          as.matrix(first)
        ),
        unname(
          as.matrix(second)
        ),
        check.attributes = FALSE
      )
    )
  }

  is_polynomial <- same_contrast_matrix(
    contrast_matrix,
    stats::contr.poly(n_levels)
  )

  # Treatment contrasts contain one positive weight and no negative
  # weight in coding-matrix form. Conversion to the corresponding
  # hypothesis matrix makes the reference-level comparison explicit.
  is_treatment <- same_contrast_matrix(
    contrast_matrix,
    stats::contr.treatment(n_levels)
  )

  is_sum <- same_contrast_matrix(
    contrast_matrix,
    stats::contr.sum(n_levels)
  )

  is_sliding_difference <- FALSE

  if (requireNamespace(
    "MASS",
    quietly = TRUE
  )) {

    is_sliding_difference <- same_contrast_matrix(
      contrast_matrix,
      MASS::contr.sdif(n_levels)
    )
  }

  if (is_treatment) {

    if (!requireNamespace(
      "hypr",
      quietly = TRUE
    )) {

      stop(
        paste0(
          "The 'hypr' package is required to interpret ",
          "treatment contrasts."
        ),
        call. = FALSE
      )
    }

    hypothesis_object <- hypr::hypr()

    hypr::cmat(hypothesis_object,add_intercept = TRUE) <- stats::contr.treatment(n_levels)

    contrast_matrix <- as.matrix(t(hypr::hmat(hypothesis_object))[, -1, drop = FALSE])
  }

  if (is_sliding_difference) {

    if (!requireNamespace(
      "hypr",
      quietly = TRUE
    )) {

      stop(
        paste0(
          "The 'hypr' package is required to interpret ",
          "sliding-difference contrasts."
        ),
        call. = FALSE
      )
    }

    hypothesis_object <- hypr::hypr()

    hypr::cmat(hypothesis_object,add_intercept = TRUE) <- original_contrast_matrix

    contrast_matrix <- as.matrix(t(hypr::hmat(hypothesis_object))[, -1, drop = FALSE])
  }

  # ============================================================
  # 7. MODEL INFORMATION
  # ============================================================

  fixed_effects <- .extract_fixef(model)

  intercept <- as.numeric(fixed_effects[1L])

  fixed_effect_names <- names(fixed_effects)

  focal_coefficient_indices <- startsWith(
    fixed_effect_names,
    predictor
  ) &
    !grepl(
      ":",
      fixed_effect_names,
      fixed = TRUE
    )

  contrast_coefficients <- fixed_effects[focal_coefficient_indices]

  p_information <- get_predictor_p_values(
    model = model,
    predictor = predictor,
    n_contrasts = n_contrasts
  )

  p_values <- p_information$p

  mixed_model <- .is_mixed_model(model)

  generalized_model <- .is_generalized_model(model)

  # ============================================================
  # 8. LEVEL-SPECIFIC ESTIMATES
  # ============================================================

  level_estimates <- get_level_estimates(
    model = model,
    predictor = predictor,
    outcome = outcome,
    data = data,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    level = level
  )

  level_estimate_values <- level_estimates$Estimate
  level_se_values <- level_estimates$SE
  level_ci_low_values <- level_estimates$CI_low
  level_ci_high_values <- level_estimates$CI_high

  # ============================================================
  # 9. SPECIAL CASE: POLYNOMIAL CONTRASTS
  # ============================================================

  if (is_polynomial) {

    if (length(contrast_coefficients) <
        ncol(contrast_matrix)) {

      stop(
        paste0(
          "The polynomial contrast coefficients could not be fully ",
          "identified in the model."
        ),
        call. = FALSE
      )
    }

    results <- data.frame(
      Contrast = character(),
      Levels = character(),
      Estimate = numeric(),
      SE = numeric(),
      CI_low = numeric(),
      CI_high = numeric(),
      Weight = numeric(),
      Polynomial_fit = numeric(),
      Significance = character(),
      stringsAsFactors = FALSE
    )

    polynomial_names <- colnames(stats::contr.poly(n_levels))

    for (contrast_index in seq_len(
      ncol(contrast_matrix)
    )) {

      polynomial_weights <- as.numeric(contrast_matrix[, contrast_index])

      polynomial_coefficient <- as.numeric(contrast_coefficients[contrast_index])

      polynomial_fit <- intercept +
        polynomial_coefficient *
        polynomial_weights

      significance_label <- if (
        is.na(p_values[contrast_index])
      ) {

        NA_character_

      } else if (
        p_values[contrast_index] < sig.level
      ) {

        "sig"

      } else {

        "notsig"
      }

      temporary_results <- data.frame(
        Contrast = paste0(
          "Contrast_",
          contrast_index,
          "_",
          polynomial_names[contrast_index]
        ),
        Levels = level_names,
        Estimate = level_estimate_values,
        SE = level_se_values,
        CI_low = level_ci_low_values,
        CI_high = level_ci_high_values,
        Weight = polynomial_weights,
        Polynomial_fit = polynomial_fit,
        Significance = significance_label,
        stringsAsFactors = FALSE
      )

      results <- rbind(
        results,
        temporary_results
      )
    }

    disclaimer <- build_disclaimer(
      estimate_type = estimate_type,
      prediction_scale = prediction_scale,
      p_source = p_information$source,
      is_mixed = mixed_model,
      is_generalized = generalized_model,
      is_polynomial = TRUE,
      is_sum = FALSE
    )

    attr(
      results,
      "disclaimer"
    ) <- disclaimer

    message(disclaimer)

    return(results)
  }

  # ============================================================
  # 10. GENERAL CASE: POSITIVE VS NEGATIVE CONTRAST GROUPS
  # ============================================================

  results <- data.frame(
    Contrast = character(),
    Levels = character(),
    Estimate = numeric(),
    SE = numeric(),
    CI_low = numeric(),
    CI_high = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
  )

  sum_grand_mean <- NULL

  if (is_sum) {

    sum_grand_mean <- get_sum_grand_mean(
      level_estimates = level_estimates,
      estimate_type = estimate_type,
      prediction_scale = prediction_scale,
      level = level
    )
  }

  for (contrast_index in seq_len(
    ncol(contrast_matrix)
  )) {

    contrast_weights <- as.numeric(contrast_matrix[, contrast_index])

    positive_levels <- which(contrast_weights > 0)

    negative_levels <- which(contrast_weights < 0)

    if (length(positive_levels) == 0L ||
        length(negative_levels) == 0L) {

      message(
        paste0(
          "Contrast_",
          contrast_index,
          ": no interpretable comparison."
        )
      )

      next
    }

    positive_group <- paste(
      level_names[positive_levels],
      collapse = ", "
    )

    negative_group <- paste(
      level_names[negative_levels],
      collapse = ", "
    )

    # Group estimates are unweighted arithmetic means of the
    # corresponding level-specific estimates.
    positive_estimate <- mean(
      level_estimate_values[positive_levels],
      na.rm = TRUE
    )

    negative_estimate <- mean(
      level_estimate_values[negative_levels],
      na.rm = TRUE
    )

    # These group-level SEs treat the displayed level estimates as
    # independent. The sum-contrast Grand Mean is handled separately
    # with the complete fixed-effect covariance matrix when predicted
    # estimates are requested.
    if (estimate_type == "predicted" && !is.null(attr(level_estimates, "model_matrix")) && !is.null(attr(level_estimates, "vcov_beta"))) {
      mm <- attr(level_estimates, "model_matrix")
      vcov_b <- attr(level_estimates, "vcov_beta")
      fam <- attr(level_estimates, "family")
      eta_val <- attr(level_estimates, "eta")

      # Positive side
      p_rows <- matrix(colMeans(mm[positive_levels, , drop = FALSE]), nrow = 1L)
      if (prediction_scale == "response" && !is.null(fam)) {
        p_rows <- matrix(colMeans(mm[positive_levels, , drop = FALSE] * fam$mu.eta(eta_val[positive_levels])), nrow = 1L)
      }
      positive_se <- sqrt(max(as.numeric(p_rows %*% vcov_b %*% t(p_rows)), 0))

      # Negative side
      n_rows <- matrix(colMeans(mm[negative_levels, , drop = FALSE]), nrow = 1L)
      if (prediction_scale == "response" && !is.null(fam)) {
        n_rows <- matrix(colMeans(mm[negative_levels, , drop = FALSE] * fam$mu.eta(eta_val[negative_levels])), nrow = 1L)
      }
      negative_se <- sqrt(max(as.numeric(n_rows %*% vcov_b %*% t(n_rows)), 0))
    } else {
      positive_se <- sqrt(sum(level_se_values[positive_levels]^2, na.rm = TRUE)) / length(positive_levels)
      negative_se <- sqrt(sum(level_se_values[negative_levels]^2, na.rm = TRUE)) / length(negative_levels)
    }

    positive_ci_low <- positive_estimate - z_value * positive_se

    positive_ci_high <- positive_estimate + z_value * positive_se

    negative_ci_low <- negative_estimate - z_value * negative_se

    negative_ci_high <- negative_estimate + z_value * negative_se

    # Sum contrasts are displayed as the focal level versus the
    # unweighted Grand Mean of all factor levels.
    if (is_sum) {

      negative_group <- "Grand Mean"

      negative_estimate <- sum_grand_mean$Estimate

      negative_se <- sum_grand_mean$SE

      negative_ci_low <- sum_grand_mean$CI_low

      negative_ci_high <- sum_grand_mean$CI_high
    }

    significance_label <- if (
      is.na(p_values[contrast_index])
    ) {

      NA_character_

    } else if (
      p_values[contrast_index] < sig.level
    ) {

      "sig"

    } else {

      "notsig"
    }

    temporary_results <- data.frame(
      Contrast = paste0(
        "Contrast_",
        contrast_index
      ),
      Levels = c(
        positive_group,
        negative_group
      ),
      Estimate = c(
        positive_estimate,
        negative_estimate
      ),
      SE = c(
        positive_se,
        negative_se
      ),
      CI_low = c(
        positive_ci_low,
        negative_ci_low
      ),
      CI_high = c(
        positive_ci_high,
        negative_ci_high
      ),
      Significance = significance_label,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    results <- rbind(
      results,
      temporary_results
    )
  }

  # ============================================================
  # 11. FINAL DISCLAIMER AND RETURN
  # ============================================================

  disclaimer <- build_disclaimer(
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    p_source = p_information$source,
    is_mixed = mixed_model,
    is_generalized = generalized_model,
    is_polynomial = FALSE,
    is_sum = is_sum
  )

  family_object <- tryCatch(
    stats::family(model),
    error = function(e) NULL
  )

  family_name <- if (!is.null(family_object)) {
    family_object$family
  } else {
    NA_character_
  }

  link_name <- if (!is.null(family_object)) {
    family_object$link
  } else {
    NA_character_
  }



  results$Display_Levels <- .format_display_levels(
    results$Levels
  )

  attr(results,"disclaimer") <- disclaimer
  attr(results, "family") <- family_name
  attr(results, "link") <- link_name
  attr(results, "estimate_type") <- estimate_type
  attr(results, "prediction_scale") <- prediction_scale

  message(disclaimer)

  results
}
