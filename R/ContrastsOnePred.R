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
#' The function retrieves the contrast (coding) matrix used for predictor and
#' derives from it the hypothesis matrix, which states the combination of
#' level-specific means estimated by each coefficient. The levels contributing
#' positively and negatively to each contrast, and their weights, are taken
#' from the hypothesis matrix (see the section "Coding matrix and hypothesis
#' matrix").
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
#' Confidence intervals of observed means respect the natural range of the
#' response. For proportions (binary outcomes or binomial models), the Wilson
#' score interval is used; for averages of several levels or cells it is
#' computed with the effective sample size \eqn{\hat{p}(1 - \hat{p}) /
#' \mathrm{SE}^2}, which reduces to the classical Wilson interval for a single
#' level. For positive responses analyzed with count or positive-continuous
#' families (for example Poisson, negative binomial, or Gamma), the interval
#' is computed on the log scale, \eqn{\exp(\log m \pm z\,\mathrm{SE}/m)}; a
#' mean of zero receives the exact Poisson upper bound
#' \eqn{-\log(\alpha/2)/n}. Other responses use symmetric intervals.
#'
#' For estimate_type = "predicted", estimates are computed from the
#' fixed-effect coefficients and their variance-covariance matrix.
#'
#' In mixed-effects models, random effects are excluded from model-based
#' predictions. Predicted estimates are therefore population-level or
#' fixed-effects predictions.
#'
#' Model-based estimates are computed on a reference grid. Other numeric
#' predictors are fixed at their observed means. Each level-specific estimate
#' averages over the levels of the other categorical predictors with equal
#' weights (a balanced average, independent of the observed cell sizes). A
#' categorical predictor that interacts with `predictor` is instead held at the
#' combination of its levels at which the coefficients of `predictor` are
#' defined, given by the first row of its hypothesis matrix: this is again the
#' balanced average for sum, Helmert, and sliding-difference coding, and the
#' reference level for treatment coding. On the linear-predictor scale, the
#' displayed comparisons therefore reconstruct the coefficients exactly; if
#' `predictor` interacts with a numeric predictor, its coefficients refer to the
#' value zero of that predictor and the function issues a warning.
#'
#' On the response scale, level estimates that average over other factors are
#' means of the response-scale predictions, and their confidence intervals are
#' obtained with the multivariate delta method.
#'
#' For generalized models with prediction_scale = "response", estimates are
#' transformed using the inverse link function. Standard errors on the response
#' scale are approximated using the delta method.
#'
#' On the response scale, the confidence interval of a prediction or of an
#' average of predictions \eqn{\hat{\theta}} is computed on the link scale and
#' transformed back,
#'
#' \deqn{g^{-1}\left(g(\hat{\theta}) \pm z \,
#'   \frac{\mathrm{SE}(\hat{\theta})}{\left|g^{-1\prime}(g(\hat{\theta}))\right|}\right),}
#'
#' so that its bounds remain within the range of the response (for example
#' between 0 and 1 for probabilities). For a single prediction this is the
#' interval computed on the linear predictor and transformed with the inverse
#' link. The interval is therefore generally asymmetric around the estimate,
#' whereas the reported standard error refers to the response scale.
#'
#' For non-polynomial contrasts involving more than one level on the positive
#' or negative side, the reported group estimate is the weighted mean of the
#' corresponding level-specific estimates, with weights taken from the
#' hypothesis matrix. For all standard codings these weights are equal, so the
#' group estimate is the unweighted arithmetic mean. Groups with unequal
#' weights, which arise with non-orthogonal custom codings, are labelled
#' `W(A .20, C .80)`.
#'
#' Standard errors of groups of observed means treat the level means as
#' independent. Standard errors of groups of model-based predictions use the
#' complete variance-covariance matrix of the fixed effects; on the response
#' scale they are obtained with the multivariate delta method, and the
#' confidence intervals are computed on the link scale as described above.
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
#' @section Coding matrix and hypothesis matrix:
#' The contrast matrix \eqn{C} describes how the factor is coded in the
#' model. It does not, in general, describe what the coefficients compare.
#' For a factor with \eqn{k} levels and \eqn{k - 1} contrasts, the
#' level-specific means satisfy
#'
#' \deqn{\mu = X_C \beta, \qquad X_C = [1, C],}
#'
#' so that
#'
#' \deqn{\beta = H \mu, \qquad H = X_C^{-1}.}
#'
#' Row \eqn{h_j} of \eqn{H} states the combination of level means estimated
#' by coefficient \eqn{\beta_j}. The function always derives the compared
#' levels from \eqn{H}. Reading them from the signs of \eqn{C} is correct only
#' for orthogonal codings; for non-orthogonal custom codings it can display a
#' comparison different from the one estimated by the coefficient. If
#' \eqn{X_C} is not full rank, the function stops.
#'
#' Every non-intercept row of \eqn{H} sums to zero and is written as
#'
#' \deqn{h_j = s_j (w^+_j - w^-_j),}
#'
#' where \eqn{w^+_j} and \eqn{w^-_j} are non-negative weights summing to one
#' (the two displayed groups) and \eqn{s_j > 0} is the scale factor returned
#' in `Scale_factor`. Consequently,
#'
#' \deqn{\beta_j = s_j \left(\sum_i w^+_{ji} \mu_i - \sum_i w^-_{ji} \mu_i\right).}
#'
#' On the linear-predictor scale this identity holds exactly for model-based
#' predictions of an additive model, whatever values the other predictors
#' are fixed at, because each row of \eqn{H} sums to zero. The function
#' verifies it and issues a warning when it fails, which usually indicates
#' that the predictor is involved in an interaction. For observed means, the
#' identity holds only when the model contains the focal predictor alone; the
#' disclaimer reports when it does not.
#'
#' For non-orthogonal custom codings, the function also compares the levels
#' suggested by the signs of each coding column with the levels involved in
#' the corresponding hypothesis. When they differ, a message describes the
#' comparison actually estimated by the coefficient.
#'
#' @section Contrast types:
#' The contrast type is identified from the coding matrix and, for matrices
#' that do not coincide with a standard generator, from the structure of the
#' hypothesis matrix. It determines only the presentation of the comparisons;
#' their content always comes from the hypothesis matrix. The function
#' recognizes:
#' \itemize{
#' \item Treatment contrasts, with any reference level (for example
#'   [stats::contr.treatment()] with any `base`, or [stats::contr.SAS()]):
#'   each level versus the reference level, \eqn{s = 1}.
#' \item Sum contrasts, with any omitted level (for example
#'   [stats::contr.sum()]): each level versus the Grand Mean, \eqn{s = 1}.
#' \item Helmert contrasts generated by [stats::contr.helmert()]: each level
#'   versus the mean of the preceding levels, \eqn{s_j = 1 / (j + 1)}.
#' \item Reverse Helmert contrasts generated by [reverse_helmert()]: each level
#'   versus the mean of the subsequent levels, \eqn{s_j = 1 / (k - j + 1)}.
#' \item Sliding-difference contrasts generated by [MASS::contr.sdif()]: each
#'   level versus the preceding level, \eqn{s = 1}.
#' \item Polynomial contrasts generated by [stats::contr.poly()].
#' \item Other user-defined (custom) contrast matrices, displayed as the two
#'   weighted groups of levels implied by the hypothesis matrix, with an
#'   arbitrary scale factor.
#' }
#'
#' Rescaled versions of these codings are recognized as well; rescaling a
#' coding column changes only the scale factor.
#'
#' Sum contrasts are displayed as one level versus the Grand Mean of all
#' levels, which includes the focal level. This presentation is equivalent to
#' the split of the hypothesis row into disjoint groups and is used because it
#' matches the usual interpretation of sum coding.
#'
#' Polynomial contrasts are handled separately because all factor levels
#' contribute to each polynomial component through their polynomial weights,
#' so that a split into positive and negative groups is not meaningful. For
#' each component, the function returns one row per factor level and reports
#' both the polynomial weight and the fitted curve implied by that component.
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
#' For non-polynomial contrasts, the output contains two rows per contrast,
#' the positive and the negative side, and includes:
#'
#' \describe{
#' \item{Contrast}{
#' Name or sequential identifier of the contrast.
#' }
#' \item{Levels}{
#' Factor level or group of levels represented by the contrast: a level
#' name, `"Mean(A, B)"` for equally weighted groups, `"W(A .20, C .80)"`
#' for groups with unequal weights, or `"Grand Mean"` for sum contrasts.
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
#' \item{Scale_factor}{
#' Scale factor \eqn{s_j} such that the model coefficient equals
#' \eqn{s_j} multiplied by the difference between the two displayed
#' quantities on the linear-predictor scale. It equals 1 for treatment, sum,
#' and sliding-difference contrasts.
#' }
#' \item{Display_Levels}{
#' Labels used for display; currently identical to `Levels`.
#' }
#' }
#'
#' For polynomial contrasts, the output contains one row per factor level for
#' each polynomial component, has no `Scale_factor` column, and additionally
#' includes:
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
#' The returned object also contains the following attributes:
#' \describe{
#' \item{disclaimer}{
#' Description of how estimates, standard errors, and confidence intervals
#' were computed, including any note on non-orthogonal custom codings. The
#' same information is printed as a message when the function is executed.
#' }
#' \item{contrast_type}{
#' The identified contrast type.
#' }
#' \item{hypothesis_matrix}{
#' The hypothesis matrix \eqn{H} without the intercept row, with one row per
#' contrast and one column per factor level.
#' }
#' \item{family, link, estimate_type, prediction_scale}{
#' Model family and link (NA when unavailable) and the settings used.
#' }
#' }
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

      # Intervals respect the natural range of the response, for
      # example [0, 1] for proportions (see .observed_interval()).
      level_intervals <- .observed_interval(
        estimate = as.numeric(level_means),
        se = level_ses,
        n = as.numeric(level_ns),
        range_type = response_range_type,
        conf_level = level
      )

      return(
        data.frame(
          Levels = level_names,
          Estimate = as.numeric(level_means),
          SE = level_ses,
          CI_low = level_intervals$ci_low,
          CI_high = level_intervals$ci_high,
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
    # Predictions are computed on a reference grid: numeric
    # predictors are fixed at their mean, and each level estimate
    # averages over the levels of the other factors (balanced, or at
    # the levels at which the focal coefficients are defined when a
    # factor interacts with the focal predictor). See
    # .reference_grid().
    #
    # Contrast matrices are not assigned manually here. They are
    # recovered from the fitted model and applied centrally by
    # .build_fixed_model_matrix().
    # ----------------------------------------------------------

    grid_information <- .reference_grid(
      model = model,
      data = data,
      focal_variables = predictor
    )

    predictions <- .compute_fixed_predictions(
      model = model,
      newdata = grid_information$grid,
      prediction_scale = prediction_scale,
      conf_level = level
    )

    # One row per level: weights of the grid rows forming that level.
    level_weights <- t(
      vapply(
        level_names,
        function(level_name) {
          .grid_cell_weights(
            grid_information,
            stats::setNames(list(level_name), predictor)
          )
        },
        numeric(nrow(grid_information$grid))
      )
    )

    level_summaries <- do.call(
      rbind,
      lapply(
        seq_len(n_levels),
        function(level_index) {
          .summarize_grid_predictions(
            predictions = predictions,
            weights = level_weights[level_index, ],
            prediction_scale = prediction_scale,
            conf_level = level
          )
        }
      )
    )

    estimates <- data.frame(
      Levels = level_names,
      Estimate = level_summaries$estimate,
      SE = level_summaries$se,
      CI_low = level_summaries$ci_low,
      CI_high = level_summaries$ci_high,
      SD = NA_real_,
      N = NA_real_,
      stringsAsFactors = FALSE
    )

    # Preserve the quantities needed to combine levels into groups
    # and to check the coefficient identity on the link scale.
    attr(estimates, "predictions") <- predictions

    attr(estimates, "level_weights") <- level_weights

    attr(estimates, "eta") <- as.numeric(level_weights %*% attr(predictions, "eta"))

    attr(estimates, "grid_information") <- grid_information

    estimates
  }

  # ============================================================
  # 4. INTERNAL FUNCTION: WEIGHTED GROUP OF LEVELS
  # ============================================================
  #
  # Each side of a contrast is a weighted mean of the level-specific
  # estimates, with weights taken from the hypothesis matrix and
  # summing to one. A single level reuses its own estimate and
  # confidence interval, so that intervals transformed from the link
  # scale are preserved.
  # ============================================================

  get_group_estimate <- function(level_estimates,
                                 weights,
                                 estimate_type,
                                 prediction_scale,
                                 level) {

    active_levels <- which(weights > 0)

    if (length(active_levels) == 1L) {

      return(
        list(
          Estimate = level_estimates$Estimate[active_levels],
          SE = level_estimates$SE[active_levels],
          CI_low = level_estimates$CI_low[active_levels],
          CI_high = level_estimates$CI_high[active_levels]
        )
      )
    }

    if (estimate_type == "observed") {

      # Observed level means are based on disjoint groups and are
      # treated as independent for variance propagation.
      group_estimate <- sum(weights * level_estimates$Estimate)

      group_se <- sqrt(sum(weights^2 * level_estimates$SE^2))

      group_interval <- .observed_interval(
        estimate = group_estimate,
        se = group_se,
        n = sum(level_estimates$N[weights > 0]),
        range_type = response_range_type,
        conf_level = level,
        single_cell = FALSE
      )

      return(
        list(
          Estimate = group_estimate,
          SE = group_se,
          CI_low = group_interval$ci_low,
          CI_high = group_interval$ci_high
        )
      )
    }

    # Model-based predictions are correlated; the shared helper uses
    # the complete fixed-effect covariance matrix and, on the
    # response scale, the multivariate delta method. The group weights
    # over levels are translated into weights over the grid rows.
    combined <- .summarize_grid_predictions(
      predictions = attr(level_estimates, "predictions"),
      weights = as.numeric(crossprod(attr(level_estimates, "level_weights"), weights)),
      prediction_scale = prediction_scale,
      conf_level = level
    )

    list(
      Estimate = combined$estimate,
      SE = combined$se,
      CI_low = combined$ci_low,
      CI_high = combined$ci_high
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
                               contrast_type,
                               scale_factors,
                               identity_holds,
                               divergence_notes,
                               range_type) {

    disclaimer_text <- character()

    if (estimate_type == "observed") {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "Estimates are observed/raw means computed directly ",
          "from the data."
        ),
        switch(
          range_type,
          "proportion" = paste0(
            "Confidence intervals of proportions are Wilson score intervals, ",
            "which remain within [0, 1]."
          ),
          "positive" = paste0(
            "Confidence intervals of positive means are computed on the log ",
            "scale, so that they remain positive."
          ),
          NULL
        )
      )

    } else {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "Estimates are model-based predictions computed from ",
          "the fixed effects and their variance-covariance matrix, on a ",
          "reference grid in which numeric predictors are fixed at their ",
          "mean and the other factors are averaged with equal weights, or ",
          "held at the levels at which the coefficients are defined when ",
          "they interact with the predictor."
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
            "response scale. Confidence intervals are computed on the ",
            "link scale and transformed using the inverse link, so that ",
            "they remain within the range of the response; standard errors ",
            "of averages of predictions use the multivariate delta method."
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

    if (identical(contrast_type, "Polynomial")) {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "For polynomial contrasts, Weight reports the polynomial ",
          "weight and Polynomial_fit reports the curve implied by each ",
          "single polynomial contrast, computed as intercept plus the ",
          "contrast coefficient multiplied by the polynomial weight."
        )
      )

    } else {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "Each contrast compares the weighted groups of levels implied ",
          "by the hypothesis matrix of the coding (the inverse of the ",
          "coding matrix with the intercept column); groups of levels ",
          "are means of the displayed level-specific estimates."
        )
      )

      if (identical(contrast_type, "Sum")) {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "For sum contrasts, Grand Mean is the unweighted arithmetic ",
            "mean of the displayed level-specific estimates on the ",
            "reported scale."
          )
        )
      }

      if (any(abs(scale_factors - 1) > 1e-8)) {

        disclaimer_text <- c(
          disclaimer_text,
          paste0(
            "Each model coefficient equals Scale_factor multiplied by the ",
            "difference between the two displayed quantities on the ",
            "linear-predictor scale; for this coding, Scale_factor differs ",
            "from 1 for at least one contrast."
          )
        )
      }
    }

    if (isFALSE(identity_holds) &&
        estimate_type == "observed") {

      disclaimer_text <- c(
        disclaimer_text,
        paste0(
          "The observed means do not reconstruct the model coefficients ",
          "exactly, because the model contains other terms, random ",
          "effects, or levels with all-zero or all-one outcomes; use ",
          "estimate_type = \"predicted\" for displayed differences that ",
          "correspond to the coefficients."
        )
      )
    }

    disclaimer_text <- c(
      disclaimer_text,
      divergence_notes
    )

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
  # 6. CONTRAST HYPOTHESES
  # ============================================================
  #
  # The coding matrix describes how the factor enters the model;
  # the hypothesis matrix describes what each coefficient compares.
  # The latter is the source of truth for the decomposition. The
  # contrast type only selects a readable presentation.
  # ============================================================

  coding_matrix <- stats::contrasts(data[[predictor]])

  if (is.null(coding_matrix)) {
    stop(
      "The predictor has no contrast matrix.",
      call. = FALSE
    )
  }

  level_names <- levels(data[[predictor]])

  hypotheses <- .contrast_hypotheses(
    contrast_matrix = coding_matrix,
    level_names = level_names
  )

  contrast_type <- hypotheses$type

  n_contrasts <- ncol(coding_matrix)

  # ============================================================
  # 7. MODEL INFORMATION
  # ============================================================

  fixed_effects <- .extract_fixef(model)

  focal_coefficient_names <- .focal_coefficient_names(
    model = model,
    predictor = predictor
  )

  if (length(focal_coefficient_names) != n_contrasts) {

    stop(
      paste0(
        "The model contains ",
        length(focal_coefficient_names),
        " coefficients for '",
        predictor,
        "', but the contrast matrix in 'data' has ",
        n_contrasts,
        " columns. The contrasts assigned to data[[predictor]] must be ",
        "those used to fit the model."
      ),
      call. = FALSE
    )
  }

  contrast_coefficients <- fixed_effects[focal_coefficient_names]

  p_information <- .extract_p_values(
    model = model,
    coefficient_names = focal_coefficient_names,
    alpha = sig.level
  )

  p_values <- p_information$table$p_value

  significance_labels <- ifelse(
    is.na(p_values),
    NA_character_,
    ifelse(
      p_values < sig.level,
      "sig",
      "notsig"
    )
  )

  mixed_model <- .is_mixed_model(model)

  generalized_model <- .is_generalized_model(model)

  family_object <- .extract_family(model)

  # Natural range of the response, used for the intervals of observed
  # means.
  response_range_type <- .response_range_type(
    family_object = family_object,
    values = data[[outcome]]
  )

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

  # ============================================================
  # 9. IDENTITY BETWEEN DISPLAYED COMPARISONS AND COEFFICIENTS
  # ============================================================
  #
  # On the linear-predictor scale, each coefficient must equal its
  # scale factor multiplied by the displayed difference. For
  # model-based predictions of an additive model this holds exactly,
  # whatever values the other predictors are fixed at. A failure
  # indicates, for example, that the focal predictor is involved in
  # an interaction. For observed means the identity holds only when
  # the model contains the focal predictor alone.
  # ============================================================

  if (estimate_type == "predicted") {

    link_level_values <- attr(level_estimates, "eta")

  } else if (!is.null(family_object) &&
             generalized_model) {

    link_level_values <- suppressWarnings(
      family_object$linkfun(level_estimates$Estimate)
    )

  } else {

    link_level_values <- level_estimates$Estimate
  }

  identity_holds <- all(
    is.finite(link_level_values)
  ) &&
    all(
      .check_hypothesis_identity(
        hypotheses = hypotheses,
        coefficients = contrast_coefficients,
        level_values = link_level_values
      )
    )

  if (!identity_holds &&
      estimate_type == "predicted") {

    warning(
      paste0(
        "The displayed comparisons do not reconstruct the model ",
        "coefficients of '",
        predictor,
        "'. This usually means that the predictor interacts with a numeric ",
        "predictor: its coefficients refer to the value zero of that ",
        "predictor, whereas the predictions fix it at its mean."
      ),
      call. = FALSE
    )
  }

  # ============================================================
  # 10. SPECIAL CASE: POLYNOMIAL CONTRASTS
  # ============================================================

  if (identical(contrast_type, "Polynomial")) {

    intercept <- as.numeric(fixed_effects[1L])

    polynomial_names <- colnames(stats::contr.poly(length(level_names)))

    results <- do.call(
      rbind,
      lapply(
        seq_len(n_contrasts),
        function(contrast_index) {

          polynomial_weights <- as.numeric(coding_matrix[, contrast_index])

          data.frame(
            Contrast = paste0(
              "Contrast_",
              contrast_index,
              "_",
              polynomial_names[contrast_index]
            ),
            Levels = level_names,
            Estimate = level_estimates$Estimate,
            SE = level_estimates$SE,
            CI_low = level_estimates$CI_low,
            CI_high = level_estimates$CI_high,
            Weight = polynomial_weights,
            Polynomial_fit = intercept +
              as.numeric(contrast_coefficients[contrast_index]) *
              polynomial_weights,
            Significance = significance_labels[contrast_index],
            stringsAsFactors = FALSE
          )
        }
      )
    )

    scale_factors <- NA_real_
    divergence_notes <- character()

  } else {

    # ============================================================
    # 11. GENERAL CASE: WEIGHTED POSITIVE VS NEGATIVE GROUPS
    # ============================================================

    # The same group (for example the Grand Mean of sum contrasts)
    # may appear in several contrasts; it is computed only once.
    group_cache <- list()

    cached_group_estimate <- function(label,
                                      weights) {

      if (is.null(group_cache[[label]])) {

        group_cache[[label]] <<- get_group_estimate(
          level_estimates = level_estimates,
          weights = weights,
          estimate_type = estimate_type,
          prediction_scale = prediction_scale,
          level = level
        )
      }

      group_cache[[label]]
    }

    results <- do.call(
      rbind,
      lapply(
        seq_len(n_contrasts),
        function(contrast_index) {

          contrast <- hypotheses$contrasts[[contrast_index]]

          positive <- cached_group_estimate(
            contrast$positive_label,
            contrast$positive_weights
          )

          negative <- cached_group_estimate(
            contrast$negative_label,
            contrast$negative_weights
          )

          data.frame(
            Contrast = paste0("Contrast_", contrast_index),
            Levels = c(contrast$positive_label, contrast$negative_label),
            Estimate = c(positive$Estimate, negative$Estimate),
            SE = c(positive$SE, negative$SE),
            CI_low = c(positive$CI_low, negative$CI_low),
            CI_high = c(positive$CI_high, negative$CI_high),
            Significance = significance_labels[contrast_index],
            Scale_factor = contrast$scale,
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
      )
    )

    scale_factors <- vapply(
      hypotheses$contrasts,
      `[[`,
      numeric(1),
      "scale"
    )

    divergence_notes <- stats::na.omit(
      vapply(
        hypotheses$contrasts,
        `[[`,
        character(1),
        "note"
      )
    )

    for (note in divergence_notes) {
      message(note)
    }
  }

  rownames(results) <- NULL

  # ============================================================
  # 12. FINAL DISCLAIMER AND RETURN
  # ============================================================

  disclaimer <- build_disclaimer(
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    p_source = p_information$source,
    is_mixed = mixed_model,
    is_generalized = generalized_model,
    contrast_type = contrast_type,
    scale_factors = scale_factors,
    identity_holds = identity_holds,
    divergence_notes = as.character(divergence_notes),
    range_type = response_range_type
  )

  # Levels holds structural labels: a factor level, "Mean(A, B)",
  # "W(A .20, C .80)", or "Grand Mean". They are already display-ready.
  results$Display_Levels <- results$Levels

  attr(results, "disclaimer") <- disclaimer
  attr(results, "family") <- if (is.null(family_object)) NA_character_ else family_object$family
  attr(results, "link") <- if (is.null(family_object)) NA_character_ else family_object$link
  attr(results, "estimate_type") <- estimate_type
  attr(results, "prediction_scale") <- prediction_scale
  attr(results, "contrast_type") <- contrast_type
  attr(results, "hypothesis_matrix") <- hypotheses$hypothesis_matrix

  message(disclaimer)

  results
}
