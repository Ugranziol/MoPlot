#' Decompose Interaction Contrasts into Interpretable Configurations
#'
#' @description
#' `ContrastMultPred()` decomposes model coefficients associated with
#' interactions among categorical predictors into interpretable
#' configurations of factor levels.
#'
#' For each selected interaction coefficient, the function identifies the
#' positive and negative level groups implied by the contrast matrices used
#' to fit the model. It then returns observed/raw means or model-based
#' predictions for the resulting configurations.
#'
#' One variable involved in each interaction is selected as the
#' \emph{split variable}. The contrasts of the remaining target variable or
#' variables are then displayed within the levels or contrast-defined groups
#' of the split variable.
#'
#' The function supports two-way and higher-order interactions.
#'
#' @param model A fitted model object. Supported model classes include objects
#'   returned by [stats::lm()], [stats::glm()], [lme4::lmer()],
#'   [lme4::glmer()], [nlme::lme()], and [glmmTMB::glmmTMB()].
#'
#' @param data Optional data frame containing the variables used to fit
#'   `model`. If `NULL`, the function attempts to recover the analysis data
#'   from the fitted model using `stats::model.frame()`. Supplying `data`
#'   explicitly remains recommended when the model does not retain a usable
#'   model frame or when a specific analysis data set should be used.
#'
#' @param response Character string indicating the name of the response
#'   variable. By default, the response is extracted from the model formula.
#'
#' @param interaction Optional character string identifying the interaction
#'   term to decompose, for example `"A:B"` or `"A:B:C"`. When supplied,
#'   the interaction is resolved using the fitted model matrix and the
#'   coefficients associated with that interaction are selected automatically.
#'   If `coefficients` is also supplied, those coefficients are validated
#'   against the requested interaction. When both `interaction` and
#'   `coefficients` are `NULL`, all non-intercept coefficients containing
#'   `":"` are analyzed, preserving the default behavior of this function.
#'
#' @param coefficients Optional character vector containing the interaction
#'   coefficients to decompose. If `interaction` is `NULL` and `coefficients`
#'   is `NULL`, all non-intercept coefficients containing `":"` are analyzed.
#'   When `interaction` is supplied, `coefficients` may be used to select a
#'   subset of the coefficients belonging to that interaction.
#'
#' @param split_rule Character string indicating how the split variable is
#'   selected for each interaction coefficient. Available options are:
#'   \describe{
#'     \item{`"last_variable"`}{
#'       Selects the last variable appearing in the interaction coefficient.
#'     }
#'     \item{`"first_variable"`}{
#'       Selects the first variable appearing in the interaction coefficient.
#'     }
#'     \item{`"lowest_complexity"`}{
#'       Selects the variable with the lowest contrast-complexity score.
#'     }
#'     \item{`"highest_complexity"`}{
#'       Selects the variable with the highest contrast-complexity score.
#'     }
#'   }
#'
#' @param top_split Optional character string identifying a variable that must
#'   be used as the split variable. When supplied, `top_split`
#'   overrides `split_rule`.
#'
#' @param conf_level Numeric value indicating the confidence level used for
#'   confidence intervals. The default is `0.95`.
#'
#' @param alpha_sig Numeric value indicating the alpha threshold used to
#'   classify the original interaction coefficients as significant or
#'   non-significant. The default is `0.05`.
#'
#' @param configuration Character string controlling which configuration
#'   labels are included in the output. Available options are:
#'   \describe{
#'     \item{`"none"`}{
#'       Does not return configuration-label columns.
#'     }
#'     \item{`"all"`}{
#'       Returns both explicit and semantic configuration labels.
#'     }
#'     \item{`"explicit"`}{
#'       Returns the complete mathematical composition of each configuration.
#'     }
#'     \item{`"semantic"`}{
#'       Returns a compact description of the levels or groups represented by
#'       each configuration.
#'     }
#'   }
#'
#' @param estimate_type Character string indicating which estimates are
#'   returned. Available options are:
#'   \describe{
#'     \item{`"raw"`}{
#'       Returns observed/raw aggregated means computed directly from `data`.
#'     }
#'     \item{`"predicted"`}{
#'       Returns model-based predictions computed from the fixed-effect
#'       coefficients and their variance-covariance matrix.
#'     }
#'   }
#'
#' @param prediction_scale Character string indicating the scale on which
#'   model-based predictions are returned. Available options are:
#'   \describe{
#'     \item{`"response"`}{
#'       Returns predictions on the response scale. For generalized models,
#'       the inverse link function is applied.
#'     }
#'     \item{`"link"`}{
#'       Returns predictions on the linear-predictor scale.
#'     }
#'   }
#'   This argument is ignored when `estimate_type = "raw"`.
#'
#' @param methodological_note Logical. Should a methodological note describing
#'   the settings relevant to the current analysis be printed? The default is
#'   `TRUE`. The note is stored as an attribute of the returned object
#'   regardless of this setting.
#'
#' @details
#' Each interaction coefficient is separated into its contrast-coded
#' components. The function then identifies:
#'
#' \itemize{
#'   \item the categorical predictor represented by each component;
#'   \item the corresponding contrast column;
#'   \item the factor levels receiving positive, negative, or zero weights.
#' }
#'
#' One predictor is selected as the split variable. The remaining
#' predictor or predictors are treated as target variables whose contrasts are
#' interpreted within the levels or contrast-defined groups of the
#' split variable.
#'
#' For example, suppose that `FactorA` has three levels and uses Helmert
#' contrasts, whereas `FactorB` has two levels and uses sliding-difference
#' contrasts. If `FactorB` is selected as the split variable, the
#' function displays the Helmert contrasts of `FactorA` separately within the
#' levels or contrast-defined groups of `FactorB`.
#'
#' Therefore, the split variable is not necessarily the variable of primary
#' interpretive interest. It defines the conditions under which the
#' contrasts of the remaining target variable or variables are displayed.
#'
#' The split variable is selected separately for each interaction coefficient.
#' Consequently, when all model interactions are analyzed simultaneously,
#' more than one value may appear in the `split_var` column. For example, in a
#' model containing `A:B`, `A:C`, `B:C`, and `A:B:C`, the variable selected by
#' `"highest_complexity"` may differ across coefficients because only the
#' variables involved in each coefficient are considered.
#'
#' Random-effect grouping variables are excluded when identifying the
#' categorical predictors involved in the fixed-effect interactions.
#'
#' @section Contrast complexity and split selection:
#' The options `"lowest_complexity"` and `"highest_complexity"` select the
#' split variable according to the contrast coding assigned to each
#' predictor.
#'
#' Contrast systems are ranked using the following approximate scores:
#'
#' \itemize{
#'   \item Treatment contrasts: `0`.
#'   \item Sum contrasts: `1`.
#'   \item Sliding-difference contrasts: `1`.
#'   \item Helmert contrasts: `2`.
#'   \item Reverse Helmert contrasts: `2`.
#'   \item Polynomial contrasts: `3`.
#'   \item Custom or unidentified contrast matrices: `4`.
#' }
#'
#' These scores affect only the choice of split variable. They do not
#' determine the order in which variables are displayed in the configuration
#' labels.
#'
#' If several variables have the same complexity score, the rightmost tied
#' variable in the interaction coefficient is selected. A warning is issued
#' to document the choice.
#'
#' The argument `top_split` bypasses the automatic rule and forces the
#' specified variable to be used as the split variable.
#'
#' @section Ordering of configuration labels:
#' Both explicit and semantic configuration labels use the following order:
#'
#' \preformatted{
#' target variable(s) -> split variable
#' }
#'
#' The target variables are displayed first, following their stable order
#' among the categorical predictors in the fixed-effects model. The
#' split variable selected by `split_rule` or `top_split` is always
#' displayed last.
#'
#' For a two-way interaction between `A` and `B`, if `A` is the split
#' variable, labels follow the order:
#'
#' \preformatted{
#' B_A
#' }
#'
#' If `B` is the split variable, labels follow the order:
#'
#' \preformatted{
#' A_B
#' }
#'
#' For a three-way interaction among `A`, `B`, and `C`, if `C` is the
#' split variable, labels follow the order:
#'
#' \preformatted{
#' A_B_C
#' }
#'
#' If `A` is the split variable, labels follow the order:
#'
#' \preformatted{
#' B_C_A
#' }
#'
#' This ordering rule affects only the presentation and interpretation of the
#' labels. It does not modify the contrast-complexity scores, the selected
#' split variable, or the numerical estimates.
#'
#' @section Explicit and semantic configurations:
#' `configuration_explicit` reports the individual cells entering an aggregate.
#' It is intended to show the exact mathematical composition of the displayed
#' estimate.
#'
#' For example:
#'
#' \preformatted{
#' (B2_A1 + B2_A2)/2
#' }
#'
#' indicates that the estimate is the arithmetic mean of cells `B2_A1` and
#' `B2_A2`.
#'
#' `configuration_semantic` reports the corresponding interpretation:
#'
#' \preformatted{
#' B2_Mean(A1,A2)
#' }
#'
#' Both labels use the same variable order. In this example, `B` is the target
#' variable and `A` is the split variable.
#'
#' With `B` as the split variable, the corresponding labels could be:
#'
#' \preformatted{
#' (A1_B2 + A2_B2)/2
#' Mean(A1,A2)_B2
#' }
#'
#' The explicit label answers the question, "Which cells were combined?",
#' whereas the semantic label answers the question, "What does the resulting
#' estimate represent?"
#'
#' @section Raw estimates:
#' For `estimate_type = "raw"`, estimates are computed directly from the
#' observed response values.
#'
#' At the individual-cell level, the function returns the observed mean,
#' standard error, and confidence interval.
#'
#' When several disjoint cells are combined, their means are aggregated using
#' normalized contrast weights. Standard errors are combined through classical
#' variance propagation.
#'
#' Confidence intervals for aggregated raw estimates use a
#' Satterthwaite-like approximation to the degrees of freedom when sufficient
#' information is available.
#'
#' @section Model-based predictions:
#' For `estimate_type = "predicted"`, predictions are computed from the
#' fixed-effect model:
#'
#' \deqn{
#'   \hat{\eta} = X\hat{\beta}
#' }
#'
#' where \eqn{X} is the fixed-effect model matrix and \eqn{\hat{\beta}} is the
#' vector of fixed-effect coefficients.
#'
#' The covariance matrix of the linear predictors is obtained from:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\eta})
#'   =
#'   X V_{\hat{\beta}} X^\mathsf{T}
#' }
#'
#' where \eqn{V_{\hat{\beta}}} is the complete fixed-effect
#' variance-covariance matrix.
#'
#' In mixed-effects models, random effects are excluded. The returned
#' predictions are therefore population-level or fixed-effect predictions.
#'
#' For `glmmTMB` models, only the conditional component is analyzed.
#' Zero-inflation and dispersion components are not currently decomposed.
#'
#' @section Response-scale predictions:
#' For generalized models with `prediction_scale = "response"`, each
#' linear-predictor value is first transformed through the inverse link
#' function. The transformed predictions are then averaged.
#'
#' Therefore, the function computes:
#'
#' \deqn{
#'   \frac{1}{n}
#'   \sum_{i=1}^{n}
#'   g^{-1}(\hat{\eta}_i)
#' }
#'
#' rather than:
#'
#' \deqn{
#'   g^{-1}\left(
#'     \frac{1}{n}
#'     \sum_{i=1}^{n}
#'     \hat{\eta}_i
#'   \right)
#' }
#'
#' These two quantities are generally different when the inverse link is
#' nonlinear.
#'
#' @section Aggregation of predicted estimates:
#' Model-based predictions for different configurations are generally
#' correlated because they depend on the same fixed-effect coefficient vector.
#'
#' `ContrastMultPred()` therefore does not combine predicted standard errors
#' under an assumption of independence.
#'
#' Instead, the prediction rows underlying each child configuration and their
#' internal weights are retained throughout the recursive decomposition.
#' When configurations are aggregated, all underlying prediction rows are
#' recombined and the aggregate is calculated directly from the complete
#' fixed-effect variance-covariance matrix.
#'
#' On the link scale, an aggregated prediction with normalized weights
#' \eqn{a_i} is:
#'
#' \deqn{
#'   \hat{\theta}
#'   =
#'   \sum_i a_i X_i\hat{\beta}
#' }
#'
#' Defining:
#'
#' \deqn{
#'   L = \sum_i a_i X_i
#' }
#'
#' its variance is:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\theta})
#'   =
#'   L V_{\hat{\beta}} L^\mathsf{T}
#' }
#'
#' On the response scale, the aggregated prediction is:
#'
#' \deqn{
#'   \hat{\theta}
#'   =
#'   \sum_i
#'   a_i
#'   g^{-1}(X_i\hat{\beta})
#' }
#'
#' Its variance is approximated using the multivariate delta method. The
#' gradient is:
#'
#' \deqn{
#'   \nabla\hat{\theta}
#'   =
#'   \sum_i
#'   a_i
#'   g^{-1\prime}(\hat{\eta}_i)
#'   X_i
#' }
#'
#' and:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\theta})
#'   \approx
#'   \nabla\hat{\theta}
#'   V_{\hat{\beta}}
#'   \nabla\hat{\theta}^{\mathsf{T}}
#' }
#'
#' This procedure preserves the covariance among the predictions entering
#' each aggregate.
#'
#' @section Supported contrast types:
#' The function supports:
#'
#' \itemize{
#'   \item Treatment contrasts generated by [stats::contr.treatment()].
#'   \item Sum contrasts generated by [stats::contr.sum()].
#'   \item Helmert contrasts generated by [stats::contr.helmert()].
#'   \item Reverse Helmert contrasts generated by [reverse_helmert()].
#'   \item Sliding-difference contrasts generated by [MASS::contr.sdif()].
#'   \item Custom contrast matrices that can be interpreted using positive,
#'     negative, and optionally zero weights.
#' }
#'
#' For treatment and sliding-difference contrasts, the coding matrix is
#' converted to its corresponding hypothesis-matrix representation using the
#' `hypr` package. This allows the positive and negative groups to be
#' interpreted directly.
#'
#' Sum contrasts are interpreted as a focal level or group relative to the
#' unweighted mean across all levels of the corresponding factor.
#'
#' Polynomial contrasts are not currently supported for interaction
#' decomposition because all levels generally contribute to each polynomial
#' component through nonzero weights and require a different interpretive
#' strategy.
#'
#' @section P-values:
#' P-values refer to the original interaction coefficients being decomposed.
#' They are not direct tests of the individual configuration-specific
#' estimates displayed in the output.
#'
#' The function recognizes coefficient-table columns named:
#'
#' \itemize{
#'   \item `"Pr(>|t|)"`;
#'   \item `"Pr(>|z|)"`;
#'   \item `"p-value"`.
#' }
#'
#' Models fitted with [lme4::lmer()] do not ordinarily provide p-values.
#' When `lmerTest` is installed, the model is converted internally using
#' [lmerTest::as_lmerModLmerTest()], and approximate p-values are extracted
#' from the resulting coefficient table.
#'
#' If no p-values are available, `p_value` is set to `NA` and
#' `significance` is set to `"NoPValue"`.
#'
#' @return
#' A data frame of class `"plot_ready_contrasts_df"` with one row per
#' decomposed configuration.
#'
#' The returned columns include:
#'
#' \describe{
#'   \item{`coefficient`}{
#'     Name of the original interaction coefficient being decomposed.
#'   }
#'   \item{`contrast_id`}{
#'     Sequential identifier assigned to the selected coefficient.
#'   }
#'   \item{`significance`}{
#'     Significance classification of the original coefficient:
#'     `"Sig"`, `"NotSig"`, or `"NoPValue"`.
#'   }
#'   \item{`p_value`}{
#'     P-value associated with the original interaction coefficient, when
#'     available.
#'   }
#'   \item{`split_var`}{
#'     Split variable selected for the coefficient.
#'   }
#'   \item{factor columns}{
#'     One column for each categorical predictor involved in the analyzed
#'     interactions. Values identify specific levels or aggregated groups.
#'   }
#'   \item{`n`}{
#'     Number of observed rows underlying the displayed configuration.
#'   }
#'   \item{`estimate`}{
#'     Observed/raw aggregate or model-based aggregate.
#'   }
#'   \item{`se`}{
#'     Standard error of the displayed estimate.
#'   }
#'   \item{`ci_low`}{
#'     Lower confidence-interval boundary.
#'   }
#'   \item{`ci_high`}{
#'     Upper confidence-interval boundary.
#'   }
#' }
#'
#' Depending on `configuration`, the output may also include:
#'
#' \describe{
#'   \item{`configuration_explicit`}{
#'     Mathematical description of the cells entering the aggregate.
#'   }
#'   \item{`configuration_semantic`}{
#'     Compact interpretation of the specific and averaged levels.
#'   }
#' }
#'
#' The returned object contains the following attributes:
#'
#' \describe{
#'   \item{`note`}{
#'     Methodological description of the estimation and aggregation
#'     procedures.
#'   }
#'   \item{`estimate_type`}{
#'     Whether the output contains `"raw"` or `"predicted"` estimates.
#'   }
#'   \item{`prediction_scale`}{
#'     Prediction scale used for model-based estimates. This is `NA` for raw
#'     estimates.
#'   }
#'   \item{`estimate_info`}{
#'     Human-readable description of the `estimate` column.
#'   }
#'   \item{`pvalue_source`}{
#'     Source of the p-values, such as `"model"`, `"lmerTest"`, or
#'     `"not_available"`.
#'   }
#'   \item{`split_rule`}{
#'     Rule used to select the split variable.
#'   }
#'   \item{`top_split`}{
#'     Split variable explicitly requested through `top_split`, when
#'     supplied.
#'   }
#' }
#'
#' @examples
#' ###############################################################
#' ## Two-way interaction
#' ###############################################################
#' \dontrun{
#'
#' library(lme4)
#' library(MASS)
#'
#' set.seed(123)
#'
#' dat <- expand.grid(
#'   ID = factor(paste0("S", 1:40)),
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2", "B3", "B4")),
#'   Rep = 1:2
#' )
#'
#' contrasts(dat$A) <- stats::contr.helmert(3)
#' contrasts(dat$B) <- MASS::contr.sdif(4)
#'
#' random_intercept <- rnorm(
#'   nlevels(dat$ID),
#'   mean = 0,
#'   sd = 1
#' )
#'
#' effect_A <- c(
#'   A1 = 0,
#'   A2 = 1,
#'   A3 = 2
#' )
#'
#' effect_B <- c(
#'   B1 = 0,
#'   B2 = 0.5,
#'   B3 = 1.5,
#'   B4 = 2.5
#' )
#'
#' dat$Score <- 10 +
#'   effect_A[dat$A] +
#'   effect_B[dat$B] +
#'   0.7 * effect_A[dat$A] * effect_B[dat$B] +
#'   random_intercept[dat$ID] +
#'   rnorm(nrow(dat), mean = 0, sd = 1.5)
#'
#' mod <- lme4::lmer(
#'   Score ~ A * B + (1 | ID),
#'   data = dat
#' )
#'
#' # A is selected as the conditioning variable.
#' # Labels follow the order B -> A.
#' out_high <- ContrastMultPred(
#'   model = mod,
#'   response = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "highest_complexity",
#'   configuration = "all"
#' )
#'
#' # B is selected as the conditioning variable.
#' # Labels follow the order A -> B.
#' out_low <- ContrastMultPred(
#'   model = mod,
#'   response = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "lowest_complexity",
#'   configuration = "all"
#' )
#'
#'
#' ###############################################################
#' ## Three-way interaction
#' ###############################################################
#'
#' dat3 <- expand.grid(
#'   ID = factor(paste0("S", 1:30)),
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2", "B3", "B4")),
#'   C = factor(c("C1", "C2")),
#'   Rep = 1:2
#' )
#'
#' contrasts(dat3$A) <- stats::contr.helmert(3)
#' contrasts(dat3$B) <- MASS::contr.sdif(4)
#' contrasts(dat3$C) <- stats::contr.treatment(2)
#'
#' dat3$Score <- 10 +
#'   effect_A[dat3$A] +
#'   effect_B[dat3$B] +
#'   1.2 * (dat3$C == "C2") +
#'   0.4 * effect_A[dat3$A] * effect_B[dat3$B] +
#'   rnorm(nrow(dat3), mean = 0, sd = 1.5)
#'
#' mod3 <- stats::lm(
#'   Score ~ A * B * C,
#'   data = dat3
#' )
#'
#' out3 <- ContrastMultPred(
#'   model = mod3,
#'   response = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "highest_complexity",
#'   configuration = "all"
#' )
#'
#'
#' ###############################################################
#' ## glmmTMB model
#' ###############################################################
#'
#' eta_count <- 1.2 +
#'   0.2 * effect_A[dat$A] +
#'   0.15 * effect_B[dat$B]
#'
#' dat$count <- rpois(
#'   nrow(dat),
#'   lambda = exp(eta_count)
#' )
#'
#' mod_tmb <- glmmTMB::glmmTMB(
#'   count ~ A * B + (1 | ID),
#'   family = poisson,
#'   data = dat
#' )
#'
#' out_tmb <- ContrastMultPred(
#'   model = mod_tmb,
#'   response = "count",
#'   estimate_type = "predicted",
#'   prediction_scale = "response",
#'   split_rule = "highest_complexity",
#'   configuration = "semantic"
#' )
#' }
#'
#' @export


ContrastMultPred <- function(model,
                             data=NULL,
                             response = all.vars(formula(model))[1],
                             interaction = NULL,
                             coefficients = NULL,
                             split_rule = c(
                               "last_variable",
                               "first_variable",
                               "lowest_complexity",
                               "highest_complexity"
                             ),
                             top_split = NULL,
                             conf_level = 0.95,
                             alpha_sig = 0.05,
                             configuration = c("none", "all", "explicit", "semantic"),
                             estimate_type = c("raw", "predicted"),
                             prediction_scale = c("response", "link"),
                             methodological_note = TRUE) {

  # ============================================================
  # 1. Argument matching and basic checks
  # ============================================================
  #
  # The function supports:
  #   - lm
  #   - glm
  #   - lmer/glmer models from lme4
  #
  # It is designed for models that include at least one interaction
  # involving categorical predictors.
  # ============================================================

  split_rule <- match.arg(split_rule)
  configuration <- match.arg(configuration)
  estimate_type <- match.arg(estimate_type)
  prediction_scale <- match.arg(prediction_scale)

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

  if (is.null(data)) {

    data <- .extract_analysis_data(model)

    if (is.null(data)) {
      stop(
        paste0(
          "`data` was not supplied and the analysis data could not be ",
          "recovered from `model`. Please provide the data frame used to fit ",
          "the model explicitly through `data`."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.",call. = FALSE)
  }

  if (!response %in% names(data)) {
    stop("The response variable is not present in 'data'.",call. = FALSE)
  }

  # ============================================================
  # 2. Extract model information
  # ============================================================
  #
  # This function extracts:
  #   - the model matrix;
  #   - coefficient names;
  #   - p-values, when available;
  #   - significance labels.
  #
  # For lmer models fitted with lme4, p-values are usually absent.
  # If lmerTest is installed, the function attempts to recover them.
  # ============================================================

  extract_model_info <- function(model,data, alpha_sig = 0.05) {

    fixed_effects <- .extract_fixef(model)

    coefficient_names <- setdiff(names(fixed_effects),"(Intercept)")

    p_information <- .extract_p_values(model = model,coefficient_names = coefficient_names, alpha = alpha_sig)

    list(
      coef_names = coefficient_names,
      sig_lookup = p_information$table,
      pvalue_source = p_information$source
    )
  }

  model_info <- extract_model_info(model, data, alpha_sig = alpha_sig)
  coef_names <- model_info$coef_names

  # Select the interaction coefficients to analyze.
  #
  # Backward compatibility is preserved:
  #   - interaction = NULL, coefficients = NULL:
  #       analyze all interaction coefficients;
  #   - interaction = NULL, coefficients supplied:
  #       analyze exactly the requested coefficients.
  #
  # When 'interaction' is supplied, reuse the same model-matrix-based
  # selection logic used by MoPlot_multpred(). This maps an interaction
  # term (e.g., "A:B:C") to all coefficient columns belonging to that
  # term and, when 'coefficients' is also supplied, validates that the
  # requested coefficients belong to the selected interaction.
  if (!is.null(interaction)) {

    selection <- .resolve_moplot_interaction_selection(
      model = model,
      interaction = interaction,
      coefficients = coefficients,
      interaction_info = FALSE
    )

    interaction <- selection$interaction
    coefficients <- selection$coefficients

  } else if (is.null(coefficients)) {

    coefficients <- coef_names[grepl(":", coef_names, fixed = TRUE)]
  }

  if (length(coefficients) == 0) {
    stop("No interaction coefficient to analyze.")
  }

  # Identify factor predictors involved in the model.
  fixed_formula <- .extract_fixed_formula(model)

  model_vars <- all.vars(fixed_formula)

  predictor_vars <- setdiff(model_vars,response)

  predictor_vars <- predictor_vars[predictor_vars %in% names(data)]

  factor_vars <- predictor_vars[vapply(data[predictor_vars],is.factor,logical(1))]

  if (length(factor_vars) == 0) {stop("There are no factor predictors in the model.")}

  # ============================================================
  # 3. Matrix comparison helper
  # ============================================================
  #
  # Used to identify which contrast system is being used by
  # comparing the observed contrast matrix with known reference
  # matrices, ignoring dimnames.
  # ============================================================

  same_matrix <- function(a, b, tol = 1e-12) {
    if (is.null(dim(a)) || is.null(dim(b))) return(FALSE)
    if (!all(dim(a) == dim(b))) return(FALSE)
    isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
  }

  # ============================================================
  # 4. Polynomial contrasts are not yet supported
  # ============================================================
  #
  # Interactions involving polynomial contrasts require a different
  # decomposition strategy because all levels contribute to each
  # contrast through polynomial weights.
  # ============================================================

  uses_contr_poly <- function(x) {
    cm <- contrasts(x)
    if (is.null(cm)) return(FALSE)
    ref_poly <- stats::contr.poly(nlevels(x))
    isTRUE(all.equal(unname(cm), unname(ref_poly), tolerance = 1e-12))
  }

  poly_vars <- factor_vars[sapply(data[factor_vars], uses_contr_poly)]

  if (length(poly_vars) > 0) {
    stop(
      paste0(
        "Work in progress: this function does not yet support categorical variables ",
        "coded with polynomial contrasts. Variables involved: ",
        paste(poly_vars, collapse = ", ")
      )
    )
  }

  # ============================================================
  # 5. Contrast complexity score
  # ============================================================
  #
  # Used when split_rule is:
  #   - "lowest_complexity"
  #   - "highest_complexity"
  #
  # Lower values indicate simpler contrast codings.
  # Higher values indicate more complex or custom contrast codings.
  # ============================================================

  contrast_complexity <- function(var) {

    contrast_matrix <- stats::contrasts(data[[var]])

    number_levels <- nlevels(data[[var]])

    if (is.null(contrast_matrix)) { return(Inf) }

    contrast_matrix_unname <- contrast_matrix
    dimnames(contrast_matrix_unname) <- NULL

    if (same_matrix(
      contrast_matrix_unname,
      stats::contr.treatment(number_levels)
    )) {
      return(0)
    }

    if (requireNamespace("MASS", quietly = TRUE)) {

      if (same_matrix(
        contrast_matrix_unname,
        MASS::contr.sdif(number_levels)
      )) {
        return(1)
      }
    }

    if (same_matrix(
      contrast_matrix_unname,
      stats::contr.sum(number_levels)
    )) {
      return(1)
    }

    if (same_matrix(
      contrast_matrix_unname,
      stats::contr.helmert(number_levels)
    )) {
      return(2)
    }

    if (same_matrix(
      contrast_matrix_unname,
      reverse_helmert(number_levels)
    )) {
      return(2)
    }

    if (same_matrix(
      contrast_matrix_unname,
      stats::contr.poly(number_levels)
    )) {
      return(3)
    }

    # Custom or unidentified contrast matrices are assigned the
    # highest complexity score.
    4
  }

  # ============================================================
  # 6. Methodological note
  # ============================================================

  disclaimer_note <- .build_multpred_methodological_note(
    model = model,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    split_rule = split_rule,
    top_split = top_split,
    configuration = configuration
  )

  if (isTRUE(methodological_note)) { message(disclaimer_note) }

  # ============================================================
  # 7. Build component lookup
  # ============================================================
  #
  # Each model coefficient is made of components.
  #
  # Example:
  #   A1:B2-1
  #
  # has two components:
  #   A1
  #   B2-1
  #
  # This lookup maps each component to:
  #   - the original variable;
  #   - the local contrast column.
  # ============================================================

  build_component_lookup <- function(var,
                                     data) {

    variable_formula <- stats::reformulate(var)

    terms_object <- stats::delete.response(stats::terms(variable_formula))

    variable_frame <- stats::model.frame(
      terms_object,
      data = data,
      drop.unused.levels = FALSE,
      na.action = stats::na.pass
    )

    # The contrast matrix is already stored in the factor contained
    # in 'data'. Therefore, it is not passed again through
    # contrasts.arg, avoiding duplicated contrast assignments and
    # the associated "contrasts dropped from factor" warning.
    model_matrix <- stats::model.matrix(terms_object,data = variable_frame)

    model_matrix <- model_matrix[ ,colnames(model_matrix) != "(Intercept)", drop = FALSE]

    if (ncol(model_matrix) == 0L) { return(NULL) }

    data.frame(
      component = colnames(model_matrix),
      variable = var,
      local_col = seq_len(ncol(model_matrix)),
      stringsAsFactors = FALSE
    )
  }

  lookup_list <- lapply(factor_vars, build_component_lookup, data = data)
  lookup <- do.call(rbind, lookup_list)

  if (is.null(lookup) || nrow(lookup) == 0) {
    stop("It was not possible to build the contrast-component lookup table.")
  }

  # ============================================================
  # 8. Data subsetting helper
  # ============================================================
  #
  # Subsets the original data according to a named list of
  # conditions, e.g. list(A = "A1", B = "B2").
  # ============================================================

  subset_by_conditions <- function(data, conditions) {

    out <- data

    if (length(conditions) == 0) return(out)

    for (nm in names(conditions)) {
      out <- out[out[[nm]] %in% conditions[[nm]], , drop = FALSE]
    }

    out
  }

  # ============================================================
  # 9. Label helpers
  # ============================================================
  #
  # These functions create readable labels describing which
  # configurations of levels are being averaged.
  # ============================================================

  collapse_source <- function(labels) {
    paste(unique(labels), collapse = "|||")
  }

  make_configuration_explicit <- function(labels) {

    labs <- unique(labels)

    if (length(labs) == 1) return(labs)

    paste0("(", paste(labs, collapse = " + "), ")/", length(labs))
  }

  parse_named_label <- function(x) {

    parts <- strsplit(x, "\\|")[[1]]
    vars <- sub("=.*$", "", parts)
    vals <- sub("^.*=", "", parts)

    stats::setNames(vals, vars)
  }

  # named_to_display <- function(x) {
  #
  #   vals <- parse_named_label(x)
  #   vals <- vals[factor_vars[factor_vars %in% names(vals)]]
  #
  #   paste(unname(vals), collapse = "_")
  # }

  # ============================================================
  # Configuration-variable ordering
  # ============================================================
  #
  # This helper defines a stable display order for configuration
  # labels.
  #
  # All target contrast variables are displayed first, following
  # their order among the categorical predictors in the model.
  #
  # The split variable is the conditioning variable and is therefore
  # always displayed last.
  # ============================================================

  order_configuration_variables <- function(variables,
                                            split_var = NULL) {

    variables <- unique(variables)

    # Preserve the stable order of categorical predictors in the
    # fixed-effects model.
    ordered_variables <- factor_vars[factor_vars %in% variables]

    # Retain any unexpected variable that is not included in
    # factor_vars.
    remaining_variables <- setdiff(variables,ordered_variables)

    ordered_variables <- c(ordered_variables,remaining_variables)

    # The split variable is the conditioning variable and is placed
    # after all target contrast variables.
    if (!is.null(split_var) &&
        split_var %in% ordered_variables) {
          ordered_variables <- c(ordered_variables[ordered_variables != split_var],split_var)
    }

    ordered_variables
  }

  make_configuration_semantic <- function(named_labels,
                                          split_var = NULL) {

    labels <- unique(named_labels)

    parsed_labels <- lapply(labels,parse_named_label)

    variables_present <- unique(unlist(lapply(parsed_labels,names)))

    variable_order <- order_configuration_variables(
      variables = variables_present,
      split_var = split_var
    )

    same_variables <- all(
      vapply(
        parsed_labels,
        function(values) {
          all(
            variable_order %in% names(values)
          )
        },
        logical(1)
      )
    )

    if (!same_variables) {

      return(
        paste0(
          "Mean(",
          length(labels),
          " cells)"
        )
      )
    }

    value_matrix <- do.call(
      rbind,
      lapply(
        parsed_labels,
        function(values) {
          unname(
            values[variable_order]
          )
        }
      )
    )

    colnames(value_matrix) <- variable_order

    varying_variables <- apply(
      value_matrix,
      2L,
      function(values) {
        length(unique(values)) > 1L
      }
    )

    # No variable varies across the configurations.
    if (!any(varying_variables)) {

      return(
        paste(
          value_matrix[1L, ],
          collapse = "_"
        )
      )
    }

    semantic_tokens <- vapply(
      variable_order,
      function(variable_name) {

        values <- unique(
          value_matrix[, variable_name]
        )

        if (length(values) == 1L) {

          return(values)
        }

        paste0(
          "Mean(",
          paste(
            values,
            collapse = ","
          ),
          ")"
        )
      },
      character(1)
    )

    paste(
      semantic_tokens,
      collapse = "_"
    )
  }




  make_variable_columns <- function(named_labels, all_vars) {

    labs <- unique(named_labels)
    parsed <- lapply(labs, parse_named_label)

    out <- stats::setNames(as.list(rep("", length(all_vars))), all_vars)

    vars_present <- unique(unlist(lapply(parsed, names)))
    vars_present <- vars_present[vars_present %in% all_vars]

    for (v in vars_present) {

      vals <- unique(vapply(parsed, function(z) z[[v]], character(1)))

      if (length(vals) == 1) {
        out[[v]] <- vals
      } else {
        out[[v]] <- paste0("Mean(", paste(vals, collapse = ","), ")")
      }
    }

    as.data.frame(out, stringsAsFactors = FALSE)
  }




  make_ordered_label_info <- function(conditions,
                                      split_var = NULL) {

    if (length(conditions) == 0L) {

      return(
        list(
          display = "ALL",
          named = "ALL=ALL"
        )
      )
    }

    variables <- names(conditions)

    variable_order <- order_configuration_variables(
      variables = variables,
      split_var = split_var
    )

    values <- vapply(
      variable_order,
      function(variable_name) {
        as.character(
          conditions[[variable_name]]
        )
      },
      character(1)
    )

    list(
      display = paste(
        values,
        collapse = "_"
      ),
      named = paste(
        paste0(
          variable_order,
          "=",
          values
        ),
        collapse = "|"
      )
    )
  }

  # ============================================================
  # 10. Raw cell statistics
  # ============================================================
  #
  # Computes observed means and classical standard errors within
  # each cell/configuration.
  # ============================================================

  compute_raw_cell_stats <- function(sub) {

    x <- sub[[response]]
    x <- x[!is.na(x)]
    n <- length(x)

    if (n == 0) {
      return(data.frame(
        n = 0,
        estimate = NA_real_,
        sd = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    }

    m <- mean(x)
    sdv <- if (n > 1) stats::sd(x) else NA_real_
    se <- if (n > 1) sdv / sqrt(n) else NA_real_

    alpha <- 1 - conf_level
    crit <- if (n > 1) stats::qt(1 - alpha / 2, df = n - 1) else NA_real_
    moe <- if (n > 1) crit * se else NA_real_

    data.frame(
      n = n,
      estimate = m,
      sd = sdv,
      se = se,
      ci_low = if (n > 1) m - moe else NA_real_,
      ci_high = if (n > 1) m + moe else NA_real_
    )
  }


  # ============================================================
  # 13. Predicted cell statistics
  # ============================================================
  #
  # Computes model-based estimates for a cell/configuration.
  #
  # For lm/lmer:
  #   estimate = mean(X beta)
  #
  # For glm/glmer on link scale:
  #   estimate = mean(eta)
  #
  # For glm/glmer on response scale:
  #   estimate = mean(linkinv(eta))
  #
  # SEs are obtained from the fixed-effect variance-covariance
  # matrix using the delta method when needed.
  # ============================================================

  compute_predicted_cell_stats <- function(sub) {

    number_observations <- nrow(sub)

    if (number_observations == 0L) {

      return(
        data.frame(
          n = 0,
          estimate = NA_real_,
          sd = NA_real_,
          se = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    # Compute the average population-level prediction for all rows
    # contributing to the current interaction cell.
    #
    # On the link scale, this corresponds to the average linear
    # predictor. On the response scale, it corresponds to the average
    # of the inverse-link-transformed predictions.
    average_prediction <- .compute_average_prediction(
      model = model,
      newdata = sub,
      prediction_scale = prediction_scale,
      conf_level = conf_level
    )

    data.frame(
      n = number_observations,
      estimate = average_prediction$estimate,
      sd = NA_real_,
      se = average_prediction$se,
      ci_low = average_prediction$ci_low,
      ci_high = average_prediction$ci_high,
      stringsAsFactors = FALSE
    )
  }

  # ============================================================
  # 14. Unified cell statistics dispatcher
  # ============================================================
  #
  # In addition to the displayed statistics, predicted estimates
  # retain the rows of newdata and their internal averaging weights.
  #
  # These internal objects are propagated through the recursive
  # decomposition and allow higher-level aggregates to be calculated
  # directly from the fixed-effect covariance matrix.
  # ============================================================

  compute_cell_stats <- function(sub) {

    if (estimate_type == "raw") {

      statistics <- compute_raw_cell_stats(sub)

      # Raw estimates do not require model-based prediction bundles.
      statistics$.prediction_data <- list(NULL)
      statistics$.prediction_weights <- list(NULL)

      return(statistics)
    }

    statistics <- compute_predicted_cell_stats(sub)

    number_rows <- nrow(sub)

    # Each row contributes equally to the cell-level prediction.
    # The weights are stored explicitly so that they can later be
    # multiplied by contrast weights during recursive aggregation.
    row_weights <- if (number_rows > 0L) {
      rep(1 / number_rows,number_rows)
    } else {
      numeric(0)
    }

    statistics$.prediction_data <- list(sub)
    statistics$.prediction_weights <- list(row_weights)

    statistics
  }

  # ============================================================
  # 15. Weighted aggregation of raw estimates and SEs
  # ============================================================
  #
  # This helper is used only when estimate_type = "raw".
  #
  # Raw cell means are based on disjoint groups and their standard
  # errors are combined through classical variance propagation.
  #
  # Predicted estimates are not aggregated with this helper because
  # model-based predictions share the same fixed-effect coefficients
  # and are generally correlated.
  # ============================================================

  combine_weighted_mean_stats <- function(estimates, ses, ns, weights) {

    ok <- !is.na(estimates) & !is.na(weights)

    estimates <- estimates[ok]
    ses <- ses[ok]
    ns <- ns[ok]
    weights <- weights[ok]

    if (length(estimates) == 0 || sum(weights) == 0) {
      return(data.frame(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        df = NA_real_
      ))
    }

    # Normalize weights so that they sum to 1.
    a <- weights / sum(weights)

    estimate_agg <- sum(a * estimates)

    # Variance propagation for a weighted mean.
    ses2 <- ifelse(is.na(ses), 0, ses)^2
    var_components <- (a^2) * ses2
    se_agg <- sqrt(sum(var_components, na.rm = TRUE))

    # Approximate degrees of freedom for raw observed estimates.
    valid_df <- !is.na(ns) & ns > 1 & !is.na(var_components)

    if (sum(valid_df) == 0 || sum(var_components[valid_df]) == 0) {
      df_agg <- NA_real_
    } else {
      num <- (sum(var_components[valid_df]))^2
      den <- sum((var_components[valid_df]^2) / (ns[valid_df] - 1))
      df_agg <- if (den > 0) num / den else NA_real_
    }

    alpha <- 1 - conf_level

    crit <- if (!is.na(df_agg) && estimate_type == "raw") {
      stats::qt(1 - alpha / 2, df = df_agg)
    } else {
      stats::qnorm(1 - alpha / 2)
    }

    data.frame(
      estimate = estimate_agg,
      se = se_agg,
      ci_low = estimate_agg - crit * se_agg,
      ci_high = estimate_agg + crit * se_agg,
      df = df_agg
    )
  }


  # ============================================================
  # 16. Exact aggregation of model-based prediction bundles
  # ============================================================
  #
  # Predicted cell estimates are generally correlated because they
  # are obtained from the same fitted coefficient vector.
  #
  # Consequently, their standard errors must not be combined as if
  # the cell estimates were independent.
  #
  # This helper:
  #
  #   1. receives the newdata rows underlying each child estimate;
  #   2. receives the within-child averaging weights;
  #   3. applies the higher-level contrast weights;
  #   4. combines all rows into one prediction dataset;
  #   5. calls .compute_average_prediction() once.
  #
  # This preserves the complete covariance structure implied by the
  # fixed-effect variance-covariance matrix.
  # ============================================================

  combine_prediction_bundles <- function(prediction_data,
                                         prediction_weights,
                                         weights) {

    if (!is.list(prediction_data) ||
        !is.list(prediction_weights)) {

      stop(
        paste0(
          "Internal prediction bundles must be stored as list columns."
        ),
        call. = FALSE
      )
    }

    if (length(prediction_data) != length(weights) ||
        length(prediction_weights) != length(weights)) {

      stop(
        paste0(
          "The number of prediction bundles does not match the ",
          "number of aggregation weights."
        ),
        call. = FALSE
      )
    }

    valid_bundle <- vapply(
      seq_along(prediction_data),
      function(index) {

        current_data <- prediction_data[[index]]
        current_weights <- prediction_weights[[index]]

        is.data.frame(current_data) &&
          nrow(current_data) > 0L &&
          is.numeric(current_weights) &&
          length(current_weights) == nrow(current_data) &&
          all(!is.na(current_weights)) &&
          sum(current_weights) > 0
      },
      logical(1)
    )

    valid_bundle <- valid_bundle &
      !is.na(weights) &
      weights > 0

    if (!any(valid_bundle)) {

      return(
        list(
          statistics = data.frame(
            estimate = NA_real_,
            se = NA_real_,
            ci_low = NA_real_,
            ci_high = NA_real_,
            stringsAsFactors = FALSE
          ),
          prediction_data = NULL,
          prediction_weights = NULL
        )
      )
    }

    prediction_data <- prediction_data[valid_bundle]
    prediction_weights <- prediction_weights[valid_bundle]
    weights <- weights[valid_bundle]

    # Normalize the higher-level contrast weights.
    outer_weights <- weights / sum(weights)

    weighted_data <- vector("list",length(prediction_data))

    weighted_rows <- vector("list",length(prediction_data))

    for (index in seq_along(prediction_data)) {

      current_data <- prediction_data[[index]]
      current_weights <- prediction_weights[[index]]

      # Normalize the weights within each child bundle. The child
      # estimate is therefore preserved before applying its outer
      # contrast weight.
      current_weights <- current_weights / sum(current_weights)

      weighted_data[[index]] <- current_data

      weighted_rows[[index]] <- current_weights * outer_weights[index]
    }

    combined_data <- do.call(rbind,weighted_data)

    rownames(combined_data) <- NULL

    combined_weights <- unlist(weighted_rows,use.names = FALSE)

    # Calculate the final average and its uncertainty directly from
    # the fixed-effect coefficient covariance matrix.
    statistics <- .compute_average_prediction(
      model = model,
      newdata = combined_data,
      prediction_scale = prediction_scale,
      conf_level = conf_level,
      weights = combined_weights
    )

    list(
      statistics = statistics,
      prediction_data = combined_data,
      prediction_weights = combined_weights
    )
  }

  # ============================================================
  # 16. Choose split variable
  # ============================================================
  #
  # This function decides which variable should be used as the
  # main splitting variable when decomposing an interaction.
  #
  # Available rules:
  #
  #   first_variable:
  #     use the first variable appearing in the interaction term.
  #
  #   last_variable:
  #     use the last variable appearing in the interaction term.
  #
  #   lowest_complexity:
  #     choose the variable with the simplest contrast coding.
  #
  #   highest_complexity:
  #     choose the variable with the most complex contrast coding.
  #
  # If several variables have the same complexity score, the
  # rightmost tied variable is selected.
  #
  # top_split overrides all automatic choices.
  # ============================================================

  choose_split_var <- function(component_df,
                               rule = c(
                                 "last_variable",
                                 "first_variable",
                                 "lowest_complexity",
                                 "highest_complexity"
                               ),
                               forced = NULL) {

    rule <- match.arg(rule)

    vars <- component_df$variable
    vars <- vars[!duplicated(vars)]

    if (!is.null(forced)) {
      if (!forced %in% vars) {
        stop(sprintf("The variable '%s' is not involved in the coefficient.", forced))
      }
      return(forced)
    }

    if (length(vars) == 1) return(vars)

    if (rule == "first_variable") {
      return(vars[1])
    }

    if (rule == "last_variable") {
      return(vars[length(vars)])
    }

    scores <- vapply(vars, contrast_complexity, numeric(1))

    if (rule == "lowest_complexity") {
      target_score <- min(scores, na.rm = TRUE)
    } else {
      target_score <- max(scores, na.rm = TRUE)
    }

    candidate_vars <- vars[scores == target_score]

    if (length(candidate_vars) > 1) {
      warning(
        paste0(
          "Multiple variables have the same complexity score (",
          paste(candidate_vars, collapse = ", "),
          "). The rightmost tied variable in the interaction is selected: ",
          candidate_vars[length(candidate_vars)],
          "."
        )
      )
    }

    candidate_vars[length(candidate_vars)]
  }

  # ============================================================
  # 17. Extract contrast information for one variable
  # ============================================================
  #
  # This function returns the levels involved in a specific
  # contrast column of a factor.
  #
  # It handles:
  #   - sum contrasts;
  #   - treatment contrasts;
  #   - sliding-difference contrasts;
  #   - standard/custom contrast matrices.
  #
  # Treatment and sliding-difference contrasts are converted to
  # a hypothesis-matrix representation because their raw coding
  # matrix is not always directly interpretable as positive vs
  # negative level groups.
  # ============================================================

  get_contrast_info <- function(var, local_col) {

    if (!is.factor(data[[var]])) {
      stop(sprintf("The variable '%s' is not a factor in 'data'.", var))
    }

    cm <- contrasts(data[[var]])

    if (is.null(cm)) {
      stop(sprintf("The variable '%s' has no available contrast matrix.", var))
    }

    levels_var <- levels(data[[var]])
    nlev <- nlevels(data[[var]])

    if (local_col > ncol(cm)) {
      stop(sprintf("Contrast index %d is out of range for variable '%s'.", local_col, var))
    }

    # Sum contrasts are interpreted as level vs grand mean.
    if (same_matrix(cm, stats::contr.sum(nlev))) {

      w <- cm[, local_col]
      names(w) <- rownames(cm)

      focal_levels <- names(w)[w > 0]
      focal_weights <- w[w > 0]

      return(list(
        contrast_type = "sum",
        levels_all = levels_var,
        focal_levels = focal_levels,
        focal_weights = focal_weights,
        label_negative = paste0("Mean(", paste(levels_var, collapse = ","), ")")
      ))
    }

    dm_names <- dimnames(cm)
    cm_un <- cm
    dimnames(cm_un) <- NULL

    is_treatment <- same_matrix(
      cm_un,
      stats::contr.treatment(nlev)
    )

    is_sliding_difference <- FALSE

    if (requireNamespace("MASS", quietly = TRUE)) {

      is_sliding_difference <- same_matrix(
        cm_un,
        MASS::contr.sdif(nlev)
      )
    }

    if ((is_treatment || is_sliding_difference) &&
        !requireNamespace("hypr", quietly = TRUE)) {

      stop(
        paste0(
          "The 'hypr' package is required to interpret treatment ",
          "and sliding-difference contrasts."
        ),
        call. = FALSE
      )
    }

    if (is_treatment) {

      hypothesis_object <- hypr::hypr()

      hypr::cmat(hypothesis_object,add_intercept = TRUE) <- cm_un

      cm_use <- t(hypr::hmat(hypothesis_object))[, -1, drop = FALSE]

      dimnames(cm_use) <- dm_names

    } else if (is_sliding_difference) {

      hypothesis_object <- hypr::hypr()

      hypr::cmat(hypothesis_object,add_intercept = TRUE) <- cm_un

      cm_use <- t(hypr::hmat(hypothesis_object))[, -1, drop = FALSE]

      dimnames(cm_use) <- dm_names

    } else {

      cm_use <- cm
    }

    w <- cm_use[, local_col]
    names(w) <- rownames(cm_use)

    list(
      contrast_type = "standard",
      weights = w,
      levels_positive = names(w)[w > 0],
      levels_negative = names(w)[w < 0],
      weights_positive = w[w > 0],
      weights_negative = abs(w[w < 0])
    )
  }

  # ============================================================
  # 18. Aggregate across actual levels
  # ============================================================
  #
  # This is used when a contrast group contains multiple levels.
  #
  # For each actual level, the function recursively builds the
  # inner groups for the remaining variables, then combines the
  # resulting estimates using the current contrast weights.
  # ============================================================

  aggregate_across_actual_levels <- function(actual_levels,
                                             actual_weights,
                                             component_df_remaining,
                                             current_var,
                                             conditions = list(),
                                             split_rule = "last_variable",
                                             top_split_var = NULL) {

    nested_list <- lapply(seq_along(actual_levels), function(i) {

      lev <- actual_levels[i]

      tmp <- build_inner_groups_for_one_level(
        component_df = component_df_remaining,
        conditions = c(conditions, setNames(list(lev), current_var)),
        split_rule = split_rule,
        top_split_var = top_split_var
      )

      tmp$actual_level <- lev
      tmp$actual_weight <- actual_weights[i]
      tmp
    })

    nested_all <- do.call(rbind, nested_list)
    keys <- unique(nested_all$inner_side)

    do.call(rbind, lapply(keys, function(k) {

      sub <- nested_all[nested_all$inner_side == k, , drop = FALSE]

      if (estimate_type == "raw") {

        agg <- combine_weighted_mean_stats(
          estimates = sub$estimate,
          ses = sub$se,
          ns = sub$n,
          weights = sub$actual_weight
        )

        aggregated_prediction_data <- NULL
        aggregated_prediction_weights <- NULL

      } else {

        prediction_aggregation <- combine_prediction_bundles(
          prediction_data = sub$.prediction_data,
          prediction_weights = sub$.prediction_weights,
          weights = sub$actual_weight
        )

        agg <- prediction_aggregation$statistics

        aggregated_prediction_data <-
          prediction_aggregation$prediction_data

        aggregated_prediction_weights <-
          prediction_aggregation$prediction_weights
      }

      labels_source <- unique(unlist(strsplit(sub$configuration_source, "\\|\\|\\|")))
      named_source <- unique(unlist(strsplit(sub$configuration_named_source, "\\|\\|\\|")))
      var_cols <- make_variable_columns(named_source, factor_vars)

      aggregated_result <- cbind(
        data.frame(
          inner_side = k,
          inner_variable_path = sub$inner_variable_path[1L],
          configuration_source = collapse_source(labels_source),
          configuration_named_source = collapse_source(named_source),
          configuration_explicit = make_configuration_explicit(
            sub$configuration_explicit
          ),
          configuration_semantic = make_configuration_semantic(
            named_source,
            split_var = top_split_var
          ),
          n = sum(sub$n, na.rm = TRUE),
          estimate = agg$estimate,
          sd = NA_real_,
          se = agg$se,
          ci_low = agg$ci_low,
          ci_high = agg$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )

      aggregated_result$.prediction_data <-
        list(aggregated_prediction_data)

      aggregated_result$.prediction_weights <-
        list(aggregated_prediction_weights)

      aggregated_result
    }))
  }

  # ============================================================
  # 19. Build inner groups for one level
  # ============================================================
  #
  # This is the recursive core of the function.
  #
  # It takes a set of remaining contrast components and:
  #   1. chooses the next variable to decompose;
  #   2. identifies positive and negative levels for its contrast;
  #   3. either computes cell statistics directly, or recursively
  #      decomposes the remaining variables.
  #
  # This allows the function to work with two-way and higher-order
  # interactions.
  # ============================================================

  build_inner_groups_for_one_level <- function(component_df,
                                               conditions = list(),
                                               split_rule = "last_variable",
                                               top_split_var = NULL) {

    # Base case: no more variables to decompose.
    # Compute the estimate for the current cell/configuration.
    if (nrow(component_df) == 0) {

      sub <- subset_by_conditions(data, conditions)
      stats_df <- compute_cell_stats(sub)
      lab <- make_ordered_label_info(conditions, split_var = top_split_var)
      var_cols <- make_variable_columns(lab$named, factor_vars)

      cell_result <- cbind(
        data.frame(
          inner_side = "cell",
          inner_variable_path = "",
          configuration_source = lab$display,
          configuration_named_source = lab$named,
          configuration_explicit = lab$display,
          configuration_semantic = make_configuration_semantic(
            lab$named,
            split_var = top_split_var
          ),
          n = stats_df$n,
          estimate = stats_df$estimate,
          sd = stats_df$sd,
          se = stats_df$se,
          ci_low = stats_df$ci_low,
          ci_high = stats_df$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )

      # Internal prediction information used by higher recursive levels.
      cell_result$.prediction_data <-
        stats_df$.prediction_data

      cell_result$.prediction_weights <-
        stats_df$.prediction_weights

      return(cell_result)
    }

    # Choose the next variable to decompose.
    current_var <- choose_split_var(component_df, rule = split_rule)
    row_idx <- which(component_df$variable == current_var)[1]
    local_col <- component_df$local_col[row_idx]

    contrast_info <- get_contrast_info(current_var, local_col)
    remaining <- component_df[component_df$variable != current_var, , drop = FALSE]

    # Identify positive and negative levels of the current contrast.
    if (contrast_info$contrast_type == "sum") {

      pos_levels <- contrast_info$focal_levels
      pos_weights <- as.numeric(contrast_info$focal_weights)
      neg_levels <- contrast_info$levels_all
      neg_weights <- rep(1, length(neg_levels))

    } else {

      weights <- contrast_info$weights

      pos_levels <- names(weights)[weights > 0]
      neg_levels <- names(weights)[weights < 0]
      pos_weights <- as.numeric(weights[weights > 0])
      neg_weights <- as.numeric(abs(weights[weights < 0]))
    }

    # If no variables remain, compute the positive and negative
    # contrast groups directly.
    if (nrow(remaining) == 0) {

      pos_cells <- lapply(
        seq_along(pos_levels),
        function(i) {

          lev <- pos_levels[i]

          conds <- c(conditions,setNames(list(lev),current_var))

          sub <- subset_by_conditions(data,conds)

          stats_df <- compute_cell_stats(sub)

          lab <- make_ordered_label_info(conds,split_var = top_split_var)

          positive_cell <- data.frame(
            label = lab$display,
            named_label = lab$named,
            n = stats_df$n,
            estimate = stats_df$estimate,
            se = stats_df$se,
            weight = pos_weights[i],
            stringsAsFactors = FALSE
          )

          positive_cell$.prediction_data <- stats_df$.prediction_data

          positive_cell$.prediction_weights <- stats_df$.prediction_weights

          positive_cell
        }
      )

      pos_cells <- do.call(rbind,pos_cells)

      neg_cells <- lapply(
        seq_along(neg_levels),
        function(i) {

          lev <- neg_levels[i]

          conds <- c(
            conditions,
            setNames(
              list(lev),
              current_var
            )
          )

          sub <- subset_by_conditions(data,conds)

          stats_df <- compute_cell_stats(sub)

          lab <- make_ordered_label_info(conds,split_var = top_split_var)

          negative_cell <- data.frame(
            label = lab$display,
            named_label = lab$named,
            n = stats_df$n,
            estimate = stats_df$estimate,
            se = stats_df$se,
            weight = neg_weights[i],
            stringsAsFactors = FALSE
          )

          negative_cell$.prediction_data <- stats_df$.prediction_data

          negative_cell$.prediction_weights <- stats_df$.prediction_weights

          negative_cell
        }
      )

      neg_cells <- do.call(rbind,neg_cells)

      if (estimate_type == "raw") {

        pos_stats <- combine_weighted_mean_stats(
          estimates = pos_cells$estimate,
          ses = pos_cells$se,
          ns = pos_cells$n,
          weights = pos_cells$weight
        )

        neg_stats <- combine_weighted_mean_stats(
          estimates = neg_cells$estimate,
          ses = neg_cells$se,
          ns = neg_cells$n,
          weights = neg_cells$weight
        )

        positive_prediction_data <- NULL
        positive_prediction_weights <- NULL

        negative_prediction_data <- NULL
        negative_prediction_weights <- NULL

      } else {

        positive_aggregation <- combine_prediction_bundles(
          prediction_data = pos_cells$.prediction_data,
          prediction_weights = pos_cells$.prediction_weights,
          weights = pos_cells$weight
        )

        negative_aggregation <- combine_prediction_bundles(
          prediction_data = neg_cells$.prediction_data,
          prediction_weights = neg_cells$.prediction_weights,
          weights = neg_cells$weight
        )

        pos_stats <- positive_aggregation$statistics
        neg_stats <- negative_aggregation$statistics

        positive_prediction_data <- positive_aggregation$prediction_data

        positive_prediction_weights <- positive_aggregation$prediction_weights

        negative_prediction_data <- negative_aggregation$prediction_data

        negative_prediction_weights <- negative_aggregation$prediction_weights
      }

      pos_var_cols <- make_variable_columns(pos_cells$named_label, factor_vars)
      neg_var_cols <- make_variable_columns(neg_cells$named_label, factor_vars)

      positive_result <- cbind(
        data.frame(
          inner_side = "positive",
          inner_variable_path = current_var,
          configuration_source = collapse_source(
            pos_cells$label
          ),
          configuration_named_source = collapse_source(
            pos_cells$named_label
          ),
          configuration_explicit = make_configuration_explicit(
            pos_cells$label
          ),
          configuration_semantic = make_configuration_semantic(
            pos_cells$named_label,
            split_var = top_split_var
          ),
          n = sum(pos_cells$n, na.rm = TRUE),
          estimate = pos_stats$estimate,
          sd = NA_real_,
          se = pos_stats$se,
          ci_low = pos_stats$ci_low,
          ci_high = pos_stats$ci_high,
          stringsAsFactors = FALSE
        ),
        pos_var_cols
      )

      positive_result$.prediction_data <- list(positive_prediction_data)

      positive_result$.prediction_weights <- list(positive_prediction_weights)


      negative_result <- cbind(
        data.frame(
          inner_side = "negative",
          inner_variable_path = current_var,
          configuration_source = collapse_source(
            neg_cells$label
          ),
          configuration_named_source = collapse_source(
            neg_cells$named_label
          ),
          configuration_explicit = make_configuration_explicit(
            neg_cells$label
          ),
          configuration_semantic = make_configuration_semantic(
            neg_cells$named_label,
            split_var = top_split_var
          ),
          n = sum(neg_cells$n, na.rm = TRUE),
          estimate = neg_stats$estimate,
          sd = NA_real_,
          se = neg_stats$se,
          ci_low = neg_stats$ci_low,
          ci_high = neg_stats$ci_high,
          stringsAsFactors = FALSE
        ),
        neg_var_cols
      )

      negative_result$.prediction_data <-
        list(negative_prediction_data)

      negative_result$.prediction_weights <-
        list(negative_prediction_weights)

      return(
        rbind(
          positive_result,
          negative_result
        )
      )
    }

    # Recursive case: there are still variables to decompose.
    pos_df <- aggregate_across_actual_levels(
      actual_levels = pos_levels,
      actual_weights = pos_weights,
      component_df_remaining = remaining,
      current_var = current_var,
      conditions = conditions,
      split_rule = split_rule,
      top_split_var = top_split_var
    )

    neg_df <- aggregate_across_actual_levels(
      actual_levels = neg_levels,
      actual_weights = neg_weights,
      component_df_remaining = remaining,
      current_var = current_var,
      conditions = conditions,
      split_rule = split_rule,
      top_split_var = top_split_var
    )

    pos_df$inner_side <- paste0(
      "positive",
      ifelse(pos_df$inner_side == "cell", "", paste0("_", pos_df$inner_side))
    )

    neg_df$inner_side <- paste0(
      "negative",
      ifelse(neg_df$inner_side == "cell", "", paste0("_", neg_df$inner_side))
    )

    pos_df$inner_variable_path <- ifelse(
      is.na(pos_df$inner_variable_path) |
        pos_df$inner_variable_path == "",
      current_var,
      paste0(
        current_var,
        "_",
        pos_df$inner_variable_path
      )
    )

    neg_df$inner_variable_path <- ifelse(
      is.na(neg_df$inner_variable_path) |
        neg_df$inner_variable_path == "",
      current_var,
      paste0(
        current_var,
        "_",
        neg_df$inner_variable_path
      )
    )

    rbind(pos_df, neg_df)
  }

  # ============================================================
  # 20. Build the full bundle for one side of the top split
  # ============================================================
  #
  # For the top split variable, this function builds all inner
  # configurations for each actual level and then aggregates them
  # according to the contrast weights.
  # ============================================================

  build_group_bundle <- function(actual_levels,
                                 actual_weights,
                                 coefficient,
                                 split_var,
                                 component_df_remaining,
                                 split_rule = "last_variable") {

    nested_list <- lapply(seq_along(actual_levels), function(i) {

      lev <- actual_levels[i]

      tmp <- build_inner_groups_for_one_level(
        component_df = component_df_remaining,
        conditions = setNames(list(lev), split_var),
        split_rule = split_rule,
        top_split_var = split_var
      )

      tmp$actual_level <- lev
      tmp$actual_weight <- actual_weights[i]
      tmp
    })

    nested_all <- do.call(rbind, nested_list)
    keys <- unique(nested_all$inner_side)

    do.call(rbind, lapply(keys, function(k) {

      sub <- nested_all[nested_all$inner_side == k, , drop = FALSE]

      if (estimate_type == "raw") {

        agg <- combine_weighted_mean_stats(
          estimates = sub$estimate,
          ses = sub$se,
          ns = sub$n,
          weights = sub$actual_weight
        )

        aggregated_prediction_data <- NULL
        aggregated_prediction_weights <- NULL

      } else {

        prediction_aggregation <- combine_prediction_bundles(
          prediction_data = sub$.prediction_data,
          prediction_weights = sub$.prediction_weights,
          weights = sub$actual_weight
        )

        agg <- prediction_aggregation$statistics

        aggregated_prediction_data <-
          prediction_aggregation$prediction_data

        aggregated_prediction_weights <-
          prediction_aggregation$prediction_weights
      }

      labels_source <- unique(unlist(strsplit(sub$configuration_source, "\\|\\|\\|")))
      named_source <- unique(unlist(strsplit(sub$configuration_named_source, "\\|\\|\\|")))
      var_cols <- make_variable_columns(named_source, factor_vars)

      bundle_result <- cbind(
        data.frame(
          coefficient = coefficient,
          split_var = split_var,
          inner_side = k,
          inner_variable_path = sub$inner_variable_path[1L],
          configuration_source = collapse_source(labels_source),
          configuration_named_source = collapse_source(named_source),
          configuration_explicit = make_configuration_explicit(
            sub$configuration_explicit
          ),
          configuration_semantic = make_configuration_semantic(
            named_source,
            split_var = split_var
          ),
          n = sum(sub$n, na.rm = TRUE),
          estimate = agg$estimate,
          se = agg$se,
          ci_low = agg$ci_low,
          ci_high = agg$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )

      bundle_result$.prediction_data <-
        list(aggregated_prediction_data)

      bundle_result$.prediction_weights <-
        list(aggregated_prediction_weights)

      bundle_result
    }))
  }

  # ============================================================
  # 21. Decompose one interaction coefficient
  # ============================================================
  #
  # This function takes a single model coefficient, such as:
  #
  #   A1:B2-1
  #
  # and decomposes it into interpretable positive/negative level
  # configurations according to the contrast matrices of the
  # variables involved.
  # ============================================================

  decompose_one_coefficient <- function(coef_name,
                                        split_rule = "last_variable",
                                        top_split = NULL) {

    parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]

    component_df <- lookup[lookup$component %in% parts, , drop = FALSE]
    component_df <- component_df[match(parts, component_df$component), , drop = FALSE]

    if (nrow(component_df) != length(parts)) {
      warning(sprintf("It was not possible to fully map coefficient '%s'. Skipped.", coef_name))
      return(NULL)
    }

    split_var <- choose_split_var(component_df, rule = split_rule, forced = top_split)

    split_row <- which(component_df$variable == split_var)[1]
    split_col <- component_df$local_col[split_row]

    contrast_info <- get_contrast_info(split_var, split_col)
    remaining <- component_df[component_df$variable != split_var, , drop = FALSE]

    if (contrast_info$contrast_type == "sum") {

      pos_levels <- contrast_info$focal_levels
      pos_weights <- as.numeric(contrast_info$focal_weights)
      neg_levels <- contrast_info$levels_all
      neg_weights <- rep(1, length(neg_levels))

    } else {

      weights <- contrast_info$weights

      pos_levels <- names(weights)[weights > 0]
      neg_levels <- names(weights)[weights < 0]
      pos_weights <- as.numeric(weights[weights > 0])
      neg_weights <- as.numeric(abs(weights[weights < 0]))
    }

    pos_bundle <- build_group_bundle(
      actual_levels = pos_levels,
      actual_weights = pos_weights,
      coefficient = coef_name,
      split_var = split_var,
      component_df_remaining = remaining,
      split_rule = split_rule
    )

    neg_bundle <- build_group_bundle(
      actual_levels = neg_levels,
      actual_weights = neg_weights,
      coefficient = coef_name,
      split_var = split_var,
      component_df_remaining = remaining,
      split_rule = split_rule
    )

    pos_bundle$split_side <- "positive"
    neg_bundle$split_side <- "negative"

    rbind(pos_bundle, neg_bundle)
  }

  # ============================================================
  # 22. Apply decomposition to all selected coefficients
  # ============================================================

  out_list <- lapply(
    coefficients,
    decompose_one_coefficient,
    split_rule = split_rule,
    top_split = top_split
  )

  out_list <- out_list[!sapply(out_list, is.null)]

  if (length(out_list) == 0) {
    stop("No coefficient was successfully processed.")
  }

  out <- do.call(rbind, out_list)

  # Internal prediction bundles are needed only during recursive
  # aggregation and are not part of the user-facing output.
  internal_prediction_columns <- c(
    ".prediction_data",
    ".prediction_weights"
  )

  out <- out[ ,setdiff(names(out),internal_prediction_columns),drop = FALSE]

  # ============================================================
  # 23. Add significance information
  # ============================================================

  contrast_lookup <- data.frame(
    coefficient = coefficients,
    contrast_id = seq_along(coefficients),
    stringsAsFactors = FALSE
  )

  contrast_lookup <- merge(
    contrast_lookup,
    model_info$sig_lookup,
    by = "coefficient",
    all.x = TRUE,
    sort = FALSE
  )

  out <- merge(
    out,
    contrast_lookup,
    by = "coefficient",
    all.x = TRUE,
    sort = FALSE
  )

  # ============================================================
  # 24. Preserve hierarchical decomposition
  # ============================================================

  if (!"inner_side" %in% names(out)) {
    stop(
      "Internal hierarchical decomposition information was unavailable.",
      call. = FALSE
    )
  }

  decomposition_columns <- c(
    "coefficient",
    "contrast_id",
    "significance",
    "p_value",
    "split_var",
    "split_side",
    "inner_side",
    "inner_variable_path",
    factor_vars,
    "n",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "configuration_explicit",
    "configuration_semantic"
  )

  decomposition_columns <- intersect(decomposition_columns,names(out))

  decomposition <- out[ ,decomposition_columns,drop = FALSE]

  decomposition$decomposition_path <- ifelse(
    is.na(decomposition$inner_side) |
      decomposition$inner_side == "" |
      decomposition$inner_side == "cell",
    decomposition$split_side,
    paste0(
      decomposition$split_side,
      "_",
      decomposition$inner_side
    )
  )

  decomposition$variable_path <- ifelse(
    is.na(decomposition$inner_variable_path) |
      decomposition$inner_variable_path == "",
    decomposition$split_var,
    paste0(
      decomposition$split_var,
      "_",
      decomposition$inner_variable_path
    )
  )

  variable_depth <- lengths(
    strsplit(
      decomposition$variable_path,
      "_",
      fixed = TRUE
    )
  )

  if (!all(variable_depth == decomposition$decomposition_depth)) {
    stop(
      "The variable hierarchy and decomposition hierarchy have inconsistent depths.",
      call. = FALSE
    )
  }

  decomposition$decomposition_depth <- lengths(
    strsplit(
      decomposition$decomposition_path,
      "_",
      fixed = TRUE
    )
  )

  rownames(decomposition) <- NULL

  # ============================================================
  # 24.1 Select and order final columns
  # ============================================================

  base_cols <- c(
    "coefficient", "contrast_id", "significance", "p_value",
    "split_var", "split_side",
    factor_vars,
    "n", "estimate",
    "se", "ci_low", "ci_high"
  )

  configuration_cols <- switch(
    configuration,
    none = character(0),
    all = c("configuration_explicit", "configuration_semantic"),
    explicit = "configuration_explicit",
    semantic = "configuration_semantic"
  )

  out <- out[, c(base_cols, configuration_cols)]
  rownames(out) <- NULL

  # ============================================================
  # 25. ESTIMATE INFORMATION
  # ============================================================

  if (estimate_type == "raw") {

    estimate_info <- paste(
      "Displayed estimate:",
      "observed/raw aggregated means."
    )

  } else {

    family_object <- .extract_family(model)

    if (prediction_scale == "response") {

      if (!is.null(family_object) &&
          identical(family_object$family, "binomial")) {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted probability on the response scale."
        )

      } else if (
        !is.null(family_object) &&
        family_object$family %in% c(
          "poisson",
          "Negative Binomial(1)",
          "Negative Binomial(2)"
        )
      ) {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted expected count on the response scale."
        )

      } else {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted value on the response scale."
        )
      }

    } else {

      estimate_info <- paste(
        "Displayed estimate:",
        "average model-predicted value on the linear-predictor/link scale."
      )
    }
  }

  # ============================================================
  # 26. Store metadata as attributes
  # ============================================================
  #
  # These attributes are useful for reproducibility and downstream
  # plotting functions.
  # ============================================================

  attr(out, "note") <- disclaimer_note
  attr(out, "estimate_type") <- estimate_type
  attr(out, "prediction_scale") <- if (estimate_type == "raw") NA_character_ else prediction_scale
  attr(out, "estimate_info") <- estimate_info
  attr(out, "pvalue_source") <- model_info$pvalue_source
  attr(out, "split_rule") <- split_rule
  attr(out, "top_split") <- top_split
  attr(out, "interaction") <- interaction
  attr(out, "decomposition") <- decomposition

  class(out) <- c("plot_ready_contrasts_df", class(out))

  return(out)
}
