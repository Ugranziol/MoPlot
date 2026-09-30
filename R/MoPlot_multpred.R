#' Plot contrasts for interacting categorical predictors
#'
#' @description
#' `MoPlot_multpred()` creates a graphical representation of model
#' interactions involving two or more categorical predictors.
#'
#' The function decomposes each selected interaction coefficient into the
#' elementary and higher-order contrasts implied by the contrast matrices used
#' to fit the model. The upper panel displays the observed or model-predicted
#' quantities involved in this decomposition. The optional lower panel displays
#' the corresponding model coefficient, its confidence interval, and an
#' approximate Cohen's d.
#'
#' By default, the function identifies the interaction of highest order in the
#' fitted model and displays all coefficients associated with that interaction.
#' If more than one interaction of the highest order is present, the function
#' asks the user to select one explicitly through `interaction`.
#'
#' @param model A fitted model supported by the package. Supported model
#'   classes include objects returned by [stats::lm()], [stats::glm()],
#'   [lme4::lmer()], [lme4::glmer()], [nlme::lme()], and
#'   [glmmTMB::glmmTMB()].
#'
#' @param data Optional data frame containing the variables used to fit
#'   `model`. If `NULL`, the function attempts to recover the analysis data
#'   from the fitted model using `stats::model.frame()`. Supplying `data`
#'   explicitly remains recommended when the model does not retain a usable
#'   model frame or when a specific analysis data set should be used.
#'
#' @param interaction Optional character string identifying the interaction
#'   term to display, for example `"A:B:C"`.
#'
#'   If `NULL`, the function automatically selects the unique interaction of
#'   highest order in the fitted model. If several interactions share the
#'   highest order, the function stops and asks the user to specify
#'   `interaction`.
#'
#' @param coefficients Optional character vector identifying specific model
#'   coefficients to display.
#'
#'   If `NULL`, all coefficients associated with the selected `interaction`
#'   are displayed. If `interaction = NULL` and `coefficients` is supplied,
#'   the corresponding interaction is inferred automatically, provided that
#'   all selected coefficients belong to the same interaction.
#'
#' @param interaction_info Logical. Should the function print the available
#'   interaction terms and their corresponding coefficient selections?
#'   The default is `TRUE`.
#'
#'   The printed table shows syntax that can be copied directly into
#'   `interaction` and `coefficients` in subsequent calls.
#'
#' @param split_rule Character string specifying how the top-level split
#'   variable is selected during the hierarchical decomposition. Available
#'   options are:
#'
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
#' @param top_split Optional character string identifying the predictor to
#'   force as the top-level split variable. When supplied, `top_split`
#'   overrides `split_rule`.
#'
#' @param siglevel Numeric value between 0 and 1 defining the significance
#'   threshold used for coefficient labels. The default is `0.05`.
#'
#' @param level Numeric value between 0 and 1 specifying the confidence level
#'   used for displayed confidence intervals. The default is `0.95`.
#'
#' @param estimate_type Character string specifying how quantities in the
#'   upper panel are computed. Available options are:
#'
#'   \describe{
#'     \item{`"raw"`}{
#'       Displays observed/raw aggregated means computed directly from `data`.
#'     }
#'     \item{`"predicted"`}{
#'       Displays model-based predictions computed from the fixed effects and
#'       their variance-covariance matrix.
#'     }
#'   }
#'
#'   The default is `"raw"`.
#'
#' @param prediction_scale Character string specifying the scale used for
#'   model-based predictions when `estimate_type = "predicted"`.
#'   Available options are:
#'
#'   \describe{
#'     \item{`"response"`}{
#'       Displays predictions on the response scale.
#'     }
#'     \item{`"link"`}{
#'       Displays predictions on the linear-predictor/link scale.
#'     }
#'   }
#'
#'   The default is `"response"`. This argument is ignored when
#'   `estimate_type = "raw"`.
#'
#' @param coefficient_scale Character string specifying the quantity displayed
#'   in the lower panel. Available options are:
#'
#'   \describe{
#'     \item{`"model"`}{
#'       Displays the coefficient on its original model scale.
#'     }
#'     \item{`"interpreted"`}{
#'       Displays the interpreted coefficient when an appropriate
#'       transformation is available, such as an odds ratio or incidence rate
#'       ratio.
#'     }
#'   }
#'
#'   The default is `"model"`.
#'
#' @param coef Logical. Should the lower coefficient panel be displayed?
#'   The default is `TRUE`.
#'
#' @param style Character string specifying the graphical representation of the
#'   interaction decomposition. Available options are:
#'
#'   \describe{
#'     \item{`"minimal"`}{
#'       Displays the observed or predicted quantities and the dashed segments
#'       representing the elementary component contrasts, without displaying
#'       intermediate delta formulas.
#'     }
#'     \item{`"tree"`}{
#'       Displays the decomposition as a hierarchical tree showing how
#'       elementary contrasts are recursively combined into higher-order
#'       contrasts and ultimately into the final interaction contrast.
#'     }
#'   }
#'
#'   The default is `"minimal"`.
#'
#' @param caption Controls whether a descriptive caption is printed in the
#'   console. It can be:
#'
#'   \itemize{
#'     \item `TRUE`, to generate and print an automatic caption;
#'     \item `FALSE` or `NULL`, to suppress the caption;
#'     \item a single character string to print a custom caption.
#'   }
#'
#'   The automatic caption describes the selected interaction, the graphical
#'   hierarchy, the contrast coding assigned to each interacting variable, the
#'   displayed estimate type, the selected graphical style, and the lower
#'   coefficient panel when present.
#'
#'   The caption is printed in the console and is not added to the graphical
#'   object.
#'
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in lower-panel estimate and effect-size labels. The default is `2`.
#'
#' @param base_size Positive numeric value defining the base font size used in
#'   both panels. The default is `13`.
#'
#' @param panel_heights Numeric vector of length two specifying the relative
#'   heights of the upper and lower panels. The default is `c(2, 1)`.
#'   This argument is used only when `coef = TRUE`.
#'
#' @details
#' `MoPlot_multpred()` is a graphical interface to the hierarchical
#' decomposition performed by [ContrastMultPred()].
#'
#' Each selected interaction coefficient is decomposed according to the
#' contrast matrices assigned to the categorical predictors involved in the
#' interaction. This allows treatment, sum, Helmert, reverse Helmert,
#' sliding-difference, and interpretable custom contrasts to be represented
#' directly in terms of the factor levels or level aggregates implied by the
#' fitted model.
#'
#' The graphical hierarchy is determined by `split_rule` or `top_split`.
#' The top-level split variable defines the outermost grouping structure,
#' whereas successive predictors define progressively lower levels of the
#' decomposition. The final predictor in the resulting hierarchy defines the
#' quantities represented on the x-axis.
#'
#' When an interaction contains more than two predictors, the decomposition is
#' recursive. Elementary contrasts are first constructed from pairs of
#' displayed quantities and are then combined into progressively higher-order
#' contrasts until the complete model interaction contrast is obtained.
#'
#' For `estimate_type = "predicted"`, predictions are population-level
#' predictions based on the fixed effects. For generalized models on the
#' response scale, uncertainty is propagated using the fixed-effect
#' variance-covariance matrix and the multivariate delta method where
#' appropriate. A methodological note describing only the settings relevant to
#' the current call is printed when the function is run.
#'
#' @section Interaction and coefficient selection:
#' The arguments `interaction` and `coefficients` provide two complementary
#' levels of control.
#'
#' To display all coefficients associated with one interaction, specify only
#' `interaction`, for example:
#'
#' \preformatted{
#' interaction = "A:B:C"
#' }
#'
#' To display only selected coefficients from that interaction, specify both:
#'
#' \preformatted{
#' interaction = "A:B:C"
#' coefficients = c("A1:B2:C2", "A2:B2:C2")
#' }
#'
#' If neither argument is supplied, the unique highest-order interaction is
#' selected automatically.
#'
#' When `interaction_info = TRUE`, the function prints the interactions found
#' in the model together with valid coefficient selections, providing
#' copy-ready syntax for subsequent calls.
#'
#' @section Graphical styles:
#' The three graphical styles differ in how much of the hierarchical
#' decomposition is shown.
#'
#' `"minimal"` emphasizes the observed or predicted quantities and elementary
#' contrasts.
#'
#' `"tree"` represents the same hierarchy through branching segments,
#' emphasizing the recursive structure of the decomposition.
#'
#' The lower coefficient panel, when requested, remains horizontally aligned
#' with the corresponding contrast block in the upper panel for all three
#' styles.
#'
#' @return
#' When `coef = FALSE`, a `ggplot` object containing the upper interaction
#' decomposition panel.
#'
#' When `coef = TRUE`, a `patchwork` object containing the aligned upper and
#' lower panels.
#'
#' @examples
#' \dontrun{
#'
#' set.seed(123)
#'
#' dat <- expand.grid(
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2")),
#'   C = factor(c("C1", "C2")),
#'   rep = seq_len(20)
#' )
#'
#' contrasts(dat$A) <- reverse_helmert(3)
#' contrasts(dat$B) <- MASS::contr.sdif(2)
#' contrasts(dat$C) <- stats::contr.treatment(2)
#'
#' dat$y <- with(
#'   dat,
#'   10 +
#'     2 * (A == "A2") +
#'     4 * (A == "A3") +
#'     2 * (B == "B2") +
#'     2 * (C == "C2") +
#'     2 * (A == "A2" & B == "B2" & C == "C2") +
#'     rnorm(nrow(dat))
#' )
#'
#' mod <- stats::lm(
#'   y ~ A * B * C,
#'   data = dat
#' )
#'
#' # Automatically display the unique highest-order interaction.
#' MoPlot_multpred(
#'   model = mod,
#'   data = dat
#' )
#'
#' # Display all coefficients from a specific interaction.
#' MoPlot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   style = "tree"
#' )
#'
#' # Display one specific coefficient using the minimal representation.
#' MoPlot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   coefficients = "A1:B2-1:C2",
#'   style = "minimal"
#' )
#'
#' # Use model-based predictions instead of raw means.
#' MoPlot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   estimate_type = "predicted",
#'   prediction_scale = "response"
#' )
#' }
#'
#' @importFrom patchwork wrap_plots
#' @export
MoPlot_multpred <- function(
    model,
    data=NULL,
    interaction = NULL,
    coefficients = NULL,
    interaction_info = TRUE,
    split_rule = c(
      "last_variable",
      "first_variable",
      "lowest_complexity",
      "highest_complexity"
    ),
    top_split = NULL,
    siglevel = 0.05,
    level = 0.95,
    estimate_type = c("raw", "predicted"),
    prediction_scale = c("response", "link"),
    coefficient_scale = c("model", "interpreted"),
    coef = TRUE,
    style = c("minimal", "tree"),
    caption = TRUE,
    digits = 2L,
    base_size = 13,
    panel_heights = c(2, 1)
) {

  # ============================================================
  # 1. ARGUMENT MATCHING
  # ============================================================

  split_rule <- match.arg(split_rule)

  estimate_type <- match.arg(estimate_type)

  prediction_scale <- match.arg(prediction_scale)

  coefficient_scale <- match.arg(coefficient_scale)

  style <- match.arg(style)



  # ============================================================
  # 1.1 RESOLVE INTERACTION AND COEFFICIENT SELECTION
  # ============================================================

  selection <- .resolve_moplot_interaction_selection(
    model = model,
    interaction = interaction,
    coefficients = coefficients,
    interaction_info = interaction_info
  )

  interaction <- selection$interaction
  coefficients <- selection$coefficients

  # ============================================================
  # 2. COMPUTE INTERACTION CONTRASTS
  # ============================================================

  contrast_output <- ContrastMultPred(
    model = model,
    data = data,
    coefficients = coefficients,
    split_rule = split_rule,
    top_split = top_split,
    conf_level = level,
    alpha_sig = siglevel,
    configuration = "all",
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    methodological_note = FALSE
  )

  # ============================================================
  # 2.1 METHODOLOGICAL NOTE
  # ============================================================

  plot_methodological_note <- .build_multpred_methodological_note(
    model = model,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    split_rule = split_rule,
    top_split = top_split,
    configuration = NULL
  )

  message(plot_methodological_note)

  # ============================================================
  # 3. PREPARE UPPER-PANEL DATA
  # ============================================================

  prep <- .prepare_moplot_multpred_data(contrast_output)

  label_key <- .collect_moplot_label_key(prep = prep, max_chars = 5L)

  # ============================================================
  # 3.1 CAPTION
  # ============================================================

  if (isTRUE(caption)) {

    caption_text <- MoPlotStele_multpred(
      model = model,
      contrast_output = contrast_output,
      prep = prep,
      interaction = interaction,
      style = style,
      siglevel = siglevel,
      coef = coef,
      coefficient_scale = coefficient_scale,
      estimate_type = estimate_type,
      prediction_scale = prediction_scale,
      label_key=label_key
    )

    cat(
      "Plot caption:\n\n",
      paste0(
        strwrap(
          caption_text,
          width = 100
        ),
        collapse = "\n"
      ),
      "\n\n",
      sep = ""
    )

  } else if (
    is.character(caption) &&
    length(caption) == 1L &&
    !is.na(caption)
  ) {

    cat(
      "Plot caption:\n\n",
      paste0(
        strwrap(
          caption,
          width = 100
        ),
        collapse = "\n"
      ),
      "\n\n",
      sep = ""
    )
  }

  # ============================================================
  # 4. BUILD UPPER PANEL
  # ============================================================

  upper_plot <- .plot_moplot_multpred_upper(
    prep = prep,
    style = style,
    base_size = base_size,
    label_key=label_key
  )


  export_size <- .moplot_multpred_export_size(
    prep = prep,
    style = style,
    coef = coef
  )

  # ============================================================
  # 5. RETURN UPPER PANEL ONLY
  # ============================================================

  if (!isTRUE(coef)) {

    attr(
      upper_plot,
      "moplot_export_size"
    ) <- export_size

    return(
      upper_plot
    )
  }

  # ============================================================
  # 6. PREPARE MODEL PARAMETERS
  # ============================================================

  parameters <- .prepare_moplot_multpred_parameters(
    model = model,
    coefficient_data = prep$coefficient_data,
    level = level,
    coefficient_scale = coefficient_scale
  )

  # ============================================================
  # 7. BUILD LOWER PANEL
  # ============================================================

  lower_plot <- .plot_moplot_multpred_lower(
    prep = prep,
    parameters = parameters,
    style = style,
    digits = digits,
    base_size = base_size
  )

  # ============================================================
  # 8. ASSEMBLE PANELS
  # ============================================================

  final_plot <- .assemble_moplot_multpred(
    upper_plot = upper_plot,
    lower_plot = lower_plot,
    panel_heights = panel_heights
  )

  attr(
    final_plot,
    "moplot_export_size"
  ) <- export_size

  final_plot
}
