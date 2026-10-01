#' Plot contrasts for one categorical predictor
#'
#' @description
#' `MoPlot_onepred()` creates a graphical representation of the contrasts
#' associated with one categorical predictor in a fitted model.
#'
#' The upper panel displays the quantities compared by each contrast, together
#' with their confidence intervals. The optional lower panel displays the
#' corresponding model parameter, its confidence interval, and an approximate
#' Cohen's d.
#'
#' The horizontal position of each parameter in the lower panel is inherited
#' from the midpoint of the corresponding contrast in the upper panel.
#'
#' @param model A fitted model supported by the package.
#'
#' @param predictor A single character string identifying the categorical
#'   predictor to display.
#'
#' @param siglevel A numeric value between 0 and 1 defining the significance
#'   threshold used for the contrast labels. The default is `0.05`.
#'
#' @param level A numeric value between 0 and 1 specifying the confidence
#'   level. The default is `0.95`.
#'
#' @param estimate_type Character string specifying how the quantities in the
#'   upper panel should be computed. Available options are:
#'
#'   \describe{
#'     \item{`"auto"`}{
#'       Uses observed means when the fixed part of the model contains only
#'       the focal predictor, the model has no random effects, and, for
#'       generalized models, `prediction_scale = "response"`. Model-based
#'       predictions are used otherwise: with other predictors or
#'       interactions, with random effects (observed means would ignore the
#'       clustering of the observations), and for estimates on the link scale.
#'     }
#'     \item{`"observed"`}{
#'       Uses observed means computed directly from the model data.
#'     }
#'     \item{`"predicted"`}{
#'       Uses model-based fixed-effect predictions.
#'     }
#'   }
#'
#' @param prediction_scale Character string specifying the scale used for
#'   model-based predictions. Available options are `"response"` and `"link"`.
#'   The default is `"response"`.
#'
#' @param coefficient_scale Character string specifying the estimate displayed
#'   in the lower panel. Available options are:
#'
#'   \describe{
#'     \item{`"model"`}{
#'       Displays the model coefficient on its original model scale.
#'     }
#'     \item{`"interpreted"`}{
#'       Displays the interpreted estimate returned by `extract_param()`,
#'       such as an odds ratio or incidence rate ratio when appropriate.
#'     }
#'   }
#'
#' @param per_difference Logical. Should the lower-panel coefficients be
#'   divided by their scale factor? The default is `FALSE`.
#'
#'   For some codings, such as Helmert and reverse Helmert contrasts or
#'   custom matrices, each model coefficient equals a scale factor `s`
#'   multiplied by the difference displayed in the upper panel (for example
#'   `s = 1/3` for the second Helmert contrast). When `per_difference = TRUE`,
#'   the lower panel reports `beta / s` and its confidence interval, so that
#'   both panels refer to the same difference; with
#'   `coefficient_scale = "interpreted"`, ratios are computed as
#'   `exp(beta / s)`, the ratio for the displayed comparison. Test statistics,
#'   p-values, and Cohen's d are unaffected. For generalized models, `beta / s`
#'   is a difference on the link scale. Codings with `s = 1` for every
#'   contrast (treatment, sum, sliding differences) are unaffected.
#'
#' @param coef Logical. Should the lower coefficient panel be displayed?
#'   The default is `TRUE`.
#'
#' @param show_data Logical. Should raw observations be displayed in the upper
#'   panel using horizontal jittering? The default is `FALSE`.
#'
#' @param style Character string specifying the graphical style. Available
#'   options are:
#'
#'   \describe{
#'     \item{`"minimal"`}{
#'       Uses a restrained grayscale palette.
#'     }
#'     \item{`"publication"`}{
#'       Uses a colorblind-safe palette and stronger visual emphasis.
#'     }
#'   }
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
#'   The caption is not added to the graphical object. Each sentence is printed
#'   on a separate console line so that the text can be copied and used
#'   independently of the plot.
#'
#' @param digits Number of decimal places used in the lower-panel labels.
#'   The default is `2`.
#'
#' @param base_size Base font size used in both panels. The default is `13`.
#'
#' @param panel_heights A numeric vector of length two specifying the relative
#'   heights of the upper and lower panels. The default is `c(2, 1)`.
#'   This argument is used only when `coef = TRUE`.
#'
#' @details
#' The function displays one focal categorical predictor at a time, but the
#' fitted model may contain additional predictors.
#'
#' Models in which the focal predictor is involved in an interaction are not
#' currently supported.
#'
#' Treatment and Sum contrasts are displayed as independent pairs to improve
#' readability. Other supported contrast types reuse repeated displayed
#' quantities where appropriate.
#'
#' Polynomial contrasts are not currently supported by this plotting function.
#'
#' No model-intercept line is displayed.
#'
#' @return
#' When `coef = FALSE`, a `ggplot` object containing the upper panel.
#'
#' When `coef = TRUE`, a `patchwork` object containing the aligned upper and
#' lower panels.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' dat <- data.frame(
#'   group = factor(
#'     rep(c("A", "B", "C"), each = 40),
#'     levels = c("A", "B", "C")
#'   ),
#'   y = c(
#'     rnorm(40, 0.0, 0.9),
#'     rnorm(40, 0.4, 0.9),
#'     rnorm(40, 0.9, 0.9)
#'   )
#' )
#'
#' contrasts(dat$group) <- contr.helmert(3)
#'
#' mod <- lm(
#'   y ~ group,
#'   data = dat
#' )
#'
#' MoPlot_onepred(
#'   model = mod,
#'   predictor = "group"
#' )
#'
#' MoPlot_onepred(
#'   model = mod,
#'   predictor = "group",
#'   show_data = TRUE,
#'   style = "publication"
#' )
#'
#' MoPlot_onepred(
#'   model = mod,
#'   predictor = "group",
#'   coef = FALSE
#' )
#' }
#'
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
MoPlot_onepred <- function(
    model,
    predictor,
    siglevel = 0.05,
    level = 0.95,
    estimate_type = c("auto", "observed", "predicted"),
    prediction_scale = c("response", "link"),
    coefficient_scale = c("model", "interpreted"),
    per_difference = FALSE,
    coef = TRUE,
    show_data = FALSE,
    style = c("minimal", "publication"),
    caption = TRUE,
    digits = 2L,
    base_size = 13,
    panel_heights = c(2, 1)
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING
  # ==========================================================================

  estimate_type <- match.arg(estimate_type)

  prediction_scale <- match.arg(prediction_scale)

  coefficient_scale <- match.arg(coefficient_scale)

  style <- match.arg(style)

  # ==========================================================================
  # 2. ARGUMENT VALIDATION
  # ==========================================================================

  if (
    !is.logical(per_difference) ||
    length(per_difference) != 1L ||
    is.na(per_difference)
  ) {
    stop(
      "`per_difference` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.logical(coef) ||
    length(coef) != 1L ||
    is.na(coef)
  ) {
    stop(
      "`coef` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.logical(show_data) ||
    length(show_data) != 1L ||
    is.na(show_data)
  ) {
    stop(
      "`show_data` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(digits) ||
    length(digits) != 1L ||
    is.na(digits) ||
    !is.finite(digits) ||
    digits < 0 ||
    digits != as.integer(digits)
  ) {
    stop(
      "`digits` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  digits <- as.integer(
    digits
  )

  if (
    !is.numeric(base_size) ||
    length(base_size) != 1L ||
    is.na(base_size) ||
    !is.finite(base_size) ||
    base_size <= 0
  ) {
    stop(
      "`base_size` must be a single positive numeric value.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(panel_heights) ||
    length(panel_heights) != 2L ||
    anyNA(panel_heights) ||
    any(!is.finite(panel_heights)) ||
    any(panel_heights <= 0)
  ) {
    stop(
      paste0(
        "`panel_heights` must contain two finite positive numeric values."
      ),
      call. = FALSE
    )
  }

  caption_is_valid <-
    isTRUE(caption) ||
    identical(caption, FALSE) ||
    is.null(caption) ||
    (
      is.character(caption) &&
        length(caption) == 1L &&
        !is.na(caption)
    )

  if (!caption_is_valid) {
    stop(
      paste0(
        "`caption` must be TRUE, FALSE, NULL, or a single character string."
      ),
      call. = FALSE
    )
  }

  # Raw observations and link-scale predictions are not directly comparable.
  if (
    isTRUE(show_data) &&
    identical(estimate_type, "predicted") &&
    identical(prediction_scale, "link") &&
    .is_generalized_model(model)
  ) {
    stop(
      paste0(
        "Raw observations cannot be displayed together with predictions ",
        "on the link scale. Use `prediction_scale = \"response\"` or ",
        "`show_data = FALSE`."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 3. PREPARE DATA
  # ==========================================================================

  prep <- .prepare_moplot_onepred(
    model = model,
    predictor = predictor,
    siglevel = siglevel,
    level = level,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    coefficient_scale = coefficient_scale,
    per_difference = per_difference
  )

  factor_levels <- levels(.extract_model_frame(model)[[predictor]])

  label_result <- .abbreviate_moplot_labels(
    labels = factor_levels,
    max_chars = 5L
  )

  label_key <- label_result$key

  # The automatic estimate selection may resolve to predicted estimates.
  if (
    isTRUE(show_data) &&
    identical(prep$estimate_type, "predicted") &&
    identical(prediction_scale, "link") &&
    .is_generalized_model(model)
  ) {
    stop(
      paste0(
        "Raw observations cannot be displayed together with predictions ",
        "on the link scale. Use `prediction_scale = \"response\"` or ",
        "`show_data = FALSE`."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 4. CAPTION
  # ==========================================================================

  caption_text <- NULL

  if (isTRUE(caption)) {

    caption_text <- MoPlotStele(
      model = model,
      predictor = predictor,
      contrast_output = prep$contrast_data,
      siglevel = siglevel,
      show_data = show_data,
      coef = coef,
      coefficient_scale = coefficient_scale,
      per_difference = per_difference,
      family = prep$family,
      link = prep$link,
      estimate_type = prep$estimate_type,
      prediction_scale = prep$prediction_scale,
      estimate_label = unique(
        prep$parameter_data$Estimate_type
      )[1L],
      label_key=label_key
    )

  } else if (
    is.character(caption) &&
    length(caption) == 1L
  ) {

    caption_text <- caption
  }

  if (!is.null(caption_text)) {

    caption_text <- trimws(caption_text)

    # Print each sentence on a separate line.
    # Decimal points, such as those in 0.050, are not modified.
    caption_console <- gsub(
      pattern = "\\.\\s+(?=[A-Z])",
      replacement = ".\n",
      x = caption_text,
      perl = TRUE
    )

    cat(
      "\nPlot caption:\n\n",
      caption_console,
      "\n\n",
      sep = ""
    )
  }

  # ==========================================================================
  # 5. UPPER PANEL
  # ==========================================================================

  upper_plot <- .plot_moplot_onepred_upper(
    prep = prep,
    show_data = show_data,
    style = style,
    base_size = base_size,
    label_key=label_key
  )

  # ==========================================================================
  # 6. RETURN UPPER PANEL ONLY
  # ==========================================================================

  if (!isTRUE(coef)) {
    return(
      upper_plot
    )
  }
  # ==========================================================================
  # 7. LOWER PANEL
  # ==========================================================================

  lower_plot <- .plot_moplot_onepred_lower(
    prep = prep,
    style = style,
    digits = digits,
    base_size = base_size
  )

  # ==========================================================================
  # 8. COMBINE PANELS
  # ==========================================================================

  combined_plot <- patchwork::wrap_plots(
    upper_plot,
    lower_plot,
    ncol = 1L,
    heights = panel_heights
  )


  combined_plot
}
