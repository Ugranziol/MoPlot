# ============================================================
# PRINT METHODS FOR MODEL CONTRAST TABLE OBJECTS
# ============================================================
#
# These S3 methods control how ModelContrastTable() results are
# displayed in the console.
#
# The underlying objects remain ordinary lists and data frames.
# Therefore, users can still extract, transform, export, or test
# the returned values using standard R operations.
# ============================================================


# ============================================================
# 1. COMPLETE MODEL CONTRAST TABLE
# ============================================================

#' Print a complete model contrast summary
#'
#' @description
#' Prints the complete output returned by
#' [ModelContrastTable()] when `type = "all"`.
#'
#' @param x An object of class `"model_contrast_table"`.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.model_contrast_table <- function(x, ...) {

  cli::cli_h1(
    "Model Contrast Summary"
  )

  cli::cli_h2(
    "Model-level information"
  )

  model_summary <- x[["Model-Level Summary"]]

  print(
    tibble::as_tibble(
      model_summary
    ),
    n = Inf,
    width = Inf
  )

  cli::cli_h2(
    "Variable-level information"
  )

  variable_summary <- x[["Variable-Level Summary"]]

  print(
    tibble::as_tibble(
      variable_summary
    ),
    n = Inf,
    width = Inf
  )

  invisible(x)
}


# ============================================================
# 2. MODEL-LEVEL SUMMARY
# ============================================================

#' Print a model-level contrast summary
#'
#' @description
#' Prints the model-level summary returned by
#' [ModelContrastTable()] when `type = "summary"`.
#'
#' @param x An object of class `"model_contrast_summary"`.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.model_contrast_summary <- function(x, ...) {

  cli::cli_h1(
    "Model-level summary"
  )

  printable_object <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )

  print(
    tibble::as_tibble(
      printable_object
    ),
    n = Inf,
    width = Inf
  )

  invisible(x)
}


# ============================================================
# 3. VARIABLE-LEVEL SUMMARY
# ============================================================

#' Print a variable-level contrast summary
#'
#' @description
#' Prints a variable-level summary returned by
#' [ModelContrastTable()] when `type` is `"x&y"`, `"x"`, or `"y"`.
#'
#' @param x An object of class `"model_contrast_variables"`.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.model_contrast_variables <- function(x, ...) {

  cli::cli_h1(
    "Variable-level summary"
  )

  printable_object <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )

  print(
    tibble::as_tibble(
      printable_object
    ),
    n = Inf,
    width = Inf
  )

  invisible(x)
}
