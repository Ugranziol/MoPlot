# ============================================================
# PRINT METHODS FOR MODEL CONTRAST TABLE OBJECTS (INTERNAL)
# ============================================================
#
# These S3 methods control how ModelContrastTable() results are
# displayed in the console internally within the package.
#
# The underlying objects remain ordinary lists and data frames.
# ============================================================


# ============================================================
# 1. COMPLETE MODEL CONTRAST TABLE
# ============================================================

#' Print a complete model contrast summary
#'
#' @param x An object of class `"model_contrast_table"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1 cli_h2
#' @importFrom tibble as_tibble
#' @noRd
print.model_contrast_table <- function(x, ...) {

  cli::cli_h1("Model Contrast Summary")

  cli::cli_h2("Model-level information")

  model_summary <- x[["Model-Level Summary"]]

  print(
    tibble::as_tibble(model_summary),
    n = Inf,
    width = Inf,
    ...
  )

  cli::cli_h2("Variable-level information")

  variable_summary <- x[["Variable-Level Summary"]]

  print(
    tibble::as_tibble(variable_summary),
    n = Inf,
    width = Inf,
    ...
  )

  invisible(x)
}


# ============================================================
# 2. MODEL-LEVEL SUMMARY
# ============================================================

#' Print a model-level contrast summary
#'
#' @param x An object of class `"model_contrast_summary"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1
#' @importFrom tibble as_tibble
#' @noRd
print.model_contrast_summary <- function(x, ...) {

  cli::cli_h1("Model-level summary")

  printable_object <- as.data.frame(x, stringsAsFactors = FALSE)

  print(
    tibble::as_tibble(printable_object),
    n = Inf,
    width = Inf,
    ...
  )

  invisible(x)
}


# ============================================================
# 3. VARIABLE-LEVEL SUMMARY
# ============================================================

#' Print a variable-level contrast summary
#'
#' @param x An object of class `"model_contrast_variables"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1
#' @importFrom tibble as_tibble
#' @noRd
print.model_contrast_variables <- function(x, ...) {

  cli::cli_h1("Variable-level summary")

  printable_object <- as.data.frame(x, stringsAsFactors = FALSE)

  print(
    tibble::as_tibble(printable_object),
    n = Inf,
    width = Inf,
    ...
  )

  invisible(x)
}
