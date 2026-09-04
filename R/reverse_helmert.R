#' Reverse Helmert Contrasts
#'
#' Creates a Reverse Helmert contrast matrix by reversing both the rows and
#' columns of a standard Helmert contrast matrix.
#'
#' @param n Number of factor levels.
#'
#' @return A numeric matrix with `n` rows and `n - 1` columns containing
#'   Reverse Helmert contrasts.
#'
#' @details
#' Reverse Helmert contrasts compare each factor level, except the last,
#' with the mean of the subsequent levels.
#'
#' The contrast matrix is constructed from [stats::contr.helmert()] by
#' reversing the row order and subsequently the column order. This preserves
#' the Helmert coding structure while reversing the direction in which the
#' sequential comparisons are defined.
#'
#' @examples
#' reverse_helmert(3)
#' reverse_helmert(4)
#'
#' dat <- data.frame(
#'   group = factor(rep(c("A", "B", "C"), each = 10)),
#'   y = rnorm(30)
#' )
#'
#' contrasts(dat$group) <- reverse_helmert(3)
#' model <- lm(y ~ group, data = dat)
#'
#' @seealso [stats::contr.helmert()]
#'
#' @export
reverse_helmert <- function(n) {
  tmp <- apply(stats::contr.helmert(n), 2L, rev)
  tmp <- tmp[, ncol(tmp):1L, drop = FALSE]
  rownames(tmp) <- seq_len(n)
  tmp
}
