#' Export a MoPlot figure using publication-ready settings
#'
#' @description
#' `MoPlot_export()` exports a plot created by `MoPlot_multpred()` using
#' dimensions recommended from the graphical complexity of the plot.
#'
#' PNG and TIFF files are exported at high resolution. TIFF files use
#' lossless LZW compression. PDF files are exported through the Cairo device
#' to improve rendering of Unicode symbols used in MoPlot figures.
#'
#' @param plot A plot object created by `MoPlot_multpred()`.
#'
#' @param filename Character string specifying the output filename. Supported
#'   extensions are `.png`, `.tif`, `.tiff`, and `.pdf`.
#'
#' @param width Optional export width in inches. If `NULL`, the recommended
#'   width stored in `plot` is used.
#'
#' @param height Optional export height in inches. If `NULL`, the recommended
#'   height stored in `plot` is used.
#'
#' @param dpi Resolution in dots per inch for raster formats. The default is
#'   `600`. This argument is ignored for PDF output.
#'
#' @param bg Background colour. The default is `"white"`.
#'
#' @return Invisibly returns the normalized path of the exported file.
#'
#' @export
MoPlot_export <- function(
    plot,
    filename,
    width = NULL,
    height = NULL,
    dpi = 600,
    bg = "white"
) {

  # ==========================================================================
  # 1. VALIDATION
  # ==========================================================================

  if (
    !is.character(filename) ||
    length(filename) != 1L ||
    is.na(filename) ||
    !nzchar(filename)
  ) {
    stop(
      "`filename` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  extension <- tolower(tools::file_ext(filename))

  supported_extensions <- c(
    "png",
    "tif",
    "tiff",
    "pdf"
  )

  if (
    !extension %in%
    supported_extensions
  ) {
    stop(
      paste0(
        "Unsupported file format '.",
        extension,
        "'. Supported formats are PNG, TIFF, and PDF."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 2. RECOMMENDED DIMENSIONS
  # ==========================================================================

  recommended_size <- attr(plot,"moplot_export_size")

  if (
    is.null(width) ||
    is.null(height)
  ) {

    if (
      is.null(recommended_size) ||
      !all(
        c(
          "width",
          "height"
        ) %in%
        names(
          recommended_size
        )
      )
    ) {
      stop(
        paste0(
          "Recommended export dimensions were not found in `plot`. ",
          "Supply both `width` and `height` explicitly."
        ),
        call. = FALSE
      )
    }

    if (is.null(width)) {
      width <- recommended_size[["width"]]
    }

    if (is.null(height)) {
      height <- recommended_size[["height"]]
    }
  }

  if (
    !is.numeric(width) ||
    length(width) != 1L ||
    is.na(width) ||
    !is.finite(width) ||
    width <= 0
  ) {
    stop(
      "`width` must be a positive finite number.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(height) ||
    length(height) != 1L ||
    is.na(height) ||
    !is.finite(height) ||
    height <= 0
  ) {
    stop(
      "`height` must be a positive finite number.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 3. EXPORT
  # ==========================================================================

  if (
    extension ==
    "pdf"
  ) {

    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = bg
    )

  } else if (
    extension %in%
    c(
      "tif",
      "tiff"
    )
  ) {

    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "tiff",
      compression = "lzw",
      bg = bg
    )

  } else {

    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "png",
      bg = bg
    )
  }

  message(
    paste0(
      "MoPlot exported to '",
      filename,
      "' at ",
      format(
        width,
        trim = TRUE
      ),
      " x ",
      format(
        height,
        trim = TRUE
      ),
      " inches",
      if (
        extension ==
        "pdf"
      ) {
        "."
      } else {
        paste0(
          " and ",
          dpi,
          " dpi."
        )
      }
    )
  )

  invisible(
    normalizePath(
      filename,
      mustWork = FALSE
    )
  )
}
