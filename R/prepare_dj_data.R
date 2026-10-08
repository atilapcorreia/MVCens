#' Prepare Quarterly Dow-Jones Data (1920-1934)
#'
#' @description
#' Provides quarterly Dow-Jones dividends and divisor observations from
#' 1920 through 1934 as a matrix-valued dataset. The observations are embedded
#' in the function; no download is required.
#'
#' @param output_dir Directory in which to save the processed dataset.
#'   Defaults to the current working directory (\code{getwd()}).
#' @param standardize Logical. If \code{TRUE} (default), also return an array
#'   standardized separately for each variable across all quarters and years.
#'
#' @details
#' The data consist of 15 annual matrices of dimension \eqn{4 \times 2}.
#' Rows identify quarters Q1-Q4; columns contain DJ dividends and DJ divisor.
#' The annual matrices are stacked into a \eqn{4 \times 2 \times 15} array,
#' with years ordered chronologically. Standardization uses the pooled mean
#' and sample standard deviation of each column variable.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{X_raw}}{Original data as a \eqn{4 \times 2 \times 15} array.}
#'   \item{\code{X_std}}{Variable-wise standardized array, or \code{NULL}
#'     when \code{standardize = FALSE}.}
#'   \item{\code{dj_data}}{Named list of the 15 original annual matrices.}
#' }
#'
#' @section Saved file:
#' Saves \code{dow_jones_1920_1934_processed.rds} in \code{output_dir},
#' overwriting an existing file of that name. Reload with \code{readRDS()}.
#'
#' @note Annual observations follow consecutive years and should not
#'   automatically be assumed independent in statistical applications.
#'
#' @seealso \code{\link{mv_fit}}, \code{\link[base]{readRDS}}
#'
#' @examples
#' \dontrun{
#' dj <- prepare_dj_data(output_dir = "~/DowJones")
#' dim(dj$X_raw)
#' dj$dj_data[["1920"]]
#' dim(dj$X_std)
#' }
#' @export
prepare_dj_data <- function(output_dir = getwd(), standardize = TRUE) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(trimws(output_dir))) {
    stop("`output_dir` must be a non-empty directory path.")
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("`standardize` must be TRUE or FALSE.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  # Each consecutive group of four values corresponds to one year.
  dividends <- c(
    31.97, 30.00, 30.00, 28.75,
    26.81, 23.81, 22.06, 22.06,
    21.13, 19.38, 19.38, 20.88,
    24.19, 26.94, 25.94, 25.94,
    30.00, 26.75, 27.20, 27.20,
    32.39, 30.70, 30.00, 27.75,
    37.13, 27.13, 29.88, 26.88,
    33.31, 29.31, 31.56, 28.81,
    31.38, 27.63, 30.63, 27.45,
    31.58, 27.08, 28.05, 28.75,
    30.05, 27.90, 27.50, 25.95,
    25.56, 22.43, 19.56, 19.93,
    15.54, 18.23, 16.51, 16.26,
    13.34, 12.94, 12.74, 12.49,
    13.90, 13.95, 14.10, 15.30
  )
  divisor <- c(
    20.00, 20.00, 20.00, 20.00,
    20.00, 20.00, 20.00, 20.00,
    20.00, 20.00, 20.00, 20.00,
    20.00, 20.00, 20.00, 20.00,
    20.00, 19.40, 19.40, 18.90,
    18.40, 18.40, 19.00, 19.00,
    17.42, 16.67, 16.67, 16.67,
    16.67, 16.67, 16.67, 16.67,
    16.67, 16.67, 16.17, 13.92,
    12.11, 10.77, 10.47, 10.47,
    10.47,  9.85, 10.38, 10.38,
    10.38, 10.38, 10.38, 10.38,
    15.46, 15.46, 15.46, 15.46,
    15.46, 15.46, 15.71, 15.71,
    15.71, 15.71, 15.74, 15.74
  )

  years <- as.character(1920:1934)
  quarters <- paste0("Q", 1:4)
  variables <- c("DJ dividends", "DJ divisor")
  stopifnot(length(dividends) == 60L, length(divisor) == 60L,
            all(is.finite(dividends)), all(is.finite(divisor)))

  X_raw <- array(
    NA_real_, dim = c(4L, 2L, 15L),
    dimnames = list(quarter = quarters, variable = variables, year = years)
  )
  for (i in seq_along(years)) {
    idx <- ((i - 1L) * 4L + 1L):(i * 4L)
    X_raw[, , i] <- cbind(dividends[idx], divisor[idx])
  }
  dj_data <- stats::setNames(
    lapply(seq_along(years), function(i) X_raw[, , i]), years
  )

  X_std <- NULL
  if (standardize) {
    center <- c(mean(dividends), mean(divisor))
    scale <- c(stats::sd(dividends), stats::sd(divisor))
    if (any(!is.finite(scale) | scale <= 0)) {
      stop("Variable-specific standard deviations must be finite and positive.")
    }
    X_std <- sweep(X_raw, 2L, center, FUN = "-")
    X_std <- sweep(X_std, 2L, scale, FUN = "/")
  }

  result <- list(X_raw = X_raw, X_std = X_std, dj_data = dj_data)
  out_file <- file.path(output_dir, "dow_jones_1920_1934_processed.rds")
  saveRDS(result, out_file)
  message("Dow-Jones data saved to: ", out_file)
  result
}
