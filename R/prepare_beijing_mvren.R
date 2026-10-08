# Column names used in dplyr data-masking expressions. Declaring these names
# prevents spurious R CMD check NOTES without modifying the calculations.
utils::globalVariables(c(
  "year", "month", "day", "period", "station", "hour",
  "pollutant", "concentration", "mean_concentration",
  "all_12_cells_observed", "obs_id", "center"
))

#' Prepare Beijing Air Quality Data for Matrix-Variate Models
#'
#' @description
#' Download the Beijing Multi-Site Air Quality dataset from the UCI Machine
#' Learning Repository and organize hourly measurements as a sample of
#' \eqn{3 \times 4} matrices. The function returns both original-scale and
#' standardized observations for matrix-variate analysis.
#'
#' @param output_dir Directory for the downloaded source files and processed
#'   results. Defaults to \code{getwd()}. The directory is created if needed.
#'   Supply a single, non-empty character string identifying a writable path.
#'
#' @details
#' Each observation represents one monitoring station on one calendar day.
#' Matrix rows correspond to PM2.5, NO2, and O3, while columns correspond to
#' four six-hour periods: 00:00--05:59, 06:00--11:59, 12:00--17:59, and
#' 18:00--23:59. Observations are arranged by station and date.
#'
#' @section Data preparation:
#' \enumerate{
#'   \item \strong{Six-hour aggregation.} For each station, day, pollutant,
#'     and period, compute the mean concentration only when at least four
#'     hourly measurements are available.
#'   \item \strong{Completeness filter.} Retain a station-day only when all
#'     12 pollutant-by-period means are available. Each retained station-day
#'     becomes one matrix in the array.
#'   \item \strong{Pollutant-wise standardization.} For each pollutant,
#'     calculate the mean and sample standard deviation across all periods
#'     and retained station-days. Use these values to standardize that
#'     pollutant's entries throughout the sample.
#' }
#' Thus, \code{X_std} is centered and scaled separately for each pollutant,
#' pooling its entries over all periods and observations.
#'
#' @return
#' A named list containing:
#' \describe{
#'   \item{\code{X_raw}}{Numeric \eqn{3 \times 4 \times n} array of
#'     six-hour mean concentrations on the original scale.}
#'   \item{\code{X_std}}{Numeric array of the same dimensions, with
#'     pollutant-wise standardized concentrations. Intended for model fitting.}
#'   \item{\code{obs_info}}{Table linking each array slice to a station,
#'     date, observation index, and identifier.}
#'   \item{\code{block_data}}{Long-form table of accepted six-hour means,
#'     observed hourly counts, and retained station-day identifiers.}
#'   \item{\code{block_data_std}}{The same long-form data, including
#'     standardization constants and standardized concentrations.}
#'   \item{\code{scaling}}{Pollutant-specific means (\code{center}) and
#'     sample standard deviations (\code{scale}).}
#'   \item{\code{monthly_means}}{Monthly pollutant means computed from
#'     accepted six-hour blocks on the original concentration scale.}
#'   \item{\code{pollutants}}{Names of the array rows.}
#'   \item{\code{periods}}{Names of the array columns.}
#' }
#'
#' @section Download and saved files:
#' Source archives and extracted CSV files are stored under the
#' \code{beijing_air_quality_uci} subdirectory of \code{output_dir}.
#' The archive is downloaded if it is not already present; extraction and
#' preprocessing occur at every call. An initial download requires internet
#' access.
#'
#' Two processed files are saved directly in \code{output_dir}, overwriting
#' existing files with the same names:
#' \describe{
#'   \item{\code{beijing_mvren_processed.rds}}{The complete returned list;
#'     restore it using \code{readRDS()}.}
#'   \item{\code{beijing_mvren_processed.RData}}{Separate objects
#'     \code{X_raw}, \code{X_std}, \code{obs_info}, \code{complete_blocks},
#'     \code{complete_blocks_std}, \code{scaling}, and
#'     \code{monthly_means}; restore them using \code{load()}.}
#' }
#'
#' @section Notes:
#' The function only prepares data; model estimation must be performed
#' separately, for example with \code{mv_fit()}. The returned arrays contain
#' no missing values after filtering. A diagnostic compares the retained
#' sample size with the reference value of 16,133 and reports a warning if
#' they differ. Successful download and extraction depend on the external
#' repository and its archive structure.
#'
#' @seealso \code{\link{mv_fit}}, \code{\link[base]{readRDS}}
#'
#' @examples
#' \dontrun{
#' # Download and process the data
#' beijing <- prepare_beijing_mvren(output_dir = "~/Beijing")
#'
#' # Access the standardized array and observation metadata
#' X <- beijing$X_std
#' dim(X)
#' head(beijing$obs_info)
#' beijing$scaling
#'
#' # Optionally fit a matrix-variate model
#' fit <- mv_fit(model = "MVREN", X = X,
#'               precision = 1e-6, max_iter = 200)
#'
#' # Read saved results without preprocessing again
#' saved <- readRDS("~/Beijing/beijing_mvren_processed.rds")
#' }
#'
#' @export
prepare_beijing_mvren <- function(output_dir = getwd()) {
  # ---- Settings ----------------------------------------------------------------
  uci_url <- paste0(
    "https://archive.ics.uci.edu/static/public/501/",
    "beijing%2Bmulti%2Bsite%2Bair%2Bquality%2Bdata.zip"
  )

  # Raw UCI download/extraction files are kept in a subfolder of output_dir.
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be a non-empty directory path.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  message("Output directory: ", output_dir)

  data_dir <- file.path(output_dir, "beijing_air_quality_uci")
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  outer_zip <- file.path(data_dir, "beijing_multi_site_air_quality_data.zip")

  pollutants <- c("PM2.5", "NO2", "O3")
  periods <- c("00-05", "06-11", "12-17", "18-23")

  # ---- Download and extract -----------------------------------------------------
  if (!file.exists(outer_zip)) {
    message("Downloading the Beijing Multi-Site Air Quality dataset from UCI...")
    utils::download.file(uci_url, destfile = outer_zip, mode = "wb", quiet = FALSE)
  }

  extract_dir <- file.path(data_dir, "extracted")
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)

  # The current UCI download is an outer archive that contains the original
  # PRSA2017_Data_20130301-20170228.zip archive. The code also works if the
  # station CSV files are exposed directly in the outer archive.
  utils::unzip(outer_zip, exdir = extract_dir)

  nested_zip <- list.files(
    extract_dir,
    pattern = "^PRSA2017_Data_20130301-20170228\\.zip$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(nested_zip) > 0L) {
    utils::unzip(nested_zip[1L], exdir = extract_dir)
  }

  station_files <- list.files(
    extract_dir,
    pattern = "^PRSA_Data_.*\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(station_files) == 0L) {
    stop(
      "No station-level PRSA_Data_*.csv files were found after extraction. ",
      "Please check the UCI archive structure."
    )
  }

  message("Station files found: ", length(station_files))

  # ---- Read the 12 station files ------------------------------------------------
  hourly <- station_files |>
    sort() |>
    purrr::map_dfr(
      ~ readr::read_csv(
        .x,
        na = c("NA", ""),
        show_col_types = FALSE,
        progress = FALSE
      )
    )

  required_columns <- c(
    "year", "month", "day", "hour", "station",
    "PM2.5", "NO2", "O3"
  )

  missing_columns <- setdiff(required_columns, names(hourly))
  if (length(missing_columns) > 0L) {
    stop(
      "The following required columns are missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Keep only variables needed for the construction used in the paper.
  hourly <- hourly |>
    dplyr::select(dplyr::all_of(required_columns)) |>
    dplyr::mutate(
      date = as.Date(sprintf("%04d-%02d-%02d", year, month, day)),
      period = dplyr::case_when(
        hour >= 0  & hour <= 5  ~ "00-05",
        hour >= 6  & hour <= 11 ~ "06-11",
        hour >= 12 & hour <= 17 ~ "12-17",
        hour >= 18 & hour <= 23 ~ "18-23",
        TRUE ~ NA_character_
      ),
      period = factor(period, levels = periods, ordered = TRUE)
    )

  if (anyNA(hourly$period)) {
    stop("Unexpected hour values were found outside 0,...,23.")
  }

  # ---- Six-hour aggregation -----------------------------------------------------
  # For each station-day, pollutant, and six-hour period, calculate the mean only
  # when at least four hourly pollutant measurements are observed.
  block_data <- hourly |>
    dplyr::select(station, date, hour, period, dplyr::all_of(pollutants)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(pollutants),
      names_to = "pollutant",
      values_to = "concentration"
    ) |>
    dplyr::mutate(
      pollutant = factor(pollutant, levels = pollutants, ordered = TRUE)
    ) |>
    dplyr::group_by(station, date, pollutant, period) |>
    dplyr::summarise(
      n_hourly_observed = sum(!is.na(concentration)),
      mean_concentration = if (
        sum(!is.na(concentration)) >= 4L
      ) mean(concentration, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  # ---- Completeness filter ------------------------------------------------------
  # Retain a station-day only if all 12 cells (3 pollutants x 4 periods) exist and
  # have an accepted six-hour mean.
  complete_station_days <- block_data |>
    dplyr::group_by(station, date) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      all_12_cells_observed = (dplyr::n() == 12L) && all(!is.na(mean_concentration)),
      .groups = "drop"
    ) |>
    dplyr::filter(all_12_cells_observed) |>
    dplyr::select(station, date)

  complete_blocks <- block_data |>
    dplyr::semi_join(complete_station_days, by = c("station", "date")) |>
    dplyr::mutate(
      pollutant = factor(pollutant, levels = pollutants, ordered = TRUE),
      period = factor(period, levels = periods, ordered = TRUE)
    ) |>
    dplyr::arrange(station, date, pollutant, period)

  # One observation = one retained station-day.
  obs_info <- complete_station_days |>
    dplyr::arrange(station, date) |>
    dplyr::mutate(
      obs_id = dplyr::row_number(),
      observation = paste(station, format(date, "%Y-%m-%d"), sep = "_")
    )

  complete_blocks <- complete_blocks |>
    dplyr::left_join(obs_info, by = c("station", "date")) |>
    dplyr::arrange(obs_id, pollutant, period)

  n_obs <- nrow(obs_info)

  # ---- Construct the p x q x n matrix-valued sample ----------------------------
  # Rows:    PM2.5, NO2, O3
  # Columns: 00-05, 06-11, 12-17, 18-23
  # Slices:  retained station-days
  X_raw <- array(
    NA_real_,
    dim = c(length(pollutants), length(periods), n_obs),
    dimnames = list(
      pollutant = pollutants,
      period = periods,
      observation = obs_info$observation
    )
  )

  array_index <- cbind(
    match(as.character(complete_blocks$pollutant), pollutants),
    match(as.character(complete_blocks$period), periods),
    complete_blocks$obs_id
  )

  X_raw[array_index] <- complete_blocks$mean_concentration

  stopifnot(
    identical(dim(X_raw)[1:2], c(3L, 4L)),
    !anyNA(X_raw)
  )

  # ---- Pollutant-wise global standardization -----------------------------------
  # The paper standardizes each pollutant row by its global mean and global SD.
  # Here "global" means across all four periods and all retained station-day
  # matrices for a given pollutant.
  pollutant_center <- vapply(
    seq_along(pollutants),
    function(r) mean(X_raw[r, , ]),
    numeric(1)
  )

  pollutant_scale <- vapply(
    seq_along(pollutants),
    function(r) stats::sd(as.vector(X_raw[r, , ])),
    numeric(1)
  )

  names(pollutant_center) <- pollutants
  names(pollutant_scale) <- pollutants

  if (any(!is.finite(pollutant_scale)) || any(pollutant_scale <= 0)) {
    stop("At least one pollutant has a non-positive or non-finite global SD.")
  }

  X_std <- sweep(X_raw, MARGIN = 1L, STATS = pollutant_center, FUN = "-")
  X_std <- sweep(X_std, MARGIN = 1L, STATS = pollutant_scale, FUN = "/")

  # Long-form standardized data, useful for plots and diagnostics.
  scaling <- tibble::tibble(
    pollutant = factor(pollutants, levels = pollutants, ordered = TRUE),
    center = unname(pollutant_center),
    scale = unname(pollutant_scale)
  )

  complete_blocks_std <- complete_blocks |>
    dplyr::left_join(scaling, by = "pollutant") |>
    dplyr::mutate(
      standardized_concentration = (mean_concentration - center) / scale
    )

  # Optional monthly summaries on the original concentration scale, after the
  # same completeness filter used to form the matrix-valued sample.
  monthly_means <- complete_blocks |>
    dplyr::mutate(month = as.Date(format(date, "%Y-%m-01"))) |>
    dplyr::group_by(month, pollutant) |>
    dplyr::summarise(
      mean_concentration = mean(mean_concentration),
      .groups = "drop"
    )

  # ---- Diagnostics --------------------------------------------------------------
  cat("\nData preparation summary\n")
  cat("------------------------\n")
  cat("Hourly rows read:             ", nrow(hourly), "\n", sep = "")
  cat("Monitoring stations:          ", dplyr::n_distinct(hourly$station), "\n", sep = "")
  cat("Retained station-days:        ", n_obs, "\n", sep = "")
  cat("Dimension of X_raw:           ", paste(dim(X_raw), collapse = " x "), "\n", sep = "")
  cat("Dimension of X_std:           ", paste(dim(X_std), collapse = " x "), "\n", sep = "")
  cat("Any missing values in X_raw?  ", anyNA(X_raw), "\n", sep = "")
  cat("Any missing values in X_std?  ", anyNA(X_std), "\n", sep = "")

  if (n_obs != 16133L) {
    warning(
      "The paper reports 16133 retained station-day matrices, but this run ",
      "produced ", n_obs, ". If the UCI data have not changed, inspect the ",
      "downloaded archive and preprocessing assumptions."
    )
  } else {
    message("The retained sample size matches the paper: n = 16133.")
  }

  # ---- Save processed objects ---------------------------------------------------
  processed <- list(
    X_raw = X_raw,
    X_std = X_std,
    obs_info = obs_info,
    block_data = complete_blocks,
    block_data_std = complete_blocks_std,
    scaling = scaling,
    monthly_means = monthly_means,
    pollutants = pollutants,
    periods = periods
  )

  rds_file <- file.path(output_dir, "beijing_mvren_processed.rds")
  rdata_file <- file.path(output_dir, "beijing_mvren_processed.RData")

  saveRDS(processed, file = rds_file)
  save(
    X_raw,
    X_std,
    obs_info,
    complete_blocks,
    complete_blocks_std,
    scaling,
    monthly_means,
    file = rdata_file
  )

  cat("\nObjects saved to:\n")
  cat("  ", rds_file, "\n", sep = "")
  cat("  ", rdata_file, "\n", sep = "")
  cat("\nRaw UCI files stored in:\n")
  cat("  ", data_dir, "\n", sep = "")
  cat("\nUse result$X_std as the 3 x 4 x n standardized array for model fitting.\n")
  cat("Example first raw matrix:\n")
  print(X_raw[, , 1L])
  cat("\nExample first standardized matrix:\n")
  print(X_std[, , 1L])

  return(processed)
}
