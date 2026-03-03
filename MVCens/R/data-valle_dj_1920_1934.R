#' Dow Jones dividends and divisor (quarterly), 1920--1934
#'
#' Quarterly Dow Jones (DJ) series used in the Valle application example.
#' The dataset contains two variables observed for each quarter from 1920 to 1934:
#' DJ dividends and the DJ divisor.
#'
#' @details
#' The object is a 3D array indexed as \code{[Quarter, Variable, Year]}.
#' Use \code{dimnames(valle_dj_1920_1934)} to see the labels.
#'
#' @format A numeric array with dimensions \code{4 x 2 x 15} and dimnames:
#' \describe{
#'   \item{Quarter}{\code{"Q1"}, \code{"Q2"}, \code{"Q3"}, \code{"Q4"}.}
#'   \item{Variable}{\code{"DJ dividends"} and \code{"DJ divisor"}.}
#'   \item{Year}{Character years from \code{"1920"} to \code{"1934"}.}
#' }
#'
#' @usage data(valle_dj_1920_1934)
#'
#' @examples
#' data(valle_dj_1920_1934)
#'
#' ## Basic checks
#' dim(valle_dj_1920_1934)
#' dimnames(valle_dj_1920_1934)
#'
#' ## Extract the 1929 quarterly values (4 x 2 matrix)
#' valle_dj_1920_1934[, , "1929"]
#'
#' ## Extract DJ dividends for all quarters/years (4 x 15 matrix)
#' dividends <- valle_dj_1920_1934[, "DJ dividends", ]
#' head(dividends)
#'
#' ## Application example
#' # res <- mv_ecm("MVSN",
#' #                samples = valle_dj_1920_1934,
#' #                cc = NULL, LS = NULL,
#' #                precision = 1e-7, MaxIter = 50
#' # )
"valle_dj_1920_1934"