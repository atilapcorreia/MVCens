## data-raw/valle_dj_1920_1934.R
## Run once to create data/valle_dj_1920_1934.rda

years <- 1920:1934

dj_data <- list(
  `1920` = matrix(c(31.97,20.00,30.00,20.00,30.00,20.00,28.75,20.00), ncol=2, byrow=TRUE),
  `1921` = matrix(c(26.81,20.00,23.81,20.00,22.06,20.00,22.06,20.00), ncol=2, byrow=TRUE),
  `1922` = matrix(c(21.13,20.00,19.38,20.00,19.38,20.00,20.88,20.00), ncol=2, byrow=TRUE),
  `1923` = matrix(c(24.19,20.00,26.94,20.00,25.94,20.00,25.94,20.00), ncol=2, byrow=TRUE),
  `1924` = matrix(c(30.00,20.00,26.75,19.40,27.20,19.40,27.20,18.90), ncol=2, byrow=TRUE),
  `1925` = matrix(c(32.39,18.40,30.70,18.40,30.00,19.00,27.75,19.00), ncol=2, byrow=TRUE),
  `1926` = matrix(c(37.13,17.42,27.13,16.67,29.88,16.67,26.88,16.67), ncol=2, byrow=TRUE),
  `1927` = matrix(c(33.31,16.67,29.31,16.67,31.56,16.67,28.81,16.67), ncol=2, byrow=TRUE),
  `1928` = matrix(c(31.38,16.67,27.63,16.67,30.63,16.17,27.45,13.92), ncol=2, byrow=TRUE),
  `1929` = matrix(c(31.58,12.11,27.08,10.77,28.05,10.47,28.75,10.47), ncol=2, byrow=TRUE),
  `1930` = matrix(c(30.05,10.47,27.90,9.85,27.50,10.38,25.95,10.38), ncol=2, byrow=TRUE),
  `1931` = matrix(c(25.56,10.38,22.43,10.38,19.56,10.38,19.93,10.38), ncol=2, byrow=TRUE),
  `1932` = matrix(c(15.54,15.46,18.23,15.46,16.51,15.46,16.26,15.46), ncol=2, byrow=TRUE),
  `1933` = matrix(c(13.34,15.46,12.94,15.46,12.74,15.71,12.49,15.71), ncol=2, byrow=TRUE),
  `1934` = matrix(c(13.90,15.71,13.95,15.71,14.10,15.74,15.30,15.74), ncol=2, byrow=TRUE)
)

valle_dj_1920_1934 <- array(
  NA_real_,
  dim = c(4, 2, length(years)),
  dimnames = list(
    Quarter = paste0("Q", 1:4),
    Column  = c("DJ dividends", "DJ divisor"),
    Year    = as.character(years)
  )
)

for (i in seq_along(years)) {
  valle_dj_1920_1934[, , i] <- dj_data[[as.character(years[i])]]
}

if (!requireNamespace("usethis", quietly = TRUE)) stop("install.packages('usethis')")
usethis::use_data(valle_dj_1920_1934, overwrite = TRUE)