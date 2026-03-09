# MVCens

`MVCens` is an R package for simulation and estimation in matrix-variate Normal and Skew-Normal models, with support for missing values and interval-censored observations.

The package is intended for data that are naturally represented as matrices, including spatio-temporal measurements, environmental records, multivariate longitudinal observations, and other structured multivariate datasets.

---

## Overview

`MVCens` provides tools for simulation and inference under matrix-variate models, with support for:

- Matrix-Variate Normal (`MVN`) models
- Matrix-Variate Normal with censoring or missingness (`MVNC`)
- Matrix-Variate Skew-Normal (`MVSN`) models
- Matrix-Variate Skew-Normal with censoring or missingness (`MVSNC`)
- Missing values
- Interval-censored observations
- ECM-based parameter estimation

By combining data generation and model fitting in a unified framework, the package is useful for simulation studies, methodological research, and applied analyses involving matrix-valued observations.

---

## Main functions

### `rmvsn()`

Generates random matrices from the Matrix-Variate Skew-Normal distribution using the stochastic representation

\[
\boldsymbol{\mathcal{X}} = \textbf{M} + |Z|\textbf{A} + \boldsymbol{\mathcal{V}},
\]

where \(Z \sim N(0,1)\) and \(\boldsymbol{\mathcal{V}}\) follows a matrix-variate Normal distribution with row covariance matrix \(U\) and column covariance matrix \(V\).

This function is useful for generating skewed matrix-valued samples with flexible dependence structures across rows and columns.

---

### `rmatrix_censored()`

Simulates matrix-valued observations under either a Matrix-Variate Normal or a Matrix-Variate Skew-Normal model, with optional incomplete-data mechanisms controlled by the argument `Ind`:

- `Ind = 1`: interval censoring
- `Ind = 2`: missing values
- `Ind = 3`: mixed mechanism with both missing and interval-censored entries

The function returns:

- `X.cens`: observed data array, containing exact values or lower limits,
- `cc`: indicator array (`1` for censored/missing entries and `0` for fully observed entries),
- `LS`: upper limits for censored entries, or `Inf` for missing values.

This function is particularly useful for constructing realistic datasets for simulation studies and for assessing estimation procedures under incomplete matrix-valued settings.

---

### `mv_ecm()`

Provides a unified interface for fitting the following models:

- `"MVN"`: Matrix-Variate Normal
- `"MVNC"`: Matrix-Variate Normal with censoring or missingness
- `"MVSN"`: Matrix-Variate Skew-Normal
- `"MVSNC"`: Matrix-Variate Skew-Normal with censoring or missingness

Depending on the selected model, the function may return quantities such as:

- the estimated location matrix,
- row and column covariance matrices,
- skewness parameters for skew-normal models,
- reconstructed values for incomplete observations when applicable,
- log-likelihood,
- BIC,
- number of iterations.

This unified interface provides a coherent workflow for simulation, estimation, and model-based analysis of matrix-variate data.

---

## Output

The object returned by `mv_ecm()` may include, depending on the fitted model:

- `mu`: estimated location matrix
- `Sigma`: estimated row covariance matrix
- `Psi`: estimated column covariance matrix
- `A`: estimated skewness matrix (for skew-normal models)
- `dadosPred`: reconstructed values for incomplete observations (for applicable models)
- `loglik`: final log-likelihood
- `BIC`: Bayesian Information Criterion
- `iter`: number of ECM iterations

Please refer to the function documentation for the complete list of arguments and returned components.

---

## Applications

`MVCens` is suitable for problems involving matrix-valued observations in areas such as:

- spatio-temporal modeling
- environmental monitoring
- multivariate longitudinal studies
- structured repeated measurements
- simulation studies with incomplete data
- model-based analysis of censored matrix-valued observations

---

## License

This package is distributed under the terms of the license provided in the repository.

## 📊 Example Matrix-Valued Dataset

The object `dj_data` is an example of a matrix-valued longitudinal dataset
represented as a named list of annual matrices spanning the years 1920 to
1934. Each list element corresponds to one year and contains a `4 × 2`
numeric matrix. This structure is useful for illustrating how real data can
be organized in matrix form before being analyzed with matrix-variate normal
or skew-normal models.

Because each observation is naturally stored as a matrix rather than as a
simple vector, `dj_data` provides a convenient example for testing functions
related to matrix-valued simulation and estimation

In particular, the dataset can be used to demonstrate how
repeated matrix observations over time may be incorporated into the modeling
framework implemented in `MVCens`.

More generally, this type of data layout is appropriate when rows and columns
have meaningful joint interpretation, such as temporal-by-variable,
location-by-measurement, or grouped multivariate response settings. Thus,
`dj_data` serves as a simple motivating example of the kind of structured
matrix-valued data for which the package was designed.

```r
dj_data <- list(
  `1920` = matrix(c(31.97, 20.00, 30.00, 20.00, 30.00, 20.00, 28.75, 20.00), ncol=2, byrow=TRUE),
  `1921` = matrix(c(26.81, 20.00, 23.81, 20.00, 22.06, 20.00, 22.06, 20.00), ncol=2, byrow=TRUE),
  `1922` = matrix(c(21.13, 20.00, 19.38, 20.00, 19.38, 20.00, 20.88, 20.00), ncol=2, byrow=TRUE),
  `1923` = matrix(c(24.19, 20.00, 26.94, 20.00, 25.94, 20.00, 25.94, 20.00), ncol=2, byrow=TRUE),
  `1924` = matrix(c(30.00, 20.00, 26.75, 19.40, 27.20, 19.40, 27.20, 18.90), ncol=2, byrow=TRUE),
  `1925` = matrix(c(32.39, 18.40, 30.70, 18.40, 30.00, 19.00, 27.75, 19.00), ncol=2, byrow=TRUE),
  `1926` = matrix(c(37.13, 17.42, 27.13, 16.67, 29.88, 16.67, 26.88, 16.67), ncol=2, byrow=TRUE),
  `1927` = matrix(c(33.31, 16.67, 29.31, 16.67, 31.56, 16.67, 28.81, 16.67), ncol=2, byrow=TRUE),
  `1928` = matrix(c(31.38, 16.67, 27.63, 16.67, 30.63, 16.17, 27.45, 13.92), ncol=2, byrow=TRUE),
  `1929` = matrix(c(31.58, 12.11, 27.08, 10.77, 28.05, 10.47, 28.75, 10.47), ncol=2, byrow=TRUE),
  `1930` = matrix(c(30.05, 10.47, 27.90, 9.85, 27.50, 10.38, 25.95, 10.38),  ncol=2, byrow=TRUE),
  `1931` = matrix(c(25.56, 10.38, 22.43, 10.38, 19.56, 10.38, 19.93, 10.38), ncol=2, byrow=TRUE),
  `1932` = matrix(c(15.54, 15.46, 18.23, 15.46, 16.51, 15.46, 16.26, 15.46), ncol=2, byrow=TRUE),
  `1933` = matrix(c(13.34, 15.46, 12.94, 15.46, 12.74, 15.71, 12.49, 15.71), ncol=2, byrow=TRUE),
  `1934` = matrix(c(13.90, 15.71, 13.95, 15.71, 14.10, 15.74, 15.30, 15.74), ncol=2, byrow=TRUE)
)
```r

## 📥 Installation

You can install the development version of `MVCens` directly from GitHub using `devtools`:

install.packages("devtools")
devtools::install_github("atilapcorreia/MVCens")
library(MVCens)

## Documentation

After installation, the documentation can be accessed with:

```r
?rmvsn
?rmatrix_censored
?mv_ecm
```

A complete index of available functions can be viewed with:

```r
help(package = "MVCens")
```
