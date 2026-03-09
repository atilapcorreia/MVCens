📦 MatrixCensored
Tools for Matrix-Variate Normal and Skew-Normal Models with Censoring

---

✨ Overview

mvCensored is an R package for simulating and fitting matrix-variate
statistical models, with dedicated support for:

- Matrix-Variate Normal (MVN) distributions
- Matrix-Variate Skew-Normal (MVSN) distributions
- Handling missing and interval-censored observations
- Parameter estimation using EM/ECM algorithms

The package is intended for multivariate datasets naturally arranged in
matrix form, such as environmental measurements, spatio-temporal data,
multivariate longitudinal collections, and other structured arrays.

---

🚀 Key Features

✔️ Simulation of Matrix-Variate Data

#### rmvsn()
Generates random matrices from the Matrix-Variate Skew-Normal (MVSN)
distribution using the representation:

  X = M + A * |Z| + V,

where Z ~ N(0,1) and V ~ MatrixNormal(0, U, V).
This allows simulation of skewed matrix-variate samples with flexible
covariance structures.

---

#### rmatrix_censored()
Simulates matrix-valued observations under either:

- a Matrix-Variate Normal (MVN) model, or
- a Matrix-Variate Skew-Normal (MVSN) model,

and optionally introduces:

- missing entries, and/or
- interval censoring via quantile-based cutoffs.

This function is useful for generating realistic datasets for testing
algorithms that handle incomplete or censored matrix-valued data.

---

### ✔️ Unified ECM Estimation

#### mv_ecm()
A single interface for fitting four matrix-variate models:

- "MVN"   – Matrix-Variate Normal
- "MVNC"  – MVN with interval censoring or missingness
- "MVSN"  – Matrix-Variate Skew-Normal
- "MVSNC" – MVSN with interval censoring or missingness

The function automatically calls the appropriate estimation routine and
performs:

- EM or ECM optimization, depending on the selected model
- Estimation of the mean matrix M, row covariance Σ, and column covariance Ψ
- Estimation of skewness parameters for skew-normal models
- Imputation of censored or missing values when applicable
- Computation of log-likelihood, BIC, and convergence diagnostics

This unified structure provides a consistent and flexible workflow for 
matrix-variate model fitting.

---

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
  `1930` = matrix(c(30.05, 10.47, 27.90, 9.85, 27.50, 10.38, 25.95, 10.38), ncol=2, byrow=TRUE),
  `1931` = matrix(c(25.56, 10.38, 22.43, 10.38, 19.56, 10.38, 19.93, 10.38), ncol=2, byrow=TRUE),
  `1932` = matrix(c(15.54, 15.46, 18.23, 15.46, 16.51, 15.46, 16.26, 15.46), ncol=2, byrow=TRUE),
  `1933` = matrix(c(13.34, 15.46, 12.94, 15.46, 12.74, 15.71, 12.49, 15.71), ncol=2, byrow=TRUE),
  `1934` = matrix(c(13.90, 15.71, 13.95, 15.71, 14.10, 15.74, 15.30, 15.74), ncol=2, byrow=TRUE)
)
