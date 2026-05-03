# MVCens

`MVCens` is an R package for simulation and likelihood-based estimation of matrix-variate statistical models, with particular emphasis on **censoring**, **missing values**, **asymmetry**, and **heavy-tailed behavior**.

The package is intentionally designed around two user-facing functions:

- `mv_random()` — unified random generation for matrix-variate models;
- `mv_fit()` — unified ECM-based estimation interface.

All remaining routines are internal computational helpers.

---

## Why MVCens?

Many modern datasets are naturally observed as matrices rather than vectors. Examples include:

- multivariate longitudinal measurements;
- spatio-temporal panels;
- environmental monitoring grids;
- repeated image-like measurements;
- financial panels observed across assets and time.

`MVCens` provides practical tools for generating and fitting structured probabilistic models in these settings.

---

## Installation

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("atilapcorreia/MVCens")
library(MVCens)
```

---

## Main Functions

## `mv_random()`

Generates matrix-valued random observations from the implemented models.

### Available Models

| Model | Description |
|------|-------------|
| `"MVSN"` | Matrix-Variate Skew-Normal |
| `"MVNC"` | Censored Matrix-Variate Normal |
| `"MVSNC"` | Censored Matrix-Variate Skew-Normal |
| `"MVST"` | Matrix-Variate Skew-t |
| `"MVNIG"` | Matrix-Variate Normal-Inverse Gaussian |
| `"MVVG"` | Matrix-Variate Variance-Gamma |

### Incomplete-data Mechanism (`Ind`)

| Value | Meaning |
|------|---------|
| `1` | Interval censoring |
| `2` | Missing values encoded as `(-Inf, Inf)` |
| `3` | Mixture of censoring and missingness |

---

## `mv_fit()`

Fits matrix-variate models using ECM-type algorithms.

### Supported Models

| Model | Description |
|------|-------------|
| `"MVN"` | Matrix-Variate Normal |
| `"MVNC"` | Censored Matrix-Variate Normal |
| `"MVSN"` | Matrix-Variate Skew-Normal |
| `"MVSNC"` | Censored Matrix-Variate Skew-Normal |
| `"MVST"` | Matrix-Variate Skew-t |

Returned objects may contain parameter estimates, reconstructed censored values, log-likelihood, BIC, convergence diagnostics, and iteration counts.

---

## Basic Setup

```r
library(MVCens)

set.seed(123)

M <- matrix(c(1.0, 1.5, 2.0,
              0.8, 1.2, 1.7,
              1.4, 1.9, 2.3),
            nrow = 3, byrow = TRUE)

A <- matrix(c( 0.25, -0.10,  0.35,
              -0.20,  0.30,  0.15,
               0.40,  0.05, -0.25),
            nrow = 3, byrow = TRUE)

U <- matrix(c(1.50, 0.35, 0.20,
              0.35, 1.20, 0.30,
              0.20, 0.30, 1.10),
            nrow = 3, byrow = TRUE)

V <- matrix(c(1.30, 0.25, 0.10,
              0.25, 1.10, 0.20,
              0.10, 0.20, 1.40),
            nrow = 3, byrow = TRUE)
```

---

## Random Generation Examples

### Matrix-Variate Skew-Normal

```r
x_mvsn <- mv_random("MVSN", n = 50, M = M, A = A, U = U, V = V)
dim(x_mvsn)
```

### Censored Matrix-Variate Normal

```r
x_mvnc <- mv_random("MVNC", n = 50, cens = 0.15, Ind = 1,
                    M = M, U = U, V = V)
names(x_mvnc)
```

### Matrix-Variate Skew-t

```r
x_mvst <- mv_random("MVST", n = 50, M = M, A = A, U = U, V = V, nu = 6)
dim(x_mvst)
```

### Matrix-Variate Normal-Inverse Gaussian

```r
x_mvnig <- mv_random("MVNIG", n = 50, cens = 0.10, Ind = 3,
                     M = M, A = A, U = U, V = V,
                     gamma_tilde = 2.5)
names(x_mvnig)
```

---

## Model Fitting Examples

### Matrix-Variate Normal

```r
fit_mvn <- mv_fit("MVN", samples = x_mvsn)
names(fit_mvn)
```

### Censored Matrix-Variate Normal

```r
fit_mvnc <- mv_fit("MVNC",
                   samples = x_mvnc$X.cens,
                   cc = x_mvnc$cc,
                   LS = x_mvnc$LS)
names(fit_mvnc)
```

### Matrix-Variate Skew-t

```r
fit_mvst <- mv_fit("MVST",
                   samples = x_mvst,
                   nu = 4,
                   get.nu = TRUE)
names(fit_mvst)
```

---

## Documentation

After installation:

```r
?mv_random
?mv_fit
help(package = "MVCens")
```

---

## Applications

`MVCens` is particularly useful for:

- simulation studies with matrix-valued data;
- censored or partially observed matrix data;
- skewed dependence structures;
- heavy-tailed matrix observations;
- financial return panels;
- environmental monitoring data;
- longitudinal multivariate systems.

---

## License

This package is distributed under the license provided in the repository.

## 📊 Quarterly Dow-Jones dividends and divisor, 1920-1934

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
  `1930` = matrix(c(30.05, 10.47, 27.90, 9.85, 27.50, 10.38, 25.95,  10.38), ncol=2, byrow=TRUE),
  `1931` = matrix(c(25.56, 10.38, 22.43, 10.38, 19.56, 10.38, 19.93, 10.38), ncol=2, byrow=TRUE),
  `1932` = matrix(c(15.54, 15.46, 18.23, 15.46, 16.51, 15.46, 16.26, 15.46), ncol=2, byrow=TRUE),
  `1933` = matrix(c(13.34, 15.46, 12.94, 15.46, 12.74, 15.71, 12.49, 15.71), ncol=2, byrow=TRUE),
  `1934` = matrix(c(13.90, 15.71, 13.95, 15.71, 14.10, 15.74, 15.30, 15.74), ncol=2, byrow=TRUE)
)
```r
