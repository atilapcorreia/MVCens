# MVCens

`MVCens` is an R package for simulation and likelihood-based estimation of matrix-variate statistical models,
with particular emphasis on **censoring**, **missing values**, **asymmetry**, and **heavy-tailed behavior**.

Supported models include matrix-variate normal, skew-normal, skew-(t),
normal-inverse Gaussian, variance-gamma, row skew-normal, row
exponential-normal, and censored counterparts where applicable.

The package is intentionally designed around two primary user-facing functions:

* `mv_random()` — unified random generation for matrix-variate models;
* `mv_fit()` — unified ECM-based estimation interface.

Most remaining routines are internal computational utilities that support
these high-level interfaces.

---

## Why MVCens?

Many modern datasets are naturally observed as matrices rather than vectors. Examples include:

* multivariate longitudinal measurements;
* spatio-temporal panels;
* environmental monitoring grids;
* repeated image-like measurements;
* financial panels observed across assets and time.

`MVCens` provides practical tools for simulating and fitting structured probabilistic models in these settings.

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

The package revolves around a small set of high-level functions that provide
a unified interface for simulation, estimation, and theoretical summaries of
matrix-variate models.

| Category | Functions | Purpose |
|----------|-----------|---------|
| Simulation | `mv_random()` | Generate random matrix-valued observations. |
| Estimation | `mv_fit()` | Fit matrix-variate models by maximum likelihood using ECM algorithms. |
| Densities | `dmvsn()`, `dmvst()`, `dmvrsn()`, `dmvren()` | Evaluate probability density functions. |
| Log-likelihoods | `loglik_mvn()`, `loglik_mvnc()`, `loglik_mvsn()`, `loglik_mvsnc()`, `loglik_mvst()`, `loglik_mvrsn()`, `loglik_mvren()` | Evaluate model log-likelihoods directly. |
| MVRSN summaries | `mvrsn_mean()`, `mvrsn_covariances()` | Compute theoretical MVRSN moments. |
| MVREN summaries | `mvren_mean()`, `mvren_covariances()` | Compute theoretical MVREN moments. |
| Reproducible parameter sets | `mvrsn_article_parameters()`, `mvren_article_parameters()` | Provide predefined parameter configurations for reproducible examples and simulation studies. |
| Monte Carlo simulation | `mvrsn_monte_carlo()`, `mvren_monte_carlo()` | Run reproducible Monte Carlo audits of the MVRSN and MVREN ECM estimators. |

### `mv_random()`

Generates matrix-valued random observations from the implemented models.

#### Available Models

| Model     | Description                            |
| --------- | -------------------------------------- |
| `"MVN"`   | Matrix-Variate Normal                  |
| `"MVSN"`  | Matrix-Variate Skew-Normal             |
| `"MVNC"`  | Censored Matrix-Variate Normal         |
| `"MVSNC"` | Censored Matrix-Variate Skew-Normal    |
| `"MVST"`  | Matrix-Variate Skew-t                  |
| `"MVNIG"` | Matrix-Variate Normal-Inverse Gaussian |
| `"MVVG"`  | Matrix-Variate Variance-Gamma          |
| `"MVRSN"` | Matrix-Variate Row Skew-Normal         |
| `"MVREN"` | Matrix-Variate Row Exponential-Normal  |

#### Incomplete-data Mechanism (`Ind`)

| Value | Meaning                                 |
| ----- | --------------------------------------- |
| `1`   | Interval censoring                      |
| `2`   | Missing values encoded as `(-Inf, Inf)` |
| `3`   | Mixture of censoring and missingness    |

---

### `mv_fit()`

Fits matrix-variate models using ECM-type algorithms.

#### Supported Models

| Model     | Description                         |
| --------- | ----------------------------------- |
| `"MVN"`   | Matrix-Variate Normal               |
| `"MVNC"`  | Censored Matrix-Variate Normal      |
| `"MVSN"`  | Matrix-Variate Skew-Normal          |
| `"MVSNC"` | Censored Matrix-Variate Skew-Normal |
| `"MVST"`  | Matrix-Variate Skew-t               |
| `"MVRSN"` | Matrix-Variate Row Skew-Normal      |
| `"MVREN"` | Matrix-Variate Row Exponential-Normal |

Returned objects may contain parameter estimates, reconstructed censored values, log-likelihood, BIC, convergence diagnostics, and iteration counts.

#### MVREN Rate Identification (`lambda_mode`)

For `model = "MVREN"`, the argument `lambda_mode` controls how the
row-specific exponential-rate vector is handled during estimation.

| Value | Meaning |
| ----- | ------- |
| `"fixed"` | Keeps the rate vector fixed at its initial value. |
| `"unconstrained"` | Estimates both `A` and `lambda` without a scale constraint; this parameterization is non-identifiable and therefore produces a warning. |
| `"unit_A_rows"` | Normalizes each row of `A` to unit Euclidean norm and rescales the corresponding rate parameter. |

The `q_policy` argument controls how a numerically non-positive-definite
matrix `Q` is handled: `"warn"` regularizes it with a warning, `"strict"`
stops with an error, and `"regularize"` applies the correction silently.

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

lambda_mvren <- rep(1, nrow(M))
```

---

The following examples illustrate the basic workflow for simulating and
fitting the supported matrix-variate models.

## Random Generation Examples

### Matrix-Variate Normal

```r
x_mvn <- mv_random("MVN", n = 50, M = M, U = U, V = V)
dim(x_mvn$X)
```

### Matrix-Variate Skew-Normal

```r
x_mvsn <- mv_random("MVSN", n = 50, M = M, A = A, U = U, V = V)
dim(x_mvsn$X)
```

### Censored Matrix-Variate Normal

```r
x_mvnc <- mv_random("MVNC", n = 50, cens = 0.15, Ind = 1,
                    M = M, U = U, V = V)
names(x_mvnc)
```

### Censored Matrix-Variate Skew-Normal

```r
x_mvsnc <- mv_random("MVSNC", n = 50, cens = 0.15, Ind = 1,
                     M = M, A = A, U = U, V = V)
names(x_mvsnc)
```

### Matrix-Variate Skew-t

```r
x_mvst <- mv_random("MVST", n = 50, M = M, A = A, U = U, V = V, nu = 6)
dim(x_mvst$X)
```

### Matrix-Variate Normal-Inverse Gaussian

```r
x_mvnig <- mv_random("MVNIG", n = 50, M = M, A = A, U = U, V = V,
                     gamma_tilde = 2.5)
names(x_mvnig)
```

### Matrix-Variate Variance-Gamma

```r
x_mvvg <- mv_random("MVVG", n = 50, M = M, A = A, U = U, V = V,
                    rate = 1.25)
names(x_mvvg)
```

### Matrix-Variate Row Skew-Normal

```r
x_mvrsn <- mv_random("MVRSN", n = 50, M = M, A = A, U = U, V = V)
names(x_mvrsn)
```

### Matrix-Variate Row Exponential-Normal

```r
x_mvren <- mv_random("MVREN",
                     n = 50,
                     M = M,
                     A = A,
                     U = U,
                     V = V,
                     lambda = lambda_mvren)
dim(x_mvren)
```

---

## Model Fitting Examples

### Matrix-Variate Normal

```r
fit_mvn <- mv_fit("MVN", samples = x_mvn$X)
names(fit_mvn)
```

### Matrix-Variate Skew-Normal

```r
fit_mvsn <- mv_fit("MVSN", samples = x_mvsn$X)
names(fit_mvsn)
```

### Censored Matrix-Variate Normal

```r
fit_mvnc <- mv_fit("MVNC",
                   samples = x_mvnc$X.cens,
                   cc = x_mvnc$cc,
                   LS = x_mvnc$LS)
names(fit_mvnc)
```

### Censored Matrix-Variate Skew-Normal

```r
fit_mvsnc <- mv_fit("MVSNC",
                    samples = x_mvsnc$X.cens,
                    cc = x_mvsnc$cc,
                    LS = x_mvsnc$LS)
names(fit_mvsnc)
```

### Matrix-Variate Skew-t

```r
fit_mvst <- mv_fit("MVST",
                   samples = x_mvst$X,
                   nu = 4,
                   get.nu = TRUE)
names(fit_mvst)
```

### Matrix-Variate Row Skew-Normal

```r
fit_mvrsn <- mv_fit("MVRSN",
                    samples = x_mvrsn$X)
names(fit_mvrsn)
```

### Matrix-Variate Row Exponential-Normal

```r
fit_mvren <- mv_fit("MVREN",
                    samples = x_mvren,
                    lambda_mode = "fixed",
                    lambda_init = lambda_mvren,
                    q_policy = "warn")
names(fit_mvren)
```

The fitted object includes `M`, `A`, `Sigma`, `Psi`, `lambda`, `loglik`,
`BIC`, the log-likelihood history, convergence information, and additional
monotonicity diagnostics.

---

## Density Evaluation

The package also exports density functions that can be used independently of
`mv_random()` and `mv_fit()`.

### Multivariate Skew-Normal

```r
dmvsn(
  y = as.vector(x_mvsn$X[, , 1]),
  mu = as.vector(M),
  Sigma = kronecker(V, U),
  lambda = as.vector(A)
)
```

### Multivariate Skew-t

```r
dmvst(
  y = as.vector(x_mvst$X[, , 1]),
  mu = as.vector(M),
  Sigma = kronecker(V, U),
  lambda = as.vector(A),
  nu = 6
)
```

### Matrix-Variate Row Skew-Normal

```r
dmvrsn(
  X = x_mvrsn$X[, , 1],
  M = M,
  A = A,
  Sigma = U,
  Psi = V
)
```

### Matrix-Variate Row Exponential-Normal

```r
dmvren(
  X = x_mvren[, , 1],
  M = M,
  A = A,
  Sigma = U,
  Psi = V,
  lambda = lambda_mvren
)
```

Set `log = TRUE` to return the corresponding log-density. The optional
`q_policy` argument provides the same numerical handling choices used by the
MVREN fitting interface.

---

## Log-Likelihood Evaluation

Observed-data log-likelihoods can also be evaluated directly without fitting
a model.

### Matrix-Variate Normal

```r
loglik_mvn(
  samples = x_mvn$X,
  M = M,
  Sigma = U,
  Psi = V
)
```

### Matrix-Variate Skew-Normal

```r
loglik_mvsn(
  dados = x_mvsn$X,
  muM = M,
  AM = A,
  SigmaM = U,
  PsiM = V
)
```

### Censored Matrix-Variate Normal

```r
loglik_mvnc(
  cc = x_mvnc$cc,
  LS = x_mvnc$LS,
  samples = x_mvnc$X.cens,
  M = M,
  Sigma = U,
  Psi = V
)
```

### Censored Matrix-Variate Skew-Normal

```r
loglik_mvsnc(
  cc = x_mvsnc$cc,
  LS = x_mvsnc$LS,
  dados = x_mvsnc$X.cens,
  muM = M,
  SigmaM = U,
  PsiM = V,
  lambdaM = A
)
```

### Matrix-Variate Skew-t

```r
loglik_mvst(
  nu = 6,
  dados = x_mvst$X,
  muM = M,
  AM = A,
  SigmaM = U,
  PsiM = V
)
```

### Matrix-Variate Row Skew-Normal

```r
loglik_mvrsn(
  X_array = x_mvrsn$X,
  M = M,
  A = A,
  Sigma = U,
  Psi = V
)
```

### Matrix-Variate Row Exponential-Normal

```r
loglik_mvren(
  X_array = x_mvren,
  M = M,
  A = A,
  Sigma = U,
  Psi = V,
  lambda = lambda_mvren
)
```

---

## Additional MVRSN Utilities

The package also provides exported utility functions for the Matrix-Variate
Row Skew-Normal (MVRSN) distribution. These functions compute theoretical
moment summaries implied by the model.

### Theoretical Mean

```r
mvrsn_mean(
  M = M,
  A = A
)
```

### Theoretical Covariance Summaries

```r
covs_mvrsn <- mvrsn_covariances(
  M = M,
  A = A,
  Sigma = U,
  Psi = V
)

covs_mvrsn$row
covs_mvrsn$column
```

### Predefined Simulation Parameters

```r
mvrsn_parameters <- mvrsn_article_parameters()
names(mvrsn_parameters)
```

The function `mvrsn_article_parameters()` returns the predefined matrices
`M`, `A`, `Sigma`, and `Psi` used in the MVRSN simulation study. By default,
`Psi` is normalized to have determinant one, consistently with the
identifiability convention adopted by the MVRSN ECM estimator.

---

## Additional MVREN Utilities

The Matrix-Variate Row Exponential-Normal (MVREN) model uses the stochastic
representation

```text
X = M + diag(W) A + V,
```

where the components of `W` are mutually independent exponential random
variables with row-specific rates `lambda`, and `V` follows a matrix-variate
normal distribution with row covariance `Sigma` and column covariance `Psi`.

### Theoretical Mean

```r
mvren_mean(
  M = M,
  A = A,
  lambda = lambda_mvren
)
```

### Theoretical Covariance Summaries

```r
covs_mvren <- mvren_covariances(
  M = M,
  A = A,
  Sigma = U,
  Psi = V,
  lambda = lambda_mvren
)

covs_mvren$row
covs_mvren$column
```

### Predefined Simulation Parameters

```r
mvren_parameters <- mvren_article_parameters()
names(mvren_parameters)
```

The function `mvren_article_parameters()` returns the predefined matrices
`M`, `A`, `Sigma`, and `Psi`, together with the row-specific exponential-rate
vector `lambda`. By default, `Psi` is normalized to have determinant one and
the reciprocal scale adjustment is applied to `Sigma`.

---

## Monte Carlo Simulation Utilities

The package provides parallel Monte Carlo audit functions for assessing the
finite-sample behavior of the MVRSN and MVREN ECM estimators. Their default
data-generating parameter configurations are supplied by
`mvrsn_article_parameters()` and `mvren_article_parameters()`, respectively.

### MVRSN Monte Carlo Audit

```r
mvrsn_mc <- mvrsn_monte_carlo(
  sample_sizes = c(50, 100),
  replications = 10,
  truth = mvrsn_article_parameters(),
  max_iter = 200,
  tol = 1e-6,
  seed = 123
)

names(mvrsn_mc)
mvrsn_mc$summary
```

For each sample size, `mvrsn_monte_carlo()` repeatedly generates MVRSN
observations, fits the model with the ECM algorithm, and records convergence
and monotonicity diagnostics, log-likelihood and BIC values, iteration counts,
and relative Frobenius errors for `M`, `A`, `Sigma`, and `Psi`.

The returned object contains `results`, with one row per replication;
`summary`, with results aggregated by sample size; and `truth`, containing the
parameters used to generate the observations.

### MVREN Monte Carlo Audit

```r
mvren_mc <- mvren_monte_carlo(
  sample_sizes = c(50, 100),
  replications = 10,
  truth = mvren_article_parameters(),
  lambda_mode = "fixed",
  max_iter = 200,
  tol = 1e-6,
  seed = 123
)

names(mvren_mc)
mvren_mc$summary
```

For each sample size, `mvren_monte_carlo()` repeatedly generates MVREN
observations, fits the model with the ECM algorithm, and records analogous
estimation and convergence diagnostics, including the relative error of the
rate vector `lambda`.

The returned object contains replication-level `results`, a sample-size
`summary`, `generating_truth`, and `comparison_truth`. The latter accounts for
the identifiability transformation used when
`lambda_mode = "unit_A_rows"`.

---

## Documentation

After installation:

```r
?mv_random
?mv_fit

## Density functions
?dmvsn
?dmvst
?dmvrsn
?dmvren

## Log-likelihood functions
?loglik_mvn
?loglik_mvnc
?loglik_mvsn
?loglik_mvsnc
?loglik_mvst
?loglik_mvrsn
?loglik_mvren

## MVRSN utilities
?mvrsn_mean
?mvrsn_covariances
?mvrsn_article_parameters
?mvrsn_monte_carlo

## MVREN utilities
?mvren_mean
?mvren_covariances
?mvren_article_parameters
?mvren_monte_carlo

help(package = "MVCens")
```

---

## Applications

`MVCens` is particularly useful for:

* simulation studies with matrix-valued data;
* censored or partially observed matrix data;
* skewed dependence structures;
* row-specific exponential asymmetry;
* heavy-tailed matrix observations;
* financial return panels;
* environmental monitoring data;
* longitudinal multivariate systems.

---

## 📊 Quarterly Dow-Jones dividends and divisor, 1920-1934

The object `dj_data` is an example of a matrix-valued longitudinal dataset
represented as a named list of annual matrices spanning the years 1920 to
1934. Each list element corresponds to one year and contains a `4 × 2`
numeric matrix. This structure is useful for illustrating how real data can
be organized in matrix form before being analyzed with matrix-variate normal
or skew-normal models.

Because each observation is naturally stored as a matrix rather than as a
simple vector, `dj_data` provides a convenient example for testing functions
related to matrix-valued simulation and estimation.

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
```

## License

This package is distributed under the license provided in the repository.
