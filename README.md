# MVCens

`MVCens` is an R package for simulation and likelihood-based estimation of matrix-variate statistical models,
with particular emphasis on **censoring**, **missing values**, **asymmetry**, and **heavy-tailed behavior**.

Supported models include matrix-variate normal, skew-normal, skew-(t),
normal-inverse Gaussian, variance-gamma, row skew-normal, and censored
counterparts where applicable.

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
| Densities | `dmvsn()`, `dmvst()`, `dmvrsn()` | Evaluate probability density functions. |
| Log-likelihoods | `loglik_mvn()`, `loglik_mvnc()`, `loglik_mvsn()`, `loglik_mvsnc()`, `loglik_mvst()`, `loglik_mvrsn()` | Evaluate model log-likelihoods directly. |
| MVRSN summaries | `mvrsn_mean()`, `mvrsn_covariances()` | Compute theoretical MVRSN moments. |

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
covs <- mvrsn_covariances(
  M = M,
  A = A,
  Sigma = U,
  Psi = V
)

covs$row
covs$column
```

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

## Log-likelihood functions
?loglik_mvn
?loglik_mvnc
?loglik_mvsn
?loglik_mvsnc
?loglik_mvst
?loglik_mvrsn

## MVRSN utilities
?mvrsn_mean
?mvrsn_covariances

help(package = "MVCens")
```

---

## Applications

`MVCens` is particularly useful for:

* simulation studies with matrix-valued data;
* censored or partially observed matrix data;
* skewed dependence structures;
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

## References

The methodology implemented in **MVCens** builds upon the following references.

### Matrix-Variate Distributions

- Gupta, A. K., & Nagar, D. K. (2000). *Matrix Variate Distributions*. Chapman & Hall/CRC.

### Matrix-Variate Skew-Normal Models

- Azzalini, A., & Capitanio, A. (2014). *The Skew-Normal and Related Families*. Cambridge University Press.

- Chen, J. T., & Gupta, A. K. (2005). On matrix variate skew-normal distributions. *Journal of Multivariate Analysis*, **94**(1), 114–126.

- Harrar, S. W., & Gupta, A. K. (2008). On matrix variate skew-normal distributions. *Journal of Multivariate Analysis*, **99**(9), 1900–1914.

- Ning, W., & Gupta, A. K. (2023). Matrix variate extended skew-normal distribution. *Journal of Multivariate Analysis*, **194**, 105122.

- Correia, A. P., Diniz, C. A. R., & Lachos, V. H. (2026). Matrix-Variate Skew Normal Distribution: Properties and Estimation. *Statistical Analysis and Data Mining: The ASA Data Science Journal*, **19**(3), e70090. https://doi.org/10.1002/sam.70090

### Heavy-Tailed Matrix-Variate Models

- Gallaugher, M. P. B., & McNicholas, P. D. (2018). Three skewed matrix variate distributions. *Pattern Recognition*, **80**, 83–93.

### ECM Algorithm

- Meng, X.-L., & Rubin, D. B. (1993). Maximum likelihood estimation via the ECM algorithm: A general framework. *Biometrika*, **80**(2), 267–278.

---

## License

This package is distributed under the license provided in the repository.
