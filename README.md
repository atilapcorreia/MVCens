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

## References

The methodology implemented in **MVCens** builds upon the following references.

### Matrix-Variate Distributions and Regression

- Allen, G. I., & Tibshirani, R. (2010). Transposable regularized covariance models with an application to missing data imputation. *The Annals of Applied Statistics*, **4**(2), 764.
- Anderlucci, L., Montanari, A., & Viroli, C. (2014). A matrix-variate regression model with canonical states: An application to elderly danish twins. *Statistica*, **74**(4), 367–381.
- Anderson, T. W. (2003). *An introduction to multivariate statistical analysis* (3rd ed.). Wiley.
- Dutilleul, P. (1999). The MLE algorithm for the matrix normal distribution. *Journal of Statistical Computation and Simulation*, **64**(2), 105–123.
- Glanz, H., & Carvalho, L. (2018). An Expectation–Maximization algorithm for the matrix normal distribution with an application in remote sensing. *Journal of Multivariate Analysis*, **167**, 31–48.
- Gupta, A. K., & Nagar, D. K. (1999). *Matrix variate distributions*. Chapman & Hall/CRC. https://doi.org/10.1201/9781482294702
- Gupta, A. K., & Nagar, D. K. (2004). *Matrix variate distributions* (Vol. 104). Chapman & Hall/CRC. https://doi.org/10.1201/9780203489247
- Kroonenberg, P. M. (2008). *Applied multiway data analysis*. Wiley-Interscience.
- Nguyen, T. T. (1997). A note on matrix variate normal distribution. *Journal of Multivariate Analysis*, **60**(1), 148–153.
- Viroli, C. (2012). On matrix-variate regression analysis. *Journal of Multivariate Analysis*, **111**, 296–309.

### Matrix-Variate Skew-Normal Models

- Chen, J. T., & Gupta, A. K. (2005). Matrix variate skew normal distributions. *Statistics*, **39**(3), 247–253.
- Correia, A. P., Diniz, C. A. R., & Lachos, V. H. (2026). Matrix-variate skew normal distribution: Properties and estimation. *Statistical Analysis and Data Mining: An ASA Data Science Journal*, e70090. https://doi.org/10.1002/sam.70090
- Ning, W., & Gupta, A. K. (2012). Matrix variate extended skew normal distributions. *Random Operators and Stochastic Equations*, **20**(4), 299–310.

### Heavy-Tailed and Asymmetric Matrix-Variate Models

- Gallaugher, M. P., & McNicholas, P. D. (2017). A matrix variate skew-t distribution. *Stat*, **6**(1), 160–170.
- Gallaugher, M. P., & McNicholas, P. D. (2019). Three skewed matrix variate distributions. *Statistics & Probability Letters*, **145**, 103–109.
- Kozubowski, T. J., Mazur, S., & Podgórski, K. (2023). Matrix variate generalized asymmetric Laplace distributions. *Theor. Probability and Math. Statist.*, **109**(109), 45–67.
- Naderi, M., Bekker, A., Arashi, M., & Jamalizadeh, A. (2020). A theoretical framework for landsat data modeling based on the matrix variate mean-mixture of normal model. *Plos One*, **15**(4), e0230773.
- Pereira, J. S., Diniz, C., & Lachos, V. (2025). *Matrix-variate shifted generalized asymmetric Laplace distribution* (Manuscript submitted for publication). Federal University of Western Pará, Federal University of Sao Carlos, University of Connecticut.
- Tomarchio, S. D. (2024). Matrix-variate normal mean-variance birnbaum–saunders distributions and related mixture models. *Computational Statistics*, **39**(2), 405–432.
- Tomarchio, S. D., Punzo, A., & Bagnato, L. (2022). On the use of the matrix-variate tail-inflated normal distribution for parsimonious mixture modeling. In N. Salvati, C. Perna, S. Marchetti, & R. Chambers (Eds.), *Studies in theoretical and applied statistics* (pp. 407–423). Springer.
- Yurchenko, Y. (2021). Matrix variate and tensor variate Laplace distributions. *arXiv Preprint arXiv:2104.05669*. https://arxiv.org/abs/2104.05669

### Multivariate Skew-Normal and Heavy-Tailed Foundations

- Ahsanullah, M. (1985). Some characterizations of the bivariate normal distribution. *Metrika*, **32**(1), 215–218.
- Alodat, M. T., & Shakhatreh, M. K. (2020). Gaussian process regression with skewed errors. *Journal of Computational and Applied Mathematics*, **370**, 112665. https://doi.org/10.1016/j.cam.2019.112665
- Andrews, D. F., & Mallows, C. L. (1974). Scale mixtures of normal distributions. *Journal of the Royal Statistical Society: Series B (Methodological)*, **36**(1), 99–102.
- Arellano-Valle, R. B., & Genton, M. G. (2010). Multivariate extended skew-t distributions and related families. *Metron*, **LXVIII**, 201–234.
- Azzalini, A. (1985). A class of distributions which includes the normal ones. *Scandinavian Journal of Statistics*, **12**(2), 171–178.
- Azzalini, A. (1986). Further results on a class of distributions which includes the normal ones. *Statistica*, **46**(2), 199–208.
- Azzalini, A., & Capitanio, A. (1999). Statistical applications of the multivariate skew normal distribution. *Journal of the Royal Statistical Society: Series B*, **61**(3), 579–602.
- Azzalini, A., & Valle, A. D. (1996). The multivariate skew-normal distribution. *Biometrika*, **83**(4), 715–726.
- Cruz, R. de la, Arellano-Valle, R. B., & Genton, M. G. (2016). The multivariate extended skew-normal distribution: Properties and a hierarchical representation. *Statistics & Probability Letters*, **118**, 54–60.
- Gupta, A. K., Varga, T., & Bodnar, T. (2013). *Elliptically contoured models in statistics and portfolio theory*. Springer.
- Kotz, S., & Nadarajah, S. (2004). *Multivariate - t distributions and their applications*. Cambridge University Press.
- Lange, K. L., Little, R. J. A., & Taylor, J. M. G. (1989). Robust statistical modeling using the t distribution. *Journal of the American Statistical Association*, **84**(408), 881–896.
- Lin, T.-I. (2016). Extending the skew-normal distribution: A review and some new developments. *Symmetry*, **8**(12), 123.
- Sahu, S. K., Dey, D. K., & Branco, M. D. (2003). A new class of multivariate skew distributions with applications to Bayesian regression models. *Canadian Journal of Statistics*, **31**(2), 129–150.

### Censored, Truncated, Missing-Data, and Mixed-Effects Models

- Alencar, F. H., Galarza, C. E., Matos, L. A., & Lachos, V. H. (2022). Finite mixture modeling of censored and missing data using the multivariate skew-normal distribution. *Advances in Data Analysis and Classification*, **16**(3), 521–557. https://doi.org/10.1007/s11634-021-00448-5
- Bandyopadhyay, D., Lachos, V. H., Castro, L. M., & Dey, D. K. (2012). Skew-normal/independent linear mixed models for censored responses with applications to HIV viral loads. *Biometrical Journal*, **54**(3), 405–425. https://doi.org/10.1002/bimj.201000173
- Bolfarine, H., Montenegro, L. C., & Lachos, V. H. (2007). Influence diagnostics for skew-normal linear mixed models. *Sankhya: The Indian Journal of Statistics*, **69**(4), 648–670.
- Galarza, C. E., Matos, L. A., Dey, D. K., & Lachos, V. H. (2022). On moments of folded and doubly truncated multivariate extended skew-normal distributions. *Journal of Computational and Graphical Statistics*, **31**(2), 455–465.
- Galarza, C. E., Matos, L. A., & Lachos, V. H. (2022). An EM algorithm for estimating the parameters of the multivariate skew-normal distribution with censored responses. *METRON*, **80**, 231–253. https://doi.org/10.1007/s40300-021-00227-4
- Klein, J. P., & Moeschberger, M. L. (2003). *Survival analysis: Techniques for censored and truncated data*. Springer.
- Lachos, V. H., Dey, D. K., & Andrade, M. G. (2013). Linear and nonlinear mixed-effects models for censored HIV viral loads using normal/independent distributions. *Statistics in Medicine*, **32**(9), 1596–1612. https://doi.org/10.1002/sim.5587
- Lachos, V. H., Ghosh, P., & Arellano-Valle, R. B. (2007). Likelihood based inference for skew-normal independent linear mixed models. *Statistica Sinica*, **17**(1), 303–322. https://www.jstor.org/stable/24307557
- Lachos, V. H., Moreno, E. J. L., Chen, K., & Cabral, C. R. B. (2017). Finite mixture modeling of censored data using the multivariate student-t distribution. *Journal of Multivariate Analysis*, **159**, 151–167.
- Lachos, V. H., Tomarchio, S. D., Punzo, A., & Ingrassia, S. (2025). An EM algorithm for fitting matrix-variate normal distributions on interval-censored and missing data. *Statistics and Computing*, **35**(2), 1–11.
- Lin, T. I., Ho, H. J., & Chen, C. L. (2009). Analysis of multivariate skew normal models with incomplete data. *Journal of Multivariate Analysis*, **100**, 2337–2351.
- Lin, T.-I., & Wang, W.-L. (2025). Multivariate contaminated normal linear mixed models applied to Alzheimer’s disease study with censored and missing data. *Statistical Methods in Medical Research*, **34**(3), 490–507.
- Matos, L. A., Prates, M. O., Chen, M.-H., & Lachos, V. H. (2013). Likelihood-based inference for mixed-effects models with censored response using the multivariate-t distribution. *Statistica Sinica*, **23**, 1323–1342.
- Meeker, W. Q., & Escobar, L. A. (1998). *Statistical methods for reliability data*. Wiley.
- Valeriano, K. A. L., Galarza, C. E., Matos, L. A., & Lachos, V. H. (2023). Likelihood-based inference for the multivariate skew-t regression with censored or missing responses. *Journal of Multivariate Analysis*, **196**, 105174. https://doi.org/10.1016/j.jmva.2023.105174
- Wang, W.-L. (2023). Multivariate contaminated normal censored regression model: Properties and maximum likelihood inference. *Journal of Computational and Graphical Statistics*, **32**(4), 1671–1684. https://doi.org/10.1080/10618600.2023.2184375
- Wu, L. (2009). *Mixed effects models for complex data*. Chapman & Hall/CRC.
- Zhong, K., Olivari, R. C., Garay, A. M., & Lachos, V. H. (2024). Linear mixed-effects models for censored data with serial correlation errors using the multivariate student’s t-distribution. *The New England Journal of Statistics in Data Science*, **2**(1), 1–16. https://doi.org/10.51387/24-NEJSDS68

### Mixture Modeling, Clustering, and Three-Way Data

- Basford, K. E., & McLachlan, G. J. (1985). The mixture method of clustering applied to three-way data. *J. Classif.*, **2**, 109–125.
- Browne, R. P., & McNicholas, P. D. (2015). A mixture of generalized hyperbolic distributions. *Canadian Journal of Statistics*, **43**(2), 176–198.
- Cabral, C. R. B., Lachos, V. H., & Prates, M. O. (2012). Multivariate mixture modeling using skew-normal independent distributions. *Computational Statistics & Data Analysis*, **56**(1), 126–142.
- Chang, W. C. (1983). On using principal components before separating a mixture of two multivariate normal distributions. *Appl. Stat.*, **32**, 267–275.
- Doğru, F. Z., Bulut, Y. M., & Arslan, O. (2016). Finite mixtures of matrix variate t distributions. *Gazi University Journal of Science*, **29**(2), 335–341.
- Fernández, C., & Green, J. (2002). Modelling spatially correlated data via mixtures: A bayesian approach. *Journal of the Royal Statistical Society*, **64**(4), 805–826.
- Franczak, B. C., Browne, R. P., & McNicholas, P. D. (2014). Mixtures of shifted asymmetric Laplace distributions. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, **36**(6), 1149–1157. https://doi.org/10.1109/TPAMI.2013.216
- Gallaugher, M. P. B., & McNicholas, P. D. (2018). Finite mixtures of skewed matrix variate distributions. *Pattern Recognition*, **80**, 83–93.
- Gordon, A. D., & Vichi, M. (1998). Partitions of partitions. *J. Classif.*, **15**, 265–285.
- Hunt, L. A., & Basford, K. E. (1999). Fitting a mixture model to three-mode three-way data with categorical and continuous variables. *J. Classif.*, **16**, 283–296.
- Jones, M. C., & Sibson, R. (1987). What is projection pursuit? (With discussion). *J. R. Stat. Soc. A*, **150**, 1–38.
- Lozano, A. C., Li, H., Niculescu-Mizil, A., Liu, Y., Perlich, C., Hosking, J., & Abe, N. (2009). KDD09. *15th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 587–596.
- Melnykov, V., & Zhu, X. (2018). On model-based clustering of skewed matrix data. *Journal of Multivariate Analysis*, **167**, 181–194.
- Naderi, M., Tamandi, M., Mirfarah, E., Wang, W.-L., & Lin, T.-I. (2024). Three-way data clustering based on the mean-mixture of matrix-variate normal distributions. *Computational Statistics and Data Analysis*, **199**, 108016. https://doi.org/10.1016/j.csda.2024.108016
- Peel, D., & McLachlan, G. J. (2000). Robust mixture modeling using the t distribution. *Statistics and Computing*, **10**(4), 339–348.
- Sarkar, S., Zhu, X., Melnykov, V., & Ingrassia, S. (2020). On parsimonious models for modeling matrix data. *Computational Statistics & Data Analysis*, **142**, 106822.
- Silva, A., Qin, X., Rothstein, S. J., & McNicholas, P. D. (2023). Finite mixtures of matrix variate poisson-log normal distributions for three-way count data. *Bioinformatics*, **39**(5), 1–8. https://doi.org/10.1093/bioinformatics/btad167
- Stephens, M. (2000). Bayesian analysis of mixtures models with an unknown number of components - an alternative to reversible jump methods. *Annals of Statistics*, **28**(1), 40–74.
- Vermunt, J. K. (2007). A hierarchical mixture model for clustering three-way data sets. *Comput. Stat. Data Anal.*, **51**, 5368–5376.
- Vichi, M. (1999). One mode classification of a three-way data set. *J. Classif.*, **16**, 27–44.
- Vichi, M., Rocci, R., & Kiers, A. L. (2007). Simultaneous component and clustering models for three-way data: Within and between approaches. *J. Classif.*, **24**, 71–98.
- Viroli, C. (2011a). Finite mixtures of matrix normal distributions for classifying three-way data. *Statistics and Computing*, **21**, 511–522.
- Viroli, C. (2011b). Model based clustering for three-way data structures. *Bayesian Analysis*, **6**(4), 573–602.
- Wolfe, J. H. (1965). A computer program for the maximum likelihood analysis of types. *Technical Bulletin*, **65**, 15.

### EM, ECM, Optimization, and Numerical Methods

- Bilmes, J. A. et al. (1998). A gentle tutorial of the EM algorithm and its application to parameter estimation for gaussian mixture and hidden markov models. *International Computer Science Institute*, **4**(510), 126.
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge University Press.
- Brent, R. P. (2013). *Algorithms for minimization without derivatives*. Courier Corporation.
- Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society: Series B*, **39**(1), 1–22.
- Guney, Y., Arslan, O., & Yavuz, F. G. (2022). Robust estimation in multivariate heteroscedastic regression models with autoregressive covariance structures using EM algorithm. *Journal of Multivariate Analysis*, **191**, 105026. https://doi.org/10.1016/j.jmva.2022.105026
- Higham, N. J. (1988). Computing a nearest symmetric positive semidefinite matrix. *Linear Algebra and Its Applications*, **103**, 103–118.
- Liu, C., & Rubin, D. B. (1994). The ECME algorithm: A simple extension of EM and ECM with faster monotone convergence. *Biometrika*, **81**, 633–648.
- McLachlan, G. J., & Krishnan, T. (2008). *The EM algorithm and extensions*. John Wiley & Sons.
- Meng, X.-L., & Rubin, D. B. (1993). Maximum likelihood estimation via the ECM algorithm: A general framework. *Biometrika*, **80**(2), 267–278.

### Software Packages and Data Resources

- Correia, A. P., Diniz, C. A. R., & Lachos, V. H. (2026). *MVCens: Matrix-variate censored normal and skew-normal models* [Computer software]. https://github.com/atilapcorreia/MVCens
- Galarza, C. E., Kan, R., & Lachos, V. H. (2024). *MomTrunc: Moments of folded and doubly truncated multivariate distributions*. https://CRAN.R-project.org/package=MomTrunc
- Murphy, R. R., Perry, E., Keisman, J., Harcum, J., & Leppo, E. W. (2023). *Baytrends: Long term water quality trend analysis*. https://CRAN.R-project.org/package=baytrends
- Thompson, G., Ripley, B., & Venables, W. (2024). *MixMatrix: Classification with matrix variate normal and t distributions*. https://CRAN.R-project.org/package=MixMatrix

---

## License

This package is distributed under the license provided in the repository.
