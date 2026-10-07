# MVCens

`MVCens` is an R package for simulation, likelihood-based estimation,
density and log-likelihood evaluation, and reproducible Monte Carlo studies
for matrix-variate statistical models. It places particular emphasis on
**censoring**, **missing values**, **asymmetry**, and **heavy-tailed behavior**.

Matrix-valued observations arise naturally in multivariate longitudinal
studies, spatio-temporal panels, environmental monitoring grids, repeated
image-like measurements, and financial panels. `MVCens` provides a unified
framework for analyzing such data without discarding their row-and-column
structure.

## Installation

Install the development version of `MVCens` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("atilapcorreia/MVCens")
```

Then load the package with:

```r
library(MVCens)
```

## Quick start

A basic complete-data workflow consists of generating matrix-valued
observations and fitting the corresponding model.

```r
library(MVCens)

p <- 2
q <- 3
n <- 50

M <- matrix(0, p, q)
Sigma <- diag(p)
Psi <- diag(q)

X <- mv_random(
  model = "MVN",
  n = n,
  M = M,
  Sigma = Sigma,
  Psi = Psi
)

dim(X)

fit <- mv_fit(
  model = "MVN",
  X = X
)

fit$M
fit$Sigma
fit$Psi
```

For censored data, `mv_random()` returns the observed censored sample together
with the censoring indicators and censoring limits required by `mv_fit()`.

```r
sim <- mv_random(
  model = "MVNC",
  n = n,
  M = M,
  Sigma = Sigma,
  Psi = Psi,
  cens = 0.15,
  Ind = 1
)

fit_cens <- mv_fit(
  model = "MVNC",
  X = sim$X.cens,
  cc = sim$cc,
  LS = sim$LS
)
```

## Main functions

The package is organized around a compact high-level interface:

- `mv_random()` generates random matrix-valued observations;
- `mv_fit()` fits supported models by maximum likelihood using ECM-type
  algorithms;
- `mv_monte_carlo()` runs reproducible simulation campaigns and summarizes
  estimator performance.

| Category | Functions | Purpose |
|:---------|:----------|:--------|
| Simulation | `mv_random()` | Generate random matrix-valued observations. |
| Estimation | `mv_fit()` | Fit supported models by maximum likelihood using ECM-type algorithms. |
| Monte Carlo | `mv_monte_carlo()` | Run reproducible simulation campaigns and summarize estimator performance. |
| Densities | `dmvsn()`, `dmvst()`, `dmvrsn()`, `dmvren()` | Evaluate probability density functions. |
| Log-likelihoods | `loglik_mvn()`, `loglik_mvnc()`, `loglik_mvsn()`, `loglik_mvsnc()`, `loglik_mvst()`, `loglik_mvrsn()`, `loglik_mvren()` | Evaluate model log-likelihoods directly. |
| MVRSN summaries | `mvrsn_mean()`, `mvrsn_covariances()` | Compute theoretical MVRSN moments. |
| MVREN summaries | `mvren_mean()`, `mvren_covariances()` | Compute theoretical MVREN moments. |
| Parameter sets | `mvrsn_article_parameters()`, `mvren_article_parameters()` | Return predefined configurations for examples and simulation studies. |

The common notation uses `M` for location, `A` for skewness or a latent-effect
matrix, `Sigma` for the row covariance or scale, and `Psi` for the column
covariance or scale. Samples are typically represented as numeric arrays with
dimensions `p x q x n`.

Model-specific ECM and Monte Carlo engines are internal dispatch targets. In
ordinary workflows, users should call the unified interfaces `mv_random()`,
`mv_fit()`, and `mv_monte_carlo()`.

## Supported models

| Model | Description | Generation | Fitting | Monte Carlo |
|:------|:------------|:----------:|:-------:|:-----------:|
| `MVN` | Matrix-Variate Normal | yes | yes | yes |
| `MVNC` | Censored Matrix-Variate Normal | yes | yes | yes |
| `MVSN` | Matrix-Variate Skew-Normal | yes | yes | yes |
| `MVSNC` | Censored Matrix-Variate Skew-Normal | yes | yes | yes |
| `MVST` | Matrix-Variate Skew-t | yes | yes | yes |
| `MVRSN` | Matrix-Variate Row Skew-Normal | yes | yes | yes |
| `MVREN` | Matrix-Variate Row Exponential-Normal | yes | yes | yes |
| `MVNIG` | Matrix-Variate Normal-Inverse Gaussian | yes | no | no |
| `MVVG` | Matrix-Variate Variance-Gamma | yes | no | no |

`MVNIG` and `MVVG` are generation-only models and are rejected by the fitting
and Monte Carlo interfaces.

### Return values from `mv_random()`

The return structure depends on the selected model and options.

For ordinary complete-data models, `mv_random()` returns a numeric array with
dimensions `p x q x n`.

```r
x_mvn <- mv_random(
  model = "MVN",
  n = 10,
  M = M,
  Sigma = Sigma,
  Psi = Psi
)

dim(x_mvn)
```

For censored models such as `MVNC` and `MVSNC`, the returned object is a list
containing:

- `X.cens`: observed censored data;
- `cc`: censoring indicators;
- `LS`: censoring limits.

```r
sim_mvnc <- mv_random(
  model = "MVNC",
  n = 10,
  M = M,
  Sigma = Sigma,
  Psi = Psi,
  cens = 0.15,
  Ind = 1
)

names(sim_mvnc)
```

For `MVRSN` and `MVREN`, setting `return_latent = TRUE` returns a list
containing the generated observations and the latent variables.

```r
A <- matrix(0.25, p, q)

sim_mvren <- mv_random(
  model = "MVREN",
  n = 10,
  M = M,
  A = A,
  Sigma = Sigma,
  Psi = Psi,
  lambda = rep(1, p),
  return_latent = TRUE
)

names(sim_mvren)
```

## Incomplete-data mechanisms

For censored generators, the `Ind` argument selects the incomplete-data
mechanism:

| `Ind` | Meaning |
|:-----:|:--------|
| `1` | Interval censoring |
| `2` | Missing values represented by `(-Inf, Inf)` |
| `3` | Mixture of censoring and missingness |

The censoring indicator array `cc` and upper-limit array `LS` have the same
dimensions as the observed sample and are passed to `mv_fit()` when fitting a
censored model.

## MVREN identification and numerical policy

The high-level MVREN interface adopts unit exponential rates as its
identifiability convention. Accordingly, `mv_random("MVREN", ...)` accepts an
omitted `lambda` or a vector containing only ones, and `mv_fit("MVREN", ...)`
estimates the identified model without a `lambda_mode` argument.

The lower-level MVREN density and moment functions retain an optional
row-specific `lambda` argument so that an explicitly supplied parameterization
can be evaluated directly.

The `q_policy` argument determines how a numerically non-positive-definite
matrix `Q` is handled:

- `"warn"` regularizes the matrix and issues a warning;
- `"strict"` stops with an error;
- `"regularize"` applies the correction silently.

## Densities and log-likelihoods

The density helpers can be used independently of simulation and estimation:

```r
dmvsn(y, mu, Sigma, lambda)
dmvst(y, mu, Sigma, lambda, nu)
dmvrsn(X, M, A, Sigma, Psi)
dmvren(X, M, A, Sigma, Psi, lambda = NULL)
```

Setting `log = TRUE` returns the corresponding log-density. The optional MVREN
argument `q_policy` provides the same numerical safeguards used by the fitting
interface.

Model log-likelihoods can likewise be evaluated directly:

```r
loglik_mvn(samples, M, Sigma, Psi)
loglik_mvnc(cc, LS, samples, M, Sigma, Psi)
loglik_mvsn(X, muM, AM, SigmaM, PsiM)
loglik_mvsnc(cc, LS, X, muM, AM, SigmaM, PsiM, lambdaM)
loglik_mvst(nu, X, muM, AM, SigmaM, PsiM)
loglik_mvrsn(X_array, M, A, Sigma, Psi)
loglik_mvren(X_array, M, A, Sigma, Psi, lambda = NULL)
```

These functions evaluate supplied parameters; they do not fit the model or
modify the input arrays.

## MVRSN and MVREN utilities

The package includes theoretical summaries and reproducible parameter sets for
the row skew-normal and row exponential-normal models:

```r
mvrsn_mean(M, A)
mvrsn_covariances(M, A, Sigma, Psi)
mvrsn_article_parameters()

mvren_mean(M, A, lambda = NULL)
mvren_covariances(M, A, Sigma, Psi, lambda = NULL)
mvren_article_parameters()
```

The MVREN model is defined through

```text
X = M + diag(W) A + V,
```

where the components of `W` are mutually independent exponential random
variables and `V` follows a matrix-variate normal distribution with row
covariance `Sigma` and column covariance `Psi`.

The general theoretical parameterization allows row-specific rates `lambda`;
the identified high-level generation and fitting interfaces use
`lambda_i = 1` for every row.

The article-parameter helpers provide predefined configurations for
reproducible examples and simulation studies. The covariance factors follow the
package's identifiability convention, with `Psi` normalized to have determinant
one and the reciprocal scale adjustment applied to `Sigma` where necessary.

## Monte Carlo studies

`mv_monte_carlo()` runs independent replications while preserving the
sequential update order within each ECM fit. A fixed seed produces reproducible
random streams independently of the number of workers. Every entry in
`sample_sizes` must be a positive integer.

A basic example is:

```r
truth <- list(
  M = M,
  Sigma = Sigma,
  Psi = Psi
)

study <- mv_monte_carlo(
  model = "MVN",
  sample_sizes = c(50, 100),
  replications = 10,
  truth = truth,
  workers = 1,
  seed = 123
)

study$summary
```

The returned object contains replication-level `results`, a `summary`
aggregated by sample size, and the generating parameters. Depending on the
model, recorded quantities include parameter errors, log-likelihood, BIC,
iteration count, convergence and monotonicity diagnostics, and error messages.

When `workers = NULL`, all logical processors detected on the host are used.
Parallelization is restricted to independent replications and does not change
the update sequence of an individual ECM run.

## Applications

`MVCens` is particularly useful for:

- simulation studies with matrix-valued data;
- censored or partially observed matrix data;
- asymmetric and row-specific latent-effect structures;
- heavy-tailed matrix observations;
- environmental and spatio-temporal monitoring data;
- longitudinal multivariate systems; and
- financial panels observed across assets and time.

## Quarterly Dow-Jones dividends and divisor, 1920-1934

The object `dj_data` is a matrix-valued longitudinal example represented as a
named list of annual matrices spanning 1920-1934. Each year contains a
`4 x 2` numeric matrix of quarterly Dow-Jones dividends and divisor values.

The example illustrates how repeated observations can retain meaningful row
and column dimensions instead of being flattened into vectors. This layout is
appropriate for temporal-by-variable, location-by-measurement, and related
structured multivariate settings.

```r
data("dj_data", package = "MVCens")

names(dj_data)
dim(dj_data[[1]])
dj_data[[1]]
```

## License

This package is distributed under the MIT License. See `LICENSE` for details.
