# MVCens

`MVCens` is an R package for simulation, likelihood-based estimation, density and log-likelihood evaluation, and reproducible Monte Carlo studies for matrix-variate statistical models. It emphasizes **censoring**, **missing values**, **asymmetry**, and **heavy-tailed behavior**.

Matrix-valued observations arise in multivariate longitudinal studies, spatio-temporal data, environmental monitoring, and financial panels. `MVCens` retains their row-and-column structure and provides a common interface for generating, fitting, and evaluating several model families.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("atilapcorreia/MVCens")
library(MVCens)
```

## Quick start

The main workflow uses `mv_random()` to generate observations and `mv_fit()` to estimate model parameters. Samples are numeric arrays of dimensions `p x q x n`, where each `p x q` slice is one matrix-valued observation.

```r
set.seed(123)
p <- 2
q <- 3
n <- 50

M <- matrix(0, p, q)
Sigma <- diag(p)
Psi <- diag(q)

X <- mv_random("MVN", n = n, M = M, Sigma = Sigma, Psi = Psi)
dim(X)

fit <- mv_fit("MVN", X = X)
fit$M
fit$Sigma
fit$Psi
```

For censored or missing data, the generator supplies the observations, censoring indicators (`cc`), and censoring limits (`LS`) needed by the corresponding fit:

```r
sim <- mv_random(
  "MVNC", n = n, M = M, Sigma = Sigma, Psi = Psi,
  cens = 0.15, Ind = 1
)

fit_cens <- mv_fit(
  "MVNC", X = sim$X.cens, cc = sim$cc, LS = sim$LS
)
```

## Main functions

| Task | Functions | Description |
|:-----|:----------|:------------|
| Simulation | `mv_random()` | Generate complete or censored matrix-valued samples. |
| Estimation | `mv_fit()` | Estimate supported models using likelihood-based EM-type procedures. |
| Monte Carlo | `mv_monte_carlo()` | Run reproducible generate-and-fit studies and summarize estimation performance. |
| Densities | `dmvsn()`, `dmvst()`, `dmvrsn()`, `dmvren()` | Evaluate distribution densities; the first two are vector-based kernels. |
| Log-likelihoods | `loglik_mvn()`, `loglik_mvnc()`, `loglik_mvsn()`, `loglik_mvsnc()`, `loglik_mvst()`, `loglik_mvrsn()`, `loglik_mvren()` | Evaluate likelihoods at supplied parameters. |
| Moments | `mvrsn_mean()`, `mvrsn_covariances()`, `mvren_mean()`, `mvren_covariances()` | Obtain theoretical mean and covariance summaries. |
| Reference parameters | `mvrsn_article_parameters()`, `mvren_article_parameters()` | Retrieve predefined simulation configurations. |
| Data preparation | `prepare_beijing_mvren()`, `prepare_dj_data()` | Construct matrix-valued datasets for analysis. |

The shared notation uses `M` for location, `A` for skewness or row-specific latent effects, `Sigma` for the row covariance/scale matrix, and `Psi` for the column covariance/scale matrix. The applicable arguments and returned objects vary by model; consult `?mv_random`, `?mv_fit`, and `?mv_monte_carlo` for details.

## Supported models

| Model | Description | Generation | Fitting | Monte Carlo |
|:------|:------------|:----------:|:-------:|:-----------:|
| `MVN` | Matrix-Variate Normal | Yes | Yes | Yes |
| `MVNC` | Censored Matrix-Variate Normal | Yes | Yes | Yes |
| `MVSN` | Matrix-Variate Skew-Normal | Yes | Yes | Yes |
| `MVSNC` | Censored Matrix-Variate Skew-Normal | Yes | Yes | Yes |
| `MVST` | Matrix-Variate Skew-t | Yes | Yes | Yes |
| `MVRSN` | Matrix-Variate Row Skew-Normal | Yes | Yes | Yes |
| `MVREN` | Matrix-Variate Row Exponential-Normal | Yes | Yes | Yes |
| `MVNIG` | Matrix-Variate Normal-Inverse Gaussian | Yes | No | No |
| `MVVG` | Matrix-Variate Variance-Gamma | Yes | No | No |

`MVNIG` and `MVVG` are **generation-only** models; they are not supported by `mv_fit()` or `mv_monte_carlo()`.

## Data generation and fitting

For complete-data models, `mv_random()` usually returns a `p x q x n` array. For `MVNC` and `MVSNC`, it returns a list that includes `X.cens`, `cc`, and `LS`. The `Ind` argument specifies the incomplete-data mechanism:

| `Ind` | Mechanism |
|:----:|:----------|
| `1` | Interval censoring |
| `2` | Missing values |
| `3` | A mixture of censoring and missingness |

For `MVRSN` and `MVREN`, `return_latent = TRUE` can additionally expose generated latent variables.

`mv_fit()` dispatches to the appropriate model-specific estimation routine. Common controls include `precision` (convergence tolerance) and `max_iter` (iteration limit). Certain models also allow user-supplied starting matrices and numerical controls. The returned fit object is model-specific and may contain estimated parameters, a likelihood path, convergence information, BIC, and diagnostics. Check convergence before interpreting estimates; a small iteration limit used for a demonstration is not sufficient evidence of convergence.

## MVREN and MVRSN

The row-specific models introduce independent latent effects for different matrix rows. In the MVREN formulation,

```text
X = M + diag(W) A + V,
```

where the entries of `W` are independent exponential random variables and `V` has a matrix-normal distribution with row scale `Sigma` and column scale `Psi`. MVRSN instead uses independent half-normal row effects.

The high-level MVREN generation and fitting interfaces use **unit exponential rates** for identifiability (`lambda_i = 1`). Lower-level MVREN density, likelihood, and moment functions accept optional non-unit rates to evaluate explicitly supplied parameterizations. For numerical handling of the MVREN `Q` matrix, `q_policy` supports `"warn"`, `"strict"`, and `"regularize"`; see `?dmvren`.

The functions `mvren_mean()`, `mvren_covariances()`, `mvrsn_mean()`, and `mvrsn_covariances()` calculate theoretical summaries rather than estimates. Reference configurations are available through `mvren_article_parameters()` and `mvrsn_article_parameters()`.

## Monte Carlo studies

`mv_monte_carlo()` generates datasets from specified parameters, fits a model to each replication, and records performance summaries. Independent replications can run in parallel; iterations within an individual EM-type fit remain sequential. A fixed seed provides reproducible random streams across worker counts.

```r
truth <- list(M = M, Sigma = Sigma, Psi = Psi)

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

The output includes replication-level `results`, an aggregated `summary`, and the generating parameters. Available measures depend on the model and may include parameter errors, log-likelihood, BIC, iteration counts, convergence indicators, and numerical diagnostics. The summaries are descriptive measures of simulation performance.

## Example datasets

### Beijing Multi-Site Air Quality

`prepare_beijing_mvren()` downloads hourly observations from the UCI repository and forms `3 x 4` matrices: rows are PM2.5, NO2, and O3; columns are four six-hour periods. A period mean requires at least four observed hourly measurements, and a station-day is retained only when all 12 matrix entries are available. Standardization is performed separately for each pollutant across retained observations.

```r
beijing <- prepare_beijing_mvren(output_dir = "~/Beijing")
X_beijing <- beijing$X_std

dim(X_beijing)
head(beijing$obs_info)
```

The function returns original and standardized arrays together with supporting tables, and saves processed `.rds` and `.RData` files. Its first download requires an internet connection. It **prepares** observations but does not fit a model.

### Quarterly Dow-Jones dividends and divisor, 1920-1934

`prepare_dj_data()` uses the embedded historical quarterly values; no download is required. It organizes the series as 15 annual `4 x 2` matrices (quarters by dividends/divisor), stacked into a `4 x 2 x 15` array. By default, the two variables are standardized separately using all quarters and years.

```r
dj <- prepare_dj_data(output_dir = "~/DowJones")

X_dj <- dj$X_std
dim(X_dj)
dj$dj_data[["1920"]]

# Original, unstandardized data
X_original <- dj$X_raw
```

As documented in the supplied manual, the function returns `X_raw`, `X_std`, and a named annual-matrix list `dj_data`. It saves `dow_jones_1920_1934_processed.rds` in `output_dir`; use `readRDS()` to reload it. Because the 15 matrices represent consecutive years, their independence should not be assumed automatically.

## Documentation

Use R's built-in help pages for full arguments, model-specific options, and return-value details:

```r
help(package = "MVCens")
?mv_random
?mv_fit
?mv_monte_carlo
?prepare_beijing_mvren
?prepare_dj_data
```

## License

This package is distributed under the MIT License. See `LICENSE` for details.
