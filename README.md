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

## Main functions

The package is organized around a compact high-level interface:

```r
mv_random(model, n, M, A = NULL, Sigma, Psi, ...)
mv_fit(model, X, cc = NULL, LS = NULL,
       precision = 1e-6, max_iter = 500L, ...)
mv_monte_carlo(model, sample_sizes, replications, truth,
               workers = NULL, ...)
```

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

The common notation uses `M` for location, `A` for skewness or a
latent-effect matrix, `Sigma` for the row covariance or scale, and `Psi`
for the column covariance or scale. Samples are represented as numeric arrays
with dimensions `p x q x n`.

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

`MVNIG` and `MVVG` are generation-only models and are rejected by the
fitting and Monte Carlo interfaces. For censored models, `mv_random()`
returns the censored observations together with the indicators and limits
required by `mv_fit()`. Complete-data models return matrix-valued samples in
the format appropriate to the selected model.

Fitted objects may include parameter estimates, reconstructed censored values,
the observed-data log-likelihood, BIC, iteration counts, convergence status,
the log-likelihood history, and monotonicity diagnostics.

## Incomplete-data mechanisms

For censored generators, the `Ind` argument selects the incomplete-data
mechanism:

| `Ind` | Meaning |
|:-----:|:--------|
| `1` | Interval censoring |
| `2` | Missing values represented by `(-Inf, Inf)` |
| `3` | Mixture of censoring and missingness |

The censoring indicator array `cc` and upper-limit array `LS` have the same
dimensions as the observed sample and are passed to `mv_fit()` when fitting
a censored model.

## MVREN identification and numerical policy

The high-level MVREN interface adopts unit exponential rates as its
identifiability convention. Accordingly, `mv_random("MVREN", ...)` accepts
an omitted `lambda` or a vector containing only ones, and `mv_fit("MVREN",
...)` estimates the identified model without a `lambda_mode` argument.

The lower-level MVREN density and moment functions retain an optional
row-specific `lambda` argument so that an explicitly supplied
parameterization can be evaluated directly.

The `q_policy` argument determines how a numerically
non-positive-definite matrix `Q` is handled. The option `"warn"`
regularizes the matrix and issues a warning, `"strict"` stops with an error,
and `"regularize"` applies the correction silently.

## Densities and log-likelihoods

The density helpers can be used independently of simulation and estimation:

```r
dmvsn(y, mu, Sigma, lambda)
dmvst(y, mu, Sigma, lambda, nu)
dmvrsn(X, M, A, Sigma, Psi)
dmvren(X, M, A, Sigma, Psi, lambda = NULL)
```

Setting `log = TRUE` returns the corresponding log-density. The optional
MVREN argument `q_policy` provides the same numerical safeguards used by the
fitting interface.

Model log-likelihoods can likewise be evaluated directly:

```r
loglik_mvn(samples, M, Sigma, Psi)
loglik_mvnc(cc, LS, samples, M, Sigma, Psi)
loglik_mvsn(X, muM, AM, SigmaM, PsiM)
loglik_mvsnc(cc, LS, X, muM, SigmaM, PsiM, lambdaM)
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
covariance `Sigma` and column covariance `Psi`. The general theoretical
parameterization allows row-specific rates `lambda`; the identified
high-level generation and fitting interfaces use `lambda_i = 1` for every
row.

The article-parameter helpers provide predefined configurations for
reproducible examples and simulation studies. The covariance factors follow
the package's identifiability convention, with `Psi` normalized to have
determinant one and the reciprocal scale adjustment applied to `Sigma` where
necessary.

## Monte Carlo studies

`mv_monte_carlo()` runs independent replications while preserving the
sequential update order within each ECM fit. A fixed seed produces
reproducible random streams independently of the number of workers. Every
entry in `sample_sizes` must be a positive integer.

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

## Quarterly Dow-Jones dividends and divisor, 1920–1934

The object `dj_data` is a matrix-valued longitudinal example represented as
a named list of annual matrices spanning 1920–1934. Each year contains a
`4 × 2` numeric matrix of quarterly Dow-Jones dividends and divisor values.

The example illustrates how repeated observations can retain meaningful row
and column dimensions instead of being flattened into vectors. This layout is
appropriate for temporal-by-variable, location-by-measurement, and related
structured multivariate settings.

```r
dj_data <- list(
  `1920` = matrix(c(31.97, 20.00, 30.00, 20.00, 30.00, 20.00, 28.75, 20.00), ncol = 2, byrow = TRUE),
  `1921` = matrix(c(26.81, 20.00, 23.81, 20.00, 22.06, 20.00, 22.06, 20.00), ncol = 2, byrow = TRUE),
  `1922` = matrix(c(21.13, 20.00, 19.38, 20.00, 19.38, 20.00, 20.88, 20.00), ncol = 2, byrow = TRUE),
  `1923` = matrix(c(24.19, 20.00, 26.94, 20.00, 25.94, 20.00, 25.94, 20.00), ncol = 2, byrow = TRUE),
  `1924` = matrix(c(30.00, 20.00, 26.75, 19.40, 27.20, 19.40, 27.20, 18.90), ncol = 2, byrow = TRUE),
  `1925` = matrix(c(32.39, 18.40, 30.70, 18.40, 30.00, 19.00, 27.75, 19.00), ncol = 2, byrow = TRUE),
  `1926` = matrix(c(37.13, 17.42, 27.13, 16.67, 29.88, 16.67, 26.88, 16.67), ncol = 2, byrow = TRUE),
  `1927` = matrix(c(33.31, 16.67, 29.31, 16.67, 31.56, 16.67, 28.81, 16.67), ncol = 2, byrow = TRUE),
  `1928` = matrix(c(31.38, 16.67, 27.63, 16.67, 30.63, 16.17, 27.45, 13.92), ncol = 2, byrow = TRUE),
  `1929` = matrix(c(31.58, 12.11, 27.08, 10.77, 28.05, 10.47, 28.75, 10.47), ncol = 2, byrow = TRUE),
  `1930` = matrix(c(30.05, 10.47, 27.90, 9.85,  27.50, 10.38, 25.95, 10.38), ncol = 2, byrow = TRUE),
  `1931` = matrix(c(25.56, 10.38, 22.43, 10.38, 19.56, 10.38, 19.93, 10.38), ncol = 2, byrow = TRUE),
  `1932` = matrix(c(15.54, 15.46, 18.23, 15.46, 16.51, 15.46, 16.26, 15.46), ncol = 2, byrow = TRUE),
  `1933` = matrix(c(13.34, 15.46, 12.94, 15.46, 12.74, 15.71, 12.49, 15.71), ncol = 2, byrow = TRUE),
  `1934` = matrix(c(13.90, 15.71, 13.95, 15.71, 14.10, 15.74, 15.30, 15.74), ncol = 2, byrow = TRUE)
)
```

## License

This package is distributed under the MIT License. See `LICENSE` for details.
