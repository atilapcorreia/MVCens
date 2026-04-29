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
