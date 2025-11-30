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

 ## 📁 Installation

 devtools::install_github("atilapcorreia/MVCens")

 ---
