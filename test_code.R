library(MVCens)

M <- matrix(
  c(
    0.5, 1.0, 0.5, 1.0,
    1.0, 1.0, 0.5, 0.5,
    1.5, 1.5, 1.0, 1.5
  ),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

A <- matrix(
  c(
    0.18, 1.62, 1.16, 0.02,
    0.25, 0.73, 0.83, 0.13,
    0.71, 1.42, 0.22, 0.81
  ),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

Sigma <- matrix(
  c(
    0.100, 0.040, 0.016,
    0.040, 0.100, 0.040,
    0.016, 0.040, 0.100
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
)

Psi <- matrix(
  c(
    2.151657, 1.721326, 1.377061, 1.101649,
    1.721326, 2.151657, 1.721326, 1.377061,
    1.377061, 1.721326, 2.151657, 1.721326,
    1.101649, 1.377061, 1.721326, 2.151657
  ),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)


samples <- rmvsn(n = 100, M, A, Sigma, Psi)

samples <- rmatrix_censored(n = 10, cens = 0.05, Ind = 1, M = M, U = Sigma, V = Psi, A = A, dist = "SN")

output <- mv_ecm(dist = "MVSNC", samples = samples$X.cens, cc = samples$cc, LS = samples$LS)

print(output)

sqrt(sum((output$Psi - Psi)^2))

#########################################################
########## Data importation test ########################
#########################################################

data(valle_dj_1920_1934)

print(valle_dj_1920_1934)









