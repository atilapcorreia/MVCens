make_posdef <- function(X, epsilon = 1e-8) {

  # ensure symmetry

  X <- (X + t(X)) / 2

  eig <- eigen(X, symmetric = TRUE)

  lambda_min <- min(eig$values)

  if (lambda_min > epsilon) {

    return(X)   # already PD

  }

  # amount needed to shift all eigenvalues above delta

  shift <- epsilon - lambda_min

  X + shift * diag(nrow(X))

}

mvn_loglik <- function(samples, M, Sigma, Psi){

  p <- dim(samples)[1]

  q <- dim(samples)[2]

  n <- dim(samples)[3]

  auxiliary_variable <- 0

  for(j in 1:n){

    auxiliary_variable <- auxiliary_variable + matrixNormal::dmatnorm(samples[ , , j], M, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)

  }

  return(auxiliary_variable)
}

mvnc_loglik <- function(cc, LS, samples, M, Sigma, Psi){

  p <- dim(samples)[1]

  q <- dim(samples)[2]

  n <- dim(samples)[3]

  ver <- matrix(0, n, 1)

  Vari <- kronecker(Psi, Sigma)

  for (j in 1:n){

    cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

    LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

    y1  <- as.vector(matrixNormal::vec(samples[ , , j]))

    mu1 <- as.vector(matrixNormal::vec(M))

    if(sum(cc1) == 0){

      ver[j] <- matrixNormal::dmatnorm(samples[ , , j], M, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)

    }

    if(sum(cc1) >= 1){

      if(sum(cc1) == p * q){

        ver[j] <- matrixNormal::pmatnorm(Lower = samples[ , , j], Upper = LS[ , , j], M, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)

      }

      else{

        muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5) %*% (y1[cc1==0]-mu1[cc1==0])

        Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5)%*%Vari[cc1==0,cc1==1]

        Sc  <- (Sc + t(Sc)) / 2

        auxM <- Vari[cc1==0, cc1==0]

        auxM <-(auxM + t(auxM)) / 2

        auxcdf <- mnormt::sadmvn(lower=as.vector(y1[cc1==1]), upper=as.vector(LS1[cc1==1]), as.vector(muc), Sc, abseps=.Machine$double.eps^0.5)

        if(auxcdf==0)  auxcdf <- .Machine$double.xmin

          ver[j] <- mnormt::dmnorm(as.vector(y1[cc1==0]), as.vector(mu1[cc1==0]), auxM, log = TRUE) + log(auxcdf)

      }
    }
  }


  if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin

  return(sum(ver))

}

somaL3 <- function(L1, Sigma, Psi){

  n <- ncol(L1)

  p <- ncol(Sigma)

  q <- ncol(Psi)

  suma2 <- matrix(0, q, q)

  suma3 <- matrix(0, p, p)

  for(j in 1:n){

      aux <- matrix(L1[,j], p, q)

      suma2 <- suma2 + t(aux) %*% solve(Sigma) %*% aux

      suma3 <- suma3 + aux %*% solve(Psi) %*% t(aux)

  }

  return(list(sPsi = suma2, sSigma = suma3))

}

#' @keywords internal

mvnc_ecm <- function(samples, cc, LS, precision = 0.0000001, MaxIter = 50){

  p <- dim(samples)[1]

  q <- dim(samples)[2]

  n <- dim(samples)[3]

  loglik <- numeric()

  criterio <- 1

  count <- 0

  # initial values

  dados1 <- samples

  dados2 <- samples

  dados2[cc == 1] <- mean(samples[cc == 0])

  mu <- apply(dados2, c(1, 2), mean)

  Psi <- diag(q)

  Sigma <- diag(p)

  Vari <- kronecker(Psi, Sigma)

  while(criterio > precision){

    count <- count + 1

    suma0 <- matrix(0, p, q)

    suma2 <- matrix(0, q, q)

    suma3 <- matrix(0, p, p)

    for (j in 1:n){

      cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

      LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

      y1 <- as.vector(matrixNormal::vec(samples[ , , j]))

      mu1 <- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){

        omega <- matrix(y1, p, q)

        suma0 <- suma0 + omega

      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p * q) {

          muc  <- mu1

          Sc   <- Vari

          aux  <- MomTrunc::meanvarTMD(y1, LS1, muc, Sc, dist = "normal")

          tuy  <- aux$mean

          dados1[ , , j] <- matrix(tuy, p, q)

        } else {

          muc <- mu1[cc1==1] + Vari[cc1==1,cc1==0] %*% solve(Vari[cc1==0,cc1==0]) %*% (y1[cc1==0] - mu1[cc1==0])

          Sc  <- Vari[cc1==1,cc1==1] - Vari[cc1==1,cc1==0] %*% solve(Vari[cc1==0, cc1==0]) %*% Vari[cc1==0, cc1==1]

          Sc  <- (Sc + t(Sc)) / 2

          aux <- MomTrunc::meanvarTMD(y1[cc1==1], LS1[cc1==1], muc, Sc, dist = "normal")

          tuy <- matrix(y1, p * q, 1)

          tuy[cc1==1] <- aux$mean

          dados1[ , , j]<- matrix(tuy, p, q)

        }

        suma0 <- suma0 + matrix(tuy, p, q)

      }

    }

    mu <- suma0 / n

    for (j in 1:n){

      cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

      LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

      y1 <- as.vector(matrixNormal::vec(samples[ , , j]))

      mu1 <- as.vector(matrixNormal::vec(mu))

      if(sum(cc1) == 0){

        omega <- matrix(y1, p, q)

        suma2 <- suma2 + t(omega - mu) %*% solve(Sigma) %*% (omega - mu)

      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p * q) {

          muc  <- mu1

          Sc   <- Vari

          aux  <- MomTrunc::meanvarTMD(y1, LS1, muc, Sc, dist = "normal")

          tuy  <- aux$mean

          tuyy <- aux$EYY

          dados1[ , , j] <- matrix(tuy, p, q)

        } else {

          muc <- mu1[cc1==1] + Vari[cc1==1, cc1==0] %*% solve(Vari[cc1==0, cc1==0]) %*% (y1[cc1==0] - mu1[cc1==0])

          Sc  <- Vari[cc1==1, cc1==1] - Vari[cc1==1,cc1==0] %*% solve(Vari[cc1==0,cc1==0]) %*% Vari[cc1==0,cc1==1]

          Sc  <- (Sc + t(Sc)) / 2

          aux <- MomTrunc::meanvarTMD(y1[cc1==1], LS1[cc1==1], muc, Sc, dist="normal")

          tuy <- matrix(y1, p * q, 1)

          tuy[cc1==1] <- aux$mean

          tuyy <- matrix(0, p * q, p * q)

          tuyy[cc1==1, cc1==1] <- aux$varcov

          tuyy <- tuyy + tuy %*% t(tuy)

          dados1[ , , j] <- matrix(tuy, p, q)

        }

        omega1 <- tuyy - (tuy) %*% t(mu1) - (mu1) %*% t(tuy) + mu1 %*% t(mu1)

        omega <- as.matrix(Matrix::nearPD(omega1)$mat)

        omega <- make_posdef(omega, epsilon = 1e-8)

        # omega <- (omega + t(omega)) / 2

        L1 <- t(chol(omega))

        aux <- somaL3(L1, Sigma, Psi)

        suma2 <- suma2 + aux$sPsi

      }

    }

    Psi <- suma2 / det(suma2)^(1 / q)

    Vari <- kronecker(Psi, Sigma)

    for (j in 1:n){

      cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

      LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

      y1 <- as.vector(matrixNormal::vec(samples[ , , j]))

      mu1 <- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){

        omega <- matrix(y1, p, q)

        suma3 <- suma3 + (omega- mu) %*% solve(Psi) %*% t(omega- mu)

      }

      if(sum(cc1)>= 1){

        if(sum(cc1) == p * q) {

          muc  <- mu1

          Sc   <- Vari

          aux  <- MomTrunc::meanvarTMD(y1, LS1, muc, Sc, dist = "normal")

          tuy  <- aux$mean

          tuyy <- aux$EYY

          dados1[ , , j] <- matrix(tuy, p, q)

        } else {

          muc <- mu1[cc1==1] + Vari[cc1==1, cc1==0] %*% solve(Vari[cc1==0, cc1==0]) %*% (y1[cc1==0] - mu1[cc1==0])

          Sc  <- Vari[cc1==1, cc1==1] - Vari[cc1==1, cc1==0] %*% solve(Vari[cc1==0, cc1==0]) %*% Vari[cc1==0, cc1==1]

          Sc  <- (Sc + t(Sc)) / 2

          aux <- MomTrunc::meanvarTMD(y1[cc1==1], LS1[cc1==1], muc, Sc, dist = "normal")

          tuy <- matrix(y1, p * q, 1)

          tuy[cc1==1] <- aux$mean

          tuyy <- matrix(0, p * q, p * q)

          tuyy[cc1==1, cc1==1] <- aux$varcov

          tuyy <- tuyy + tuy %*% t(tuy)

          dados1[ , , j] <- matrix(tuy, p, q)

        }

        omega1 <- tuyy - (tuy) %*% t(mu1) - (mu1) %*% t(tuy) + mu1 %*% t(mu1)

        omega <- as.matrix(Matrix::nearPD(omega1)$mat)

        omega <- make_posdef(omega, epsilon = 1e-8)

        # omega <- (omega + t(omega)) / 2

        L1 <- t(chol(omega))

        aux <- somaL3(L1, Sigma, Psi)

        suma3 <- suma3 + aux$sSigma

      }

    }

    Sigma <- suma3 / (q * n)

    Vari <- kronecker(Psi, Sigma)

    loglik[count] <- mvnc_loglik(cc, LS, samples, mu, Sigma, Psi)

    if (count > 2){

      at <- (loglik[count] - loglik[count - 1]) / (loglik[count - 1] - loglik[count - 2])

      criterio <- (loglik[count] - loglik[count - 1]) / (1 - at)

    }

    if (count==MaxIter){

      criterio <- 0.000000000000001

    }


  }

  npar <- (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1

  BIC   <- -2 * loglik[count] + npar * log(n)    # to be minimized

  obj.out <- list(mu = mu, Sigma = Sigma, Psi = Psi, dadosPred = dados1, loglik = loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixN"
  return(obj.out)

}

#' @keywords internal

mvn_ecm <- function(samples, precision = 0.0000001, MaxIter = 50){

  p <- dim(samples)[1]

  q <- dim(samples)[2]

  n <- dim(samples)[3]

  loglik <- numeric()

  criterio <- 1

  count <- 0

  mu <- apply(samples, c(1, 2), mean)

  # initial values

  Psi <- diag(q)

  Sigma <- diag(p)

  while(criterio > precision){

    count <- count + 1

    suma2 <- matrix(0, q, q)

    suma3 <- matrix(0, p, p)

    for (j in 1:n){

      suma2 <- suma2 + t(samples[ , , j] - mu) %*% solve(Sigma) %*% (samples[ , , j] - mu)

    }

    Psi <- suma2 / det(suma2)^(1 / q)

    for (j in 1:n){

      suma3 <- suma3 + (samples[ , , j] - mu) %*% solve(Psi) %*% t(samples[ , , j] - mu)

    }

    Sigma <- suma3 / (q * n)

    loglik[count] <- mvn_loglik(samples, mu, Sigma, Psi)

    if (count > 2){

      at <- (loglik[count] - loglik[count - 1]) / (loglik[count - 1] - loglik[count - 2])

      criterio <-(loglik[count] - loglik[count - 1]) / (1 - at)

    }

    if (count == MaxIter){

      criterio <- 0.000000000000001

    }

  }

  npar <- (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1

  BIC   <- -2 * loglik[count] + npar * log(n)    # to be minimized

  obj.out <- list(mu = mu, Sigma = Sigma, Psi = Psi, loglik = loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixN"

  return(obj.out)

}

