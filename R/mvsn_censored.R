sqrtm <- function(A) {

  if(length(A)==1)

    Asqrt = sqrt(A)

  else{

    sva <- svd(A)

    if (min(sva$d) >= 0){

      Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)

    } else{

      stop("Matrix square root is not defined")

    }
  }
  return(as.matrix(Asqrt))
}


invmills = function(x,mu=0,sd=1){

  z = (x-mu)/sd

  if(z < -1e4){

    return(-z/sd)

  } else{

    return(exp(dnorm(x,mu,sd,log = TRUE) - pnorm(q = x,mean = mu,sd = sd,log.p = TRUE)))

  }
}


# This function uses the best (in terms of time/precision) function from the packages above depending of the dimension of the vector,
# degrees of freedom and even if one method collapses

prob_opt <- function(lower, upper, mean, sigma, nu = NULL, uselog2 = FALSE) {

    # Coerce inputs

    lower <- c(lower)

    upper <- c(upper)

    mean  <- c(mean)

    sigma <- as.matrix(sigma)

    p     <- ncol(sigma)

    # Input validation

    if (length(mean)  != p) stop("mean length must match ncol(sigma)")

    if (length(lower) != p) stop("lower bound length must match ncol(sigma)")

    if (length(upper) != p) stop("upper bound length must match ncol(sigma)")

    # ---- Normal case ----

    if (is.null(nu)) {

      if (p < 10) {

        prob <- tryCatch(mvtnorm::pmvnorm(lower, upper, mean, sigma)[1], error = function(e) NA_real_)

    } else {

        prob <- NA_real_

    }

      if (is.na(prob)) {

        prob <- tlrmvnmvt::pmvn(lower, upper, mean, sigma, uselog2 = FALSE)[[1]]

      }

    }

    # ---- Student-t case ----

    else {

      if (p < 10 && nu %% 1 == 0) {

      # mvtnorm::pmvt doesn’t support mean → shift bounds

          prob <- tryCatch(mvtnorm::pmvt(lower - mean, upper - mean, sigma = sigma, df = nu)[1], error = function(e) NA_real_)

      } else {

          prob <- NA_real_

      }

      if (is.na(prob)) {

        prob <- tlrmvnmvt::pmvt(lower, upper, mean = mean, sigma = sigma, df = nu, uselog2 = FALSE)[[1]]

      }

    }

    # Numerical safeguard

    prob <- max(min(prob, 1), 0)

    # Return

    if (uselog2) return(log2(prob)) else return(prob)

  }

dmvSN <- function(y, mu, Sigma, lambda){

  SigmaM <- Sigma + lambda %*% t(lambda)

  Aux1 <- t(lambda) %*% solve(Sigma) %*% (y - mu)

  Aux2<- 1 + t(lambda) %*% solve(Sigma) %*% (lambda)

  dens <- 2 * mvtnorm::dmvnorm(y, mu, SigmaM) * pnorm(Aux1 / sqrt(Aux2))

  return(dens)

}

loglikelSN <- function(dados, muM, AM, SigmaM, PsiM){

  p <- dim(dados)[1]

  q <- dim(dados)[2]

  n <- dim(dados)[3]

  suma1 <- 0

  for(j in 1:n){

    suma1 <- suma1 + log(dmvSN(vec(dados[ , , j]), vec(muM), kronecker(PsiM, SigmaM), vec(AM)))

  }

  return(suma1)

}

##### Log-likelihood Matrix skew-normal censored

logLikCensSN <- function(cc, LS, dados, muM, SigmaM, PsiM, lambdaM){

  p <- dim(dados)[1]

  q <- dim(dados)[2]

  m <- dim(dados)[3]

  ver <- matrix(0, m, 1)

  Vari <- kronecker(PsiM, SigmaM)

  ## Azzalini's version of the parameters

  lambdaMaux <- matrix(matrixNormal::vec(lambdaM), p * q,1)

  Sigma <- Vari + lambdaMaux %*% t(lambdaMaux)

  lambda <- solve(sqrtm(Sigma)) %*% lambdaMaux / as.numeric(sqrt(1 - t(lambdaMaux) %*% solve(Sigma) %*% lambdaMaux))

  mu <- as.vector(matrixNormal::vec(muM))

  for (j in 1:m){

    cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

    LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

    y1 <- as.vector(matrixNormal::vec(dados[ , , j]))


    if(sum(cc1)==0){

      ver[j] <- MomTrunc::dmvESN(x = as.vector(y1), mu = as.vector(mu), Sigma = Sigma, lambda = lambda, tau = 0)

    }

    if(sum(cc1) > 0){

      if(sum(cc1) == p * q){

        ver[j] <- MomTrunc::pmvESN(lower = c(y1), upper = c(LS1), mu = as.vector(mu), Sigma = Sigma, lambda = lambda, tau = 0)

      }

      else{

        muc <- mu[cc1 == 1] + Sigma[cc1 == 1, cc1 == 0] %*% solve(Sigma[cc1 == 0, cc1 == 0]) %*% (y1[cc1 == 0] - mu[cc1 == 0])

        Sc <- Sigma[cc1 == 1, cc1 == 1] - Sigma[cc1 == 1, cc1 == 0] %*% solve(Sigma[cc1 == 0, cc1 == 0]) %*% Sigma[cc1 == 0,cc1 == 1]

        varphi = solve(sqrtm(Sigma)) %*% lambda

        tilvarphi.o = varphi[cc1 == 0] + solve(Sigma[cc1 == 0, cc1 == 0]) %*% Sigma[cc1 == 0, cc1 == 1] %*% varphi[cc1 == 1]

        c.oc = as.numeric(1 / sqrt(1 + t(varphi[cc1 == 1]) %*% Sc %*% varphi[cc1 == 1]))

        tau.co = as.numeric(t(tilvarphi.o) %*% (y1[cc1 == 0] - mu[cc1 == 0]))

        ver[j] <- MomTrunc::dmvESN(x = as.vector(y1[cc1 == 0]), mu = as.vector(mu[cc1 == 0]), Sigma = Sigma[cc1 == 0, cc1 == 0], lambda = c.oc*sqrtm(Sigma[cc1==0,cc1==0])%*%tilvarphi.o,tau = 0)*(MomTrunc::pmvESN(lower = c(y1[cc1==1]), upper = c(LS1[cc1==1]), mu = as.vector(muc),Sigma = Sc,lambda = sqrtm(Sc)%*%varphi[cc1==1],tau = tau.co))#(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))

      }
    }
  }

  if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin

  return(sum(log(ver)))

}

somaL3 <- function(L1, Sigma, Psi){

  n <- ncol(L1)

  p <- ncol(Sigma)

  q <- ncol(Psi)

  suma2 <- matrix(0, q, q)

  suma3 <- matrix(0, p, p)

  for(j in 1:n){

    aux <- matrix(L1[ , j], p, q)

    suma2 <- suma2 + t(aux) %*% solve(Sigma) %*% aux

    suma3 <- suma3 + aux %*% solve(Psi) %*% t(aux)

  }

  return(list(sPsi = suma2, sSigma = suma3))

}

#' Generate Matrix Variate Normal or Skew Normal Data with Censoring
#'
#' This function generates matrix-valued observations from either a Matrix Variate Normal
#' (MVN) or Matrix Variate Skew Normal (MVSN) distribution, optionally applying
#' missing or interval-censoring based on a quantile cutoff.
#'
#' @param n Integer. Number of generated matrices.
#' @param cens Numeric between 0 and 1. Censoring proportion used as the
#'   quantile cutoff. For example, cens = 0.15 means censoring values below
#'   the 15 percent quantile.
#' @param Ind Integer indicating the type of censoring:
#'   1 = symmetric interval censoring (left and right),
#'   2 = missing,
#'   3 = mixed censoring (half missing values and half censored values).
#' @param M Mean matrix with dimensions p by q.
#' @param U Row covariance matrix with dimensions p by p.
#' @param V Column covariance matrix with dimensions q by q.
#' @param A Skewness matrix with dimensions p by q. Used only when
#'   dist = "SN".
#' @param dist Character string indicating the distribution. Use "Normal"
#'   for MVN or "SN" for MVSN.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{X.cens}: Censored data as an array of size p by q by n.
#'   \item \code{cc}: Censoring indicator matrix (1 means censored, 0 means observed).
#'   \item \code{LS}: Upper limits associated with censoring.
#' }
#'
#' @details
#' Matrix-variate data are generated using the function
#' LaplacesDemon::rmatrixnorm. When dist = "SN", skewness is introduced by
#' adding A multiplied by the absolute value of a standard normal draw.
#'
#' Censoring is implemented using a quantile-based threshold controlled by
#' the argument cens. The mechanism is specified by Ind.
#'
#' @importFrom stats quantile sd rnorm
#'
#' @export
#'
#' @importFrom LaplacesDemon rmatrixnorm
#'
#' @examples
#' # Example: Generate censored matrix-variate data from
#' #          Normal and Skew-Normal models
#'
#' set.seed(123)
#'
#' # Dimensions: matrices of size p x q, sample size n
#' p <- 2
#' q <- 2
#' n <- 20
#'
#' # Model parameters
#' M <- matrix(0, p, q)        # Location matrix
#' U <- diag(p)                # Row covariance
#' V <- diag(q)                # Column covariance
#' A <- matrix(1, p, q)        # Skewness matrix
#'
#' # censoring proportion
#' cens <- 0.15
#'
#' # ---------------------------------------------------------
#' # 1) Generate censored MVN (Normal) matrix-variate data
#' # ---------------------------------------------------------
#' out_normal <- rmatrix_censored(
#'   n = n,
#'   cens = cens,
#'   Ind = 1,      # interval-censored
#'   M = M,
#'   U = U,
#'   V = V,
#'   A = A,
#'   dist = "Normal"
#' )
#'
#' # Extract components
#' X_normal  <- out_normal$X.cens   # true observations or lower censoring limits
#' cc_normal <- out_normal$cc       # censoring indicators
#' LS_normal <- out_normal$LS       # upper censoring limits
#'
#'
#' # ---------------------------------------------------------
#' # 2) Generate censored MVSN (Skew-Normal) matrix data
#' # ---------------------------------------------------------
#' out_skew <- rmatrix_censored(
#'   n = n,
#'   cens = cens,
#'   Ind = 1,      # interval-censored
#'   M = M,
#'   U = U,
#'   V = V,
#'   A = A,
#'   dist = "SN"
#' )
#'
#' # Extract components
#' X_skew  <- out_skew$X.cens
#' cc_skew <- out_skew$cc
#' LS_skew <- out_skew$LS
#'
#' # Display the first generated matrix
#' X_skew[ , , 1]

rmatrix_censored <- function(n = n, cens = cens, Ind = 1, M = M, U = U, V = V, A = A, dist = "SN"){

  p <- ncol(U)

  q <- ncol(V)

  X.or <- array(data = NA, dim = c(p, q, n))

  if(dist == "Normal"){

    for (i in 1:n) {

        X.or[,,i] <- LaplacesDemon::rmatrixnorm(M = M,U = U, V = V)

    }

  }


  if(dist == "SN"){

  for (i in 1:n){

      Xo <- rnorm(1)

      Aux <- LaplacesDemon::rmatrixnorm(M = M, U = U, V = V)

      X.or[,,i] <- A*abs(Xo)+Aux

      }

  }

  if(Ind==1){

  LS <- X.cens <- X.or

  cutoff <- stats::quantile(X.cens, prob = cens)

  cc<-(X.cens < cutoff) + 0

  LS[cc==1] <- X.cens[cc==1] + 2 * sd(X.or[cc==1]) # LS (when censoring) contains the upper limit points (for missing use +Inf)

  LS[cc==0] <- 0   ## LS, when observed (cc==0) can be any value for data points

  X.cens[cc==1] <- X.cens[cc==1] - 2 * sd(X.or[cc==1]) # allocate the lower limit in the  data (for missing use -Inf)

}

  if(Ind==2){

    LS <- X.cens <- X.or

    cutoff <- stats::quantile(X.cens, prob=cens)

    cc<-(X.cens< cutoff)+0

    LS[cc==1]<- +Inf

    LS[cc==0]<- 0  ## LS, when observed (cc==0) can be any value for data points

    X.cens[cc==1]<- -Inf

    }


  if(Ind==3){

    LS <- X.cens <- X.or

    cutoff <- stats::quantile(X.cens, prob = cens)

    cc <- (X.cens < cutoff) + 0


    # Consider 50% censored and 50% missing of the total 15%


    n1 <- round(0.5 * sum(cc))

    n2 <- round(0.5 * sum(cc))

    LS[cc==1][1:n1] <- +Inf

    LS[cc==1][(n1+1):(n1+n2)] <- X.cens[cc==1][(n1+1):(n1+n2)]+2*sd(X.or[cc==1][(n1+1):(n1+n2)])

    LS[cc==0] <- 0  ## LS, when observed (cc==0) can be any value for data points

    X.cens[cc==1][1:n1] <- -Inf

    X.cens[cc==1][(n1+1):(n1+n2)] <- X.cens[cc==1][(n1+1):(n1+n2)]-2*sd(X.or[cc==1][(n1+1):(n1+n2)]) # allocate the lower limit in the  data (for missing use -Inf)

    }

  return(list(X.cens = X.cens, cc = cc, LS = LS))

}

#' @keywords internal

mvsn_ecm <- function(dados, precision=0.00000001, MaxIter=50){

  p <- dim(dados)[1]

  q <- dim(dados)[2]

  n <- dim(dados)[3]

  loglik <- numeric()

  criterio <- 1

  count <- 0

  # initial values

  mu <- apply(dados, c(1,2), mean)

  Psi <- diag(q)

  Sigma <- diag(p)

  Vari <- kronecker(Psi, Sigma)

  A <- apply(dados, c(1,2), mean)

  while(criterio > precision){

    count <- count + 1

    print(count)

    suma00 <- 0

    suma0 <- matrix(0, p, q)

    suma1 <- matrix(0, p, q)

    suma2 <- matrix(0, q, q)

    suma3 <- matrix(0, p, p)

    mu1 <- as.vector(matrixNormal::vec(mu))

    A1 <- as.vector(matrixNormal::vec(A))

    for (j in 1:n){

          y1 <- as.vector(matrixNormal::vec(dados[ , , j]))

          Mtij2 <- as.numeric(1 / (1 + t(A1) %*% solve(Vari) %*% A1))

          Mtij <- sqrt(Mtij2)

          mutij <- Mtij2 * t(A1) %*% solve(Vari) %*% (y1 - mu1)

          dj <- mutij / Mtij

          prob <- pnorm(dj)

          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin

          E = dnorm(dj) / prob

          t1 <- as.numeric(mutij + Mtij * E)

          t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

          suma00 <- suma00 + t2

          suma0 <- suma0 + (matrix(y1, p, q) - t1 * A)

          suma1 <- suma1 + t1 * (matrix(y1, p, q) - mu)

    }

    mu <- suma0 / n

    mu1 <- as.vector(matrixNormal::vec(mu))

    A <- suma1 / as.numeric(suma00)

    A1 <- as.vector(matrixNormal::vec(A))

    for (j in 1:n){

          y1 <- as.vector(matrixNormal::vec(dados[ , , j]))

          Mtij2 <- as.numeric(1 / (1 + t(A1) %*% solve(Vari) %*% A1))

          Mtij <- sqrt(Mtij2)

          mutij <- Mtij2*t(A1) %*% solve(Vari) %*% (y1 - mu1)

          dj <- mutij / Mtij

          prob <- pnorm(dj)

          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin

          E = dnorm(dj) / prob

          t1 <- as.numeric(mutij + Mtij * E)

          t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

          omega1 <- (y1 - mu1 - t1 * A1) %*% t(y1 - mu1 - t1 * A1) + (t2 - t1^2) * A1 %*% t(A1)

          omega <- as.matrix(Matrix::nearPD(omega1)$mat)

          omega <- (omega + t(omega)) / 2

          L1 <- t(Matrix::chol(omega))

          aux <- somaL3(L1, Sigma, Psi)

          suma2 <- suma2 + aux$sPsi

      }

    Psi <- suma2 / det(suma2)^(1/q)

    Vari <- kronecker(Psi, Sigma)

    for (j in 1:n){

          y1 <- as.vector(matrixNormal::vec(dados[ , , j]))

          Mtij2 <- as.numeric(1 / (1 + t(A1) %*% solve(Vari) %*% A1))

          Mtij <- sqrt(Mtij2)

          mutij <- Mtij2*t(A1)%*%solve(Vari) %*% (y1 - mu1)

          dj <- mutij / Mtij

          prob <- pnorm(dj)

          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin

          E = dnorm(dj) / prob

          t1 <- as.numeric(mutij + Mtij * E)

          t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

          omega1 <- (y1 - mu1 - t1 * A1) %*% t(y1 - mu1 - t1 * A1) + (t2 - t1^2) * A1 %*% t(A1)

          omega <- as.matrix(Matrix::nearPD(omega1)$mat)

          omega <- (omega + t(omega)) / 2

          L1 <- t(Matrix::chol(omega))

          aux <- somaL3(L1, Sigma, Psi)

          suma3<-suma3+aux$sSigma

      }

    Sigma <- suma3 / (q * n)

    Vari <- kronecker(Psi, Sigma)

    loglik[count] <- loglikelSN(dados, mu, A, Sigma, Psi)

    if (count > 2){

      at <- (loglik[count] - loglik[count - 1]) / (loglik[count - 1] - loglik[count - 2])

      criterio <- (loglik[count] - loglik[count - 1]) / (1 - at)

    }

    if (count == MaxIter){

      criterio <- 0.000000000000001

    }

 }


  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1

  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(mu = mu, A = A, Sigma = Sigma, Psi = Psi, loglik = loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixSN"

  return(obj.out)

}

#' @keywords internal

mvsnc_ecm <- function(dados, cc, LS, precision=0.0000001, MaxIter = 50){

  ## dados: if censoring it should contain the lower limit, if missing use 0 (array)
  ## cc: indicator matrix 1: censoring 0: observed    (array)
  ## LS: superior Limit of the interval censoring, if missing then use +Inf   (array)
  ## precision:  is the precision used in the convergence
  ## MaxIter: Maximum number of iterations

  p <- dim(dados)[1]

  q <- dim(dados)[2]

  n <- dim(dados)[3]

  loglik <- numeric()

  criterio <- 1

  count <- 0

  # initial values

  dados1 <- dados

  dados2 <- dados

  dados2[cc==1] <- mean(dados[cc==0])

  muM <- apply(dados2, c(1,2), mean)

  DeltaM <- apply(dados2, c(1,2), mean)

  PsiM <- diag(q)

  SigmaM <- diag(p)

  Gamma <- kronecker(PsiM, SigmaM)

  Delta <- matrix(matrixNormal::vec(DeltaM), p * q, 1)

  Sigma <- Gamma + Delta %*% t(Delta)

  shape <- solve(sqrtm(Sigma)) %*% Delta / as.numeric(sqrt(1 - t(Delta) %*% solve(Sigma) %*% Delta))

  delta <- shape / as.numeric(sqrt(1 + t(shape) %*% shape))

  mu <- as.vector(matrixNormal::vec(muM))

  while(criterio > precision){

    count <- count + 1

    print(count)

    sumaw2 <- 0

    suma0 <- matrix(0, p, q)

    suma1 <- matrix(0, p, q)

    suma2 <- matrix(0, q, q)

    suma3 <- matrix(0, p, p)

        SSj <- sqrtm(Sigma)

        iSSj <- solve(SSj)

        M2 <- 1 / (1 + sum(shape^2))

        varphi <- iSSj %*% shape


    for (j in 1:n){

      cc1 <- as.vector(matrixNormal::vec(cc[ , , j]))

      LS1 <- as.vector(matrixNormal::vec(LS[ , , j]))

      y1 <- as.vector(matrixNormal::vec(dados[ , , j]))



      if(sum(cc1)==0){

            #aux values

            aux1 = t(varphi) %*% (y1 - mu)

            aux2 = M2 * t(Delta) %*% solve(Gamma) %*% (y1 - mu)

            #output

            T1.y = as.numeric(aux2 + sqrt(M2) * invmills(aux1))

            T2.y = as.numeric(aux2^2 + sqrt(M2) * aux2 * invmills(aux1) + M2)

            y.hat  = matrix(y1, p * q, 1)

            yy.hat =  y1 %*% t(y1)

            t1.hat = T1.y

            t2.hat = T2.y

            ty.hat = T1.y * y.hat

          } else{

            if(sum(cc1)==p*q){

              eta = sqrt(2 / pi) / sqrt(1 + sum(shape^2))

              Sigma = (Sigma + t(Sigma)) / 2

              aux3 = MomTrunc::meanvarTMD(lower = c(y1), upper = c(LS1), mu = mu, Sigma = Sigma, lambda = shape, dist = "SN")

              y.hat  = aux3$mean

              yy.hat =  aux3$EYY

              w0.hat = MomTrunc::onlymeanTMD(lower = c(y1), upper = c(LS1), mu = c(mu), Sigma = Gamma, dist = "normal")

              Ltemp  = prob_opt(lower = c(y1), upper = c(LS1), mean = c(mu), sigma = Gamma, nu = NULL, uselog2 = TRUE)

              LLtemp = MomTrunc::pmvSN(lower = c(y1), upper = c(LS1), mu = as.vector(mu), Sigma = Sigma, lambda = shape, log2 = TRUE)

              gamma  = eta*2^(Ltemp - LLtemp)

              val       = c(t(shape) %*% iSSj %*% (y.hat - mu))

              gamma.ap  = invmills(val)

              boolgamma = is.na(gamma) |  gamma <= 0 | is.infinite(gamma) | abs(gamma) > 100*abs(gamma.ap)

              gamma    = ifelse(boolgamma,gamma.ap,gamma)

              #output

              t1.hat = as.numeric(M2*t(Delta)%*%solve(Gamma)%*%(y.hat-mu) + sqrt(M2)*gamma)

              t2.hat = as.numeric(M2^2*t(Delta)%*%solve(Gamma)%*%(yy.hat - 2*y.hat%*%t(mu) + mu%*%t(mu))%*%solve(Gamma)%*%Delta + M2 +
                                  gamma*M2^(3/2)*t(Delta)%*%solve(Gamma)%*%(w0.hat - mu))

              ty.hat = as.numeric(M2*(yy.hat - y.hat%*%t(mu))%*%solve(Gamma)%*%Delta + gamma*sqrt(M2)*w0.hat)

            }else{

              #new conditional parameters

              muc         = mu[cc1==1] + Sigma[cc1==1,cc1==0] %*% solve(Sigma[cc1==0,cc1==0]) %*% (y1[cc1==0] - mu[cc1==0])

              Sc          = Sigma[cc1==1,cc1==1] - Sigma[cc1==1,cc1==0] %*% solve(Sigma[cc1==0,cc1==0]) %*% Sigma[cc1==0,cc1==1]

              Sc          = (Sc + t(Sc)) / 2

              tilvarphi.o = varphi[cc1==0] + solve(Sigma[cc1==0,cc1==0]) %*% Sigma[cc1==0,cc1==1] %*% varphi[cc1==1]

              tau.co      = as.numeric(t(tilvarphi.o) %*% (y1[cc1==0] - mu[cc1==0]))

              lambda.co   = sqrtm(Sc) %*% varphi[cc1==1]

              tautil.cc   = tau.co / sqrt(1 + sum(lambda.co^2))

              #auxiliar conditional parameters

              SS.cc     = sqrtm(Sc)

              iSS.cc    = solve(SS.cc)

              varphi.cc = iSS.cc %*% lambda.co

              Delta.cc  = SS.cc %*% lambda.co / sqrt(1 + sum(lambda.co^2))

              Gamma.cc  = Sc - Delta.cc %*% t(Delta.cc)

              mub.cc    = tautil.cc * Delta.cc

              eta.cc    = invmills(tau.co, 0, sqrt(1 + sum(lambda.co^2)))

              #first and second conditional moments

              aux3   = MomTrunc::meanvarTMD(lower = c(y1[cc1==1]),upper = c(LS1[cc1==1]),mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,dist = "ESN")

              w1.hat = aux3$mean

              w2.hat =  aux3$EYY

              Ltemp  = prob_opt(lower = c(y1[cc1==1]),upper = c(LS1[cc1==1]),c(muc - mub.cc),sigma = Gamma.cc, nu = NULL, uselog2 = TRUE)

              LLtemp = MomTrunc::pmvESN(lower = y1[cc1==1], upper = LS1[cc1==1],mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,log2 = TRUE)

              gamma.cc  = eta.cc*2^(Ltemp - LLtemp)

              # ratio approximation

              val         = c(tau.co + t(lambda.co)%*%iSS.cc%*%(w1.hat - muc))

              gamma.cc.ap = invmills(val)

              boolgamma   = is.na(gamma.cc) |  gamma.cc <= 0 | is.infinite(gamma.cc) | abs(gamma.cc) > 100*abs(gamma.cc.ap)

              gamma.cc    = ifelse(boolgamma,gamma.cc.ap,gamma.cc)

              w0.hat = MomTrunc::onlymeanTMD(lower = y1[cc1==1],upper = LS1[cc1==1],mu = c(muc - mub.cc),Sigma = Gamma.cc,dist = "normal")

              y0.hat = matrix(y1,p*q,1)

              y0.hat[cc1==1] = w0.hat

              aux4 = gamma.cc*M2^(3/2)*t(Delta)%*%solve(Gamma)%*%(y0.hat - mu)

              aux5 = gamma.cc*sqrt(M2)*y0.hat

              y.hat  = matrix(y1,p*q,1)

              y.hat[cc1==1] = w1.hat

              yy.hat =  y1%*%t(y1)

              yy.hat[cc1==0,cc1==1] = y1[cc1==0]%*%t(w1.hat)

              yy.hat[cc1==1,cc1==0] = w1.hat%*%t(y1[cc1==0])

              yy.hat[cc1==1,cc1==1] = w2.hat

              t1.hat = as.numeric(M2*t(Delta)%*%solve(Gamma)%*%(y.hat-mu) + sqrt(M2)*gamma.cc)

              t2.hat = as.numeric(M2^2*t(Delta)%*%solve(Gamma)%*%(yy.hat - 2*y.hat%*%t(mu) + mu%*%t(mu))%*%solve(Gamma)%*%Delta + M2 + aux4)

              ty.hat = as.numeric(M2*(yy.hat - y.hat%*%t(mu))%*%solve(Gamma)%*%Delta + aux5)

            }

          }


     suma0 <- suma0 + matrix(y.hat, p, q) - t1.hat * DeltaM

     suma1 <- suma1 + matrix(ty.hat, p, q) - t1.hat * muM

     sumaw2 <- sumaw2 + t2.hat

     omega1 <- yy.hat-y.hat %*% t(mu) - ty.hat %*% t(Delta) - mu %*% t(y.hat) + mu %*% t(mu)+t1.hat*(mu%*%t(Delta)) -Delta%*%t(ty.hat) +t1.hat*(Delta%*%t(mu))+ t2.hat*(Delta%*%t(Delta))


    omega <- as.matrix(Matrix::nearPD(omega1)$mat)

    omega <- make_posdef(omega, epsilon = 1e-8)

    # omega <- (omega + t(omega)) / 2

    L1 <- t(Matrix::chol(omega))

    aux <- somaL3(L1, SigmaM, PsiM)

    suma2 <- suma2 + aux$sPsi

    suma3 <- suma3 + aux$sSigma

   }

   muM <- suma0 / n

   DeltaM <- suma1 / sumaw2

   SigmaM <- suma3 / (q * n)

   SigmaM <- (SigmaM + t(SigmaM)) / 2

   PsiM <- suma2 / det(suma2)^(1/q)

   PsiM<- (PsiM + t(PsiM)) / 2

   Gamma <- kronecker(PsiM, SigmaM)

   Delta <- matrix(matrixNormal::vec(DeltaM), p * q, 1)

   Sigma <- Gamma + Delta %*% t(Delta)

   shape <- solve(sqrtm(Sigma)) %*% Delta / as.numeric(sqrt(1 - t(Delta) %*% solve(Sigma) %*% Delta))

   delta <- shape / as.numeric(sqrt(1 + t(shape) %*% shape))

   mu <- as.vector(matrixNormal::vec(muM))



    loglik[count] <- logLikCensSN(cc, LS, dados, muM, SigmaM, PsiM, DeltaM)

    if (count>2){

      at <- (loglik[count] - loglik[count - 1]) / (loglik[count - 1] - loglik[count - 2])

      criterio <- (loglik[count] - loglik[count - 1]) / (1 - at)

    }

    if (count==MaxIter){

      criterio <- 0.000000000000001

    }


  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2)

  BIC   <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(mu = muM, Sigma = SigmaM, Psi = PsiM, A = DeltaM, loglik = loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixSN.cens"
  return(obj.out)

}
