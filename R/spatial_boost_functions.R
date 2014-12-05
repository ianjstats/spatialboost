#### Header ####
# File: spatial_boost_functions.R
# Author: Ian Johnston
# Last modified: 29NOV2014
# Purpose: Compile a list of all of the functions needed to
#          use the Spatial Boost model on binary traits
#### Set global parameters and load libraries ####
library(car) # for wcrossprod
library(corpcor) # for fast.svd
library(snow) # for clusters (gibbs sampling)
library(rlecuyer) # for RNG
#### Load in Polya-Gamma code ####
M_PI <- pi
M_1_PI <- 1 / pi
M_2_PI <- 2 / pi
M_PI_2 <- pi / 2
M_4_PI <- (2.0 * M_2_PI)
M_2_PI2 <- (M_2_PI * M_1_PI)
TRUNC <- (0.64)
#' Polya-Gamma A1N function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param x 
#' @keywords A1N
#' @export
#' @examples
#' A1N()
A1N <- function(x) {
  (K * exp(-0.5 * K * K * (x)))
}
#' Polya-Gamma A2N function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param x 
#' @keywords A2N
#' @export
#' @examples
#' A2N()
A2N <- function(x) {
  (K * exp(-1.5 * log(M_PI_2 * (x)) - M_2_PI2 * K * K / (x)))
}
#' Polya-Gamma AN function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param x 
#' @keywords AN
#' @export
#' @examples
#' AN()
AN <- function(x) {
  if(x > TRUNC)
  {
    return(A1N(x))
  }else
  {
    return(A2N(x))
  }
}
#' Polya-Gamma lpnorm function
#'
#' This function returns the log of the area to the 
#' left of the given value under the standard normal
#' probability distribution function.
#' @param z z-score
#' @keywords lpnorm
#' @export
#' @examples
#' lpnorm()
lpnorm <- function(z) {
  pnorm(z, log = TRUE)
}
#' Polya-Gamma probexp function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param z
#' @param fz 
#' @keywords probexp
#' @export
#' @examples
#' probexp()
probexp <- function(z, fz) {
  b <- sqrt(1.0 / TRUNC) * (TRUNC * z - 1)
  a <- -sqrt(1.0 / TRUNC) * (TRUNC * z + 1)
  x0 <- log(fz) + fz * TRUNC
  xb <- x0 - z + lpnorm(b)
  xa <- x0 + z + lpnorm(a)
  return (1.0 / (1.0 + M_4_PI * (exp(xb) + exp(xa))))
}
#' Polya-Gamma rtigauss function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param z 
#' @keywords rtigauss
#' @export
#' @examples
#' rtigauss()
rtigauss <- function(z) {
  x <- TRUNC + 1.0
  if (z * TRUNC < 1.0) {
    alpha <- 0.0
    while (runif(1) > alpha) {
      E1 <- rexp(1)
      E2 <- rexp(1)
      while (E1 * E1 > 2 * E2 / TRUNC) {
        E1 <- rexp(1)
        E2 <- rexp(1)
      }
      x <- 1 + E1 * TRUNC
      x <- TRUNC / (x * x)
      alpha <- exp(-0.5 * z * z * x)
    }
  }else {
    mu <- 1.0 / z
    while (x > TRUNC) {
      y <- rnorm(1)
      muy2 <- 0.5 * mu * y
      muy2 <- muy2 * muy2
      x <- mu + 2 * muy2 * (1 - sqrt(1 + mu / muy2))
      if (runif(1) > mu / (mu + x)) {
        x <- mu * mu / x
      }
    }
  }
  return(x)
}
#' Polya-Gamma rpg function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param z
#' @param fz
#' @param pq 
#' @keywords rpg
#' @export
#' @examples
#' rpg()
rpg <- function(z, fz, pq) {
  while(1) {
    x <- numeric(1)
    s <- numeric(1)
    y <- numeric(1)
    n <- 0 # iteration counter
    within <- 1 # flag for y <= s
    assign("K", M_PI_2, env = globalenv()) #n = 0, K = pi / 2
    if(runif(1) < pq) { # truncated Exp(1)?
      x <- TRUNC + rexp(1) / fz
    }else {
      x <- rtigauss(z)
    } # truncated inverse Gaussian
    s <- AN(x)
    y <- runif(1) * s
    while (within) {
      n <- n + 1
      assign("K", K + M_PI, env = globalenv())
      if (n %% 2 == 1) {
        s <- s - AN(x)
        if (y <= s) {
          return(0.25 * x)
        }
      }else {
        s <- s + AN(x)
        within <- 1*(y <= s)
      }
    }
  }
}
#' Polya-Gamma genpg function
#'
#' This function is one of several that is used to simulate a 
#' Polya-Gamma random variable.
#' @param n an integer (first parameter of the distribution)
#' @param Z a floating point number (second parameter of the distribution)
#' @keywords genpg
#' @export
#' @examples
#' genpg()
genpg <- function(n, Z) {
  i <- numeric(1)
  z <- 0.5 * abs(Z)
  fz <- 0.5 * (M_PI_2 * M_PI_2 + z * z)
  pq <- probexp(z, fz)
  r <- 0
  for(i in (0:(n - 1))) {
    r <- r + rpg(z, fz, pq)
  }
  return(r)
}
#### Numerically precise logit related functions ####
logLowL <- -20
logUppL <- 20
logepsilon <- log(5 ^ (-304))
#' l function
#'
#' This function is one of several that is used to compute
#' the log of a sum of exponents in a way that avoids overflow.
#' @param l1 a floating point number (the first exponent in the sum)
#' @param l2 a floating point number (the second exponent in the sum)
#' @keywords l
#' @export
#' @examples
#' l()
l <- function(l1, l2)
{
  mi <- min(l1, l2); ma <- max(l1, l2)
  if (mi == -Inf || mi - ma < logepsilon) return(ma)
  return(ma + log(1 + exp(mi - ma)))
}
#' lse function
#'
#' This function is one of several that is used to compute
#' the log of a sum of exponents in a way that avoids overflow.
#' @param x a vector of exponents for which to compute log(sum(exp(x))).
#' @keywords lse
#' @export
#' @examples lse(c(0,0)) # log(exp(0) + exp(0)) = log(2) = 0.6931472
#' lse()
lse <- function(x) Reduce(l, x)
#' invlogit function
#'
#' This function computes the inverse of the logit function.
#' @param x 
#' @keywords invlogit
#' @export
#' @examples invlogit(0) # 0.5
#' invlogit()
invlogit <- function(x)
{
  y <- x
  ind <- which(x < logUppL)
  if(length(ind) == 0)
  {
    y <- rep(1, length(x))
  }else
  {
    y[ind] <- exp(x[ind]) / (1 + exp(x[ind]))
    y[-ind] <- 1
  }
  return(y)
}
#' logit function
#'
#' This function computes the logit of a given value.
#' @param x a number between 0 and 1
#' @keywords logit
#' @export
#' @examples logit(0.5) # 0
#' logit()
logit <- function(x)
{
  y <- x
  ind0 <- which(y==0)
  ind1 <- which(y==1)
  indx <- intersect(which(y > 0), which(y < 1))
  if(length(ind0) > 0) {
    y[ind0] <- logLowL
  }
  if(length(ind1) > 0) {
    y[ind1] <- logUppL
  }  
  if(length(indx) > 0) {
    y[indx] <- log(y[indx]) - log(1 - y[indx])
  }
  return(y)
}
#### Functions to simulate GWAS data ####
#' Simulate GWAS Data function
#'
#' This function can be used to simulate a simple 
#' GWAS data set with a binary trait response variable.
#' @param n sample size, defaults to 50
#' @param p number of markers, defaults to 10
#' @param alpha the prior probability of association, defaults to 0.2
#' @param KAPPA value of kappa in the prior on beta, defaults to 100
#' @param sig2 value of sigma squared in the prior on beta, defaults to 1
#' @keywords simulate data
#' @export
#' @examples DATA <- simulate_gwas_data(n = 100,
#'                                      p = 1000,
#'                                      alpha = 0.01,
#'                                      KAPPA = 1e6,
#'                                      sig2 = 1e-4)
#'           which(DATA$THETA == 1) # list of causal markers
#'           summary(colMeans(DATA$X) / 2) # summary of minor allele frequencies
#' simulate_gwas_data()
simulate_gwas_data <- function(n = 50,
                               p = 10,
                               alpha = 0.2,
                               KAPPA = 1e2,
                               sig2 = 1) {
  sig <- sqrt(sig2)
  MAF <- runif(p, 0.05, 0.5)
  X <- matrix(0, n, p)
  for(j in 1:p) {
    X[,j] <- rbinom(n, 2, MAF[j])
  }
  THETA <- rbinom(p, 1, alpha)
  BETA <- rnorm(p, 0, sig * (sqrt(KAPPA) * THETA +
                               1 - THETA))
  BETA0 <- rnorm(1, 0, sig * sqrt(KAPPA))
  Y <- rbinom(n, 1, invlogit(BETA0 + X %*% BETA))
    return(list(X = X,
                Y = Y,
                THETA = THETA,
                BETA = c(BETA0, BETA)))
}
#' Truncated SVD approximation to X
#'
#' This function computes the truncated singular value
#' decomposition to a given set of genotypes, X, and a 
#' given set of additional covariates, Z, for a desired
#' tolerance on the mean squared error.
#' @param X a matrix of genotypes
#' @param Z a matrix of additional covariates, defaults to NULL
#' @param MSEtol the desired tolerance on the MSE, defaults to 1\%
#' @keywords SVD
#' @export
#' @examples DATA <- simulate_gwas_data()
#'           SVDX <- SB_optimal_SVD(DATA$X)
#' SB_optimal_SVD()
SB_optimal_SVD <- function(X,
                           Z = NULL,
                           MSEtol = 0.01) {
  ZX <- cbind(Z, X)
  print("Computing the singular value decomposition of the design matrix...")
  SVDX <- fast.svd(ZX)
  L <- length(SVDX$d)
  X_approx <- NULL
  print("Building approximation to the design matrix...")
  if(MSEtol > 0) {
    X_approx <- matrix(0, nrow(ZX), ncol(ZX))
    for(i in 1:L) {
      print(paste0("Checking the approximation of the design matrix using the top ",
                   i,
                   " singular vectors"))
      X_approx <- X_approx + SVDX$d[i] * outer(SVDX$u[, i],
                                               SVDX$v[, i])
      MSE <- mean((ZX - X_approx) ^ 2)
      if(MSE <= MSEtol) {
        print(paste0("Achieved the desired tolerance on MSE with the top ",
                     i, " singular vectors"))
        return(list(U = SVDX$u[, 1:i],
                    D = SVDX$d[1:i],
                    V = SVDX$v[, 1:i]))
      }
    }
    return(list(U = SVDX$u[, 1:i],
                D = SVDX$d[1:i],
                V = SVDX$v[, 1:i]))
  } else {
      return(list(U = SVDX$u[, 1:L],
                  D = SVDX$d[1:L],
                  V = SVDX$v[, 1:L]))
  }
}
#### Expectation-Maximization functions for qualitative traits ####
#' EM E-step on theta
#'
#' This function performs the E-step of our EM algorithm for the 
#' spatial boost model and computes the expected value of each 
#' indicator variable in our model (the theta's)
#' @param KAPPA the value of kappa in the prior on beta
#' @param GENEBOOST the negative gene boost of each marker
#' @param BETA the current value of beta 
#' @param SIG2_0 the current value of sigma squared
#' @keywords Expected value of theta
#' @export
#' @examples
#' EM_theta_qual_expectation()
EM_theta_qual_expectation <- function(KAPPA,
                                      GENEBOOST,
                                      BETA,
                                      SIG2_0)
{
  BETACONST <- 0.5*(1 / SIG2_0) * (1 / KAPPA - 1)
  logS <- GENEBOOST + 0.5 * log(KAPPA) + BETACONST * BETA ^ 2
  return(exp(-1 * apply(cbind(0, logS), 1, lse)))
}
#' EM M-step on sigma squared
#'
#' This function performs part of the M-step of our EM algorithm 
#' and computes the optimal value of sigma squared
#' @param KAPPA the value of kappa in the prior on beta
#' @param V a hyper-parameter for the prior on sigma squared (nu)
#' @param LAM a hyper-parameter for the prior sigma squared (lambda)
#' @param THETA the current expected values of theta
#' @param BETA the current values of beta
#' @keywords Optimized value of sigma squared
#' @export
#' @examples
#' EM_sigma_qual_maximize()
EM_sigma_qual_maximize <- function(KAPPA,
                                   V,
                                   LAM,
                                   THETA,
                                   BETA)
{
  p <- length(THETA)
  ETHETA_SUM <- THETA / KAPPA + (1 - THETA)
  LAM_STAR <- LAM + 0.5 * sum(BETA ^ 2 * ETHETA_SUM)
  V_STAR <- V + p / 2 + 1
  SIG2_0 <- LAM_STAR / (V_STAR + 1)
  return(SIG2_0)
}
#' EM M-step on beta utilizing SVD approximation
#'
#' This function performs part of the M-step of our EM algorithm 
#' and computes the optimal value of beta utilizing a computational
#' speed-up by using the truncated SVD approximation to the design 
#' matrix. 
#' @param U 
#' @param V
#' @param SIGMA
#' @param W
#' @param Z
#' @keywords Optimized value of beta
#' @export
#' @examples
#' EM_beta_qual_delta()
EM_beta_qual_delta <- function(U,
                               V,
                               SIGMA,
                               W,
                               Z)
{
  k <- ncol(U)
  hSIGMA <- sqrt(SIGMA)
  C <- chol(wcrossprod(U, U, W))
  S1 <- sweep(tcrossprod(C, V), 2, hSIGMA, "*")
  Cs <- chol(diag(k) + tcrossprod(S1))
  S2 <- backsolve(Cs, S1, transpose = TRUE)
  Z1 <- Z * hSIGMA
  Z2 <- crossprod(S2, S2 %*% Z1)
  return(hSIGMA * (Z1 - Z2))
}
#' EM M-step on beta without utilizing SVD approximation
#'
#' This function performs part of the M-step of our EM algorithm 
#' and computes the optimal value of beta without utilizing any 
#' computational speed-up.
#' @param Y the vector of response variables
#' @param X the full design matrix (additional covariates in front)
#' @param BETA the current value of beta
#' @param THETA the current expected values of theta
#' @param SIG2_0 the current value of sigma squared
#' @param KAPPA the value of kappa in the prior on beta
#' @keywords Optimized value of beta
#' @export
#' @examples
#' EM_beta_qual_maximize()
EM_beta_qual_maximize <- function(Y,
                                  X,
                                  BETA,
                                  THETA,
                                  SIG2_0,
                                  KAPPA)
{
  n <- nrow(X)
  MU <- invlogit(X %*% BETA)
  LIKPRCM <- diag(as.numeric(MU) * as.numeric(1 - MU), n)
  XTYM <- crossprod(X, Y - MU)
  XTWX <- crossprod(X, LIKPRCM %*% X)
  p <- length(THETA)
  PRIORPRCM <- THETA / KAPPA + (1 - THETA) / SIG2_0
  BETA <- solve(XTWX + diag(PRIORPRCM, p), XTWX %*% BETA + XTYM)
  return(BETA)
}
#' EM algorithm
#'
#' This function runs our EM algorithm once on a given GWAS data 
#' set, assuming a binary trait for the response variable.
#' @param X the matrix of genotypes
#' @param Z the matrix of additional covariates (optional)
#' @param Y the vector of response variables
#' @param weights the normalized sums of gene weights * relevances
#' @param XI0 the hyper-parameter xi_0
#' @param XI1 the hyper-parameter xi_1
#' @param KAPPA the value of kappa in the prior on beta
#' @param V a hyper-parameter in the prior on sigma squared (nu)
#' @param LAM a hyper-parameter in the prior on sigma squared (lambda)
#' @param maxiter maximum number of iterations for the EM algorithm
#' @param tol tolerance for convergence of the EM algorithm
#' @param print boolean flag for whether or not to print information
#' @param k number of additional covariates, defaults to 0
#' @param SVD_X list containing the SVD of cbind(Z, X)
#' @keywords EM algorithm
#' @export
#' @examples
#' EM_qual()
EM_qual <- function(X = NULL,
                    Z = NULL,
                    Y,
                    weights = NULL,
                    XI0 = -2,
                    XI1 = 1,
                    KAPPA = 1e3,
                    V = NULL,
                    LAM = NULL,
                    maxiter = 20,
                    tol = 1e-4,
                    print = TRUE,
                    k = 0,
                    SVD_X = NULL)
{
  n <- NULL
  p <- NULL
  ZX <- NULL
  ppl <- NULL
  if(!is.null(Z)) {
    k <- ncol(Z)
    ZX <- cbind(Z, X)
  }
  l <- NULL
  U <- NULL
  UMAT <- NULL
  VMAT <- NULL
  LAMVEC <- NULL
  if(is.null(SVD_X)) {

    n <- nrow(X)
    p <- ncol(X)
  } else {
    UMAT <- SVD_X$U
    VMAT <- SVD_X$V
    LAMVEC <- SVD_X$D
    U <- sweep(UMAT, 2, LAMVEC, "*")
    n <- nrow(UMAT)
    p <- nrow(VMAT) - k
    l <- length(LAMVEC)
  }
  GENEBOOST <- rep(0, p + k)
  if(is.null(weights)) {
    GENEBOOST <- rep(-1 * XI0, p + k)
  } else {
    if(k > 0) {
        GENEBOOST[1:k] <- -Inf
        GENEBOOST[(k + 1):(p + k)] <- -XI0 - XI1 * weights
    } else {
        GENEBOOST <- XI0 - XI1 * weights
    }
  }
  EMNORM <- 1e4
  EMITER <- 1
  BETA <- rep(0, p + k)
  THETA <- rep(0, p + k)
  if(k > 0) {
    THETA[1:k] <- 1
  }
  if(is.null(V) | is.null(LAM)) {
    V <- 3
    LAM <- 1e-1
  }
  SIG2_0 <- LAM/(V + 1)
  if(print) {
    print("Beginning EM Filter on X...")
  }
  while((EMNORM > tol) & (EMITER < maxiter))
  {
    CBETA <- BETA
    CTHETA <- THETA
    CSIG2_0 <- SIG2_0
    SIG2_0 <- EM_sigma_qual_maximize(KAPPA,
                                     V,
                                     LAM,
                                     THETA,
                                     BETA)
    THETA <- EM_theta_qual_expectation(KAPPA,
                                       GENEBOOST,
                                       BETA,
                                       SIG2_0)
    if(k > 0) {
      THETA[1:k] <- 1
    }
    SIGMA = SIG2_0 / (THETA / KAPPA + (1 - THETA))
    if(is.null(U)) {
      BETA <- EM_beta_qual_maximize(Y,
                                    ZX,
                                    CBETA,
                                    THETA,
                                    SIG2_0,
                                    KAPPA)
      MU <- invlogit(ZX %*% CBETA)
      ppl <- sum((Y - MU) ^ 2 + MU * (1 - MU))
    } else {
      MU <- invlogit(U %*% crossprod(VMAT, CBETA))
      W <- as.numeric(MU * (1 - MU))
      delta <- EM_beta_qual_delta(U,
                                  VMAT,
                                  SIGMA,
                                  W,
                                  VMAT %*% crossprod(U, Y - MU) - CBETA / SIGMA)
      BETA <- CBETA + delta
      MU <- invlogit(U %*% crossprod(VMAT, BETA))
      ppl <- sum((Y - MU) ^ 2 + MU * (1 - MU))
    }
    EMNORM <- sqrt(sum((THETA - CTHETA) ^ 2,
                       (BETA - CBETA) ^ 2,
                       (SIG2_0 - CSIG2_0) ^ 2))
    if(print)
    {
      print(paste("EM iteration ",
                  EMITER,
                  "; norm = ",
                  round(EMNORM, 4),
                  sep = ""))
    }
    EMITER <- EMITER + 1
  }
  return(list(ETHETA = THETA,
              BETA = BETA,
              SIG2 = SIG2_0,
              PPL = ppl))
}
#### Gibbs sampler functions for qualitative traits ####
#' Gibbs sampler for theta
#'
#' This function samples new values for the thetas from their 
#' conditional posterior distributions given the current values 
#' of the other parameters in the model.
#' @param KAPPA the value of kappa in the prior on beta
#' @param GENEBOOST the negative of the gene boost term for each marker
#' @param BETA the current value of beta
#' @param SIG2_0 the current value of sigma squared
#' @keywords gibbs sampler for theta
#' @export
#' @examples
#' SB_qual_thetasampler()
SB_qual_thetasampler <- function(KAPPA,
                                 GENEBOOST,
                                 BETA,
                                 SIG2_0)
{
  BETACONST <- 0.5 / SIG2_0 * (1 / KAPPA - 1)
  logS <- GENEBOOST + 0.5 * log(KAPPA) + BETACONST * BETA ^ 2
  logP <- -1 * apply(cbind(0, logS), 1, lse)
  THETA <- rbinom(length(BETA), 1, exp(logP))
  return(THETA)
}
#' Gibbs sampler for sigma squared
#'
#' This function samples new values for sigma squared from its
#' conditional posterior distribution given the current values 
#' of the other parameters in the model.
#' @param KAPPA the value of kappa in the prior on beta
#' @param V a hyper-parameter for the prior on sigma squared (nu)
#' @param LAM a hyper-parameter for the prior on sigma squared (lambda)
#' @param THETA the current values of theta
#' @param BETA the current values of beta
#' @keywords gibbs sampler for sigma squared
#' @export
#' @examples
#' SB_qual_sigma2sampler()
SB_qual_sigma2sampler = function(KAPPA,
                                 V,
                                 LAM,
                                 THETA,
                                 BETA)
{
  p <- length(THETA)
  THETA_SUM <- THETA * KAPPA + (1 - THETA)
  LAM_STAR <- LAM + 0.5 * sum(BETA ^ 2 / THETA_SUM)
  V_STAR <- V + p / 2
  return(1 / rgamma(1, shape = V_STAR, scale = 1 / LAM_STAR))
}
#' Gibbs sampler for beta
#'
#' This function samples new values for beta from its
#' conditional posterior distribution given the current values 
#' of the other parameters in the model.
#' @param KAPPA the value of kappa in the prior on beta
#' @param OMEGA the current values of the latent polya gamma variables
#' @param X the design matrix
#' @param Y the response variables
#' @param THETA the current values of theta
#' @param SIG2_0 the current values of sigma squared
#' @keywords gibbs sampler for beta
#' @export
#' @examples
#' SB_qual_betasampler()
SB_qual_betasampler = function(KAPPA,
                               OMEGA,
                               X,
                               Y,
                               THETA,
                               SIG2_0)
{
  U <- (Y - 0.5) / OMEGA
  p <- length(THETA)
  X1 <- sweep(X, 1, sqrt(OMEGA), "*")
  V <- crossprod(X1)
  diag(V) <- diag(V) + 1 / (SIG2_0 *
                              (THETA * KAPPA + (1 - THETA)))
  w <- crossprod(X1, U * sqrt(OMEGA))
  C <- chol(V)
  Z <- rnorm(p)
  BETA <- backsolve(C, Z + backsolve(C, w, transpose = TRUE))
  return(BETA)
}
#' Gibbs sampler
#'
#' This function runs our overall gibbs sampler for the 
#' spatial boost model making use of the data-augmentation 
#' strategy that exploits polya-gamma latent variables
#' @param GX the design matrix
#' @param Y the response variables
#' @param weights the normalized sums of the gene weights * relevances
#' @param SIG2 the initial value of sigma squared
#' @param KAPPA the value of kappa in the prior on beta
#' @param NSAMPS the desired number of samples
#' @param NCLUST the number of clusters to use in the parallel sampling
#' @param XI0 the hyper-parameter xi_0
#' @param XI1 the hyper-parameter xi_1
#' @param V a hyper-parameter in the prior on sigma squared (nu)
#' @param LAM a hyper-parameter in the prior on sigma squared (lambda)
#' @param k the number of additional covariates, defaults to 0
#' @param print boolean flag for whether or not to print information
#' @keywords gibbs sampler
#' @export
#' @examples
#' SB_qual_gibbssampler()
SB_qual_gibbssampler <- function(GX,
                                 Y,
                                 weights = NULL,
                                 SIG2 = 1e-2,
                                 KAPPA = 1e2,
                                 NSAMPS = 1e3,
                                 NCLUST = 1,
                                 XI0 = -2,
                                 XI1 = 0,
                                 V = NULL,
                                 LAM = NULL,
                                 k = 0,
                                 print = TRUE) {
  n <- nrow(GX)
  p <- ncol(GX) - k
  GGBOOST <- rep(0, p + k)
  if(is.null(weights)) {
    GGBOOST <- rep(-1 * XI0, p + k)
  } else {
    if(k > 0) {
      GGBOOST[1:k] <- -Inf
      GGBOOST[(k + 1):(p + k)] = -XI0 - XI1 * weights
    } else {
      GGBOOST <- -XI0 - XI1 * weights
    }
  }
  if(is.null(V) | is.null(LAM)) {
    V <- 3
    LAM <- 1e-1
  }
  if(print) {
    print("Starting clusters...")
  }
  cl = makeCluster(rep("localhost", NCLUST),type = "SOCK")
  if(print) {
    print("Initializing RNGs...")
  }
  clusterSetupRNG(cl)
  cluster.setup = function()
  {
    #### Load all necessary code and libraries ####
    library(spatialboost) # FIXME: this should be whatever the package is called
  }
  if(print) {
    print("Initializing cluster functions...")
  }
  clusterCall(cl, cluster.setup)
  GBETA <- matrix(0, nrow = NSAMPS, ncol = p + k)
  GTHETA <- matrix(0, nrow = NSAMPS, ncol = p + k)
  if(k > 0) {
    GTHETA[1, 1:k] <- 1
  }
  GSIG2_0 <- rep(SIG2, NSAMPS) # input (usually from EM)
  OMEGA <- matrix(1, nrow = NSAMPS, ncol = n)
  if(print) {
    print("Starting the Gibbs sampler for X...")
  }
  for(j in 2:NSAMPS)
  {
    GTHETA[j,] <- SB_qual_thetasampler(KAPPA,
                                       GGBOOST,
                                       GBETA[j - 1,],
                                       GSIG2_0[j - 1])
    if(k > 0) {
      GTHETA[j, 1:k] <- 1
    }
    GBETA[j,] <- SB_qual_betasampler(KAPPA,
                                     OMEGA[j - 1,],
                                     GX,
                                     Y,
                                     GTHETA[j,],
                                     GSIG2_0[j - 1])
    GSIG2_0[j] <- SB_qual_sigma2sampler(KAPPA,
                                        V,
                                        LAM,
                                        GTHETA[j,],
                                        GBETA[j,])
    OMEGA[j,] <- par.omega.final(GX, GBETA[j,], cl)
    if(j %% floor(0.1*NSAMPS) == 0)
    {
      if(print) {
        print(paste("Finished iteration ", j, sep = ""))
      }
    }
  }
  stopCluster(cl)
  return(list(THETA = GTHETA,
              BETA = GBETA,
              SIG2 = GSIG2_0,
              OMEGA = OMEGA))
}
#' Polya-Gamma variable (omega) sampler
#'
#' This function is one of two used to sample new values of 
#' the polya-gamma latent variables (omega) in parallel.
#' @param INPUT
#' @keywords polya-gamma sampler
#' @export
#' @examples
#' par.omega.sampler()
par.omega.sampler <- function(INPUT)
{
  return(genpg(INPUT[1], INPUT[2]))
}
#' Polya-Gamma variable (omega) final sampler
#'
#' This function is one of two used to sample new values of 
#' the polya-gamma latent variables (omega) in parallel.
#' @param X the design matrix
#' @param BETA the current values of beta
#' @param cl the cluster on which to run the sampler
#' @keywords polya-gamma sampler
#' @export
#' @examples
#' par.omega.final()
par.omega.final <- function(X, BETA, cl)
{
  input <- cbind(1, X %*% BETA)
  OMEGA <- parRapply(cl,
                    x = input,
                    fun = par.omega.sampler)
  return(OMEGA)
}
#' BFDR
#'
#' This function, given a vector of posterior probabilities, 
#' computes the Bayesian False Discovery Rate
#' @param ppa the vector of posterior probabilities
#' @keywords BFDR
#' @export
#' @examples
#' BFDR()
BFDR <- function(ppa) {
  N <- length(ppa)
  th.res <- rep(0, N)
  th.seq <- sort(ppa)
  for(i in 1:N)
  {
    ppa_t <- 1*(ppa >= th.seq[i])
    num <- sum((1 - ppa) * ppa_t)
    denom <- sum(ppa_t)
    if(denom > 0)
    {
      th.res[i] <- num / denom
    }
  }
  return(th.res)
}
#' Find optimal choice of kappa based on EMBFDR
#'
#' This function, given a matrix of posterior probabilities (PPA_MAT)
#' obtained when using a variety of kappa values (KAPPA_VEC), finds 
#' the optimal choice of kappa based on our BFDR criterion
#' @param PPA_MAT the matrix of conditional posterior probabilities
#' @param KAPPA_VEC the vector of kappa values used in the fitting
#' @param bfdr_threshold the maximum desired EMBFDR
#' @param prt_threshold the maximum desired percent of markers selected
#' @keywords EMBFDR
#' @export
#' @examples
#' BFDR_optimal()
BFDR_optimal <- function(PPA_MAT,
                         KAPPA_VEC,
                         bfdr_threshold,
                         prt_threshold) {
  B <- ncol(PPA_MAT)
  P <- nrow(PPA_MAT)
  BFDR <- apply(PPA_MAT, 2, BFDR)
  percent_retained <- function(Bcol) {
    bestrow <- which.min((Bcol - bfdr_threshold) ^ 2)
    return(1 - bestrow / P)
  }
  prts <- apply(BFDR, 2, percent_retained)
  optkap <- which.min((prts - prt_threshold) ^ 2)
  return(KAPPA_VEC[optkap])
}
#' Expectation-Maximization filter
#'
#' This function runs through the EM filter we propose and reduces 
#' the number of markers in an initial data set to a desired number. 
#' It returns a list containing the fitted model parameters and the 
#' posterior predictive loss (PPL) at each iteration of the filter.
#' We recommend stopping the EM filter when the PPL is minimized.
#' @param X the matrix of genotypes
#' @param Z the matrix of additional covariates (optional)
#' @param Y the vector of response variables
#' @param weights the normalized sums of gene weights * relevances
#' @param XI0 the hyper-parameter xi_0
#' @param XI1 the hyper-parameter xi_1
#' @param XI1_CONST the desired ratio of XI0 / XI1 to keep
#' @param KAPPA a vector of values of kappa to use in the fitting
#' @param V a hyper-parameter in the prior on sigma squared (nu)
#' @param LAM a hyper-parameter in the prior on sigma squared (lambda)
#' @param maxiter the maximum number of iterations for the EM algorithm
#' @param tol the tolerance to achieve convergence for the EM algorithm
#' @param print a boolean flag for whether or not to print information
#' @param k the number of additional covariates, defaults to 0
#' @param SVD_X a list object containing the SVD approximation to cbind(Z, X)
#' @param per2rm percent of markers to remove at each iteration of the filter
#' @param desiredP the final desired number of markers
#' @param bfdr_threshold the threshold on EMBFDR
#' @param ns_threshold the threshold on the number of markers to be selected
#' @keywords EM filter
#' @export
#' @examples
#' EM_qual_filter()
EM_qual_filter <- function(X = NULL,
                           Z = NULL,
                           Y,
                           weights = NULL,
                           XI0 = -2,
                           XI1 = 1,
                           XI1_CONST = NULL,
                           KAPPA = 1e3,
                           V = NULL,
                           LAM = NULL,
                           maxiter = 20,
                           tol = 1e-4,
                           print = TRUE,
                           k = 0,
                           SVD_X = NULL,
                           per2rm = 0.25,
                           desiredP = 1e2,
                           bfdr_threshold = 0.05,
                           ns_threshold = 10) {
  p0 <- NULL
  NKAP <- length(KAPPA)
  SVD_flag <- FALSE
  if(is.null(X) & !is.null(SVD_X)) {
    SVD_flag <- TRUE
    p0 <- nrow(SVD_X$V) - k
  }else {
    p0 <- ncol(X)
  }
  prt_th <- ns_threshold / p0
  NFILTERS <- ceiling((log(desiredP) - log(p0)) / log(1 - per2rm))
  if(NFILTERS < 1) {
    print("Error: Set a smaller number of desired variables to retain!")
    return(NULL)
  }
  if(print) {
    print(paste0("Running the EM filter 1 of ", NFILTERS, " times"))
  }
  if(is.null(XI1_CONST)) {
    XI1_CONST <- XI1 / XI0
  }
  EMRES <- list()
  EMRES[[1]] <- list()
  TEMPRES <- vector("list", NKAP)
  PPA_MAT <- matrix(0, nrow = p0, ncol = NKAP)
  for(kp in 1:NKAP) {
    if(print) {
      print(paste0("Fitting model using KAPPA = ", KAPPA[kp]))
    }
    TEMPRES[[kp]] <- EM_qual(X,
                             Z,
                             Y,
                             weights,
                             XI0,
                             XI1,
                             KAPPA[kp],
                             V,
                             LAM,
                             maxiter,
                             tol,
                             print,
                             k,
                             SVD_X)
    if(k == 0) {
      PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA
    } else {
      PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA[-(1:k)]
    }
  }
  K2USE <- BFDR_optimal(PPA_MAT,
                        KAPPA,
                        bfdr_threshold,
                        prt_th)
  KIND <- which(KAPPA == K2USE)
  EMRES[[1]] <- TEMPRES[[KIND]]
  EMRES[[1]]$KIND <- KIND
  EMRES[[1]]$UIND <- 1:p0
  EMRES[[1]]$GWS <- weights
  if(NFILTERS == 1) {
    return(EMRES)
  } else {
    C <- p0 * invlogit(XI0)
    for(i in 2:NFILTERS) {
      if(print) {
        print(paste0("Running the EM filter ", i, " of ", NFILTERS, " times"))
      }
      EMRES[[i]] <- list()
      p1 <- NULL
      if(k > 0) {
        p1 <- floor((1 - per2rm) * (length(EMRES[[i - 1]]$ETHETA) - k))
        PPA_MAT <- matrix(0, nrow = p1, ncol = NKAP)
        prt_th <- ns_threshold / p1
        XI0_new <- min(-1, logit(C / p1))
        XI1_new <- XI0_new * XI1_CONST
        sETHETA <- sort(EMRES[[i-1]]$ETHETA[-(1:k)],
                        decreasing = TRUE,
                        index.return = TRUE)$ix
        ind1 <- EMRES[[i-1]]$UIND[sETHETA[1:p1]]
        GWS <- quantile(weights, rank(weights[ind1])/length(ind1))
        SVD_new <- NULL
        if(SVD_flag) {
          SVD_new <- SVD_X
          SVD_new$V <- SVD_new$V[c(1:k, ind1 + k),]
          for(kp in 1:NKAP) {
            if(print) {
              print(paste0("Fitting model using KAPPA = ", KAPPA[kp]))
            }
            TEMPRES[[kp]] <- EM_qual(X,
                                     Z,
                                     Y,
                                     GWS,
                                     XI0 = XI0_new,
                                     XI1 = XI1_new,
                                     KAPPA[kp],
                                     V,
                                     LAM,
                                     maxiter,
                                     tol,
                                     print,
                                     k,
                                     SVD_X = SVD_new)
            PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA[-(1:k)]
          }
        } else {
          for(kp in 1:NKAP) {
            if(print) {
              print(paste0("Fitting model using KAPPA = ", KAPPA[kp]))
            }
            TEMPRES[[kp]] <- EM_qual(X[,ind1],
                                     Z,
                                     Y,
                                     GWS,
                                     XI0 = XI0_new,
                                     XI1 = XI1_new,
                                     KAPPA[kp],
                                     V,
                                     LAM,
                                     maxiter,
                                     tol,
                                     print,
                                     k,
                                     SVD_X = NULL)
            PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA[-(1:k)]
          }
        }
        K2USE <- BFDR_optimal(PPA_MAT,
                              KAPPA,
                              bfdr_threshold,
                              prt_th)
        KIND <- which(KAPPA == K2USE)
        EMRES[[i]] <- TEMPRES[[KIND]]
        EMRES[[i]]$KIND <- KIND
        EMRES[[i]]$UIND <- ind1 # always keeping track of SNP column
        EMRES[[i]]$GWS <- GWS
      } else {
        p1 <- floor((1 - per2rm) * length(EMRES[[i - 1]]$ETHETA))
        prt_th <- ns_threshold / p1
        PPA_MAT <- matrix(0, nrow = p1, ncol = NKAP)
        XI0_new <- min(-1, logit(C / p1))
        XI1_new <- XI0_new * XI1_CONST
        sETHETA <- sort(EMRES[[i-1]]$ETHETA,
                        decreasing = TRUE,
                        index.return = TRUE)$ix
        ind1 <- EMRES[[i-1]]$UIND[sETHETA[1:p1]]
        GWS <- quantile(weights, rank(weights[ind1])/length(ind1))
        SVD_new <- NULL
        if(SVD_flag) {
          SVD_new <- SVD_X
          SVD_new$V <- SVD_new$V[ind1,]
          for(kp in 1:NKAP) {
            if(print) {
              print(paste0("Fitting model using KAPPA = ", KAPPA[kp]))
            }
            TEMPRES[[kp]] <- EM_qual(X,
                                     Z,
                                     Y,
                                     GWS,
                                     XI0 = XI0_new,
                                     XI1 = XI1_new,
                                     KAPPA[kp],
                                     V,
                                     LAM,
                                     maxiter,
                                     tol,
                                     print,
                                     k,
                                     SVD_X = SVD_new)
            PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA
          }
        } else {
          for(kp in 1:NKAP) {
            if(print) {
              print(paste0("Fitting model using KAPPA = ", KAPPA[kp]))
            }
            TEMPRES[[kp]] <- EM_qual(X[,ind1],
                                     Y,
                                     GWS,
                                     XI0 = XI0_new,
                                     XI1 = XI1_new,
                                     KAPPA[kp],
                                     V,
                                     LAM,
                                     maxiter,
                                     tol,
                                     print,
                                     k,
                                     SVD_X = NULL)
            PPA_MAT[,kp] <- TEMPRES[[kp]]$ETHETA
          }
        }
        K2USE <- BFDR_optimal(PPA_MAT,
              KAPPA,
              bfdr_threshold,
              prt_th)
        KIND <- which(KAPPA == K2USE)
        EMRES[[i]] <- TEMPRES[[KIND]]
        EMRES[[i]]$KIND <- KIND
        EMRES[[i]]$UIND <- ind1
        EMRES[[i]]$GWS <- GWS
      }
    }
  }
  return(EMRES)
}
#### Range parameter (phi) selecting function ####
#' Range parameter (phi) selecting function
#'
#' This function breaks a chromosome into sub-regions and, for 
#' each block, determines the optimal value of the range parameter, 
#' phi, from a given set of values.
#' @param X the matrix of genotypes
#' @param SX the vector of marker positions
#' @param breakpoint the minimum distance between adjacent blocks
#' @param phigrid a vector of possible values for phi 
#' @keywords phi
#' @export
#' @examples
#' phi_optimal()
phi_optimal <- function(X,
                        SX,
                        breakpoint = 3e4,
                        phigrid = c(1e4,2.5e4,5e4,7.5e4,1e5)) {
  p <- ncol(X)
  diffvec <- diff(SX)
  optphivec <- rep(0, p)
  nphi <- length(phigrid)
  bps <- which(diffvec >= breakpoint)
  nbps <- length(bps)
  makeB <- function(s1, s2, phi) {
    2 * (1 - pnorm(abs(s2 - s1) / phi))
  }
  if(nbps > 0) {
    if(bps[nbps] < p) {
      bps <- c(bps, p)
      nbps <- nbps + 1
    }
    for(b in 1:nbps) {
      MSEvec <- rep(0, nphi)
      print(paste0("Fitting phi on block ", b, " of ", nbps))
      block <- NULL
      if(b == 1) {
        block <- 1:bps[b]
      } else {
        block <- (bps[b - 1] + 1):bps[b]
      }
      A <- abs(cor(X[, block, drop = FALSE]))
      for(i in 1:nphi) {
        print(paste0("Fitting phi = ", phigrid[i]))
        B <- outer(SX[block],
                   SX[block],
                   makeB,
                   phigrid[i])
        MSEvec[i] <- mean((A - B) ^ 2)
      }
      optphivec[block] <- rep(phigrid[which.min(MSEvec)],
                              length(block))
    }
  } else {
    MSEvec <- rep(0, nphi)
    print("No breakpoints detected; fitting phi on the full X")
    A <- abs(cor(X))
    for(i in 1:nphi) {
      print(paste0("Fitting phi = ", phigrid[i]))
      makeB <- function(s1, s2, phi) {
        2 * (1 - pnorm(abs(s2 - s1) / phi))
      }
      B <- outer(SX,
                 SX,
                 makeB,
                 phigrid[i])
      MSEvec[i] <- mean((A - B) ^ 2)
    }
    optphivec <- rep(phigrid[which.min(MSEvec)], p)
  }
  return(optphivec)
}
#' Complete Analysis Using the Spatial Boost Model
#'
#' This function runs through the EM filter we propose and reduces 
#' the number of markers in an initial data set to a desired number
#' and then runs the Gibbs sampler on the final set of retained markers. 
#' It returns a list containing the fitted model parameters and the 
#' posterior predictive loss (PPL) at each iteration of the filter, 
#' as well as a list containing the samples drawn from the gibbs sampler.
#' @param X the matrix of genotypes
#' @param Z the matrix of additional covariates (optional)
#' @param Y the vector of response variables
#' @param weights the normalized sums of gene weights * relevances
#' @param XI0 the hyper-parameter xi_0
#' @param XI1 the hyper-parameter xi_1
#' @param XI1_CONST the desired ratio of XI0 / XI1 to keep
#' @param KAPPA a vector of values of kappa to use in the fitting
#' @param V a hyper-parameter in the prior on sigma squared (nu)
#' @param LAM a hyper-parameter in the prior on sigma squared (lambda)
#' @param maxiter the maximum number of iterations for the EM algorithm
#' @param tol the tolerance to achieve convergence for the EM algorithm
#' @param print a boolean flag for whether or not to print information
#' @param k the number of additional covariates, defaults to 0
#' @param SVD_X a list object containing the SVD approximation to cbind(Z, X)
#' @param per2rm percent of markers to remove at each iteration of the filter
#' @param desiredP the final desired number of markers
#' @param bfdr_threshold the threshold on EMBFDR
#' @param ns_threshold the threshold on the number of markers to be selected
#' @param NCLUST the number of clusters to use in the parallel sampling of omega
#' @param NSAMPS the desired number of samples
#' @param FINDEX the index of the EM filter to use
#' @keywords EM filter
#' @export
#' @examples
#' SB_qual_complete_analysis()
SB_qual_complete_analysis <- function(X = NULL,
                                      Z = NULL,
                                      Y = NULL,
                                      weights = NULL,
                                      EMXI0 = -2,
                                      EMXI1 = 0,
                                      XI1_CONST = NULL,
                                      GSXI0 = -2,
                                      GSXI1 = 0,
                                      KAPPA = 1e4,
                                      V = 3,
                                      LAM = 1e-1,
                                      maxiter = 1e2,
                                      tol = 1e-4,
                                      print = TRUE,
                                      k = 0,
                                      SVD_X = NULL,
                                      per2rm = 0.25,
                                      desiredP = 1e2,
                                      bfdr_threshold = 0.05,
                                      ns_threshold = 10,
                                      NCLUST = 1,
                                      NSAMPS = 1e3,
                                      FINDEX = NULL) {
  if((k > 0) & is.null(Z)) {
    print(paste0("Error: You have specified ", k, " additional covariates, but Z is NULL!"))
    return(NULL)
  } 
  if(!is.null(Z)) {
    if(k != ncol(Z)) {
      print(paste0("Error: You have specified ", k, " additional covariates, but Z has ", ncol(Z), " columns!"))
      return(NULL)
    }
    if(nrow(X) != nrow(Z)) {
      print(paste0("Error: X has ", nrow(X), " rows, but Z has ", nrow(Z), " rows!"))
    }
  }
  SVD_flag <- FALSE
  if(!is.null(SVD_X)) {
    SVD_flag <- TRUE
    if(is.null(X)) {
      print("Please specify X so that the Gibbs sampler can run correctly!")
      return(NULL)
    } else {
      if(nrow(SVD_X$V) != (ncol(X) + k)) {
        print("Error: Size of SVD does not match total number of markers and additional covariates!")
        return(NULL)
      }
    }
  }
  EMRES <- NULL
  if(SVD_flag) {
    EMRES <- EM_qual_filter(X = NULL,
                            Z,
                            Y,
                            weights,
                            XI0 = EMXI0,
                            XI1 = EMXI1,
                            XI1_CONST,
                            KAPPA,
                            V,
                            LAM,
                            maxiter,
                            tol,
                            print,
                            k,
                            SVD_X,
                            per2rm,
                            desiredP,
                            bfdr_threshold,
                            ns_threshold)
  } else {
    EMRES <- EM_qual_filter(X,
                            Z,
                            Y,
                            weights,
                            XI0 = EMXI0,
                            XI1 = EMXI1,
                            XI1_CONST,
                            KAPPA,
                            V,
                            LAM,
                            maxiter,
                            tol,
                            print,
                            k,
                            SVD_X,
                            per2rm,
                            desiredP,
                            bfdr_threshold,
                            ns_threshold)
  }
  NFILT <- length(EMRES)
  PPLVEC <- rep(0, NFILT)
  for(j in 1:NFILT) {
    PPLVEC[j] <- EMRES[[j]]$PPL
  }
  F2USE <- NULL
  if(is.null(FINDEX)) {
    F2USE <- which.min(PPLVEC)
  } else {
    F2USE <- FINDEX
  }
  FGWS <- quantile(weights,
                   rank(weights[EMRES[[F2USE]]$UIND]) / 
                     length(EMRES[[F2USE]]$UIND))
  GSRES <- NULL
  if(NSAMPS > 0) {
    if(k > 0) {
      GSRES <- SB_qual_gibbssampler(cbind(Z, X[,EMRES[[F2USE]]$UIND]),
                                    Y,
                                    FGWS,
                                    SIG2 = EMRES[[F2USE]]$SIG2,
                                    KAPPA[EMRES[[F2USE]]$KIND],
                                    NSAMPS,
                                    NCLUST,
                                    XI0 = GSXI0,
                                    XI1 = GSXI1,
                                    V,
                                    LAM,
                                    k,
                                    print)
    } else {
      GSRES <- SB_qual_gibbssampler(X[,EMRES[[F2USE]]$UIND],
                                    Y,
                                    FGWS,
                                    SIG2 = EMRES[[F2USE]]$SIG2,
                                    KAPPA[EMRES[[F2USE]]$KIND],
                                    NSAMPS,
                                    NCLUST,
                                    XI0 = GSXI0,
                                    XI1 = GSXI1,
                                    V,
                                    LAM,
                                    k,
                                    print)
    }
  }
  return(list(EMRES = EMRES,
              GSRES = GSRES))
}
#' Spatial Boost Term Calculation 
#'
#' This function calculates the normalized spatial boost term in our prior 
#' on theta, crossprod(wj(phi), r), given a set of marker positions 
#' and their chromosomes, as well as the starting and ending positions of 
#' features of interest (e.g. genes) and their relevances.
#' @param SX a vector of marker positions
#' @param CHRS a vector of marker chromosomes
#' @param RELBLOCKS a list of matrices containing, for each chromosome, 
#'                  a matrix of starting and ending positions, and relevances
#' @param phi the range parameter to use for each marker, defaults to 10,000
#' @keywords spatial boost
#' @export
#' @examples 
#' SX <- sort(runif(25, 0, 10))
#' CHRS <- rep(1, 25)
#' SITES <- list()
#' SITES[[1]] <- matrix(0, nrow = 5, ncol = 2)
#' SITES[[1]][1,] <- c(1, 3)
#' SITES[[1]][2,] <- c(2, 3)
#' SITES[[1]][3,] <- c(5, 6)
#' SITES[[1]][4,] <- c(7, 9)
#' SITES[[1]][5,] <- c(9.75, 10)
#' rel <- list()
#' rel[[1]] <- c(1, 1, 1, 1, 1)
#' RELBLOCKS <- list()
#' RELBLOCKS[[1]] <- relblock(SITES[[1]], rel[[1]])
#' spatialterm <- calculate_spatial_boost(SX,
#'                                        CHRS,
#'                                        RELBLOCKS,
#'                                        phi = 0.05)
#' plot(SX,
#'      spatialterm,
#'      xlab = "Marker position",
#'      ylab = "Spatial term",
#'      main = expression("Simple example with " * phi * " = .05"),
#'      pch = 20,
#'      col = "black",
#'      xlim = c(0, 10),
#'      ylim = c(-1, 1))
#' for(i in c(1,3,4,5)) {
#'   lines(x = SITES[[1]][i, 1:2],
#'         y = c(-.5, -.5),
#'         lwd = 2,
#'         col = "purple")
#' }
#' lines(x = SITES[[1]][2, 1:2],
#'       y = c(-.7, -.7),
#'       lwd = 2,
#'       col = "purple")
#' calculate_spatial_boost()
calculate_spatial_boost <- function(SX,
                                    CHRS = NULL,
                                    RELBLOCKS = NULL,
                                    phi = 1e4) {
  UCHRS <- sort(unique(CHRS))
  NCHRS <- length(UCHRS)
  p <- length(SX)
  if(length(phi) == 1) {
    phi <- rep(phi, p)
  }
  rel <- list()
  weights <- rep(0, p)
  if(is.null(RELBLOCKS)) {
    for(chr in 1:NCHRS) {
      rel[[chr]] <- rep(1, nrow(RELBLOCKS[[chr]]))
    }
  } else {
    for(chr in 1:NCHRS) {
      rel[[chr]] <- RELBLOCKS[[chr]][,3]
    }
  }
  tally <- 1
  for(chr in UCHRS) {
    chrind <- which(CHRS == chr)
    if(length(chrind) > 0) {
      nexons <- nrow(RELBLOCKS[[tally]])
      get.weights = function(exon.lr, snp.x, phi) {
        nl = exon.lr[1] - snp.x
        nr = exon.lr[2] - snp.x
        return(pnorm(nr, 0, phi) - pnorm(nl, 0, phi))
      }
      for(s in chrind) {
        ind = which.min(abs(RELBLOCKS[[tally]][, 1] -
                              (SX[s] - 3 * phi[s]))):
          which.min(abs(RELBLOCKS[[tally]][, 2] -
                          (SX[s] + 3 * phi[s])))
        if(length(ind) > 0) {
          weights[s] = sum(apply(RELBLOCKS[[tally]][ind, 1:2, drop = FALSE],
                                 1,
                                 get.weights,
                                 SX[s],
                                 phi[s]) * rel[[tally]][ind])
        }
      }
    }
    tally <- tally + 1
  }
  if(max(weights) > 0) {
    return(weights / max(weights))
  }else {
    return(weights)
  }
}
#' Relevance block calculation function
#'
#' This function breaks a given set of (possibly) overlapping 
#' features on a chromosome into unique blocks and averages the 
#' relevance values for the overlapping regions.
#' @param SITES the matrix of starting and ending positions of features of interest
#' @param rel the vector of relevance values for each feature
#' @keywords relevance block
#' @export
#' @examples
#' relblock()
relblock <- function(SITES,
                     rel) {
  SIND <- sort(SITES[,1], index.return = TRUE)$ix
  GENE <- SITES[SIND, ]  
  REL <- rel[SIND]
  OBP <- unique(sort(c(GENE)))
  NBP <- length(OBP)
  FX <- NULL
  tGENE <- c(OBP[1], OBP[2])
  lIND <- which(GENE[,1] <= tGENE[1])
  rIND <- which(GENE[,2] >= tGENE[2])
  iIND <- intersect(lIND, rIND)
  if(length(iIND) > 0) {
    FX <- rbind(FX, c(tGENE[1], tGENE[2], mean(REL[iIND])))
  }
  for(b in 2:(NBP-1)) {
    tGENE <- c(OBP[b], OBP[b+1])
    lIND <- which(GENE[,1] <= tGENE[1])
    rIND <- which(GENE[,2] >= tGENE[2])
    iIND <- intersect(lIND, rIND)
    if(length(iIND) > 0) {
      FX <- rbind(FX, c(tGENE[1], tGENE[2], mean(REL[iIND])))
    }
  }
  return(FX)
}
