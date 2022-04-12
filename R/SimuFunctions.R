## based on R's  stats::rnbinom() function.
## nus=list(nu1, nu2, nu3): vectors of mean counts in the three groups;
## ns=c(n1, n2, n3): sample sizes in the three groups
## kappa, a are two shape parameters in the NBP model.
#' First type of simulation scheme designed in `superdelta2` package
#' to generate random number of read count data from a Negative Binomial Poisson (NBP)
#' model with 3 groups, given mean counts parameter.
#'
#' This function is developed based on R's  stats::rnbinom() function.
#' `nus = list(nu1, nu2, nu3)` are the vectors of mean counts in the three groups.
#' `ns = c(n1, n2, n3)` are the sample sizes in the three groups.
#' `kappa` and `a` are two shape parameters in the NBP model.
#' `l` and `u` are the lower and upper bounds (uniform distribution) of simulated sample specific noise.
#'
#' @param nus `nus = list(nu1, nu2, nu3)` are the mean parameter vectors of the th$$ree groups respectively.
#' @param ns `ns = c(n1, n2, n3)` are the sample sizes of the three groups respectively.
#' @param kappa The first shape parameter.
#' @param a The second shape parameter.
#' @param l Lower bound of uniform distribution of simulated sample specific noise alpha_j.
#' @param u Upper bound of uniform distribution of simulated sample specific noise alpha_j.
#'
#' @details An NBP distribution is an integer-valued distribution with three parameters,
#'   the location parameter `mu`, and two shape parameters `kappa` and `a`.
#'   The mean and variance of X ~ NBP(mu, kappa, a) are:
#'   $E(X) = mu$
#'   $Var(X) = mu + mu*kappa^a$
#'   Note that the NBP distribution describes a nonlinear relationship between the mean and variance of genes.
#'   Other technical details of the NBP distribution, such as the probability density function
#'   and its relationship with the negative binomial (NB) distribution are not covered here.
#'
#' @return This function returns a count matrix `Y` simulated by the scheme discussed above.
#'
#' @export
#'
#' @keywords datagen distribution
#'
#' @examples
#' ngenes <- 5000; n1 <- n2 <- n3 <- 50; ns <- c(n1,n2,n3)
#' nu1 <- nu2 <- nu3 <- rep(100, ngenes)
#' nu2[seq_len(600] <- 150; nu3[401:1000] <- 75)
#' nus <- list(nu1, nu2, nu3)
#' set.seed(2020)
#' SIM1 <- Sim1(nus, ns, kappa = 0.06, a = 2.2, l = 12, u = 30)
Sim1 <- function(nus, ns, kappa, a, l, u) {
  mm <- sapply(nus, length)
  G <- length(mm)
  if (stats::var(mm) != 0) {
    stop("nus must be a list of three vectors of mean counts with the same lengths.")
  } else {
    ngenes <- mm[1]
  }
  X <- NULL
  for (g in seq_len(G)) {
    gg <- nus[[g]]^(2 - a) / kappa
    pg <- gg / (nus[[g]] + gg)
    Xg <- t(sapply(seq_len(ngenes), function(i) stats::rnbinom(ns[1], size = gg[i], prob = pg[i])))
    X <- cbind(X, Xg)
  }
  ## add sample specific noise
  alphas <- stats::runif(sum(ns), min = l, max = u)
  Y <- round(sweep(X, 2, alphas, "*"))
  ## return the counts
  return(Y)
}

#' @export
Sim2 <- function(mus, ns, kappa, a, l, u) {
  mm <- sapply(mus, length)
  G <- length(mm)
  if (stats::var(mm) != 0) {
    stop("mus must be a list of three vectors of mean log-counts with the same lengths.")
  } else {
    ngenes <- mm[1]
  }
  Y <- NULL
  alphas <- c()
  for (g in seq_len(G)) {
    ## generate alpha.j
    alphas.g <-  stats::runif(ns[g], min = l, max = u)
    alphas <- c(alphas, alphas.g)
    nus.g <- 2^(mus[[g]]) %*% t(alphas.g)
    ## turn nus.g, kappa, and a into gamma_g and p_g in NB model
    gg <- nus.g^(2 - a) / kappa
    pg <- gg / (nus.g + gg)
    Yg <- matrix(mapply(stats::rnbinom, size = gg, prob = pg, n = 1), nrow = ngenes)
    Y <- cbind(Y, Yg)
  }
  return(list(counts = Y, alphas = alphas))
}
