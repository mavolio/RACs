#' Sample a Dirichlet Multinomial distribution over a random subset of integers 
#' ranging from one to a maximum (gamma) diversity to simulate a local assemblage.
#'
#' @n Number of replicates for each of `m` samples.
#' @m Number of samples, potentially representing `n` repeated measures of 
#' abundance within the same community.
#' @size Number of individuals in a sample.
#' @gamma "Global" species richness (maximum richness in pooled samples)
#' @alpha "Local" species richness (richness in each sample)
#' @beta Variability among the alpha species chosen from a pool of size gamma
#' @theta Similarity of expected relative abundances across alpha species
#' @sigma Variability among sampled relative abundances across `n` samples
#' @shift Adjust abundance of all species after dirichlet multinomial sample; the default shift
#' adds one to enforce richness equal to `alpha`
#'
#' @import slam
#' 
#' @example 
#' # Generate 5 samples with a fixed parameterization
#' samp <- rDM(n=25, m=1000, alpha=10, beta=1, gamma=20, theta=0, sigma=1)
#' as.matrix(samp)
#' 
#' @example 
#' # Compare the relationship between evenness and richness
#' # across a range of parameter settings.
#' library(ggplot2)
#' library(vegan)
#' param <- expand.grid(
#'   alpha=2^(1:7),                          # Local richness
#'   beta=c(0.1, 1, 10),                     # Variability in species labels
#'   theta=c(0.01, 0.1, 1, 10, 100, 1000),   #
#'   sigma=c(1, 10, 100, 1000)               #
#' )
#' sims <- mapply(rDM, n=25, m=1000, gamma=max(param$alpha),
#'                alpha=param$alpha,
#'                beta=param$beta,
#'                theta=param$theta,
#'                sigma=param$sigma,
#'                SIMPLIFY=FALSE
#' )
#' param$avg_richness <- sapply(sims,
#'                              function(sim) mean(specnumber(as.matrix(sim)))
#' )
#' param$avg_evenness <- sapply(sims,
#'                              function(sim) {
#'                                x <- as.matrix(sim)
#'                                mean((diversity(x, 'invsimpson') - 1) / (specnumber(x) - 1))
#'                              })
#' ggplot(param, aes(x=avg_richness, y=avg_evenness, color=log10(theta))) +
#'   facet_grid(sigma ~ beta) +
#'   geom_point() +
#'   scale_y_continuous(limits=c(0, 1)) +
#'   scale_x_continuous(trans='log2')
#' 
#' @example 
#' # Compare pairwise similarities for samples within and across function calls
#' library(ggplot2)
#' library(vegan)
#' library(dplyr)
#' param <- expand.grid(
#'   rep=1:3,
#'   beta=c(0.1, 1, 10),
#'   sigma=c(1, 10, 100, 1000)
#' )
#' sims <- mapply(rDM, n=3, m=1000, alpha=4, gamma=8, theta=1000,
#'                beta=param$beta,
#'                sigma=param$sigma,
#'                SIMPLIFY=FALSE)
#' sims <- sapply(sims, as.matrix, simplify='array')
#' param <- param[rep(1:nrow(param), each=dim(sims)[[1]]), ]
#' param$id <- rownames(param)
#' rownames(param) <- NULL
#' sims <- aperm(sims, c(1, 3, 2))
#' dim(sims) <- c(nrow(param), dim(sims)[[3]])
#' param <- cbind(param, sims)
#' param <- select(param, -rep) %>%
#'   group_by(beta, sigma) %>%
#'   do({
#'     x <- select(., -beta, -sigma, -id)
#'     y <- combn(sub('\\..*', '', .$id), 2)
#'     data.frame(
#'       within=(y[1,]==y[2,]),
#'       bray=as.numeric(vegdist(x)))
#'   })
#' ggplot(param, aes(x=within, y=bray)) +
#'   facet_grid(sigma ~ beta) +
#'   geom_point() +
#'   xlab('Repeated measure from one assemblage') +
#'   ylab('Bray-Curtis Dissimilarity')
#'   
rabundance <- function(n, size, sites, iterations, alpha, gamma, gamma_rho, site_rho, sigma, theta, shift=TRUE) {
  require(MASS)
#  require(tidyr)
#  require(slam)

  jj <- gamma * sites
  kk <- iterations + 1
  
  # initial parameters
  a <- diag(sigma, jj)
  bdf <- expand.grid(species=1:gamma, site=1:sites)
  bdf <- lapply(bdf, factor)
  b <- model.matrix( ~ -1 + species, data=bdf)
  b <- cbind(b, model.matrix( ~ -1 + site, data=bdf))
  
  # species correlation
  gamma_mu <- matrix(0, nrow=gamma)
  gamma_Sigma <- matrix(gamma_rho, nrow=gamma, ncol=gamma)
  diag(gamma_Sigma) <- 1
  gamma_Sigma <- t(chol(gamma_Sigma))
  
  # site correlation
  site_mu <- matrix(0, nrow=sites)
  site_Sigma <- matrix(site_rho, nrow=sites, ncol=sites)
  diag(site_Sigma) <- 1
  site_Sigma <- t(chol(site_Sigma))
  
  # update b with correlations
  zero <- matrix(0, nrow=gamma, ncol=sites)
  Sigma <- rbind(
    cbind(gamma_Sigma, zero),
    cbind(t(zero), site_Sigma))
  b <- b %*% Sigma
  nn <- gamma + sites

  # initialize at long run expectation
  # (assuming -1 < eigs(a) < 1)
  Sigma <- solve(diag(1, jj) - a %*% t(a), b %*% t(b))
  x <- matrix(0, nrow=jj, ncol=kk)
  x[, 1] <- mvrnorm(1, mu=rep(0, jj), Sigma=Sigma)
  
  # simulate
  for (i in 2:kk) {
    x[, i] <- a %*% x[, i-1] + b %*% rnorm(nn, 0, 1)
  }
  x <- x[, 2:kk]
  dim(x) <- c(gamma, sites, iterations)
  
  # exp transform
  x <- exp(x)
  
  # select alpha species from x
  I <- apply(x, c(2, 3), function(x) sample.int(gamma, alpha, prob = x))
  ijk <- which(I > 0, arr.ind=TRUE)
  ijk[, 1] <- c(I)
  x <- x[ijk]
  dim(x) <- c(alpha, sites, iterations)
  
  # apply evenness power transform
  x <- x^(1/theta)
  
  # normalize to probabilities
  p <- x / rep(apply(x, c(2, 3), sum), each=alpha)
  
  # sample n multinomials for each site and iteration
  if (shift) {
    s <- as.integer(shift)
    abund <- apply(p, c(2, 3), rmultinom, n=n, size=(size - s * alpha)) + s
  } else {
    abund <- apply(p, c(2, 3), rmultinom, n=n, size=size)
  }
  dim(abund) <- c(alpha, n, sites, iterations)

  # assign to species number (out of 1:gamma)
  ijkl <- which(abund > 0, arr.ind=TRUE)
  Ijkl <- ijkl
  Ijkl[ , 1] <- I[matrix(ijkl[ , c(1, 3, 4)], ncol=3)]
  
  # return result as data frame
  result <- data.frame(Ijkl)
  result <- as.data.frame(lapply(result, factor))
  colnames(result) <- c('species', 'rep', 'site', 'iteration')
  result$abundance <- abund[ijkl]
  return(result)
}

