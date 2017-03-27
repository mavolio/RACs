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
rDM <- function(n, m, size, gamma, alpha, beta, lambda, theta, sigma, shift=TRUE) {

  # Choose sampling probability for each of gamma species, where
  # beta increases species turnover between function calls
  # FIXME not exactly sure what beta will do
  prob <- rgamma(gamma, lambda)
  prob <- prob / sum(prob)
  pos <- sample.int(gamma, gamma, prob=(rank(prob)/gamma) ^ (1/beta))
  prob <- prob[pos]

  # Choose a subset of alpha (maximum local species richness) species
  # according to prob, for each of n communities
  J <- replicate(m, sample.int(gamma, alpha, prob=prob))
  J <- t(J)

  # Set the mean of a Dirichlet distribution as
  # a power of the probs of the alpha species within each 
  # sample, and the sum of the Dirichlet params to sigma
  prob <- prob^(1/theta)
  prob <- matrix(prob[J], nrow=m)
  prob <- prob / rowSums(prob)
  a <- sigma * prob

  # Sample one Dirichlet for each community
  p <- matrix(rgamma(m * alpha, a), nrow=m)
  p <- p / rowSums(p)
  
  # Sample n Multinomials for each community
  if (shift) {
    s <- as.integer(shift)
    abund <- apply(p, 1, rmultinom, n=n, size=(size - s*alpha)) + s
  } else {
    abund <- apply(p, 1, rmultinom, n=n, size=size)
  }
  abund <- array(abund, dim=c(alpha, n, m))
  abund <- aperm(abund, c(2,1,3))

  # Assign to species number (out of 1:gamma), by
  # populating a sparse array of abundances
  ijk <- which(abund > 0, arr.ind=TRUE)
  iJk <- ijk
  iJk[, 2] <- J[ijk[, c(3, 2)]]
  simple_sparse_array(i=iJk, v=abund[ijk], dim=c(n, gamma, m))
}

