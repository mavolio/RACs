#' Sample a Dirichlet Multinomial distribution over a random subset of integers 
#' ranging from one to a maximum (gamma) diversity to simulate a local assemblage.
#'
#' @n Number of samples, potentially representing `n` repeated measures of 
#' abundance within the same community.
#' @m Number of individuals in a sample.
#' @gamma "Global" species richness (maximum richness in pooled samples)
#' @alpha "Local" species richness (richness in each sample)
#' @beta Variability among the alpha species chosen from a pool of size gamma
#' @theta Similarity of expected relative abundances across alpha species
#' @sigma Variability among sampled relative abundances across `n` samples
#' @shift Adjust abundance of all species after dirichlet multinomial sample; the default shift
#' adds one to enforce richness equal to `alpha`
#'
#' @import Matrix
#' @examples
#' # Generate 5 samples with a fixed parameterization
#' rDM(n=5, m=1000, alpha=10, beta=10, gamma=20, theta=10, sigma=1000)
#' 
#' 
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
#' sims <- mapply(rDM, n=1, m=1000, gamma=max(param$alpha),
#'                alpha=param$alpha,
#'                beta=param$beta, 
#'                theta=param$theta,
#'                sigma=param$sigma,
#'                SIMPLIFY=FALSE
#' )
#' param$richness <- sapply(sims, function(sim) length(sim@x))
#' param$evenness <- sapply(sims, function(sim) (diversity(sim@x, 'invsimpson') - 1) / (length(sim@x) - 1))
#' ggplot(param, aes(x=richness, y=evenness, color=log10(theta))) +
#'   facet_grid(sigma ~ beta) +
#'   geom_point() +
#'   scale_y_continuous(limits=c(0, 1)) +
#'   scale_x_continuous(trans='log2')
#' 
#'     
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
rDM <- function(n, m, gamma, alpha, beta, theta, sigma, shift=TRUE) {

  ## does succession allow species turnover too?
  
  # Choose sampling probability for each of gamma species, where
  # beta increases species turnover between function calls
  prob <- rgamma(gamma, beta)
  prob <- prob / sum(prob)
  prob <- sort(prob, decreasing=TRUE)
  
  # Choose a subset of alpha (maximum local species richness) species
  # according to prob
  j <- sample(1:gamma, alpha, prob=prob)
  
  # Choose a mean relative abundance (on the alpha-simplex), where
  # theta increases evenness by centering the mean
  x <- rgamma(alpha, theta)
  x <- x / sum(x)
  
  # Choose n dirichlet multinomial samples, where
  # sigma decreases turnover between samples
  a <- t(matrix(rgamma(n * alpha, x * sigma), nrow=alpha))
  a <- a / rowSums(a)
  if (shift) {
    s <- as.integer(shift)
    abund <- apply(a, 1, rmultinom, n=1, size=(m-s*alpha)) + s
  } else {
    abund <- apply(a, 1, rmultinom, n=1, size=m)
  }
  
  # Assign to species number (out of 1:gamma), by
  # populating a sparse array of abundances
  ji <- which(abund > 0, arr.ind=TRUE)
  sparseMatrix(i=ji[, 2], j=j[ji[, 1]], x=abund[ji], dims=c(n, gamma))
}

