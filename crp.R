#' Generate a random integer partition through the Chinese
#' Restaurant Process (CRP).
#'
#' In the literature on the Unified Neutral Theory of Biodiveristy,
#' a special case of the CRP makes an appearance as the Ewens' Sampling Formula, a 
#' probability distribution on vectors of species' abundances. \code{rCRP}
#' includes simulating from this family of distributions as the default parameterization,
#' but allows the full range of parameters for non UNTB-based models of a
#' species abundance distribution (SAD).
#'
#' @n The integer being partitioned.
#' @theta The "concentration" parameter, $\theta > -\alpha$.
#' @alpha The "discount" parameter, $0 <= \alpha < 1$.
#' @kappa An alternative parameterization setting $\alpha = -\kappa < 0$.
#' @m A positive integer in the alternative parameterization.
#' @zeros Preserve zeros.
#' @return A vector of integers that sum to n. If zeros == TRUE and m is
#' not null, the vector is of length m.
#' @references Notation as in section 3.2 of \doi{10.1007/b11601500}.
#'
#' @examples
#' # Simulate a SAD with at most 42 species
#' rCRP(n = 100, kappa = 3.1, m = 42)
#' 
#' # Show effect of (alternative) parameters on richness and evenness.
#' rep <- 1:10
#' kappa <- c(0.5, 1, 10)
#' m <- c(30, 60, 90)
#' df <- expand.grid(rep = rep, kappa = kappa, m = m)
#' SAD <- mapply(rCRP, kappa = df$kappa, m = df$m, MoreArgs = list(n = 150))
#' df$richness <- vapply(SAD, vegan::specnumber, 0)
#' df$evenness <- vapply(SAD,
#'   function(x) {(1 - vegan::diversity(x, index = "inv")) / (1 - vegan::specnumber(x))},
#'   0)
#' plot(df$richness, df$evenness, col = df$m, pch = df$kappa)
#' 
rCRP = function(n, theta, alpha = 0, kappa = NULL, m = NULL, zeros = TRUE) {
  if (!is.null(kappa) & !is.null(m)) {
    if (m == round(m) & m > 0 & kappa > 0) {
      alpha <- -kappa
      theta <- kappa * m
    } else {
      stop("Parameter m must be a positive integer and kappa must be non-negative.")
    }
  } else {
    if (alpha < 0 | 1 <= alpha | -alpha < theta) { 
      stop("Without kappa or m, parameters must satisfy 0 <= alpha < 1 & theta > -alpha.")
    }
  }
  if (!is.null(kappa) & !is.null(m) & zeros) {
    # Sample so as to preserve zero abundances, given integer m
    result <- rep(0, m)
    extant <- which(result != 0)
    j <- m - length(extant)
    p <- rep(NA, m)
    idx <- sample(m, 1)
    for (k in 1:n) {
      result[[idx]] <- result[[idx]] + 1
      if (result[[idx]] == 1) {
        # respond to new class
        extant <- c(extant, idx)
        j <- j + 1
        p[-extant] <- (theta + j * alpha) / (m - j)
        p[[idx]] <- 1 - alpha
      } else {
        # step up the sampling prob for extant class
        p[[idx]] <- p[[idx]] + 1
      }
      # sample a class for the next object
      idx <- sample(m, 1, prob = p / (theta + k))
    }
  } else if (zeros) {
    stop("Preserving zeros not implemented for non-null kappa and m")
  } else {
    # Initialize with a single instance of the first class.
    result <- c(1)
    k <- length(result)
    # Iterate according to sample size
    for (j in 2:n) {
      if (runif(1, 0, 1) < (theta + k * alpha) / (theta + j)) {
        # Add a new class.
        result <- c(result, 1)
        k <- k + 1
      } else {
        # Add to an existing class,
        # with probability related to current class abundance
        i <- sample(k, size = 1, prob = result - alpha)
        result[[i]] <- result[[i]] + 1
      }
    }
  }
  return(result)
}