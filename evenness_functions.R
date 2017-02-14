# very early work to develop new functions to calculate evenness

# implement the calculation of E_var from Smith and
# Wilson, 1996 (http://www.jstor.org/stable/3545749)
# 
# notes from Ian:
#   - no idea whether this is a unbiased estimator
#   - apparent that 0 abundances will be a problem

#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
#' @H the Shannon-Weiner diversity index
E_var <- function(x, S = length(x), N = sum(x), p = x / N, H = -sum(p * log(p))) {
  lnx <- log(x)
  theta <- (S - 1) / S * var(lnx)
  return(1 - 2 / pi * atan(theta))
}

# implement the calculation or __ from Maignan et
# al., 2003 (http://ssrn.com/abstract=389043)

# (a consistent estimator for the population Gini coefficient is on Wikipedia)
# there's also a gini function in the reldist package

x <- c(10, rep(10, 99))
reldist::gini(x)
E_var(x)

x <- c(10, rep(0, 99))
reldist::gini(x)
E_var(x)
