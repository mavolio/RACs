# very early work to develop new functions to calculate evenness

# implement the calculation of E_var from Smith and
# Wilson, 1996 (http://www.jstor.org/stable/3545749)
# 
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

#implement the calculation of E_q from Smith and Wilson 1996. ## email and ask about negative sign.
x <- c(80,40,20,10,1)

E_q<-function(x){ #this is how you say make a funciton, and x is the input into the function, here x is a vector of abundnaces
  r<-rank(x, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]#double bracket say take number only.
  2/pi*atan(b)
}
#.113

# implement the calculation or __ from Maignan et
# al., 2003 (http://ssrn.com/abstract=389043)

# (a consistent estimator for the population Gini coefficient is on Wikipedia)
# there's also a gini function in the reldist package

x <- c(80,40,20,10)
reldist::gini(x)
E_var(x)

x <- c(10, rep(0, 99))
reldist::gini(x)
E_var(x)

# a mult-plot example
library(dplyr)
df <- data.frame(
  plot=c(1, 1, 1, 2, 2, 2, 2),
  sp=c('a','b','c','a','b','c','d'),
  n=c(5,4,3,15,15,15,1))
plot_evenness <- group_by(df, plot) %>% summarize(E_q=E_q(n))
plot_evenness

