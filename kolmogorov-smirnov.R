library(ggplot2)

## Kolmogorov-Smirnov Distance
# 
# I didn't use species' rank as the x-axis in this example. To use
# rank you have to resolve the question of whether "relative rank"
# is meant to normalize across sampling effort, whether your metric is 
# intended to capture differences in richness, and maybe others.
# In addition, the KS distance is defined in the context of continuous
# distributions, but rank is discrete.
#
# So this example uses frequency as the x-axis, and pretends each species
# is an independent sample from a continuous distribution on (0, 1). Then
# The "two-sample" KS distance is calculated using `stat::ks.test`.


# this is just a temporary function to get a random community
rcom <- function(a, n, m) {
  p <- rgamma(n, a)
  p <- p / sum(p)
  x <- rmultinom(1, m, p)
  return(x[x > 0,1])
}

# calculate sampled frequencies
x_1 <- rcom(0.01, 500, 100)
q_1 <- x_1 / sum(x_1)
x_2 <- rcom(0.4, 100, 100)
q_2 <- x_2 / sum(x_2)
df <- data.frame(
  sp=c(rep('a', length(x_1)), rep('b', length(x_2))),
  q=c(q_1, q_2)
)

# plot 
ggplot(df, aes(q, col=sp)) +
  stat_ecdf(geom = "step") + 
  scale_x_continuous(limits=c(0, 1))

ks <- ks.test(q_1, q_2, exact=FALSE)
ks$statistic

## But they're not independent
# i'm thinking about the multivariate nature of the response variable.
# one could suppose each vec of frequency is a sample from a dirichlet distribution
# then the stats test is whether they are from the same distribution or not?
# but the diversity is unknown in advance.
# so what do you use for infinite mixtures ... dirichlet process.
# have to think about that one.
