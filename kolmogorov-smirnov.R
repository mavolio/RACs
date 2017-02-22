library(ggplot2)
library(dplyr)

## Manual calculation of Kolmogorov-Smirnov Distance
## for cumulative abundance by relative rank curves

# FIXME should relative rank start at zero? (and normalize by n-1)

df <- data.frame(
  plot_id=1,
  year=factor(c(1996, 1996, 1996, 2012, 2012, 2012, 2012)),
  rank=c(1/3, 2/3, 1, 1/4, 2/4, 3/4, 1),
  abun=c(0.5, 0.9, 1, 0.6, 0.8, 0.95, 1)
)

ggplot(df, aes(x=rank, y=abun, group=year)) +
  geom_point() +
  geom_step(aes(color=year)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(0, 1))

cc <- df %>%
  group_by(plot_id) %>%
  do({
    y <- unique(.$year)
    df1 <- filter(., year==y[[1]])
    df2 <- filter(., year==y[[2]])
    sf1 <- stepfun(df1$rank, c(0, df1$abun))
    sf2 <- stepfun(df2$rank, c(0, df2$abun))
    r <- sort(unique(c(0, df1$rank, df2$rank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    data.frame(
      Dmax=max(h),
      Dstar=sum(w*h))
  })

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

## I also want there to exist a abundance curve
# that is independent of species richness.
# With an infinite number of individuals, the curves for a perfectly
# even community of any size would be identical
# possibly the smoothed/expected rarefaction curve, normalized to one.

# for frequencies q_i, the expected number of species after n samples is
# i dunno ...
# s_0 = 
# s_1 = 1
# Pr[s_{n+1} == k| s_n] = {p1 if k = s_n, 1-p1 if k = s_n + 1



