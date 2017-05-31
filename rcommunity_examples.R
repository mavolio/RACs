library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

## Compare the relationship between evenness and richness
## across a range of parameter settings.
# set parameter ranges
param <- expand.grid(
  rep = 1:25,
  alpha = 2^(1:6),                        
  theta = c(0.1, 0.5, 1, 2, 10),          
  sigma = c(0, 0.5, 0.99),
  beta = c(0, 0.1, 0.2, 0.4, 0.6)
)
# run simulations
sims <- mapply(rcommunity, n = 1, size = 1000,
               gamma = max(param$alpha),
               alpha = param$alpha,
               theta = param$theta,
               sigma = param$sigma,
               beta = param$beta,
               SIMPLIFY = FALSE
)
# calculate avg richness and evenness across parameter rep
mats <- lapply(sims, function(sim) {
  sim %>%
    spread('species', 'abundance') %>%
    select(-sample, -site, -iteration)
})
param$richness <- sapply(mats, function(mat) specnumber(mat))
param$evenness <- sapply(mats, function(mat) (diversity(mat, 'invsimpson') - 1) / (specnumber(mat) - 1))
param <- param %>%
  group_by(alpha, theta, sigma, beta) %>%
  summarise(avg_richness = mean(richness), avg_evenness = mean(evenness))
# visualize result
param$theta <- factor(param$theta)
ggplot(param, aes(x = avg_richness, y = avg_evenness, color = theta)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(trans = 'log2') +
  facet_grid(beta ~ sigma, labeller = label_both)
# Comments:
# The plots confirm the model design: `sigma` and `beta` have no effect on site richness or evenness., the 
# parameter `alpha` sets the richness level exactly (when `shift` takes the default value of TRUE), and
# parameter `theta` corresponds to evenness (here inverse Simpson shifted and scaled to span [0, 1]) of the abundance distribution.

## Compare pairwise similarities for comparing samples across
## sites versus across iterations.
# set parameter ranges
param <- expand.grid(
  rep = 1:25,
  theta = c(0.1, 0.5, 1, 2, 10),          
  sigma = c(0, 0.5, 0.99),
  beta = c(0, 0.1, 0.2, 0.4, 0.6)
)
# run simulations
sims <- mapply(rcommunity, n=1, size=500, sites=3, iterations=3,
               alpha=16, gamma=32,
               theta=param$theta,
               sigma=param$sigma,
               beta=param$beta,
               SIMPLIFY=FALSE
)
# calculate avg community dissimilarity across reps, separately for across-site pairs and across-iteration pairs.
sims <- lapply(sims, function(sim) {
  sim <- spread(sim, 'species', 'abundance', fill = 0, sep='_')
  dist_across_iteration <- sim %>%
    group_by(site) %>%
    do({
      x <- select(., starts_with('species'))
      bc <- as.matrix(vegdist(x))[-1,]
      data.frame(bray=diag(bc))
    }) %>%
    ungroup() %>%
    summarise(across='iteration', avg_bray = mean(bray))
  dist_across_site <- sim %>%
    group_by(iteration) %>%
    do({
      x <- select(., starts_with('species'))
      bc <- as.matrix(vegdist(x))[-1,]
      data.frame(bray=diag(bc))
    }) %>%
    ungroup() %>%
    summarise(across='site', avg_bray = mean(bray))
  rbind(dist_across_site, dist_across_iteration)
})
param <- param[rep(seq_len(nrow(param)), each = 2), ]
param <- cbind(param, do.call('rbind', sims))
# visualize result
param$theta <- factor(param$theta)
ggplot(param, aes(x=across, y=avg_bray, fill=theta)) +
  facet_grid(beta ~ sigma, labeller = label_both) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0, 1)) +
  xlab('Comparisons across ...') +
  ylab('Bray-Curtis Dissimilarity')
# Comments:
# The main effect of increasing `sigma` is to increase similarity across iterations, and it slightly increases the variability of across site disimilarity.
# The main effect of increasing `beta` is to increase dissimilarity across sites, and it does not effect similarity across iterations.
# The effect of increasing theta is context dependent, although overall it causes all pairwise dissimilarity to approach a constant value (here about 0.5). For high values of 
# theta, variation among species in sampling probability that exists between sites or iterations is washed out, as anticipated. These simulations do not show
# that `alpha` (here 16) and `gamma` (here 32) control the asymptotic (for high theta) value of disimilarity.


## Scenarios for RACs project

# targets
# richness: 5, 20, 50
# evenness (Simpsons)

param <- expand.grid(
  alpha = c(5, 20, 50)
)
  

samp <- rcommunity(n=1, size=1000, sites=2, iterations=1,
                   alpha=80, gamma=100, beta=1, sigma=0.9999, theta=0.1)
spread(samp, 'species', 'abundance', fill = 0, sep='_') %>%
  select(starts_with('species')) %>%
  as.matrix() %>%
  vegdist()

### argh
size <- 100; sites <- 2; iterations <- 2; alpha <- 64; gamma <- 128; beta <- 0; sigma <- 0.3; theta <- 1; n <- 1
df <- spread(result, 'species', 'abundance', fill = 0, sep='_')
df[, 1:5]
vegdist(select(df, -rep, -site, -iteration))

### introspection plots
gamma <- 2
sites <- 3
n <- gamma * sites
m <- 10 + 1

# parameters
x <- matrix(0, nrow=n, ncol=m)
a <- diag(0.1, n)
b <- diag(0.1, n)

# initialize at long run expectation
# (assuming -1 < eigs(a) < 1)
x[, 1] <- init

# simulate
for (i in 2:m) {
  x[, i] <- a %*% x[, i-1] + b %*% rnorm(n, 0, 1)
}

## convergence of expectation
# x_hat[, 1] = a %*% x_hat[, 0]
# x_hat[, 2] = a %*% a %*% x_hat[, 0]
# ...
# x_hat[, n] = a %^% n %*% x_hat[, 0]
# so, if -1 < eigs(a) < 1, converges to zero

## convergence of covariance
# xxT_hat[1] = a %*% xxT_hat[0] %*% aT + bbT
# xxT_hat[2] = a %*% a %*% xxT_hat[0] %*% aT %*% aT + a %*% bbT %*% aT + bbT
# ...
# xxT_hat[n] = a %^% n %*% xxT_hat[0] %*% aT %^% n + sum_{i=k}^{n-1} a %^% k %*% bbT %*% aT %^% k
# so, as above first term disapears as n -> Inf ...
# xxT_hat[Inf] = sum_{k=0}^{Inf} a %^% k %*% bbT %*% aT %^% k


# plot x samples by site
#x <- mvrnorm(100, mu=rep(0, jj), Sigma=Sigma)
#x <- t(x)
df <- data.frame(factor(1:gamma), factor(rep(1:sites, each=gamma)))
df <- cbind(df, x)
colnames(df) <- c('species', 'site', 1:ncol(x))
df <- gather(df, 'iteration', 'x', -species, -site, convert=TRUE)
df <- spread(df, site, x, sep='_')
GGally::ggpairs(df, aes(colour = species))

# plot x phase
df <- data.frame(1:gamma, factor(rep(1:sites, each=gamma)))
df <- cbind(df, x)
colnames(df) <- c('species', 'site', 1:iterations)
df <- gather(df, "iteration", "x", -species, -site, convert=TRUE)
df <- spread(df, species, x, sep='_')
df <- arrange(df, iteration)
ggplot(df, aes(species_1, species_2, color=site)) +
  geom_text(aes(label=iteration)) +
  geom_path(alpha=0.2)


# plot p phase
df <- data.frame(site=factor(1:sites))
df <- cbind(df, p[1, ,])
df <- gather(df, "time", "p", -site, convert=TRUE)
df <- spread(df, site, p, sep='_')
df <- arrange(df, time)
ggplot(df, aes(site_1, site_2)) +
  geom_point(aes(color=time)) +
  geom_path(color='black')


sim <- rcommunity(n=1, size=1000, sites=1, iterations=5,
                  alpha=128, gamma=256,
                  theta=1,
                  sigma=-0.5,
                  beta=0)
sim <- spread(sim, 'species', 'abundance', sep='_', fill=0)
vegdist(as.matrix(sim[, -(1:3)]))

