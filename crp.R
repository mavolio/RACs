library('vegan')
library('ggplot2')
library('dplyr')

# Two parameter function for sampling integer partitions, regretably known as
# the "Chinese restaurant process" (CRP).
# In the more common parameterization, 0 <= alpha < 1, theta > -alpha

crp.sample = function(size, theta, alpha = 0) {
  # Initialize with a single member of the first class.
  result <- c(1)
  k <- length(result)
  # Iterate according to sample size
  for (n in 2:size) {
    if (runif(1, 0, 1) < (theta + k * alpha) / (theta + n)) {
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
  return(result)
}

# With alpha = 0 (the default), the process samples the
# UNTB/Ewens distribution with `size` individuals observed.

J <- 1000
theta <- 5.2
x <- crp.sample(J, theta)
rad <- radfit(x)
plot(rad)

# Scatter of Evenness vs. Richness
m <- 100
J <- 1000
theta <- 10
alpha <- 0
df <- data.frame(richness = rep(NA, m), diversity = rep(NA, m))
for (i in 1:nrow(df)) {
  x <- crp.sample(J, theta = theta, alpha = alpha)
  df$diversity[i] <- diversity(x, "inv")
  df$richness[i] <- length(x)
}
df$even <- df$diversity/df$richness
ggplot(df, aes(x = richness, y = even)) +
  geom_point() +
  geom_smooth(method = 'lm')

# How do parameters relate to average richness and evenness?

## (Slowly) generate data for plot of richness and evenness surfaces
J <- 1000
rep <- 1:10
alpha <- seq(0, 0.99, 0.1)
theta = seq(0, 100, 5)
df <- expand.grid(list(rep = rep, theta = theta, alpha = alpha))
for (i in 1:nrow(df)) {
  x <- crp.sample(J, df$theta[i], df$alpha[i])
  df$diversity[i] <- diversity(x, "inv")
  df$richness[i] <- length(x)
}
avg_df <- df %>%
  filter(richness > 1) %>%
  mutate(even = (1 - diversity) / (1 - richness)) %>%
  group_by(theta, alpha) %>%
  summarize(avg_richness = mean(richness), avg_even = mean(even, na.rm=T))

## Mean Richness Surface
ggplot(avg_df, aes(x = theta, y = alpha, z = avg_richness)) +
  geom_raster(aes(fill = avg_richness))

## Mean Evenness Surface
ggplot(avg_df, aes(x = theta, y = alpha, z = avg_even)) +
  geom_raster(aes(fill = avg_even))

# The GUILDS package samples from the Etienne Sampling Formula, which
# converges to the Ewens Sampling Formula when the "immigration rate"
# is very large. With lower immigration rates, the function is sampling from
# the dispersal limited local community in Hubbell's thoery.
library(GUILDS)

## (Slowly) generate data for plot of richness and evenness surfaces
J <- 1000
rep <- 1:10
theta <- seq(1, 10, 1)
eye <- seq(1, 10 * J, 100)
df <- expand.grid(list(rep = rep, theta = theta, eye = eye))
for (i in 1:nrow(df)) {
  x <- generate.ESF(theta = df$theta[i], I = df$eye[i], J = J)
  df$diversity[i] <- diversity(x, "inv")
  df$richness[i] <- length(x)
}
avg_df <- df %>%
  filter(richness > 1) %>%
  mutate(even = (1 - diversity) / (1 - richness)) %>%
  group_by(theta, eye) %>%
  summarize(avg_richness = mean(richness), avg_even = mean(even, na.rm=TRUE))

## Mean Richness Surface
ggplot(avg_df, aes(x = theta, y = eye, z = avg_richness)) +
  geom_raster(aes(fill = avg_richness))

## Mean Evenness Surface
ggplot(avg_df, aes(x = theta, y = eye, z = avg_even)) +
  geom_raster(aes(fill = avg_even)) 

# The CRP has another, less common, parameterization, that more tightly
# contrains the species richness. It has no application (known to me) is
# ecological literature. In this parameterization of the
# CRP; alpha < 0, k in {1,2,3,...}, and theta = -k*alpha

J <- 1000
rep <- 1:10
k <- seq(2, 102, 10)
alpha <- seq(-3, -0.01, 0.01)
df <- expand.grid(list(rep = rep, k = k, alpha = alpha))
df$theta = abs(df$k * df$alpha)
for (i in 1:nrow(df)) {
  x <- crp.sample(J, df$theta[i], df$alpha[i])
  df$diversity[i] <- diversity(x, "inv")
  df$richness[i] <- length(x)
}
avg_df <- df %>%
  filter(richness > 1) %>%
  mutate(even = (1 - diversity) / (1 - richness)) %>%
  group_by(k, alpha) %>%
  summarize(avg_richness = mean(richness), avg_even = mean(even, na.rm=T))

## Mean Richness Surface
ggplot(avg_df, aes(x = k, y = alpha, z = avg_richness)) +
  geom_raster(aes(fill = avg_richness))

## Mean Evenness Surface
ggplot(avg_df, aes(x = k, y = alpha, z = avg_even)) +
  geom_raster(aes(fill = avg_even))

# Scatter of Evenness vs. Richness
m <- 100
J <- 1000
alpha <- -1
k <- 60
theta <- -k*alpha
df <- data.frame(richness = rep(NA, m), diversity = rep(NA, m))
for (i in 1:nrow(df)) {
  x <- crp.sample(J, theta = theta, alpha = alpha)
  df$diversity[i] <- diversity(x, "inv")
  df$richness[i] <- length(x)
}
df$even <- (1 - df$diversity) / (1 - df$richness)
ggplot(df, aes(x = richness, y = even)) +
  geom_point() +
  geom_smooth(method = 'lm')

---
  
## Notes
  
# 1. The evenness quantity shown is (1-D)/(1-R), where D is the Inverse Simpson index and R is richness. D ranges from 1 to R, so this scales evenness between [0, 1].
# 2. Dispersal limitation doesn't seem to have any affect on evenness, to my surprise.
# 3. In theory, density dependence should help with evenness, but there's no sampling formula for that.
  
---
  
  
