rep <- 1:10
kappa <- c(0.1, 1, 10)
m <- c(5, 60, 200) 
df <- expand.grid(rep = rep, kappa = kappa, m = m)
SAD <- mapply(rCRP, kappa = df$kappa, m = df$m, MoreArgs = list(n = 1500))
df$richness <- vapply(SAD, vegan::specnumber, 0)
df$evenness <- vapply(SAD,
  function(x) {(1 - vegan::diversity(x, index = "inv")) / (1 - vegan::specnumber(x))},
  0)
plot(df$richness, df$evenness, col = df$m, pch = df$kappa)

x<-rCRP(n=100, kappa=1, m=40)

