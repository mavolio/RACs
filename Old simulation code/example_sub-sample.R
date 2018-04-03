# generate a large community, a vector of species abundances
full.community <- rCRP(n = 1000, kappa = 0.4, m = 13) 

# the total number of individuals
population <- sum(full.community)

# set a multinomial sampling probability according to relative abundance
prob <- full.community / population

# draw `n` community with `size` individuals, as if picking individuals (with
# replacement) from the full.community
samp.communities <- rmultinom(n = 3, size = 100, prob = prob)
