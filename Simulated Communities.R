library(tidyr)
library(dplyr)


' Generate a random integer partition through the Chinese
#' Restaurant Process (CRP).

#STEP1 - Run the function to generate communities

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

#have 10 replicates for each 9 community types.
rep<-1:10
kappa <- c(0.4, 2, 10)
m <- c(5, 20, 50) 
df <- expand.grid(kappa = kappa, m = m, rep=rep)

community=data.frame(row.names=1)

for(i in 1:length(df$rep)) {
  # generate a large community, a vector of species abundances
  full.community <- rCRP(n = 15000, kappa=df$kappa[i], m=df$m[i]) 
  # the total number of individuals
  population <- sum(full.community)
  # set a multinomial sampling probability according to relative abundance
  prob <- full.community / population
  # draw `n` community with `size` individuals, as if picking individuals (with
  # replacement) from the full.community
  samp.communities <- as.data.frame(rmultinom(n = 10, size = 1000, prob = prob))
  #label species, richness, evenness, and replicate
  samp.communities$species<-seq(1:nrow(samp.communities))
  samp.communities$kappa<-df$kappa[i]
  samp.communities$m<-df$m[i]
  samp.communities$rep<-df$rep[i]
  samp.communities$ComType<-paste(samp.communities$kappa, samp.communities$m, sep="_")
  
  community=rbind(samp.communities, community)  
}

comm2<-community%>%
  gather(key=timestep, abundance, V1:V10)

write.csv(comm2, "~/Documents/SESYNC/SESYNC_RACs/R Files/SimCom.csv")

###just have a single community with no replicates
kappa <- c(0.4, 2, 10)
m <- c(5, 20, 50) 
df <- expand.grid(kappa = kappa, m = m)

community1=data.frame(row.names=1)

for(i in 1:9) {
  # generate a large community, a vector of species abundances
  full.community <- rCRP(n = 15000, kappa=df$kappa[i], m=df$m[i]) 
  # the total number of individuals
  population <- sum(full.community)
  # set a multinomial sampling probability according to relative abundance
  prob <- full.community / population
  # draw `n` community with `size` individuals, as if picking individuals (with
  # replacement) from the full.community
  samp.communities <- as.data.frame(rmultinom(n = 10, size = 1000, prob = prob))
  #label species, richness, evenness, and replicate
  samp.communities$species<-seq(1:nrow(samp.communities))
  samp.communities$kappa<-df$kappa[i]
  samp.communities$m<-df$m[i]
  samp.communities$ComType<-paste(samp.communities$kappa, samp.communities$m, sep="_")
  
  community1=rbind(samp.communities, community1)  
}

comm3<-community1%>%
  gather(key=timestep, abundance, V1:V10)

write.csv(comm3, "~/Documents/SESYNC/SESYNC_RACs/R Files/SimCom_noreps.csv")

