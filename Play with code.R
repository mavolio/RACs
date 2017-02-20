library(tidyr)
library(dplyr)

rep <- 1:10
kappa <- c(0.4, 2, 10)
m <- c(5, 20, 50, 100) 
df <- expand.grid(rep = rep, kappa = kappa, m = m)
SAD <- mapply(rCRP, kappa = df$kappa, m = df$m, MoreArgs = list(n = 15000))
df$richness <- vapply(SAD, vegan::specnumber, 0)
df$evenness <- vapply(SAD,
  function(x) {(1 - vegan::diversity(x, index = "inv")) / (1 - vegan::specnumber(x))},
  0)
plot(df$richness, df$evenness, col = df$m, pch = df$kappa)


####okay doing this for 1 community
take1<-as.data.frame(rCRP(n=1000, kappa=1, m=20))
names(take1)[1]<-paste("abundance")
sorted<-as.data.frame(take1[order(-take1$abundance),])
names(sorted)[1]<-paste("abundance")

take1a<-sorted%>%
  mutate(species=letters[seq(1:20)])%>%
  mutate(t1=abundance*rand(c(-5,5)))

community=data.frame(row.names=1)
time<-as.data.frame(seq(1:10))
names(time)[1]<-paste("time")

####
#loop to get dissimilarity for each house
  for(i in 1:length(time$time)) {
    #simulate community
    spool<-rCRP(n=500, kappa=1, m=50)
   
    community=rbind(spool, community)  
    
  }


