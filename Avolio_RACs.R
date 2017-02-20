library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)


#' Generate a random integer partition through the Chinese
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


####
#loop to get dissimilarity for each community #sppool 50
# 
# community50=data.frame(row.names=1)
# time<-as.data.frame(seq(1:10))
# names(time)[1]<-paste("time")
# 
# for(i in 1:length(time$time)) {
#   #simulate community
#   spool<-as.data.frame(rCRP(n=500, kappa=1, m=50))
#   names(spool)[1]<-paste("abundance")
#   spool$timestep<-time$time[i]
#   spool$species<-seq(1:50)
#   spool$rank=rank(-spool$abundance, ties.method="average")
#   
#   community50=rbind(spool, community50)  
# }
# 
# richness<-community50%>%
#   mutate(present=ifelse(abundance>1,1,0))%>%
#   tbl_df()%>%
#   group_by(timestep)%>%
#   summarise(richness=sum(present))
#     
# turnover<-turnover(df=community50, time.var="timestep", species.var="species", abundance.var="abundance", metric="total")
# 
# immigration<-turnover(df=community50, time.var="timestep", species.var="species", abundance.var="abundance", metric="appearance")
# 
# loss<-turnover(df=community50, time.var="timestep", species.var="species", abundance.var="abundance", metric="disappearance")
# 
# mrs<-rank_shift(df=community50, time.var = "timestep", species.var = "species", abundance.var="abundance")%>%
#   separate(year_pair, c("year1","year2"), sep="-", remove=F)%>%
#   select(-year_pair, -year1)
# 
# m<-merge(turnover, immigration, by="timestep")
# m1<-merge(m, loss, by="timestep")
# m2<-merge(m1, richness, by="timestep")
# codyn_metrics50<-merge(m2, mrs, by.x="timestep", by.y="year2")
# 
# ####do kendall rank correlation
# ##1 code all species by rank in year 1
# 
# year1<-community50%>%
#   filter(timestep==1)%>%
#   mutate(rankyr1=rank)%>%
#   select(-abundance, -timestep, -rank)
# 
# ##get matrix to correlate
# community_rank<-merge(year1, community50, by="species")%>%
#   select(-abundance)%>%
#   mutate(timestep2=paste("t", timestep, sep=""))%>%
#   select(-timestep)%>%
#   spread(timestep2, rank)
# 
# ###
# tao<-c(
# as.numeric(Kendall(community_rank$rankyr1, community_rank$t2)[1]),
# as.numeric(Kendall(community_rank$t2, community_rank$t3)[1]),
# as.numeric(Kendall(community_rank$t3, community_rank$t4)[1]),
# as.numeric(Kendall(community_rank$t4, community_rank$t5)[1]),
# as.numeric(Kendall(community_rank$t5, community_rank$t6)[1]),
# as.numeric(Kendall(community_rank$t6, community_rank$t7)[1]),
# as.numeric(Kendall(community_rank$t7, community_rank$t8)[1]),
# as.numeric(Kendall(community_rank$t8, community_rank$t9)[1]),
# as.numeric(Kendall(community_rank$t9, community_rank$t10)[1]))
# 
# codyn_metrics50$Cor<-tao
# codyn_metrics50$rich<-"high"
# 
# 
# #loop to get dissimilarity for each community #sppool 20
# 
# community20=data.frame(row.names=1)
# time<-as.data.frame(seq(1:10))
# names(time)[1]<-paste("time")
# 
# for(i in 1:length(time$time)) {
#   #simulate community
#   spool<-as.data.frame(rCRP(n=500, kappa=1, m=20))
#   names(spool)[1]<-paste("abundance")
#   spool$timestep<-time$time[i]
#   spool$species<-seq(1:20)
#   spool$rank=rank(-spool$abundance, ties.method="average")
#   
#   community20=rbind(spool, community20)  
# }
# 
# richness<-community20%>%
#   mutate(present=ifelse(abundance>1,1,0))%>%
#   tbl_df()%>%
#   group_by(timestep)%>%
#   summarise(richness=sum(present))
# 
# turnover<-turnover(df=community20, time.var="timestep", species.var="species", abundance.var="abundance", metric="total")
# 
# immigration<-turnover(df=community20, time.var="timestep", species.var="species", abundance.var="abundance", metric="appearance")
# 
# loss<-turnover(df=community20, time.var="timestep", species.var="species", abundance.var="abundance", metric="disappearance")
# 
# mrs<-rank_shift(df=community20, time.var = "timestep", species.var = "species", abundance.var="abundance")%>%
#   separate(year_pair, c("year1","year2"), sep="-", remove=F)%>%
#   select(-year_pair, -year1)
# 
# m<-merge(turnover, immigration, by="timestep")
# m1<-merge(m, loss, by="timestep")
# m2<-merge(m1, richness, by="timestep")
# codyn_metrics20<-merge(m2, mrs, by.x="timestep", by.y="year2")
# 
# ####do kendall rank correlation
# ##1 code all species by rank in year 1
# 
# year1<-community20%>%
#   filter(timestep==1)%>%
#   mutate(rankyr1=rank)%>%
#   select(-abundance, -timestep, -rank)
# 
# ##get matrix to correlate
# community_rank<-merge(year1, community20, by="species")%>%
#   select(-abundance)%>%
#   mutate(timestep2=paste("t", timestep, sep=""))%>%
#   select(-timestep)%>%
#   spread(timestep2, rank)
# 
# ###
# tao<-c(
#   as.numeric(Kendall(community_rank$rankyr1, community_rank$t2)[1]),
#   as.numeric(Kendall(community_rank$t2, community_rank$t3)[1]),
#   as.numeric(Kendall(community_rank$t3, community_rank$t4)[1]),
#   as.numeric(Kendall(community_rank$t4, community_rank$t5)[1]),
#   as.numeric(Kendall(community_rank$t5, community_rank$t6)[1]),
#   as.numeric(Kendall(community_rank$t6, community_rank$t7)[1]),
#   as.numeric(Kendall(community_rank$t7, community_rank$t8)[1]),
#   as.numeric(Kendall(community_rank$t8, community_rank$t9)[1]),
#   as.numeric(Kendall(community_rank$t9, community_rank$t10)[1]))
# 
# codyn_metrics20$Cor<-tao
# codyn_metrics20$rich<-"mid"
# 
# #loop to get dissimilarity for each community #sppool 5
# 
# community5=data.frame(row.names=1)
# time<-as.data.frame(seq(1:10))
# names(time)[1]<-paste("time")
# 
# for(i in 1:length(time$time)) {
#   #simulate community
#   spool<-as.data.frame(rCRP(n=500, kappa=1, m=5))
#   names(spool)[1]<-paste("abundance")
#   spool$timestep<-time$time[i]
#   spool$species<-seq(1:5)
#   spool$rank=rank(-spool$abundance, ties.method="average")
#   
#   community5=rbind(spool, community5)  
# }
# 
# richness<-community5%>%
#   mutate(present=ifelse(abundance>1,1,0))%>%
#   tbl_df()%>%
#   group_by(timestep)%>%
#   summarise(richness=sum(present))
# 
# turnover<-turnover(df=community5, time.var="timestep", species.var="species", abundance.var="abundance", metric="total")
# 
# immigration<-turnover(df=community5, time.var="timestep", species.var="species", abundance.var="abundance", metric="appearance")
# 
# loss<-turnover(df=community5, time.var="timestep", species.var="species", abundance.var="abundance", metric="disappearance")
# 
# mrs<-rank_shift(df=community5, time.var = "timestep", species.var = "species", abundance.var="abundance")%>%
#   separate(year_pair, c("year1","year2"), sep="-", remove=F)%>%
#   select(-year_pair, -year1)
# 
# m<-merge(turnover, immigration, by="timestep")
# m1<-merge(m, loss, by="timestep")
# m2<-merge(m1, richness, by="timestep")
# codyn_metrics5<-merge(m2, mrs, by.x="timestep", by.y="year2")
# 
# ####do kendall rank correlation
# ##1 code all species by rank in year 1
# 
# year1<-community5%>%
#   filter(timestep==1)%>%
#   mutate(rankyr1=rank)%>%
#   select(-abundance, -timestep, -rank)
# 
# ##get matrix to correlate
# community_rank<-merge(year1, community5, by="species")%>%
#   select(-abundance)%>%
#   mutate(timestep2=paste("t", timestep, sep=""))%>%
#   select(-timestep)%>%
#   spread(timestep2, rank)
# 
# ###
# tao<-c(
#   as.numeric(Kendall(community_rank$rankyr1, community_rank$t2)[1]),
#   as.numeric(Kendall(community_rank$t2, community_rank$t3)[1]),
#   as.numeric(Kendall(community_rank$t3, community_rank$t4)[1]),
#   as.numeric(Kendall(community_rank$t4, community_rank$t5)[1]),
#   as.numeric(Kendall(community_rank$t5, community_rank$t6)[1]),
#   as.numeric(Kendall(community_rank$t6, community_rank$t7)[1]),
#   as.numeric(Kendall(community_rank$t7, community_rank$t8)[1]),
#   as.numeric(Kendall(community_rank$t8, community_rank$t9)[1]),
#   as.numeric(Kendall(community_rank$t9, community_rank$t10)[1]))
# 
# codyn_metrics5$Cor<-tao
# codyn_metrics5$rich<-"low"
# 
# 
# codyn_metrics<-rbind(codyn_metrics5, codyn_metrics20, codyn_metrics50)
# pairs(codyn_metrics[5:7])


###what have we learned? MRS is highly dependant on richness, but Kendall's rank correlation is not. There is also no relationship between Kendall's rank correltion and MRS.

####OKAY. Trying to do this in a simpler way that looks at both richness and evenness

rep <- 1:10
kappa <- c(0.4, 2, 10)
m <- c(5, 20, 50) 
df <- expand.grid(rep = rep, kappa = kappa, m = m)

community=data.frame(row.names=1)

for(i in 1:length(df$rep)) {
  #simulate community
  spool<-as.data.frame(rCRP(n=15000, kappa=df$kappa[i], m=df$m[i]))
  names(spool)[1]<-paste("abundance")
  spool$timestep<-df$rep[i]
  spool$species<-seq(1:nrow(spool))
  spool$rank=rank(-spool$abundance, ties.method="average")
  spool$kappa<-df$kappa[i]
  spool$m<-df$m[i]
  
  community=rbind(spool, community)  
}
community$replicate<-paste(community$kapp, community$m, sep="_")

mrs<-rank_shift(df=community, time.var = "timestep", species.var = "species", abundance.var="abundance", replicate.var = "replicate")%>%
  separate(year_pair, c("year1","timestep"), sep="-", remove=F)%>%
  select(-year_pair, -year1)

richness<-community%>%
  mutate(present=ifelse(abundance>1,1,0))%>%
  tbl_df()%>%
  group_by(kappa, m, timestep)%>%
  summarise(Richness=sum(present))

replicate<-community%>%
  select(replicate)%>%
  unique()

evenness=data.frame(row.names=1)

for(i in 1:length(replicate$replicate)) {

  subset<-community%>%
    filter(replicate==replicate$replicate[i])%>%
  select(-rank)%>%
  spread(species, abundance, fill=0)

key<-subset[,1:4]
S<-specnumber(subset[,5:ncol(subset)])
InvD<-diversity(subset[,5:ncol(subset)],"inv")
SimpEven<-InvD/S
even<-cbind(key, SimpEven)

evenness=rbind(even, evenness)  
}

m<-merge(richness, evenness, by=c("timestep","kappa","m"))
m2<-merge(m, mrs, by=c("timestep","replicate"))

pairs(metrics[5:7])


####do kendall rank correlation
##1 code all species by rank in year 1

year1=data.frame(row.names=1)

for(i in 1:length(replicate$replicate)) {
  
  subset<-community%>%
    filter(replicate==replicate$replicate[i])%>%
  filter(timestep==1)%>%
  mutate(rankyr1=rank)%>%
  select(-abundance, -timestep, -rank)
  
  year1<-rbind(year1, subset)
}

##get matrix to correlate
community_rank<-merge(year1, community, by=c("species", "kappa","m","replicate"))%>%
  select(-abundance)%>%
  mutate(timestep2=paste("t", timestep, sep=""))%>%
  select(-timestep)%>%
  spread(timestep2, rank)

###
kendall_rank<-data.frame(row.names=1)

for(i in 1:length(replicate$replicate)) {
  
  subset<-community_rank%>%
    filter(replicate==replicate$replicate[i])
  
  tao<-as.data.frame(seq(from=2, to=10))
  names(tao)[1]<-paste("timestep")
  
tao_calc<-c(
  as.numeric(Kendall(subset$rankyr1, subset$t2)[1]),
  as.numeric(Kendall(subset$t2, subset$t3)[1]),
  as.numeric(Kendall(subset$t3, subset$t4)[1]),
  as.numeric(Kendall(subset$t4, subset$t5)[1]),
  as.numeric(Kendall(subset$t5, subset$t6)[1]),
  as.numeric(Kendall(subset$t6, subset$t7)[1]),
  as.numeric(Kendall(subset$t7, subset$t8)[1]),
  as.numeric(Kendall(subset$t8, subset$t9)[1]),
  as.numeric(Kendall(subset$t9, subset$t10)[1]))

tao$RankCor<-tao_calc
tao$replicate<-replicate$replicate[i]

kendall_rank<-rbind(kendall_rank, tao)
}

metrics<-merge(m2, kendall_rank, by=c("timestep","replicate"))%>%
  mutate(richness=as.character(m))%>%
  mutate(RelMRS=MRS/Richness)

pairs(metrics[c(5,6,7,8,10)])


theme_set(theme_bw(16))
rich_even<-ggplot(data=metrics, aes(x=Richness, y=SimpEven))+
  geom_point(size=3)

MRS_rich<-ggplot(data=metrics, aes(x=Richness, y=mrs))+
  geom_point(size=3)
MRS_even<-ggplot(data=metrics, aes(x=SimpEven, y=MRS,  color=richness))+
  geom_point(size=3)
grid.arrange(MRS_rich, MRS_even, ncol=2)

relMRS_rich<-ggplot(data=metrics, aes(x=Richness, y=RelMRS))+
  geom_point(size=3)
relMRS_even<-ggplot(data=metrics, aes(x=SimpEven, y=RelMRS,  color=richness))+
  geom_point(size=3)
grid.arrange(relMRS_rich, relMRS_even, ncol=2)

RC_rich<-ggplot(data=metrics, aes(x=Richness, y=RankCor))+
  geom_point(size=3)
RC_even<-ggplot(data=metrics, aes(x=SimpEven, y=RankCor,  color=richness))+
  geom_point(size=3)
grid.arrange(RC_rich, RC_even, ncol=2)

MRS_RC<-ggplot(data=metrics, aes(x=MRS, y=RankCor))+
  geom_point(size=3)
relMRS_RC<-ggplot(data=metrics, aes(x=RelMRS, y=RankCor))+
  geom_point(size=3)
grid.arrange(MRS_RC, relMRS_RC, ncol=2)


#####read in pplots and try to do top species for it
#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")
all<-read.csv("SpeciesRelativeAbundance_Nov2016.csv")
pplots<-all%>%
  filter(project_name=="pplots")%>%
  filter(treatment=="N2P0"|treatment=="N1P0")%>%
  select(-plot_id)
summedtrt<-pplots%>%
  tbl_df()%>%
  group_by(treatment_year, treatment, genus_species)%>%
  summarize(relcov=mean(relcov))
topspecies<-summedtrt%>%
  mutate(rank=rank(-relcov, ties.method="average"))%>%
  filter(rank<6)%>%
  tbl_df()%>%
  select(genus_species)%>%
  unique()

sp_subset<-merge(summedtrt, topspecies, by="genus_species")%>%
  spread(genus_species, relcov, fill=0)%>%
  gather(key="genus_species", value="relcov",3:16)%>%
  tbl_df()%>%
  group_by(treatment_year, treatment)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
        timestep=paste("t", treatment_year, sep=""))%>%
  tbl_df()%>%
  select(-relcov, -treatment_year)%>%
  spread(timestep, rank)

unique<-sp_subset%>%
  select(treatment)%>%
  unique()


kendall_rank<-data.frame(row.names=1)

for(i in 1:length(unique$treatment)) {
  
  subset<-sp_subset%>%
    filter(treatment==unique$treatment[i])
  
  tao<-as.data.frame(seq(from=2, to=12))
  names(tao)[1]<-paste("timestep")

#use cor on matrix of times and then extract lower diagonal
#m is the matics of cor outputs
# diag(m[2:n,1:n-1])
  tao_calc<-c(
    as.numeric(Kendall(subset$t1, subset$t2)[1]),
    as.numeric(Kendall(subset$t2, subset$t3)[1]),
    as.numeric(Kendall(subset$t3, subset$t4)[1]),
    as.numeric(Kendall(subset$t4, subset$t5)[1]),
    as.numeric(Kendall(subset$t5, subset$t6)[1]),
    as.numeric(Kendall(subset$t6, subset$t7)[1]),
    as.numeric(Kendall(subset$t7, subset$t8)[1]),
    as.numeric(Kendall(subset$t8, subset$t9)[1]),
    as.numeric(Kendall(subset$t9, subset$t10)[1]),
    as.numeric(Kendall(subset$t10, subset$t11)[1]),
    as.numeric(Kendall(subset$t11, subset$t12)[1]))
  
  tao$RankCor<-tao_calc
  tao$treatment<-unique$treatment[i]
  
  kendall_rank<-rbind(kendall_rank, tao)
}

S<-ggplot(data=kendall_rank, aes(x=timestep, y=RankCor, color=treatment))+
  geom_point(size=5)+
  geom_line()

###all species
A_summedtrt<-pplots%>%
  tbl_df()%>%
  group_by(treatment_year, treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  spread(genus_species, relcov, fill=0)%>%
  gather(key="genus_species", value="relcov",3:73)%>%
  tbl_df()%>%
  group_by(treatment_year, treatment)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         timestep=paste("t", treatment_year, sep=""))%>%
  tbl_df()%>%
  select(-relcov, -treatment_year)%>%
  spread(timestep, rank)

A_unique<-A_summedtrt%>%
  select(treatment)%>%
  unique()


A_kendall_rank<-data.frame(row.names=1)

for(i in 1:length(A_unique$treatment)) {
  
  subset<-A_summedtrt%>%
    filter(treatment==A_unique$treatment[i])
  
  tao<-as.data.frame(seq(from=2, to=12))
  names(tao)[1]<-paste("timestep")
  
  tao_calc<-c(
    as.numeric(Kendall(subset$t1, subset$t2)[1]),
    as.numeric(Kendall(subset$t2, subset$t3)[1]),
    as.numeric(Kendall(subset$t3, subset$t4)[1]),
    as.numeric(Kendall(subset$t4, subset$t5)[1]),
    as.numeric(Kendall(subset$t5, subset$t6)[1]),
    as.numeric(Kendall(subset$t6, subset$t7)[1]),
    as.numeric(Kendall(subset$t7, subset$t8)[1]),
    as.numeric(Kendall(subset$t8, subset$t9)[1]),
    as.numeric(Kendall(subset$t9, subset$t10)[1]),
    as.numeric(Kendall(subset$t10, subset$t11)[1]),
    as.numeric(Kendall(subset$t11, subset$t12)[1]))
  
  tao$RankCor<-tao_calc
  tao$treatment<-A_unique$treatment[i]
  
  A_kendall_rank<-rbind(A_kendall_rank, tao)
}

A<-ggplot(data=A_kendall_rank, aes(x=timestep, y=RankCor, color=treatment))+
  geom_point(size=5)+
  geom_line()

grid.arrange(A, S, ncol=2)
### do correlations for the experiments we are subsetting

setwd("~/Dropbox/converge_diverge")

data<-readRDS("relAbundAgg.RDS")%>%
  filter(site_code!="JSP")%>%
  select(-rrich, -control_anpp, -experiment_length, -public, -n, -site_proj_comm_trt)

siteprojcom<-data%>%
  select(site_proj_comm)%>%
  unique

ranks<-data.frame(row.names=1)

for(i in 1:length(siteprojcom$site_proj_comm)) {

subset<-data%>%
  filter(site_proj_comm==siteprojcom$site_proj_comm[i])%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, treatment_year, treatment, genus_species)%>%
  spread(genus_species, relcov_agg, fill=0)

ranked<-subset%>%
  gather(key="genus_species", value="relcov",9:ncol(subset))%>%
  group_by(site_code, project_name, community_type, treatment_year, treatment)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         ntrt=paste("n", n_treatment, sep=""))%>%
  tbl_df()%>%
  select(-relcov, -n_treatment, -treatment)%>%
  spread(ntrt, rank)%>%
  mutate(siteprojcomyear=paste(site_proj_comm, treatment_year, sep="_"))

ranks<-rbind(ranked, ranks)
}


siteprojcom_year<-ranks%>%
  mutate(siteprojcomyear=paste(site_proj_comm, treatment_year, sep="_"))%>%
  select(siteprojcomyear)%>%
  unique()

kendall_rank<-data.frame(row.names=1)

for(i in 1:length(siteprojcom_year$siteprojcomyear)) {
  
  subset<-ranks%>%
    filter(siteprojcomyear==siteprojcom_year$siteprojcomyear[i])
  
  tao_calc<-as.numeric(Kendall(subset$n0, subset$n1)[1])
  temp.df <- data.frame(siteprojcom=subset$site_proj_comm[1],
                        treatment_year=subset$treatment_year[1],
                        RankCor=tao_calc)
kendall_rank<-rbind(kendall_rank, temp.df)
}
head(kendall_rank)
ggplot(kendall_rank, aes(x=treatment_year, y=tao))+geom_point()+
  facet_grid(~siteprojcom)
head(temp.df)
head(subset)

write.csv(kendall_rank, "kendall_rank.csv")
