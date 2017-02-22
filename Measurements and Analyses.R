library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)

sim<-read.csv("~/Documents/SESYNC/SESYNC_RACs/R Files/SimCom.csv")%>%
  separate(timestep, c("v","time"), sep=1)%>%
  select(-X, -kappa, -m, -v)%>%
  mutate(time=as.numeric(time),
         id=paste(ComType, rep, sep="::"))
  
codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")

###CLEANING CODYN DATASET
#restrict to species that are present in an experiment
splist<-codyndat%>%
  group_by(site_code, project_name, community_type, species)%>%
  summarize(present=sum(abundance))%>%
  filter(present!=0)%>%
  select(-present)

#merge back and will drop species that do not exist in a dataset
codyndat_clean<-merge(codyndat, splist, by=c("site_code","project_name","community_type","species"))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))

#####CALCULATING DIVERSITY METRICS WITHIN A TIME STEP FOR EACH REPLICATE AND THEN AVERAGING LATER

#function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

#function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
E_q<-function(x){
  x1<-x[x!=0]
  if (length(x1)==1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r<-rank(x1, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x1)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]
  2/pi*atan(b)
}


#function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
E_simp<-function(x, S=length(x[x!=0]), N=sum(x[x!=0]), ps=x[x!=0]/N, p2=ps*ps ){
D<-sum(p2)
(1/D)/S
}

#calculating gini coefficeint using the gini function in the reldist package
#' @x the vector of abundances of each species
#' this tive the inverse of other measures of evenness??
Gini<-function(x){
  x1<-x[x!=0]
  1-reldist::gini(x1)
}

##need to get this working with NAs for mean calculations
codyndat_diversity <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q),
            Gini=mean(Gini),
            E_simp=mean(E_simp))

sim_diversity<-group_by(sim, ComType, time, rep)%>%
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  group_by(ComType, time)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q),
            Gini=mean(Gini),
            E_simp=mean(E_simp))

###graph this
pairs(sim_diversity[3:6])
pairs(codyndat_diversity[3:6])

#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS
#gains and losses
codyndat_loss<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="disappearance")
codyndat_gain<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="appearance")
codyndat_gains_loss<-merge(codyndat_gain, codyndat_loss, by=c("experiment_year","id"))%>%
  separate(id, c("site_project_comm", "plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

sim_loss<-turnover(df=sim, time.var="time", species.var="species", abundance.var="abundance", replicate.var="id", metric="disappearance")
sim_gain<-turnover(df=sim, time.var="time", species.var="species", abundance.var="abundance", replicate.var="id", metric="appearance")
sim_gains_loss<-merge(sim_gain, sim_loss, by=c("time","id"))%>%
  separate(id, c("ComType", "plot_id"), sep="::")%>%
  group_by(ComType, time)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

####New appraoch to Rank Shifts
###ranks - taking into account that all speices are not always present.
###Give all species with zerio abundace the S+1 rank for that year.
###includes species that are not present in year X but appear in year X+1 or are present in year X and disappear in year X+1

##Codyn dataset first

##add ranks dropping zeros
codyndat_rank_pres<-codyndat_clean%>%
  filter(abundance!=0)%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year, plot_id)%>%
  mutate(rank=rank(-abundance, ties.method = "average"))%>%
  tbl_df()

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
codyndat_zeros<-codyndat_clean%>%
  filter(abundance==0)
##get species richness for each year
codyndat_S<-group_by(codyndat_clean, site_project_comm, experiment_year, plot_id)%>%
  summarize(S=S(abundance))
##merge together make zero abundances rank S+1
codyndat_zero_rank<-merge(codyndat_zeros, codyndat_S, by=c("site_project_comm","experiment_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
codyndat_rank<-rbind(codyndat_rank_pres, codyndat_zero_rank)

##calculate reordering between time steps 3 ways, rank correlations, mean rank shifts not corrected, and mean ranks shifts corrected for the size of the speceis pool

reordering=data.frame(id=c(), experiment_year=c(), MRSc=c())#expeiment year is year of timestep2

spc_id<-unique(codyndat_rank$id)
  
for (i in 1:length(spc_id)){
  subset<-codyndat_rank%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  
  splist<-subset%>%
    select(species)%>%
    unique()
  sppool<-length(splist$species)
  
#now get all timestep within an experiment
  timestep<-sort(unique(subset$experiment_year))    

  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(experiment_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(experiment_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

codyndat_reorder<-reordering%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarise(MRSc=mean(MRSc))

##SIM dataset

##add ranks dropping zeros
sim_rank_pres<-sim%>%
  filter(abundance!=0)%>%
  tbl_df()%>%
  group_by(ComType, time, rep)%>%
  mutate(rank=rank(-abundance, ties.method = "average"))%>%
  tbl_df()

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
sim_zeros<-sim%>%
  filter(abundance==0)
##get species richness for each year
sim_S<-group_by(sim, ComType, time, rep)%>%
  summarize(S=S(abundance))
##merge together make zero abundances rank S+1
sim_zero_rank<-merge(sim_zeros, sim_S, by=c("ComType","time","rep"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
sim_rank<-rbind(sim_rank_pres, sim_zero_rank)

##calculating re-ordering

reordering=data.frame(id=c(), time=c(), MRSc=c())#expeiment year is year of timestep2

spc_id<-unique(sim_rank$id)

for (i in 1:length(spc_id)){
  subset<-sim_rank%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(time==timestep[i])
    
    subset_t2<-subset%>%
      filter(time==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
   
    metrics<-data.frame(id=id, time=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

sim_reorder<-reordering%>%
  separate(id, c("ComType","rep"), sep="::")%>%
  group_by(ComType, time)%>%
  summarise(MRSc=mean(MRSc))

#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps.
##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.
###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

#Codyn dataset first

#make a new dataframe with just the label;
site_project_comm_u<-unique(codyndat_clean$site_project_comm)

#makes an empty dataframe
bray_curtis=data.frame(site_project_comm=c(), experiment_year=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 

##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(site_project_comm_u)) {
  
  #subsets out each dataset
  subset=codyndat_clean%>%
    filter(site_project_comm==site_project_comm_u[i])%>%
    select(site_project_comm, experiment_year, species, abundance, plot_id)
  
  #get years
  experiment_years<-sort(unique(subset$experiment_year))
  
  #transpose data
  species=subset%>%
    spread(species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,4:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$experiment_year, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))

  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(site_project_comm=site_project_comm_u[i],
                           experiment_year=experiment_years[2:length(experiment_years)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(site_project_comm=site_project_comm_u[i],
                  experiment_year=species$experiment_year,
                  plot_id=species$plot_id,
                  dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm, experiment_year)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(site_project_comm=site_project_comm_u[i],
                           experiment_year=experiment_years[2:length(experiment_years)],
                           dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("site_project_comm","experiment_year"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}

codyndat_braycurtis<-bray_curtis

#Sim dataset

#make a new dataframe with just the label;
ComType_u<-unique(sim$ComType)

#makes an empty dataframe
bray_curtis=data.frame(ComType=c(), time=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 

#Calculating bc mean change and dispersion
for(i in 1:length(ComType_u)) {
  
  #subsets out each dataset
  subset=sim%>%
    filter(ComType==ComType_u[i])%>%
    select(ComType, time, species, abundance, rep)
  
  #get years
  timestep<-sort(unique(subset$time))
  
  #transpose data
  species=subset%>%
    spread(species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,4:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$time, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(ComType=ComType_u[i],
                           time=timestep[2:length(timestep)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(ComType=ComType_u[i],
                   time=species$time,
                   rep=species$rep,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(ComType, time)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(ComType=ComType_u[i],
                      time=timestep[2:length(timestep)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("ComType","time"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}
sim_bray_curtis<-bray_curtis

####Looking at the shape of the curve - cc

d_output=data.frame(id=c(), experiment_year=c(), Darea=c())#expeiment year is year of timestep2

spc_id<-unique(codyndat_clean$id)

for (i in 1:length(spc_id)){
  subset<-codyndat_clean%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  
  ranks<-subset%>%
  filter(abundance!=0)%>%
  group_by(calendar_year, treatment, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov))

result <- average_test %>%
  group_by(plot_id) %>%
  do({
    y <- unique(.$calendar_year)###assumption this is a length 2 list
    df1 <- filter(., calendar_year==y[[1]])
    df2 <- filter(., calendar_year==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- c(0, df1$relrank, df2$relrank)
    data.frame(Dmax=max(abs(sf1(r) - sf2(r))))#do has to output a dataframe
  })

d_output=data.frame(id=c(), experiment_year=c(), Dmax=c(), Darea=c())#expeiment year is year of timestep2

spc_id<-unique(codyndat_rank$id)

for (i in 1:length(spc_id)){
  subset<-codyndat_rank%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  
  splist<-subset%>%
    select(species)%>%
    unique()
  sppool<-length(splist$species)
  
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$experiment_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(experiment_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(experiment_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    MRSt<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/sppool
    MRSnc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))
    n<-nrow(subset_t12)
    MRSm<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/(n+n/2-n/4+1.5)
    
    #RC<-abs(cor(subset_t12$rank.x, subset_t12$rank.y, method="kendall"))
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRSc=MRSc, MRSnc=MRSnc, MRSt=MRSt, MRSm=MRSm)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}


####MERGING TO A SINGE DATASET
#codyn
merge1<-merge(codyndat_diversity, codyndat_gains_loss, by=c("site_project_comm","experiment_year"))
merge2<-merge(merge1, codyndat_reorder, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge2, codyndat_braycurtis, by=c("site_project_comm","experiment_year"))


#sim
merge1<-merge(sim_diversity, sim_gains_loss, by=c("ComType","time"))
merge2<-merge(merge1, sim_reorder, by=c("ComType","time"))
sim_allmetrics<-merge(merge2, sim_bray_curtis, by=c("ComType","time"))

#graphing this
pairs(codyndat_allmetrics[,c(3:4,7:11)])
pairs(sim_allmetrics[,c(3:4,7:11)])

##correlations CODYN
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$E_Q)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$MRSc)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$MRSc)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$loss)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$loss)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$loss)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$mean_change)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$mean_change)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$mean_change)
cor.test(codyndat_allmetrics$gain, codyndat_allmetrics$mean_change)
cor.test(codyndat_allmetrics$loss, codyndat_allmetrics$mean_change)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$dispersion_diff)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$dispersion_diff)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$dispersion_diff)
cor.test(codyndat_allmetrics$gain, codyndat_allmetrics$dispersion_diff)
cor.test(codyndat_allmetrics$loss, codyndat_allmetrics$dispersion_diff)

##correlations SIM
cor.test(sim_allmetrics$S, sim_allmetrics$E_Q)
cor.test(sim_allmetrics$S, sim_allmetrics$MRSc)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$MRSc)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$gain)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$loss)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$loss)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$loss)
cor.test(sim_allmetrics$S, sim_allmetrics$mean_change)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$mean_change)
cor.test(sim_allmetrics$MRSc, sim_allmetrics$mean_change)
cor.test(sim_allmetrics$gain, sim_allmetrics$mean_change)
cor.test(sim_allmetrics$loss, sim_allmetrics$mean_change)
cor.test(sim_allmetrics$S, sim_allmetrics$dispersion_diff)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$dispersion_diff)
cor.test(sim_allmetrics$MRSc, sim_allmetrics$dispersion_diff)
cor.test(sim_allmetrics$gain, sim_allmetrics$dispersion_diff)
cor.test(sim_allmetrics$loss, sim_allmetrics$dispersion_diff)

###Supplental figure of averaging.
ave_codyndat_allmetrics<-codyndat_allmetrics%>%
  group_by(site_project_comm)%>%
  summarize(S=mean(S),
            E_Q=mean(E_Q),
            gain=mean(gain),
            loss=mean(loss),
            MRSc=mean(MRSc),
            mean_change=mean(mean_change),
            dispersion_diff=mean(dispersion_diff))
pairs(ave_codyndat_allmetrics[,c(2:8)])



