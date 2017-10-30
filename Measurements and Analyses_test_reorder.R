library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

# Read in Data ------------------------------------------------------------
#home
sim<-read.csv("~/Dropbox/SESYNC/SESYNC_RACs/R Files/SimCom_Sept28.csv")%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)
  
codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")

#work
sim<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files/SimCom_Sept28.csv')%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))

codyndat<-read.csv('C:\\Users\\megha\\Dropbox\\CoDyn\\R Files\\11_06_2015_v7\\relative cover_nceas and converge_12012015_cleaned.csv')%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("C:\\Users\\megha\\Dropbox\\CoDyn\\R Files\\11_06_2015_v7\\siteinfo_key.csv")%>%
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



# Richness Evenness Metrics -----------------------------------------------


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
            E_Q=mean(E_q, na.rm=T),
            Gini=mean(Gini),
            E_simp=mean(E_simp))

sim_diversity<-group_by(sim, id, site, time)%>%
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  ungroup()%>%
  group_by(id, time)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q, na.rm=T),
            Gini=mean(Gini),
            E_simp=mean(E_simp))%>%
  ungroup()%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(S=mean(S),
            E_Q=mean(E_Q),
            Gini=mean(Gini),
            E_simp=mean(E_simp))
  
###graph this
pairs(sim_diversity[3:6])
pairs(codyndat_diversity[3:6])

#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS

codyndat_diversity_diff <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  ungroup()%>%
  group_by(site_project_comm, plot_id)%>%
  arrange(site_project_comm, plot_id, experiment_year)%>%
  mutate(S_diff=c(NA, diff(S)),
         E_diff=c(NA, diff(E_q)))%>%
  ungroup()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(S_diff=mean(S_diff, na.rm=T),
            E_diff=mean(E_diff, na.rm=T))%>%
  na.omit


sim_diversity_diff<-group_by(sim, id, site, time)%>%
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  ungroup()%>%
  group_by(id, site)%>%
  arrange(id, site, time)%>%
  mutate(S_diff=c(NA, diff(S)),
         E_diff=c(NA, diff(E_q)))%>%
  ungroup()%>%
  group_by(id, time)%>%
  summarize(S_diff=mean(S_diff, na.rm=T),
            E_diff=mean(E_diff, na.rm=T))%>%
  na.omit%>%
  ungroup()%>%
    separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(S_diff=mean(S_diff),
            E_diff=mean(E_diff))


# Gains and Losses --------------------------------------------------------


#gains and losses
codyndat_loss<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="disappearance")
codyndat_gain<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="appearance")
codyndat_gains_loss<-merge(codyndat_gain, codyndat_loss, by=c("experiment_year","id"))%>%
  separate(id, c("site_project_comm", "plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

sim_loss<-turnover(df=sim, time.var="time", species.var="species", abundance.var="abundance", replicate.var="id2", metric="disappearance")
sim_gain<-turnover(df=sim, time.var="time", species.var="species", abundance.var="abundance", replicate.var="id2", metric="appearance")
sim_gains_loss<-merge(sim_gain, sim_loss, by=c("time","id2"))%>%
  separate(id2, into=c("id", 'site'), sep="::")%>%
  group_by(id, time)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(gain=mean(gain),
            loss=mean(loss))

# Rank Shifts -------------------------------------------------------------


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

reordering=data.frame(id=c(), experiment_year=c(), MRSc=c(),SRS=c(), RSW=c(), MR=c(), SM=c())#expeiment year is year of timestep2

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
    SRS<-sum(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    RSW<-mean(abs(subset_t12$rank.x-subset_t12$rank.y)*abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    MR<-mean(abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    SM<-sum(abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRSc=MRSc, SRS=SRS, RSW=RSW, MR=MR, SM=SM)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

codyndat_reorder<-reordering%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarise(MRSc=mean(MRSc),
            SRS=mean(SRS),
            RSW=mean(RSW),
            MR=mean(MR),
            SM=mean(SM))

##SIM dataset
##add in zeros

##add ranks 
sim_rank_pres<-sim%>%
  filter(abundance!=0)%>%
  tbl_df()%>%
  group_by(id, time, site)%>%
  mutate(rank=rank(-abundance, ties.method = "average"))%>%
  tbl_df()%>%
  select(site, time, id, id2, species, abundance, rank)

#adding zeros
sim_addzero <- sim %>%
  group_by(id) %>%
  nest() %>%
  mutate(spread_df = purrr::map(data, ~spread(., key=species, value=abundance, fill=0) %>%
                                  gather(key=species, value=abundance, -site, -time, -id2))) %>%
  unnest(spread_df)

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
sim_zeros<-sim_addzero%>%
  filter(abundance==0)
##get species richness for each year
sim_S<-group_by(sim, id, time, site, id2)%>%
  summarize(S=S(abundance))
##merge together make zero abundances rank S+1
sim_zero_rank<-merge(sim_zeros, sim_S, by=c("id","time","site", "id2"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
sim_rank<-rbind(sim_rank_pres, sim_zero_rank)%>%
  mutate(abundance=as.numeric(abundance))

##calculating re-ordering

reordering=data.frame(id=c(), experiment_year=c(), MRSc=c(),SRS=c(), RSW=c(), MR=c(), SM=c())#expeiment year is year of timestep2

spc_id<-unique(sim_rank$id2)

for (i in 1:length(spc_id)){
  subset<-sim_rank%>%
    filter(id2==spc_id[i])
  id2<-spc_id[i]
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(time==timestep[i])
    
    subset_t2<-subset%>%
      filter(time==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id2"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12) 
    SRS<-sum(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    RSW<-mean(abs(subset_t12$rank.x-subset_t12$rank.y)*abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    MR<-mean(abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    SM<-sum(abs(subset_t12$abundance.x-subset_t12$abundance.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id2, time=timestep[i+1], MRSc=MRSc, SRS=SRS, RSW=RSW, MR=MR, SM=SM)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

sim_reorder<-reordering%>%
  separate(id, c("id","site"), sep="::")%>%
  group_by(id, time)%>%
  summarise(MRSc=mean(MRSc),
            SRS=mean(SRS),
            RSW=mean(RSW),
            MR=mean(MR),
            SM=mean(SM))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(MRSc=mean(MRSc),
            SRS=mean(SRS),
            RSW=mean(RSW),
            MR=mean(MR),
            SM=mean(SM))

# Mean Change and Dispersion ----------------------------------------------


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
id_u<-unique(sim$id)

#makes an empty dataframe
bray_curtis=data.frame(id=c(), time=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 

#Calculating bc mean change and dispersion
for(i in 1:length(id_u)) {
  
  #subsets out each dataset
  subset=sim%>%
    filter(id==id_u[i])%>%
    select(id, time, species, abundance, site)
  
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
  cent_dist_yrs=data.frame(id=id_u[i],
                           time=timestep[2:length(timestep)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(id=id_u[i],
                   time=species$time,
                   site=species$site,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(id, time)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(id=id_u[i],
                      time=timestep[2:length(timestep)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("id","time"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}
sim_bray_curtis<-bray_curtis%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(mean_change=mean(mean_change),
            dispersion_diff=mean(dispersion_diff))

# Curve change ------------------------------------------------------------


####Looking at the shape of the curve - cc
###this compares the areas of difference between two curves that are consequtive time steps for a plot.

#codyn dat first
d_output=data.frame(site_project_comm=c(), experiment_year=c(), plot_id=c(), Darea=c())#expeiment year is year of timestep2

spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
 
  ranks<-subset%>%
    filter(abundance!=0)%>%
    group_by(experiment_year, plot_id)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()

  spc_id2<-spc[i]
  
  timestep<-sort(unique(ranks$experiment_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(experiment_year==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(plot_id)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(experiment_year==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(plot_id)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="plot_id")
#dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="plot_id")
    
#dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(experiment_year, plot_id)%>%
      mutate(numplots=length(plot_id))%>%
      ungroup()%>%
      group_by(plot_id)%>%
      mutate(min=min(numplots))%>%
      select(plot_id, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="plot_id")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(experiment_year, plot_id)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
    group_by(plot_id) %>%
    do({
      y <- unique(.$experiment_year)###assumption this is a length 2 list
      df1 <- filter(., experiment_year==y[[1]])
      df2 <- filter(., experiment_year==y[[2]])
      sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
      sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
      r <- sort(unique(c(0, df1$relrank, df2$relrank)))
      h <- abs(sf1(r) - sf2(r))
      w <- c(diff(r), 0)
      data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

    d_output1=data.frame(site_project_comm=spc_id2, experiment_year=timestep[i+1], plot_id=result$plot_id, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }
}
codyndat_dstar<-d_output%>% 
  group_by(site_project_comm, experiment_year)%>%
  summarise(Dstar=mean(Dstar))

####Looking at the shape of the curve - cc
#sim dataset
d_output=data.frame(id=c(), time=c(), site=c(), Darea=c())#expeiment year is year of timestep2

com<-unique(sim$id)

for (i in 1:length(com)){
  subset<-sim%>%
    filter(id==com[i])
  
  ranks<-subset%>%
    filter(abundance!=0)%>%
    group_by(time, site)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  com_id2<-com[i]
  
  timestep<-sort(unique(ranks$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(time==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(site)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(time==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(site)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="site")
    #dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
    ##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="site")
    
    #dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(time, site)%>%
      mutate(numplots=length(site))%>%
      ungroup()%>%
      group_by(site)%>%
      mutate(min=min(numplots))%>%
      select(site, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="site")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(time, site)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
      group_by(site) %>%
      do({
        y <- unique(.$time)###assumption this is a length 2 list
        df1 <- filter(., time==y[[1]])
        df2 <- filter(., time==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(Dstar=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(id=com_id2, time=timestep[i+1], site=result$site, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }
}

sim_dstar<-d_output%>% 
  group_by(id, time)%>%
  summarise(Dstar=mean(Dstar))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(Dstar=mean(Dstar))
  

# looking at spatial differences, testing that scenarios work well --------


#######trying to look at spatial differences.
sim_subset<-sim%>%
  separate(id, into=c("alpha", "even", "comtype", "rep"), sep="_")%>%
  filter(time==1&rep==1)%>%
  mutate(alphaeven=paste(alpha, even, sep="_"))

sp_output=data.frame()

id<-unique(sim_subset$alphaeven)

for (i in 1:length(id)){
  species<-sim_subset%>%
    filter(alphaeven==id[i])%>%
  spread(species, abundance, fill=0)

mds<-metaMDS(species[,9:ncol(species)])

info<-species[,1:8]

scores <- data.frame(scores(mds, display="sites"))
scores2<- cbind(info, scores)

sp_output<-rbind(sp_output,scores2)
}

theme_set(theme_bw(12))
ggplot(data=sp_output, aes(x=NMDS1, y=NMDS2, color=comtype))+
  geom_point()+
  scale_color_manual(values=c("black","red","green","blue"))+
  facet_wrap(~alphaeven, ncol=3, scales="free")

# Merging all metrics to single datasets ----------------------------------


####MERGING TO A SINGE DATASET
#codyn
merge1<-merge(codyndat_diversity, codyndat_gains_loss, by=c("site_project_comm","experiment_year"))
merge2<-merge(merge1, codyndat_reorder, by=c("site_project_comm","experiment_year"))
merge3<-merge(merge2, codyndat_braycurtis, by=c("site_project_comm","experiment_year"))
merge4<-merge(merge3, codyndat_dstar, by=c("site_project_comm","experiment_year"))
merge5<-merge(merge4, codyndat_diversity_diff, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge5, codyndat_info, by="site_project_comm")

#sim
merge1<-merge(sim_diversity, sim_gains_loss, by=c("id3","time"))
merge2<-merge(merge1, sim_reorder, by=c("id3","time"))
merge3<-merge(merge2, sim_bray_curtis, by=c("id3","time"))
merge4<-merge(merge3, sim_diversity_diff, by=c("id3","time"))
sim_allmetrics<-merge(merge4, sim_dstar, by=c("id3","time"))%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")

sim_allmetrics$comtype2<-as.factor(sim_allmetrics$comtype)

write.csv(codyndat_allmetrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\codyn_allmetrics_diff_reorder.csv')
write.csv(sim_allmetrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_diff_reorder.csv')

# pair plot graphs --------------------------------------------------------


#graphing this
panel.pearson <- function(x, y, ...) {
  horizontal <- (par("usr")[1] + par("usr")[2]) / 2; 
  vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
  text(horizontal, vertical, format(cor(x,y), digits=3, cex=10)) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=2, col="red")
}



#dropoutliers in reordering
codyndat_allmetrics2<-codyndat_allmetrics%>%
  filter(SM<0.4)
pairs(codyndat_allmetrics2[,c(3,4, 17,18,7:16)], font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))

pairs(sim_allmetrics[,c(5,6, 18,19,9:17,20)], font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))

summary(lm(mean_change~gain+loss+E_diff+S_diff+MRSc, data=codyndat_allmetrics2))
        