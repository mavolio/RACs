library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

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
            E_Q=mean(E_q, na.rm=T),
            Gini=mean(Gini),
            E_simp=mean(E_simp))

sim_diversity<-group_by(sim, ComType, time, rep)%>%
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  group_by(ComType, time)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q, na.rm=T),
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
    arrange(abundance)%>%
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
d_output=data.frame(site_project_comm=c(), experiment_year=c(), plot_id=c(), Darea=c())#expeiment year is year of timestep2

com<-unique(sim$ComType)

for (i in 1:length(com)){
  subset<-sim%>%
    filter(ComType==com[i])
  
  ranks<-subset%>%
    filter(abundance!=0)%>%
    group_by(time, rep)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  com_id2<-com[i]
  
  timestep<-sort(unique(ranks$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(time==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(rep)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(time==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(rep)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="rep")
    #dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
    ##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="rep")
    
    #dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(time, rep)%>%
      mutate(numplots=length(rep))%>%
      ungroup()%>%
      group_by(rep)%>%
      mutate(min=min(numplots))%>%
      select(rep, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="rep")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(time, rep)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
      group_by(rep) %>%
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
    
    d_output1=data.frame(ComType=com_id2, time=timestep[i+1], rep=result$rep, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }
}

sim_dstar<-d_output%>% 
  group_by(ComType, time)%>%
  summarise(Dstar=mean(Dstar))

  
####MERGING TO A SINGE DATASET
#codyn
merge1<-merge(codyndat_diversity, codyndat_gains_loss, by=c("site_project_comm","experiment_year"))
merge2<-merge(merge1, codyndat_reorder, by=c("site_project_comm","experiment_year"))
merge3<-merge(merge2, codyndat_braycurtis, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge3, codyndat_dstar, by=c("site_project_comm","experiment_year"))

#sim
merge1<-merge(sim_diversity, sim_gains_loss, by=c("ComType","time"))
merge2<-merge(merge1, sim_reorder, by=c("ComType","time"))
merge3<-merge(merge2, sim_bray_curtis, by=c("ComType","time"))
sim_allmetrics<-merge(merge3, sim_dstar, by=c("ComType","time"))

#graphing this
pairs(codyndat_allmetrics[,c(3:4,7:12)])
pairs(sim_allmetrics[,c(3:4,7:12)])

##correlations CODYN
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$E_Q)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$Gini)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$E_simp)
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
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$E_Q, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$MRSc, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$gain, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$loss, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$mean_change, codyndat_allmetrics$Dstar)
cor.test(codyndat_allmetrics$dispersion_diff, codyndat_allmetrics$Dstar)

##correlations SIM
cor.test(sim_allmetrics$S, sim_allmetrics$E_Q)
cor.test(sim_allmetrics$S, sim_allmetrics$MRSc)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$MRSc)
cor.test(sim_allmetrics$S, sim_allmetrics$gain)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$gain)
cor.test(sim_allmetrics$MRSc, sim_allmetrics$gain)
cor.test(sim_allmetrics$S, sim_allmetrics$loss)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$loss)
cor.test(sim_allmetrics$MRSc, sim_allmetrics$loss)
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
cor.test(sim_allmetrics$S, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$E_Q, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$MRSc, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$gain, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$loss, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$mean_change, sim_allmetrics$Dstar)
cor.test(sim_allmetrics$dispersion_diff, sim_allmetrics$Dstar)

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


##graphing this
theme_set(theme_bw(10))

###Codyn graphs

simgraph<-merge(codyndat_allmetrics, codyndat_info, by="site_project_comm")

#richness y axis
se<-ggplot(data=codyndat_graphs, aes(x=E_Q, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,90))+
  xlab("Evenness")+
  ylab("Richness")
sg<-ggplot(data=codyndat_graphs, aes(x=gain, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,90))+
  ylab("")
sl<-ggplot(data=codyndat_graphs, aes(x=loss, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,90))+
  ylab("")
sr<-ggplot(data=codyndat_graphs, aes(x=MRSc, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,90))+
  ylab("")
sm<-ggplot(data=codyndat_graphs, aes(x=mean_change, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,90))+
  ylab("")
sd<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=S))+
  geom_point(size=0.5, aes(color=taxa))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,90))+
  ylab("")
sds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=S))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,90))+
  xlab("")+
  ylab("")

#evenness y-axis
eg<-ggplot(data=codyndat_graphs, aes(x=gain, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("Sp. Gains")+
  ylab("Evenness")
el<-ggplot(data=codyndat_graphs, aes(x=loss, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
er<-ggplot(data=codyndat_graphs, aes(x=MRSc, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
em<-ggplot(data=codyndat_graphs, aes(x=mean_change, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
ed<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
eds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=E_Q))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")

#gain y-axis
gl<-ggplot(data=codyndat_graphs, aes(x=loss, y=gain))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("Sp. Losses")+
  ylab("Sp. Gains")
gr<-ggplot(data=codyndat_graphs, aes(x=MRSc, y=gain))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
gm<-ggplot(data=codyndat_graphs, aes(x=mean_change, y=gain))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
gd<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=gain))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
gds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=gain))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")

##losses y-axis
lr<-ggplot(data=codyndat_graphs, aes(x=MRSc, y=loss))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("Reordering")+
  ylab("Sp. Losses")
lm<-ggplot(data=codyndat_graphs, aes(x=mean_change, y=loss))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
ld<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=loss))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")
lds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=loss))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")

#reorder y-axis
rm<-ggplot(data=codyndat_graphs, aes(x=mean_change, y=MRSc))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.4))+
  xlab("Changes Community")+
  ylab("Reordering")
rd<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=MRSc))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.4))+
  xlab("")+
  ylab("")
rds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=MRSc))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.4))+
  xlab("")+
  ylab("")

#mean change y-axis
md<-ggplot(data=codyndat_graphs, aes(x=dispersion_diff, y=mean_change))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("Dispersion")+
  ylab("Community Change")
mds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=mean_change))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,1))+
  xlab("")+
  ylab("")

#dispersion y axis
dds<-ggplot(data=codyndat_graphs, aes(x=Dstar, y=dispersion_diff))+
  geom_point(size=0.5, (aes(color=taxa)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=site_project_comm, color=taxa), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(-.3,.4))+
  xlab("D Star")+
  ylab("Dispersion")

legend=gtable_filter(ggplot_gtable(ggplot_build(dds)), "guide-box") 
grid.draw(legend)

#blank
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

grid.arrange(arrangeGrob(se+theme(legend.position="none"),
                         sg+theme(legend.position="none"),
                         sl+theme(legend.position="none"),
                         sr+theme(legend.position="none"),
                         sm+theme(legend.position="none"),
                         sd+theme(legend.position="none"),
                         sds+theme(legend.position="none"),
                         blankPlot,
                          eg+theme(legend.position="none"),
                         el+theme(legend.position="none"),
                         er+theme(legend.position="none"),
                         em+theme(legend.position="none"),
                         ed+theme(legend.position="none"),
                         eds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         gl+theme(legend.position="none"),
                         gr+theme(legend.position="none"),
                         gm+theme(legend.position="none"),
                         gd+theme(legend.position="none"),
                         gds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         lr+theme(legend.position="none"),
                         lm+theme(legend.position="none"),
                         ld+theme(legend.position="none"),
                         lds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         rm+theme(legend.position="none"),
                         rd+theme(legend.position="none"),
                         rds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         md+theme(legend.position="none"),
                         mds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         dds+theme(legend.position="none"),
                         ncol=7), legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)

###sim graphs
simgraph<-sim_allmetrics%>%
  separate(ComType, c("Simulated_Evenness","Simulated_Richness"), sep="_", remove = F)


#richness y axis
se<-ggplot(data=simgraph, aes(x=E_Q, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,50))+
  xlab("Evenness")+
  ylab("Richness")
sg<-ggplot(data=simgraph, aes(x=gain, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,50))+
  ylab("")
sl<-ggplot(data=simgraph, aes(x=loss, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,50))+
  ylab("")
sr<-ggplot(data=simgraph, aes(x=MRSc, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,50))+
  ylab("")
sm<-ggplot(data=simgraph, aes(x=mean_change, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,50))+
  ylab("")
sd<-ggplot(data=simgraph, aes(x=dispersion_diff, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  xlab("")+
  scale_y_continuous(limit=c(0,50))+
  ylab("")
sds<-ggplot(data=simgraph, aes(x=Dstar, y=S))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,50))+
  xlab("")+
  ylab("")

#evenness y-axis
eg<-ggplot(data=simgraph, aes(x=gain, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("Sp. Gains")+
  ylab("Evenness")
el<-ggplot(data=simgraph, aes(x=loss, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("")+
  ylab("")
er<-ggplot(data=simgraph, aes(x=MRSc, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("")+
  ylab("")
em<-ggplot(data=simgraph, aes(x=mean_change, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("")+
  ylab("")
ed<-ggplot(data=simgraph, aes(x=dispersion_diff, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("")+
  ylab("")
eds<-ggplot(data=simgraph, aes(x=Dstar, y=E_Q))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.6))+
  xlab("")+
  ylab("")

#gain y-axis
gl<-ggplot(data=simgraph, aes(x=loss, y=gain))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("Sp. Losses")+
  ylab("Sp. Gains")
gr<-ggplot(data=simgraph, aes(x=MRSc, y=gain))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")
gm<-ggplot(data=simgraph, aes(x=mean_change, y=gain))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")
gd<-ggplot(data=simgraph, aes(x=dispersion_diff, y=gain))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")
gds<-ggplot(data=simgraph, aes(x=Dstar, y=gain))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")

##losses y-axis
lr<-ggplot(data=simgraph, aes(x=MRSc, y=loss))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("Reordering")+
  ylab("Sp. Losses")
lm<-ggplot(data=simgraph, aes(x=mean_change, y=loss))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")
ld<-ggplot(data=simgraph, aes(x=dispersion_diff, y=loss))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")
lds<-ggplot(data=simgraph, aes(x=Dstar, y=loss))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")

#reorder y-axis
rm<-ggplot(data=simgraph, aes(x=mean_change, y=MRSc))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.24))+
  xlab("Changes Community")+
  ylab("Reordering")
rd<-ggplot(data=simgraph, aes(x=dispersion_diff, y=MRSc))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.24))+
  xlab("")+
  ylab("")
rds<-ggplot(data=simgraph, aes(x=Dstar, y=MRSc))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.24))+
  xlab("")+
  ylab("")

#mean change y-axis
md<-ggplot(data=simgraph, aes(x=dispersion_diff, y=mean_change))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("Dispersion")+
  ylab("Community Change")
mds<-ggplot(data=simgraph, aes(x=Dstar, y=mean_change))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(0,.1))+
  xlab("")+
  ylab("")

#dispersion y axis
dds<-ggplot(data=simgraph, aes(x=Dstar, y=dispersion_diff))+
  geom_point(size=2, aes(color=Simulated_Evenness, shape=Simulated_Richness ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(aes(group=ComType), method="lm", se=F, size=.25, alpha=0.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_y_continuous(limit=c(-.02,.02))+
  xlab("D Star")+
  ylab("Dispersion")

legend=gtable_filter(ggplot_gtable(ggplot_build(dds)), "guide-box") 
grid.draw(legend)

#blank
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

grid.arrange(arrangeGrob(se+theme(legend.position="none"),
                         sg+theme(legend.position="none"),
                         sl+theme(legend.position="none"),
                         sr+theme(legend.position="none"),
                         sm+theme(legend.position="none"),
                         sd+theme(legend.position="none"),
                         sds+theme(legend.position="none"),
                         blankPlot,
                         eg+theme(legend.position="none"),
                         el+theme(legend.position="none"),
                         er+theme(legend.position="none"),
                         em+theme(legend.position="none"),
                         ed+theme(legend.position="none"),
                         eds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         gl+theme(legend.position="none"),
                         gr+theme(legend.position="none"),
                         gm+theme(legend.position="none"),
                         gd+theme(legend.position="none"),
                         gds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         lr+theme(legend.position="none"),
                         lm+theme(legend.position="none"),
                         ld+theme(legend.position="none"),
                         lds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         rm+theme(legend.position="none"),
                         rd+theme(legend.position="none"),
                         rds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         md+theme(legend.position="none"),
                         mds+theme(legend.position="none"),
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         blankPlot,
                         dds+theme(legend.position="none"),
                         ncol=7), legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)

