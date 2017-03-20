library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)

###get corre dataset to work on and focus on familiar examples
corre<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Dec2016.csv")%>%
  filter(project_name=="pplots"|project_name=="herbdiv"|project_name=="WENNDEx")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type,sep="_"))%>%
  select(-X)

#function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

###
### Functions that work on a plot level, no comparison
###

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

###Calculating the metrics. No need to do anything different here.
corre_diversity <- group_by(corre, site_project_comm, calendar_year, treatment, plot_id) %>% 
  summarize(S=S(relcov),
            E_q=E_q(relcov))%>%
  tbl_df()%>%
  group_by(site_project_comm, treatment, calendar_year)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q))

##comparing plots through time (RACS). This is through time, do not need a treatment, can do on any dataset but a treatment column is also possible.

##use codyn package to calculate appearances and dissapearances. To do this need to incorporate year and treatment into the replicate var. There should be a better way to do this.

corre_id<-corre%>%
  mutate(id=paste(site_project_comm, treatment, plot_id, sep="::"))

loss<-turnover(df=corre_id, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="id", metric="disappearance")

gain<-turnover(df=corre_id, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="id", metric="appearance")

gains_loss<-merge(loss, gain, by=c("calendar_year","id"))%>%
  separate(id, c("site_project_comm","treatment", "plot_id"), sep="::")%>%
  group_by(site_project_comm, calendar_year, treatment)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

##reordering

##for this to work need zero abundances to be filled in for all species.

corre_zero_filled<-data.frame(site_code=c(),project_name=c(), community_type=c(), calendar_year=c(), treatment=c(), plot_id=c(), genus_species=c(), relcov=c())

explist<-unique(corre$site_project_comm)

for (i in 1:length(explist)){
  subset<-corre%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)%>%
    gather_(genus_species, relcov,colnames(10:ncol(subset)))###why can't i get this to work??
}


##add ranks dropping zeros
corre_rank_pres<-corre%>%
  filter(relcov!=0)%>%
  tbl_df()%>%
  group_by(site_project_comm, calendar_year, treatment, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method = "average"))%>%
  tbl_df()

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
corre_zeros<-corre%>%
  filter(relcov==0)
##get species richness for each year
corre_S<-group_by(corre, site_project_comm, calendar_year, plot_id)%>%
  summarize(S=S(relcov))
##merge together make zero abundances rank S+1
corre_zero_rank<-merge(corre_zeros, corre_S, by=c("site_project_comm","calendar_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
corre_rank<-rbind(codyndat_rank_pres, codyndat_zero_rank)

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