library(tidyverse)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)


#read in the data FIX THE PATH LATER

corredat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0")

corredat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat<-rbind(corredat1, azi, jrn, knz, sak)

#problems
#gvn face - only 2 years of data so will only have one point for the dataset.


plotinfo<-corredat%>%
  select(site_project_comm, calendar_year, plot_id, treatment, treatment_year)%>%
  unique()


#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS FOR EACH REPLICATE

##need to get this working with NAs for mean calculations
spc<-unique(corredat$site_project_comm)
codyn_div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  codyn_div_eq<-rbind(codyn_div_eq, out)
}


diversity <- group_by(corredat, site_project_comm, calendar_year, plot_id) %>% 
  summarize(S=S(relcov),
            E_q=E_q(relcov))%>%
  ungroup()%>%
  group_by(site_project_comm, plot_id)%>%
  arrange(site_project_comm, plot_id, calendar_year)%>%
  mutate(S_diff=c(999, diff(S)),
         E_diff=c(999, diff(E_q)))%>%
  filter(S_diff!=999)

#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS FOR EACH REPLICATE
explist<-unique(corredat$site_project_comm)

gain_loss<-data.frame()

for (i in 1:length(explist)){
  subset<-corredat%>%
    filter(site_project_comm==explist[i])
  
  loss<-turnover(df=subset, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="disappearance")
  gain<-turnover(df=subset, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="appearance")
  
  gain$site_project_comm<-explist[i]
  loss$site_project_comm<-explist[i]
  
  gl<-merge(gain, loss, by=c("site_project_comm","plot_id","calendar_year"))
  
  gain_loss<-rbind(gain_loss, gl)
}


####New appraoch to Rank Shifts
###ranks - taking into account that all speices are not always present.
###Give all species with zerio abundace the S+1 rank for that year.
###includes species that are not present in year X but appear in year X+1 or are present in year X and disappear in year X+1

##add ranks 
ranks<-corredat%>%
  filter(relcov!=0)%>%
  tbl_df()%>%
  group_by(site_project_comm, calendar_year, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method = "average"))%>%
  tbl_df()

#adding zeros
addzero<-data.frame()

for (i in 1:length(explist)){
  subset<-corredat%>%
    filter(site_project_comm==explist[i])
  
  wide<-subset%>%
    spread(key=genus_species, value=relcov, fill=0)
  
  long<-wide%>%
    gather(key=genus_species, value=relcov, 10:ncol(wide))
  
  addzero<-rbind(addzero, long)
}

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
corre_zeros<-addzero%>%
  filter(relcov==0)
##get species richness for each year
corre_S<-group_by(corredat, site_project_comm, calendar_year, plot_id)%>%
  summarize(S=S(relcov))
##merge together make zero abundances rank S+1
corre_zero_rank<-merge(corre_zeros, corre_S, by=c("site_project_comm","calendar_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
corre_rank<-rbind(ranks, corre_zero_rank)%>%
  mutate(exp_plot=paste(site_project_comm, plot_id, sep="::"))

##calculate reordering between time steps by mean ranks shifts corrected for the size of the speceis pool

reorder=data.frame(id=c(), calendar_year=c(), MRSc=c())#expeiment year is year of timestep2

exp_plot_list<-unique(corre_rank$exp_plot)

for (i in 1:length(exp_plot_list)){
  subset<-corre_rank%>%
    filter(exp_plot==exp_plot_list[i])
  id<-exp_plot_list[i]
  
  splist<-subset%>%
    select(genus_species)%>%
    unique()
  sppool<-length(splist$genus_species)
  
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$calendar_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(calendar_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(calendar_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("genus_species","site_project_comm","site_code","project_name","community_type"), all=T)%>%
      filter(relcov.x!=0|relcov.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id, calendar_year=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reorder=rbind(metrics, reorder)  
  }
}

reordering<-reorder%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")

#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps for all plots within a treatment.

##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.

###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

###list of all treats within an experiment

corredat$exptreat<-paste(corredat$site_project_comm, corredat$treatment, sep="::")

exptreatlist<-unique(corredat$exptreat)

#makes an empty dataframe
bray_curtis=data.frame(site_project_comm_treat=c(), calendar_year=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 
##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(exptreatlist)) {
  
  #subsets out each dataset
  subset=corredat%>%
    filter(exptreat==exptreatlist[i])%>%
    select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id)
  
  #get years
  experiment_years<-sort(unique(subset$calendar_year))
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$calendar_year, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(site_project_comm_treat=exptreatlist[i],
                           calendar_year=experiment_years[2:length(experiment_years)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(site_project_comm_treat=exptreatlist[i],
                   calendar_year=species$calendar_year,
                   plot_id=species$plot_id,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm_treat, calendar_year)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(site_project_comm_treat=exptreatlist[i],
                      calendar_year=experiment_years[2:length(experiment_years)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("site_project_comm_treat","calendar_year"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}

corre_braycurtis<-bray_curtis%>%
  separate(site_project_comm_treat, into=c("site_project_comm","treatment"), sep="::")

##getting the average for each treatment in a year

corre_diversity<-merge(plotinfo, diversity, by=c("site_project_comm","calendar_year","plot_id"))%>%
  group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
  summarize(S_diff=mean(S_diff), even_diff=mean(E_diff, na.rm=T))

corre_gainloss<-merge(plotinfo, gain_loss, by=c("site_project_comm","calendar_year","plot_id"))%>%
  group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
  summarize(gain=mean(appearance), loss=mean(disappearance))

corre_reordering<-merge(plotinfo, reordering, by=c("site_project_comm","calendar_year","plot_id"))%>%
  group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
  summarize(MRSc=mean(MRSc))

####MERGING TO A SINGE DATASET and exporting

merge1<-merge(corre_diversity, corre_gainloss, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
merge2<-merge(merge1, corre_reordering, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
all_metrics<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","treatment"), all=T)

write.csv(all_metrics, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allyears_2.csv")


write.csv(all_metrics, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allyears_2.csv")

#no longer necessary because everything is compared across years
# merge1<-merge(corre_diversity, corre_gainloss, by=c("site_project_comm","calendar_year","treatment_year","treatment"))
# merge2<-merge(merge1, corre_reordering, by=c("site_project_comm","calendar_year","treatment_year","treatment"))
# all_metrics2<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","treatment"))
# 
# write.csv(all_metrics2, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_compareyears.csv")
# 
# write.csv(all_metrics2, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_compareyears.csv")


###Getting b-C distnace of each plot to itself comparing t1 to t2.

corredat$expplot<-paste(corredat$site_project_comm, corredat$plot_id, sep="::")

exp_plot_list<-unique(corredat$expplot)


#makes an empty dataframe
bray_curtis_dissim=data.frame(site_project_comm_plot=c(), calendar_year=c(), bc_dissim=c()) 

##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(exp_plot_list)) {
  
  #subsets out each dataset
  subset=corredat%>%
    filter(expplot==exp_plot_list[i])%>%
    select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id)
  
  #get years
  experiment_years<-sort(unique(subset$calendar_year))
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=as.matrix(vegdist(species[,5:ncol(species)], method="bray"))
  
  ###experiment_year is year x+1
  bc_dis=data.frame(site_project_comm_plot=exp_plot_list[i],
                    calendar_year=experiment_years[2:length(experiment_years)],
                    bc_dissim=diag(bc[2:nrow(bc),1:(ncol(bc)-1)]))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis_dissim=rbind(bc_dis, bray_curtis_dissim)  
}

corre_braycurtis<-bray_curtis_dissim%>%
  separate(site_project_comm_plot, into=c("site_project_comm","plot_id"), sep="::")

###merging to a single dataset and adding treatment information
merge1<-merge(gain_loss, diversity, by=c("site_project_comm","calendar_year","plot_id"))
merge2<-merge(merge1, reordering,by=c("site_project_comm","calendar_year","plot_id")) 
merge3<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","plot_id"))
corre_all<-merge(plotinfo, merge3, by=c("site_project_comm","calendar_year","plot_id"))

write.csv(corre_all, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allReplicates.csv")

write.csv(corre_all, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allReplicates_2.csv")


#######Doing difference
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)

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

#read in the data FIX THE PATH LATER

corredat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotinfo<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, plot_mani)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))



corredat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_May2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotinfo<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, plot_mani)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

#problems
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE)
#zero-length inputs cannot be mixed with those of non-zero length 

##fill in the zeros
explist<-unique(corredat$site_project_comm)


corredat_sppool<-data.frame()

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corredat%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  ##make long and get averages of each species by treatment
  long<-subset%>%
    gather(genus_species, relcov, 10:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, treatment_year, genus_species)%>%
    summarize(relcov=mean(relcov))%>%
    ungroup
  
  corredat_sppool<-rbind(corredat_sppool, long)
}


###richness and evenness
diversity <- group_by(corredat_sppool, site_project_comm, calendar_year, treatment_year, treatment) %>% 
  summarize(S=S(relcov),
            Even=E_q(relcov))%>%
  tbl_df()

#subtract treatment from controls
control<-merge(diversity, plotinfo, by=c("site_project_comm", "calendar_year","treatment"))%>%
  filter(plot_mani==0)%>%
  mutate(controlS=S,
         controlEven=Even)%>%
  select(-S, -Even)

div_diff<-merge(control, diversity, by=c("site_project_comm","calendar_year","treatment_year"))%>%
  mutate(PCSdiff=(S-controlS)/controlS,
         PCEvendiff=(Even-controlEven)/controlEven,
         treatment=treatment.y)%>%
  filter(treatment.x!=treatment.y)%>%
  select(site_project_comm, calendar_year, treatment_year, treatment, PCSdiff, PCEvendiff)

###calculate species differences and reordering

#label the control versus treatment plots

corredat_treat_control<-merge(plotinfo, corredat,by=c("site_code","project_name","community_type","calendar_year",'treatment',"site_project_comm"))%>%
  mutate(expyear=paste(site_project_comm, calendar_year, sep="::"))


reordering_ct=data.frame(site_project_comm=c(), treatment=c(), calendar_year=c(), MRSc_diff=c(), spdiffc=c())

explist<-unique(corredat_treat_control$site_project_comm)

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corredat_treat_control%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  spc<-explist[i]
  
  ##make long and get averages of each species by treatment
  long<-subset%>%
    gather(genus_species, relcov, 12:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, plot_mani, genus_species)%>%
    summarize(relcov=mean(relcov))
  
  ##add ranks dropping zeros
  rank_pres<-long%>%
    filter(relcov!=0)%>%
    tbl_df()%>%
    group_by(site_project_comm, calendar_year, treatment, plot_mani)%>%
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-long%>%
    filter(relcov==0)
  ##get species richness for each year
  rich<-group_by(long, site_project_comm, calendar_year, treatment, plot_mani)%>%
    summarize(S=S(relcov))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year", "treatment", "plot_mani"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  timestep<-sort(unique(rank$calendar_year)) 
  
  for(i in 1:(length(timestep))){
    
    time<-rank%>%
      filter(calendar_year==timestep[i])
    
    time_id<-timestep[i]
    
    #fitler out control plots
    control<-time%>%
      filter(plot_mani==0)
    
    treat_list<-unique(subset(time, plot_mani!=0)$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("site_project_comm", "calendar_year","genus_species"), all=T)%>%
        filter(relcov.x!=0|relcov.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      
      spdiff<-subset_ct%>%
        filter(relcov.x==0|relcov.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      metrics<-data.frame(site_project_comm=spc, treatment=treat_id, calendar_year=time_id, MRSc_diff=MRSc_diff, spdiffc=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering_ct=rbind(metrics, reordering_ct)  
    }
  }
}

###mean change and dispersion

#####Calculating Bray-Curtis both comparing the mean community change between treatment and control plots in a time step

###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

###list of all treats within an experiment

exp_year<-unique(corredat_treat_control$expyear)

#makes an empty dataframe
bray_curtis=data.frame() 
##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(exp_year)) {
  
  #subsets out each dataset
  subset<-corredat_treat_control%>%
    filter(expyear==exp_year[i])%>%
    select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id, plot_mani)
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,6:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(site_project_comm_year=exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #renaming column
  colnames(cent_C_T)[3]<-"mean_change"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2=data.frame(site_project_comm_year=exp_year[i],
                   treatment=species$treatment,
                   plot_mani=species$plot_mani,
                   plot_id=species$plot_id,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm_year, treatment, plot_mani)%>%
    summarize(dispersion=mean(dist))
  
  control<-disp2$dispersion[disp2$plot_mani==0]
  
  ##subtract control from treatments
  disp_treat=disp2%>%
    mutate(disp_diff=dispersion-control)%>%
    select(-dispersion)
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_C_T, disp_treat, by=c("site_project_comm_year","treatment"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(bray_curtis, distances)  
}

corre_braycurtis_control_treat<-bray_curtis%>%
  separate(site_project_comm_year, into=c("site_project_comm","calendar_year"), sep="::")%>%
  filter(plot_mani!=0)

merge1<-merge(div_diff, reordering_ct, by=c("site_project_comm","calendar_year","treatment"))
all_Cont_Treat_Compare<-merge(merge1, corre_braycurtis_control_treat,by=c("site_project_comm","calendar_year","treatment"))

write.csv(all_Cont_Treat_Compare, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_ContTreat_Compare_OCT2017.csv")

