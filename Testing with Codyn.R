library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(gtools)


# Read in Data ------------------------------------------------------------
#home

codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="."))

#work

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
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep = "."))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))



# Richness Evenness Metrics -----------------------------------------------

#codyn dataset
spc<-unique(codyndat_clean$site_project_comm)
codyn_div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  codyn_div_eq<-rbind(codyn_div_eq, out)
}

###diversity

codyn_diversity<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  codyn_diversity<-rbind(codyn_diversity, out)
}
###Looking at RAC Changes
##Codyn Dataset

codyndat_rac_change<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- RAC_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyndat_rac_change<-rbind(codyndat_rac_change, out)  
}

pdata<-subset(codyndat_clean, site_project_comm == "KNZ.pplots.0")
pplots <- RAC_change(df = pdata, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plotid")

test <- RAC_change(df = codyndat_clean, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "id")

###Looking at abundance Changes
##Codyn Dataset

codyndat_abund_change<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- RAC_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyndat_abund_change<-rbind(codyndat_abund_change, out)  
}
# Mean Change and Dispersion ----------------------------------------------
#codyn dataset

#codyndat_mult_change <- multivariate_change(df = codyndat_clean, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "id")
#it doesnt work this way, give error: Error in rowSums(x, na.rm = TRUE) : 'x' must be numeric

codyn_multchange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- multivariate_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_multchange<-rbind(codyn_multchange, out)  
}


# Curve change ------------------------------------------------------------

####codyn first


codyn_curvechange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- curve_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_curvechange<-rbind(codyn_curvechange, out)  
}


