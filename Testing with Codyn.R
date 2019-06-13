library(tidyverse)
library(vegan)
library(devtools)

<<<<<<< HEAD
install_github("NCEAS/codyn", ref = github_pull(83))
=======
>>>>>>> d4912af31ca1b2f8bee07dbf0eeccc9a8c3ee1db
install_github("NCEAS/codyn", ref = "anderson")
library(codyn)

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
  
  out<-community_diversity(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id')
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


## test odd year values
test<-subset(codyndat_clean, site_project_comm =="OND.ZOOPS.0")
testing<-RAC_change(df = test, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "id", reference.time = 1976.429)

###Looking at abundance Changes
##Codyn Dataset

codyndat_abund_change<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- abundance_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
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

write.csv(codyn_multchange, "~/Dropbox/SESYNC/SESYNC_RACs/R Files/anderson_codyn_mult_new.csv", row.names = F )

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

#####putting it all together
codyndat_rac_change2<-codyndat_rac_change%>%
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
  summarise_at(vars(richness_change:losses), mean, na.rm = T)%>%
  ungroup()%>%
  mutate(experiment_year=as.character(experiment_year), experiment_year2=as.character(experiment_year2))

codyn_multchange2<-codyn_multchange%>%
  separate(experiment_year_pair, into = c("experiment_year","experiment_year2"), sep= "-")%>%
  left_join(codyn_dissimchange)%>%
  left_join(codyndat_rac_change2)%>%
  na.omit

pairs(codyn_multchange2[,c(3,4,6:12)])

summary(lm(BC_between_change~richness_change+evenness_change+rank_change+gains+losses, data=codyn_multchange2))
summary(lm(centroid_distance_change~richness_change+evenness_change+rank_change+gains+losses, data=codyn_multchange2))


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 1),
                   symbols = c("*", " "))
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=5, col="red")
}        
pairs(codyn_multchange2[,c(3,4,6:12)], labels = c("cent_dist", "dispersion","bc_btwn", "bc_within","rich","even","rank","gain","loss"), font.labels=0.5, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)