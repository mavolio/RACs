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
sim<-read.csv("~/Dropbox/SESYNC/SESYNC_RACs/R Files/SimCom_Sept28.csv")%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="."))
  
codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="."))

#work
sim<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files/SimCom_Sept28.csv')%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="."))

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

codyn_div_esimp<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', evenness = "SimpEven")
  out$site_project_comm<-spc[i]
  
  codyn_div_esimp<-rbind(codyn_div_esimp, out)
}

codyn_div<-merge(codyn_div_eq, codyn_div_esimp, by=c("site_project_comm",'plot_id','richness','experiment_year'))


codyndat_diversity_mean <- codyn_div%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(Sp=mean(richness),
            EQ=mean(evenness_EQ, na.rm=T),
            ESimp=mean(evenness_Simpson))


#calculating gini coefficeint using the gini function in the reldist package
#' @x the vector of abundances of each species
Gini<-function(x){
  x1<-x[x!=0]
  1-reldist::gini(x1)
}

codyndat_div_gini <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(Gini=Gini(abundance))%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(EGini=mean(Gini))

codyn_div_all<-merge(codyndat_diversity_mean, codyndat_div_gini, by=c("experiment_year",'site_project_comm'))

#sim dataset
com_rep<-unique(sim$id)

sim_div_eq<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site")
  out$id<-com_rep[i]
  
  sim_div_eq<-rbind(sim_div_eq, out)  
}

sim_div_esimp<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", evenness = "SimpEven")
  out$id<-com_rep[i]
  
  sim_div_esimp<-rbind(sim_div_esimp, out)  
}

sim_div<-merge(sim_div_eq, sim_div_esimp, by=c("id",'site','richness','time'))

sim_diversity_mean<-sim_div%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(Sp=mean(richness),
            EQ=mean(evenness_EQ, na.rm=T),
            ESimp=mean(evenness_Simpson))

sim_div_gini <- group_by(sim, id, time, site) %>% 
  summarize(Gini=Gini(abundance))%>%
  tbl_df()%>%
  group_by(id, time)%>%
  summarize(EGini=mean(Gini))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(EGini=mean(EGini))
  

sim_div_all<-merge(sim_diversity_mean, sim_div_gini, by=c("time",'id3'))


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

codyndat_rac_change_average<-codyndat_rac_change%>%
  group_by(site_project_comm, experiment_year_pair)%>%
  summarise(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))

##SIM dataset
sim_rac_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- RAC_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_rac_change<-rbind(sim_rac_change, out)  
}

sim_rac_change_mean<-sim_rac_change%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time_pair)%>%
  summarize(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))

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

#Sim dataset
sim_mult_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- multivariate_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_mult_change<-rbind(sim_mult_change, out)  
}

sim_multchange_mean<-sim_mult_change%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time_pair)%>%
  summarize(composition_change=mean(composition_change),
            dispersion_change=mean(dispersion_change))

# Curve change ------------------------------------------------------------

####codyn first

#clean the data first, drop plots that only have 1 species
codyn_clean2<-subset(codyndat_clean, abundance!=0)
codyn_clean2$present<-1
codyn_clean_spnum<-aggregate(present~experiment_year+plot_id+site_project_comm, sum, data=codyn_clean2)
codyn_clean2<-subset(codyn_clean_spnum, present>1)

codyn_mult_sp<-merge(codyndat_clean, codyn_clean2, by=c("site_project_comm","experiment_year","plot_id"))

#codyn_cc <-curve_change(df = codyn_clean2, "experiment_year", "species", "abundance", "id")
#not working


codyn_curvechange<-data.frame()
spc<-unique(codyn_mult_sp$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyn_mult_sp%>%
    filter(site_project_comm==spc[i])
  
  out <- curve_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_curvechange<-rbind(codyn_curvechange, out)  
}

codyn_curvechange_mean<-codyn_curvechange%>% 
  group_by(site_project_comm, experiment_year_pair)%>%
  summarise(curve_change=mean(curve_change))

#sim dataset
sim_curve_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- curve_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_curve_change<-rbind(sim_curve_change, out)  
}

sim_cc_ave<-sim_curve_change%>% 
  group_by(id, time_pair)%>%
  summarise(curve_change=mean(curve_change))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time_pair)%>%
  summarize(curve_change=mean(curve_change))
  

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
merge1<-merge(codyn_multchange, codyn_curvechange_mean, by=c("site_project_comm","experiment_year_pair"))
merge2<-merge(merge1, codyndat_rac_change_average, by=c("site_project_comm","experiment_year_pair"))%>%
  separate(experiment_year_pair, into=c("time1","experiment_year"), sep="-", remove=F)
merge3<-merge(merge2, codyn_div_all, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge3, codyndat_info, by="site_project_comm")

#sim
merge1<-merge(sim_multchange_mean, sim_rac_change_mean, by=c("id3","time_pair"))
merge2<-merge(merge1, sim_cc_ave, by=c("id3","time_pair"))%>%
  separate(time_pair, into=c("time1","time"), sep="-", remove=F)
sim_all_metrics<-merge(sim_div_all, merge2, by=c('time','id3'))%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")
sim_all_metrics$comtype2<-as.factor(sim_all_metrics$comtype)

write.csv(codyndat_allmetrics,'~/Dropbox/SESYNC/SESYNC_RACs/R Files/codyn_allmetrics_Jan2018.csv')
write.csv(sim_all_metrics,'~/Dropbox/SESYNC/SESYNC_RACs/R Files/sim_allmetrics_Jan2018.csv')

# pair plot graphs --------------------------------------------------------

codyndat_allmetrics<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\codyn_allmetrics_diff_corrected.csv')%>%
  select(-X)
sim_allmetrics<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_diff_corrected.csv')%>%
  select(-X)



#graphing this
panel.pearson <- function(x, y, ...) {
  horizontal <- (par("usr")[1] + par("usr")[2]) / 2; 
  vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
  text(horizontal, vertical, format(cor(x,y), digits=3, cex=10)) 
}

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

sim_allmetrics2<-sim_allmetrics%>%
  mutate(comtype3=as.factor(paste(alpha, even, sep="_")))
# #color by turnover and sucession
# pairs(sim_allmetrics2[,c(9:16)], col=sim_allmetrics$comtype2, labels=c("Richness\nChange", "Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange","Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor, oma=c(4,4,4,10))
# 
# #color by richness_evenness
# pairs(sim_allmetrics[,c(14:15,9:13,16)], col=sim_allmetrics$comtype3, labels=c("Richness \nChange", "Evenness \nChange","Species \nGains","Species \nLosses","Reordering","Mean \nChange","Dispersion \nDifferences","Curve \nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor, oma=c(4,4,4,10))
# 
# ## reodering static richness/eveness
# pairs(sim_allmetrics[,c(5,6,11)], col=sim_allmetrics$comtype3, labels=c("Richness", "Evenness","Reordering"), font.labels=2, cex.labels=2, upper.panel = panel.cor, oma=c(4,4,4,10))
# 
# ## gain loss static richness/eveness
# pairs(sim_allmetrics[,c(5,6,9:10)], col=sim_allmetrics$comtype3, labels=c("Richness", "Evenness","Gains","Losses"), font.labels=2, cex.labels=2, upper.panel = panel.cor, oma=c(4,4,4,10))
# 
# ## disp, mc, dstar static richness/eveness
# pairs(sim_allmetrics[,c(5,6,12,13,16)], col=sim_allmetrics$comtype3, labels=c("Richness", "Evenness","Mean \nChange","Dispersion \nDifference","Curve \nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor, oma=c(4,4,4,10))



##codyn graphs
pairs(codyndat_allmetrics[,c(3:6)],labels=c("Richness", "Evenness \n(EQ)","Evenness \n(Simpsons)","Evenness \n(Gini)"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)

pairs(codyndat_allmetrics[,c(7:14)], col=codyndat_allmetrics$taxa, labels=c("Richness\nChange", "Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange","Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)


#how do these correlate with experiment parameters. #remove outliers
codyndat_allmetrics2<-codyndat_allmetrics%>%
  mutate(spatialExtent=log(spatial_extent),
         plotSize=log(plot_size))

pairs(codyndat_allmetrics2[,c(10:11, 23,33, 32,36,37)], labels=c("Mean \nChange","Dispersion \nDifference","MAP","MAT","Number \nPlots","Spatial \nExtent","Plot \nSize"), font.labels=2, cex.labels=2, upper.panel = panel.cor)

pairs(codyndat_allmetrics2[,c(7:14,32,36,37)], font.labels=2, cex.labels=2, upper.panel = panel.cor)
cor.test(codyndat_allmetrics$mean_change, codyndat_allmetrics$MAP_mm)

cor(codyndat_allmetrics2[,c(7:14,32,36,37)], method ="pearson")

#evenness
cor.test(sim_allmetrics$S, sim_allmetrics$E_Q)
cor.test(sim_allmetrics$S, sim_allmetrics$Gini)
cor.test(sim_allmetrics$S, sim_allmetrics$E_simp)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$E_Q)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$Gini)
cor.test(codyndat_allmetrics$S, codyndat_allmetrics$E_simp)



###LOOKING AT STATIC RICHNESS AND EVENNESS
rich_gain<-ggplot(data=sim_allmetrics2, aes(x=Sp, y=G, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"), labels=c("High Turnover &\nHigh Spatial Variability","High Turnover &\nLow Spatial Variability","Low Turnover &\nHigh Spatial Variability","Low Turnover &\nLow Spatial Variability"))+
  xlab("Richness")+
  ylab("Species Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_loss<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=L, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Species Losses")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_reorder<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=R, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Rank Change")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_comp<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=mean_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Compositional Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_disp<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=dispersion_diff, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Dispersion Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_curve<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=Dstar, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

even_gain<-
  ggplot(data=sim_allmetrics, aes(x=E_Q, y=G, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Species Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_loss<-
  ggplot(data=sim_allmetrics, aes(x=E_Q, y=L, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Species Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_reorder<-
  ggplot(data=sim_allmetrics, aes(x=E_Q, y=R, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Reordering")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_comp<-
  ggplot(data=sim_allmetrics, aes(x=E_Q, y=mean_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Compositional Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_disp<-ggplot(data=sim_allmetrics, aes(x=E_Q, y=dispersion_diff, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Dispersion Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_curve<-
  ggplot(data=sim_allmetrics, aes(x=E_Q, y=Dstar, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank())

grid.arrange(arrangeGrob(rich_reorder+theme(legend.position="none"),
                         even_reorder+theme(legend.position="none"),
                         rich_gain+theme(legend.position="none"),
                         even_gain+theme(legend.position="none"),
                         rich_loss+theme(legend.position="none"),
                         even_loss+theme(legend.position="none"),
                         rich_comp+theme(legend.position="none"),
                         even_comp+theme(legend.position="none"),
                         rich_disp+theme(legend.position="none"),
                         even_disp+theme(legend.position="none"),
                         rich_curve+theme(legend.position="none"),
                         even_curve+theme(legend.position="none"),
                         ncol=2))


# example of a curve comparision for the paper ----------------------------


#######example in paper for curve comparision.
time<-c(1,1,1,1,1,1,2,2,2,2,2,2)
sp<-c(1,2,3,4,5,6,1,3,4,6,7,8)
relrank<-c(0.333,0.5, 0.667, 0.167, 1, 0.833, 0.167, 0.583, 0.333, 1, 0.833, 0.583)
cumabund<-c(90,110,125,50,132,131,70,130,110,163,161,150)

df<-data.frame(time, sp, relrank, cumabund)%>%
  arrange(relrank)

result <- df %>%
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

ggplot(df, aes(x=relrank, y=cumabund, group=time))+
  geom_step(color=time)+
  xlab("Relative Rank")+
  ylab("Cumulative Abundance")+
  theme_bw()
