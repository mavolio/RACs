library(tidyverse)
library(codyn)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(reldist)
#library(gtools)


# Read in Data ------------------------------------------------------------
# for the sim dataset, communites are differentiated by alpha (richness), theta (evenness) and scenario (rate of turnover and spatial heterogeniety: four scenarios: a: high turnover, high spatial heterogeniety; b: low turnover, low spatial heterogeniety; c: low turnover, high spatial heterogeniety; d: high turnover, low spatial heterogeniety"). For each richness-evennes combination (9 combinations) there are each community type. Each of these 10 community types have 10 replicates, called "sites" at a given point in time. Each community type, time, and site is then replicated 10 time.

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
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

#work
sim<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\SimCom_Sept28.csv')%>%
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
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep = "_"))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))



# Richness Evenness Metrics -----------------------------------------------

#codyn dataset
spc<-unique(codyndat_clean$site_project_comm)
codyn_div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "EQ")
  out$site_project_comm<-spc[i]
  
  codyn_div_eq<-rbind(codyn_div_eq, out)
}

codyn_div_esimp<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "SimpsonEvenness")
  out$site_project_comm<-spc[i]
  
  codyn_div_esimp<-rbind(codyn_div_esimp, out)
}

codyn_div_evar<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "Evar")
  out$site_project_comm<-spc[i]
  
  codyn_div_evar<-rbind(codyn_div_evar, out)
}


codyn_div1<-merge(codyn_div_eq, codyn_div_esimp, by=c("site_project_comm",'plot_id','richness','experiment_year'))
codyn_div<-merge(codyn_div1, codyn_div_evar, by=c("site_project_comm",'plot_id','richness','experiment_year'))

codyndat_diversity_mean <- codyn_div%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(Sp=mean(richness),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(SimpsonEvenness),
            Evar=mean(Evar, na.rm=T))


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
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "EQ")
  out$id<-com_rep[i]
  
  sim_div_eq<-rbind(sim_div_eq, out)  
}

sim_div_esimp<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "SimpsonEvenness")
  out$id<-com_rep[i]
  
  sim_div_esimp<-rbind(sim_div_esimp, out)  
}

sim_div_evar<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "Evar")
  out$id<-com_rep[i]
  
  sim_div_evar<-rbind(sim_div_evar, out)  
}

sim_div1<-merge(sim_div_eq, sim_div_esimp, by=c("id",'site','richness','time'))
sim_div<-merge(sim_div1, sim_div_evar, by=c("id",'site','richness','time'))

sim_diversity_mean<-sim_div%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, site)%>%#average over replicates
  summarize(Sp=mean(richness),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(SimpsonEvenness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(id3, time)%>%#average over sites
  summarize(Sp=mean(Sp),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(ESimp),
            Evar=mean(Evar, na.rm=T))

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


# RAC changes -------------------------------------------------------------


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
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
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
  group_by(id3, time, time2, site)%>%
  summarize(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))%>%
  ungroup()%>%
  group_by(id3, time, time2)%>%
  summarize(S=mean(S),
            E=mean(E,na.rm=T),
            R=mean(R),
            G=mean(G),
            L=mean(L))

# Mean Change and Dispersion ----------------------------------------------
#codyn dataset

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
  group_by(id3, time, time2)%>%
  summarize(composition_change=mean(composition_change),
            dispersion_change=mean(dispersion_change))

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

codyn_curvechange_mean<-codyn_curvechange%>% 
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
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
  group_by(id, time, time2)%>%
  summarise(curve_change=mean(curve_change))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2)%>%
  summarize(curve_change=mean(curve_change))



# Merging all metrics to single datasets ----------------------------------

####MERGING TO A SINGE DATASET

codyndat_allmetrics<-codyn_multchange%>%
  left_join(codyndat_rac_change_average)%>%
  left_join(codyn_curvechange_mean)%>%
  left_join(codyn_div_all)%>%
  left_join(codyndat_info)

#sim
sim_all_metrics<-sim_multchange_mean%>%
  left_join(sim_rac_change_mean)%>%
  left_join(sim_cc_ave)%>%
  left_join(sim_div_all)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))

write.csv(codyndat_allmetrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\codyn_allmetrics_April2018.csv', row.names = F)
write.csv(sim_all_metrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_April2018.csv', row.names = F)

# pair plot graphs --------------------------------------------------------

theme_set(theme_bw(12))

codyndat_allmetrics<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\codyn_allmetrics_April2018.csv')%>%
  mutate(absS = abs(S),
         absE = abs(E))

sim_allmetrics<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_Jan2018.csv')

#graphing this
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

##codyn graphs
colnames(codyndat_allmetrics)
pairs(codyndat_allmetrics[,c(6:10,3,4,11)], col=codyndat_allmetrics$taxa, labels=c("Richness\nChange", "Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange","Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)

#same figure with absolute value of S and E
pairs(codyndat_allmetrics[,c(38,39, 8:10,3,4,11)], col=codyndat_allmetrics$taxa, labels=c(" Abs. Richness\nChange", "Abs. Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange","Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)

#how do these correlate with experiment parameters. #remove outliers
codyndat_allmetrics2<-codyndat_allmetrics%>%
  mutate(spatialExtent=log(spatial_extent),
         plotSize=log(plot_size))
colnames(codyndat_allmetrics2)

pairs(codyndat_allmetrics2[,c(3,4,38,39,8:11,34,40,41)], upper.panel = panel.cor)

cor(codyndat_allmetrics2[,c(3,4,38,39,8:11,34,40,41)], method ="pearson")



#evenness
cor.test(sim_div_all$Sp, sim_div_all$EQ)
cor.test(sim_div_all$Sp, sim_div_all$EGini)
cor.test(sim_div_all$Sp, sim_div_all$ESimp)
cor.test(sim_div_all$Sp, sim_div_all$Evar)
cor.test(codyn_div_all$Sp, codyn_div_all$EQ)
cor.test(codyn_div_all$Sp, codyn_div_all$EGini)
cor.test(codyn_div_all$Sp, codyn_div_all$ESimp)
cor.test(codyn_div_all$Sp, codyn_div_all$Evar)

pairs(codyn_div_all[,4:7])
pairs(sim_div_all[,4:7])

#other correlations
cor.test(sim_allmetrics$Sp, sim_allmetrics$R)
cor.test(sim_allmetrics$Sp, sim_allmetrics$G)
cor.test(sim_allmetrics$Sp, sim_allmetrics$L)
cor.test(sim_allmetrics$Sp, sim_allmetrics$composition_change)
cor.test(sim_allmetrics$Sp, sim_allmetrics$dispersion_change)
cor.test(sim_allmetrics$Sp, sim_allmetrics$curve_change)
cor.test(sim_allmetrics$EQ, sim_allmetrics$R)
cor.test(sim_allmetrics$EQ, sim_allmetrics$G)
cor.test(sim_allmetrics$EQ, sim_allmetrics$L)
cor.test(sim_allmetrics$EQ, sim_allmetrics$composition_change)
cor.test(sim_allmetrics$EQ, sim_allmetrics$dispersion_change)
cor.test(sim_allmetrics$EQ, sim_allmetrics$curve_change)



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
  ggplot(data=sim_allmetrics, aes(x=Sp, y=composition_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Compositional Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_disp<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=dispersion_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Dispersion Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

rich_curve<-
  ggplot(data=sim_allmetrics, aes(x=Sp, y=curve_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Richness")+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

even_gain<-
  ggplot(data=sim_allmetrics, aes(x=EQ, y=G, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Species Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_loss<-
  ggplot(data=sim_allmetrics, aes(x=EQ, y=L, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Species Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_reorder<-
  ggplot(data=sim_allmetrics, aes(x=EQ, y=R, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Reordering")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_comp<-
  ggplot(data=sim_allmetrics, aes(x=EQ, y=composition_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Compositional Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_disp<-ggplot(data=sim_allmetrics, aes(x=EQ, y=dispersion_change, col=comtype2))+
  geom_point()+
  scale_color_manual(name="Community Type", values=c("black","red","green","blue"), breaks=c("a","d","c","b"))+
  xlab("Evenness")+
  ylab("Dispersion Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank())

even_curve<-
  ggplot(data=sim_allmetrics, aes(x=EQ, y=curve_change, col=comtype2))+
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


df <- df[order(df$time, df$cumabund),]

timestep2 <- unique(df$time)#assumes this is a length of 2

df1 <- df[df$time == timestep2[1],]
df2 <- df[df$time == timestep2[2],]
sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
r <- sort(unique(c(0, df1$relrank, df2$relrank)))
h <- abs(sf1(r) - sf2(r))
w <- c(diff(r), 0)
CC=sum(w*h)

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
