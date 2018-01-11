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
  summarize(Sp=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(Sp=mean(Sp),
            E_Q=mean(E_q, na.rm=T),
            Gini=mean(Gini),
            E_simp=mean(E_simp))

sim_diversity<-group_by(sim, id, site, time)%>%
  summarize(Sp=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))%>%
  ungroup()%>%
  group_by(id, time)%>%
  summarize(Sp=mean(Sp),
            E_Q=mean(E_q, na.rm=T),
            Gini=mean(Gini),
            E_simp=mean(E_simp))%>%
  ungroup()%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(Sp=mean(Sp),
            E_Q=mean(E_Q),
            Gini=mean(Gini),
            E_simp=mean(E_simp))
  
###graph this
pairs(sim_diversity[3:6])
pairs(codyndat_diversity[3:6])


###Looking at RAC Changes
##Codyn Dataset

#problems with _ in rep names
#RUNNING all by combine site_project_name with plot_id
codyn_rac_change<-RAC_change(df = codyndat_clean, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "id")

codyndat_rac_change_average<-codyn_rac_change%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year_pair)%>%
  summarise(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))

##SIM dataset
sim_rac_changes <- RAC_change(df = sim, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "id2")

sim_rac_changes_mean<-sim_rac_changes%>%
  separate(id2, c("id","site"), sep="::")%>%
  group_by(id, time_pair)%>%
  summarise(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time_pair)%>%
  summarize(S=mean(S),
            E=mean(E,na.rm=T),
            R=mean(R),
            G=mean(G),
            L=mean(L))

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

#clean the data first, drop plots that are not measured both comparision years, and plots that only have 1 species.
codyn_clean2=data.frame()

spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])%>%
    filter(abundance!=0)
  
  spc_id2<-spc[i]
  
  timestep<-sort(unique(subset$experiment_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(experiment_year==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(plot_id)%>%
      unique()
    
    subset_t2<-subset%>%
      filter(experiment_year==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(plot_id)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="plot_id")
    #dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
    ##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="plot_id")
    
    #dropping plots with only 1 species in any of the years    
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
      select(-min)%>%
      ungroup()
    
    codyn_clean2<-rbind(codyn_clean2, subset_t12_3)
  }
}
    
    
#codyn_cc <-curve_change(df = codyn_clean2, "experiment_year", "species", "abundance", "id")
#not working

codyn_curvechange<-data.frame()
spc<-unique(codyn_clean2$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyn_clean2%>%
    filter(site_project_comm==spc[i])
  
  out <- curve_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_curvechange<-rbind(codyn_curvechange, out)  
}

codyndat_dstar<-d_output%>% 
  group_by(site_project_comm, experiment_year)%>%
  summarise(Dstar=mean(Dstar))

#sim dataset

sim_cc <- curve_change(sim, "time", "species", "abundance", "id2")

sim_info<-sim%>%
  select(id, id2, id3)%>%
  unique()

sim_cc_merge<-merge(sim_info, sim_cc, by="id2")

sim_cc_ave<-sim_cc_merge%>% 
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
merge1<-merge(codyndat_diversity, codyndat_reorder, by=c("site_project_comm","experiment_year"))
merge2<-merge(merge1, codyndat_braycurtis, by=c("site_project_comm","experiment_year"))
merge3<-merge(merge2, codyndat_dstar, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge3, codyndat_info, by="site_project_comm")

#sim
merge1<-merge(sim_diversity, sim_reorder, by=c("id3","time"))
merge2<-merge(merge1, sim_bray_curtis, by=c("id3","time"))
sim_allmetrics<-merge(merge2, sim_dstar, by=c("id3","time"))%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")

sim_allmetrics$comtype2<-as.factor(sim_allmetrics$comtype)

write.csv(codyndat_allmetrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\codyn_allmetrics_diff_corrected.csv')
write.csv(sim_allmetrics,'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_diff_corrected.csv')

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
