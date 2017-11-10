##worked example with pplots

setwd("~/Dropbox/converge_diverge/datasets/Longform")
setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")


library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)

theme_set(theme_bw(20))

# read in data ------------------------------------------------------------


dat<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Oct2017.csv")

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\SpeciesRelativeAbundance_Oct2017.csv")
traits<-read.csv("C:\\Users\\megha\\Dropbox\\pplots\\site_review_2017\\traits_final.csv")

##pplots
pplots<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2011, treatment=="N1P0"|treatment=="N2P3")%>%
  tbl_df()%>%
  group_by(calendar_year, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         treatment=factor(treatment, levels=c("N1P0","N2P3")))###need to do this b/c otherwise R will remember every treatment


# NMDS --------------------------------------------------------------------


##step 1. do NMDS of pretreatment and last year of data
pplots_wide<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2011,treatment=="N1P0"|treatment=="N2P3")%>%
  select(treatment, calendar_year, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-pplots_wide[,1:3]
mds<-metaMDS(pplots_wide[,4:53], autotransform=FALSE, shrink=FALSE)
mds

adonis(pplots_wide[,4:53]~treatment, pplots_wide)

#differences in dispersion?
dist<-vegdist(pplots_wide[,4:53])
betadisp<-betadisper(dist,pplots_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run


toplot<-scores2
##plot NMDS
ggplot(subset(toplot, treatment=="N1P0"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","darkgray"), labels=c("Pre-Treatment","Experiment Year 12"))+
  scale_x_continuous(limits=c(-.8,1))+
  scale_y_continuous(limits=c(-1,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

ggplot(subset(toplot, treatment=="N2P3"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","darkgray"), labels=c("Pre-Treatment","Experiment Year 12"))+
  scale_x_continuous(limits=c(-.8,1))+
  scale_y_continuous(limits=c(-1,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


# Make RACs ---------------------------------------------------------------


##plot RACS

traits2<-traits%>%
  mutate(genus=tolower(Genus),
         genus_species=paste(genus, Species, sep=" "))

ractoplot<-merge(pplots, traits2, by="genus_species", all=T)%>%
  select(-Spnum, -Genus, -Species, -Family, -cot, -canopy, -clonality, -native, -bloom, -Seed_weight__g_, -genus)%>%
  mutate(cat=ifelse(life=="P"&form=="G"&C3_C4=="C4","G1", ifelse(life=="P"&form=="G"&C3_C4=="C3","G2", ifelse(life=="A"&form=="G", "G3", ifelse(life=="P"&form=="F"&n_fixer=="N", "F1",ifelse(life=="P"&form=="S"&n_fixer=="N", "F1", ifelse(life=="P"&form=="F"&n_fixer=="Y", "F2", ifelse(life=="P"&form=="S"&n_fixer=="Y","F2", "F3"))))))))%>%
  na.omit%>%
  mutate(colorfill=ifelse(genus_species=="andropogon gerardii","ag",ifelse(genus_species=="andropogon scoparius","as", ifelse(genus_species=="sorghastrum nutans", "sn", ifelse(genus_species=="solidago canadensis","sc",ifelse(genus_species=="solidago missouriensis", "sm", ifelse(genus_species=="oxalis stricta", "os", "other")))))))

#label top 5 species in contol plots in 2002 and everythign else gray. 

ggplot(data=subset(ractoplot, treatment=="N1P0"&calendar_year==2002) , aes(x=rank, y=relcov))+
   geom_line(aes(group=plot_id), color="black", size=1)+
   geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, treatment=="N1P0"&calendar_year==2011) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

ggplot(data=subset(ractoplot, treatment=="N2P3"&calendar_year==2002) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="black", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, treatment=="N2P3"&calendar_year==2011) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","pink","green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

# make CC examples --------------------------------------------------------


###figuring out cumulative curve
average<-ractoplot%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(sumabund=sum(relcov),
         relcov2=relcov/sumabund)%>%
  mutate(rank=rank(-relcov2, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov2))%>%
  mutate(year=as.character(calendar_year))
  
ccplot<-ractoplot%>%
  mutate(sumabund=sum(relcov),
         relcov2=relcov/sumabund)%>%
  mutate(rank=rank(-relcov2, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov2))%>%
  mutate(year=as.character(calendar_year))


ggplot(data=subset(ccplot, treatment=="N1P0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=1)+
  scale_color_manual(values=c("black","darkgray"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
    ylab("Cumulative Relative Abundance")+
    xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~plot_id, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


ggplot(data=subset(ccplot, treatment=="N2P3"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=1)+
  scale_color_manual(values=c("black","darkgray"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~plot_id, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


# doing SERGL -------------------------------------------------------------


#doing SERGL on pplots data
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

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


##rank shift
rank_pres<-pplots%>%
  filter(relcov!=0)%>%
  tbl_df()%>%
  group_by(treatment, calendar_year, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method = "average"))%>%
  tbl_df()%>%
  select(-X)

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
addzero <- pplots %>%
  select(-X, -rank)%>%
  spread(genus_species, relcov, fill=0)%>%
  gather(genus_species, relcov, 9:58)


zeros<-addzero%>%
  filter(relcov==0)
##get species richness for each year
pplots_S<-group_by(pplots, treatment, calendar_year, plot_id)%>%
  summarize(S=S(relcov))
##merge together make zero abundances rank S+1
zero_rank<-merge(zeros, pplots_S, by=c("treatment","calendar_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
rank<-rbind(rank_pres, zero_rank)

##calculate reordering between time steps 3 ways, rank correlations, mean rank shifts not corrected, and mean ranks shifts corrected for the size of the speceis pool

reordering=data.frame(id=c(), experiment_year=c(), S=c(), E=c(), R=c(), G=c(), L=c())#expeiment year is year of timestep2

spc_id<-unique(rank$plot_id)

for (i in 1:length(spc_id)){
  subset<-rank%>%
    filter(plot_id==spc_id[i])
  id<-spc_id[i]
  
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
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("genus_species","plot_id"), all=T)%>%
      filter(relcov.x!=0|relcov.y!=0)
    
    #reodering
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    #evnness and richenss
    s_t1 <- S(subset_t12$relcov.x)
    e_t1 <- E_q(subset_t12$relcov.x)
    s_t2 <- S(subset_t12$relcov.y)
    e_t2 <- E_q(subset_t12$relcov.y)
    
    sdiff<-abs(s_t1-s_t2)/nrow(subset_t12)
    ediff<-abs(e_t1-e_t2)/nrow(subset_t12)
    
    #gains and losses
    subset_t12$gain<-ifelse(subset_t12$relcov.x==0, 1, 0)
    subset_t12$loss<-ifelse(subset_t12$relcov.y==0, 1, 0)
    
    gain<-sum(subset_t12$gain)/nrow(subset_t12)
    loss<-sum(subset_t12$loss)/nrow(subset_t12)
    
    metrics<-data.frame(plot_id=id, calendar_year=timestep[i+1], S=sdiff, E=ediff, R=MRSc, G=gain, L=loss)
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }

}

##doing curve change
d_output=data.frame(calendar_year=c(), plot_id=c(), Darea=c())#expeiment year is year of timestep2

 rel_ranks<-pplots%>%
    filter(relcov!=0)%>%
    group_by(calendar_year, plot_id)%>%
    mutate(rank=rank(-relcov, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-relcov)%>%
    mutate(cumabund=cumsum(relcov))%>%
    ungroup()
  
timestep<-sort(unique(rel_ranks$calendar_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-rel_ranks%>%
      filter(calendar_year==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(plot_id)%>%
      unique()
    
    subset_t2<-rel_ranks%>%
      filter(calendar_year==timestep[i+1])
    
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
      group_by(calendar_year, plot_id)%>%
      mutate(numplots=length(plot_id))%>%
      ungroup()%>%
      group_by(plot_id)%>%
      mutate(min=min(numplots))%>%
      select(plot_id, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="plot_id")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(calendar_year, plot_id)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
      group_by(plot_id) %>%
      do({
        y <- unique(.$calendar_year)###assumption this is a length 2 list
        df1 <- filter(., calendar_year==y[[1]])
        df2 <- filter(., calendar_year==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(Dstar=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(calendar_year=timestep[i+1], plot_id=result$plot_id, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }

trts<-pplots%>%
  ungroup()%>%
  select(plot_id, treatment)%>%
  unique()

merge1<-merge(reordering, d_output, by=c("plot_id", "calendar_year"))
allmetrics<-merge(merge1, trts, by="plot_id")%>%
  gather(metric, value, S:Dstar)%>%
  group_by(treatment, calendar_year, metric)%>%
    summarize(vmean=mean(value),
              vn=length(plot_id),
              vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))
            
theme_set(theme_bw(12))
S<-ggplot(data=subset(allmetrics, metric=="S"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Richness Changes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
E<-ggplot(data=subset(allmetrics, metric=="E"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Evenness Changes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
R<-ggplot(data=subset(allmetrics, metric=="R"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Rank Changes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate('text', label="*", x=1.5, y = 0.25, size=10)
G<-ggplot(data=subset(allmetrics, metric=="G"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Speices Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
L<-ggplot(data=subset(allmetrics, metric=="L"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Species Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
C<-ggplot(data=subset(allmetrics, metric=="Dstar"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N+P"))+
  ylab("Curve Changes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(S, E,R, G, L, C, ncol=3)

allmetrics_full<-allmetrics<-merge(merge1, trts, by="plot_id")%>%
  gather(metric, value, S:Dstar)

t.test(value~treatment, data=subset(allmetrics_full, metric=="S")) # p = 0.261
t.test(value~treatment, data=subset(allmetrics_full, metric=="E"))  # p = 0.354
t.test(value~treatment, data=subset(allmetrics_full, metric=="R"))  # p < 0.001
t.test(value~treatment, data=subset(allmetrics_full, metric=="G")) # p = 0.196
t.test(value~treatment, data=subset(allmetrics_full, metric=="L")) # p = 0.088
t.test(value~treatment, data=subset(allmetrics_full, metric=="Dstar")) # p = 0.0422


# doing Curve Differences -------------------------------------------------


####doing curve change
average_test<-pplots%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov))

##compare in 2002
pretreat<-subset(average_test, calendar_year==2002)

df1 <- filter(pretreat, treatment=="N1P0")
df2 <- filter(pretreat, treatment=="N2P3")
sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
r <- sort(unique(c(0, df1$relrank, df2$relrank)))
h <- abs(sf1(r) - sf2(r))
w <- c(diff(r), 0)
Dstar=sum(w*h)

nine<-subset(average_test, calendar_year==2011)

df1 <- filter(nine, treatment=="N1P0")
df2 <- filter(nine, treatment=="N2P3")
sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
r <- sort(unique(c(0, df1$relrank, df2$relrank)))
h <- abs(sf1(r) - sf2(r))
w <- c(diff(r), 0)
Dstar=sum(w*h)

# appendix fig of change difference through time --------------------------


theme_set(theme_bw(16))
# ###read in datasets where we already have this
# 
vpplots<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\CORRE_RAC_Metrics_Nov2017_allReplicates.csv")%>%
  filter(site_project_comm=="KNZ_pplots_0")%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  select(-bc_dissim)%>%
  gather(metric, value, S:L)%>%
  group_by(treatment, calendar_year, metric)%>%
  summarize(vmean=mean(value),
            vn=length(plot_id),
            vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))
  
vpplots_bc<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\CORRE_RAC_Metrics_Nov2017.csv")%>%
  filter(site_project_comm=="KNZ_pplots_0")%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  select(-S, -E, -R, -G, -L)%>%
  gather(metric, value, mean_change:dispersion_diff)

bc<-
  ggplot(data=subset(vpplots_bc, metric=="mean_change"), aes(x=calendar_year, y=value, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","darkgray"))+
  geom_line(size=1)+
  ylab("Compositional Change")+
  xlab("Year")+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())
disp<-
  ggplot(data=subset(vpplots_bc, metric=="dispersion_diff"), aes(x=calendar_year, y=value, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","darkgray"), labels=c("Control","N+P"))+
  geom_line(size=1)+
  ylab("Dispersion Change")+
  xlab("Year")+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())


S<-
ggplot(data=subset(vpplots,metric=="S"), aes(x=calendar_year, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","darkgray"))+
  ylab("Richness Change")+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
E<-
ggplot(data=subset(vpplots,metric=="E"), aes(x=calendar_year, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","darkgray"))+
  ylab("Evenness Change")+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
R<-
  ggplot(data=subset(vpplots,metric=="R"), aes(x=calendar_year, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","darkgray"))+
  ylab("Rank Change")+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
G<-
  ggplot(data=subset(vpplots,metric=="G"), aes(x=calendar_year, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","darkgray"))+
  ylab("Species Gain")+
  xlab("Year")+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
L<-
  ggplot(data=subset(vpplots,metric=="L"), aes(x=calendar_year, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","darkgray"), name="Treatment", labels=c("Control", "N+P"))+
  ylab("Species Loss")+
  geom_line(size=1)+
  xlab("Year")+
  scale_x_continuous(breaks=c(2004, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

legend=gtable_filter(ggplot_gtable(ggplot_build(L)), "guide-box") 
grid.draw(legend)


grid.arrange(arrangeGrob(bc+theme(legend.position="none"),
                         disp+theme(legend.position="none"),
                         S+theme(legend.position="none"),
                         E+theme(legend.position="none"),
                         R+theme(legend.position="none"),
                         G+theme(legend.position="none"),
                         L+theme(legend.position="none"),
                         ncol=2),legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)


hpplots<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\CORRE_ContTreat_Compare_Nov2017.csv")%>%
  filter(site_project_comm=="KNZ_pplots_0")%>%
  filter(treatment=="N2P3"|treatment=="N2P3")%>%
  select(-plot_mani)%>%
  gather(metric, value, Sd:disp_diff)

bc_d<-
ggplot(data=subset(hpplots, metric=="mean_change"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Compositional Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
disp_d<-
  ggplot(data=subset(hpplots, metric=="disp_diff"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Dispersion Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
bc_d<-
  ggplot(data=subset(hpplots, metric=="mean_change"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Compositional Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
s_d<-
  ggplot(data=subset(hpplots, metric=="Sd"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Richness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
e_d<-
  ggplot(data=subset(hpplots, metric=="Ed"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Evenness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
r_d<-
  ggplot(data=subset(hpplots, metric=="Rd"), aes(x=calendar_year, y=value))+
  geom_point(color="black", size=3)+
  geom_line(color="black", size=1)+
  scale_x_continuous(breaks=c(2004, 2012))+
  ylab("Rank Difference")+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sp_d<-
    ggplot(data=subset(hpplots, metric=="spd"), aes(x=calendar_year, y=value))+
    geom_point(color="black", size=3)+
    geom_line(color="black", size=1)+
    scale_x_continuous(breaks=c(2004, 2012))+
    ylab("Species Differences")+
    xlab("Year")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  

grid.arrange(bc_d, disp_d, s_d, e_d, r_d, sp_d, ncol=2)


###examples for the appendix
t1=c(40,20,15,50,1,6,0,0)
E_q(t1)
t2=c(70,0,20,40,0,2,11,20)
E_q(t2)
t1_rel=c(0.30303,0.151515,0.113636,0.378788,0.007576,0.045455,0,0)
E_q(t1_rel)
