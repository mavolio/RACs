##worked example with pplots

setwd("~/Dropbox/converge_diverge/datasets/Longform")

library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)

theme_set(theme_bw(20))

dat<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Oct2017.csv")

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\SpeciesRelativeAbundance_Oct2017.csv")

##pplots
pplots<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2014, treatment=="N1P0"|treatment=="N2P3")%>%
  tbl_df()%>%
  group_by(calendar_year, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         treatment=factor(treatment, levels=c("N1P0","N2P3")))###need to do this b/c otherwise R will remember every treatment

##step 1. do NMDS of pretreatment and last year of data
pplots_wide<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2014,treatment=="N1P0"|treatment=="N2P3")%>%
  select(treatment, calendar_year, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-pplots_wide[,1:3]
mds<-metaMDS(pplots_wide[,4:54], autotransform=FALSE, shrink=FALSE)
mds

adonis(pplots_wide[,4:54]~treatment, pplots_wide)

#differences in dispersion?
dist<-vegdist(pplots_wide[,4:54])
betadisp<-betadisper(dist,pplots_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run

##facet_wrap
toplot<-scores2
##plot NMDS
ggplot(subset(toplot, treatment=="N1P0"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","lightgray"), labels=c("Pre-Treatment","Experiment Year 12"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

ggplot(subset(toplot, treatment=="N2P3"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","lightgray"), labels=c("Pre-Treatment","Experiment Year 12"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")





##plot RACS
ractoplot<-pplots

ggplot(data=ractoplot, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species), size=3)+
  geom_line(aes(group=plot_id), size=0.2)+
  facet_wrap(~treatment+calendar_year, ncol=2)


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
  mutate(year=as.factor(calendar_year))
  

ggplot(data=average, aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=1)+
  facet_wrap(~treatment, ncol=2)

###cherry picked examples
average_test<-pplots%>%
  group_by(calendar_year, treatment, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov))

result <- average_test %>%
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
    data.frame(
     # Dmax=max(h),
      Dstar=sum(w*h))
  })

trt<-pplots%>%
  tbl_df()%>%
  select(plot_id, treatment)%>%
  unique()
results_test<-merge(result, trt, by="plot_id")

plot(results_test$Dmax, results_test$Dstar)

boxplot(results_test$Dmax~results_test$treatment)
boxplot(results_test$Dstar~results_test$treatment)
anova(lm(Dstar~treatment, data=results_test))
TukeyHSD(aov(Dstar~treatment, data=results_test))

ggplot(data=average_test, aes(x=relrank, y=cumabund, group=interaction(calendar_year,plot_id)))+
  geom_point(size=3, aes(color=as.factor(calendar_year)))+
  geom_step(size=1, aes(color=as.factor(plot_id)))+##should we graph with geom_step?
  facet_wrap(~treatment)


theme_set(theme_bw(12))
###read in datasets where we already hae this
vpplots<-read.csv("CORRE_RAC_Metrics_Oct2017_allyears_2.csv")%>%
  filter(site_project_comm=="KNZ_pplots_0")%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  gather(metric, value, S_diff:dispersion_diff)

ggplot(data=vpplots, aes(x=calendar_year, y=value, color=treatment))+
  geom_point()+
  geom_line()+
  facet_wrap(~metric, ncol=3, scale="free")

hpplots<-read.csv("CORRE_ContTreat_Compare_OCT2017.csv")%>%
  filter(site_project_comm=="KNZ_pplots_0")%>%
  filter(treatment=="N2P3")%>%
  select(-plot_mani)%>%
  gather(metric, value, PCSdiff:disp_diff)

ggplot(data=hpplots, aes(x=calendar_year, y=value))+
  geom_point()+
  geom_line()+
  facet_wrap(~metric, ncol=3, scale="free")
