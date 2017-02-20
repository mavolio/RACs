setwd("~/Dropbox/converge_diverge/datasets/Longform")

library(ggplot2)
library(vegan)

dat<-read.csv("SpeciesRelativeAbundance_Dec2016.csv")

##pplots
pplots<-dat%>%
  filter(project_name=="pplots",calendar_year==2014)%>%
  tbl_df()%>%
  group_by(treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(rank=rank(-relcov, ties.method="average"))

##RACS

ggplot(data=pplots, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

pplots_common<-pplots%>%
  filter(rank<11)

ggplot(data=pplots_common, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

ggplot(data=pplots, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

###NMDS
pplots_wide<-dat%>%
  filter(project_name=="pplots",calendar_year==2014)%>%
  select(treatment, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-pplots_wide[,1:2]
mds<-metaMDS(pplots_wide[,3:53], autotransform=FALSE, shrink=FALSE)
mds

adonis(pplots_wide[,3:53]~treatment, pplots_wide)

#differences in dispersion?
dist<-vegdist(pplots_wide[,3:53])
betadisp<-betadisper(dist,pplots_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=treatment))+
  geom_point(size=5)


##wenndex
wenndex<-dat%>%
  filter(project_name=="WENNDEx", calendar_year==2012)%>%
  tbl_df()%>%
  group_by(treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(rank=rank(-relcov, ties.method="average"))

ggplot(data=wenndex, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

wenndex_common<-wenndex%>%
  filter(rank<11)

ggplot(data=wenndex_common, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

###NMDS
wenndex_wide<-dat%>%
  filter(project_name=="WENNDEx", calendar_year==2012)%>%
  select(treatment, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-wenndex_wide[,1:2]
mds<-metaMDS(wenndex_wide[,3:39], autotransform=FALSE, shrink=FALSE)
mds

adonis(wenndex_wide[,3:39]~treatment, wenndex_wide)

#differences in dispersion?
dist<-vegdist(wenndex_wide[,3:39])
betadisp<-betadisper(dist,wenndex_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=treatment))+
  geom_point(size=5)


##herbdiv
herbdiv<-dat%>%
  filter(project_name=="herbdiv", calendar_year==2015)%>%
  tbl_df()%>%
  group_by(treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(rank=rank(-relcov, ties.method="average"))

ggplot(data=herbdiv, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

herbdiv_common<-herbdiv%>%
  filter(rank<11)

ggplot(data=herbdiv_common, aes(x=rank, y=relcov))+
  geom_point(aes(color=genus_species))+
  facet_wrap(~treatment)

###NMDS
herbdiv_wide<-dat%>%
  filter(project_name=="herbdiv", calendar_year==2015)%>%
  select(treatment, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-herbdiv_wide[,1:2]
mds<-metaMDS(herbdiv_wide[,3:85], autotransform=FALSE, shrink=FALSE)
mds

adonis(herbdiv_wide[,3:85]~treatment, herbdiv_wide)

#differences in dispersion?
dist<-vegdist(herbdiv_wide[,3:39])
betadisp<-betadisper(dist,herbdiv_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=treatment))+
  geom_point(size=5)


