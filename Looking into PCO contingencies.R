setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")

library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(devtools)
install_github("mavolio/codyn", ref = "RACs_cleaner")
library(codyn)

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)
traits<-read.csv("C:\\Users\\megha\\Dropbox\\pplots\\site_review_2017\\traits_final.csv")


##pplots

pplots_allyears<-dat%>%
  filter(project_name=="pplots", treatment=="N1P0"|treatment=="N2P3"|treatment=="N2P0")

pplots_wide<-pplots_allyears%>%
  select(treatment, calendar_year, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)%>%
  mutate(trt_year = paste(treatment, calendar_year, sep = '_'))

#get centroid distances all in one step
dist<-vegdist(pplots_wide[,4:86])
betadisp<-betadisper(dist,pplots_wide$trt_year,type="centroid")
cent_dist <- as.matrix(vegdist(betadisp$centroids, method = "euclidean"))
#extracting distances
names<-row.names(cent_dist)
names2<-names[2:39]
timestep <- sort(unique(pplots_wide$calendar_year))
cent_dist_yrs <- data.frame(
    trt_year = names2,
    composition_change_all = diag(
    as.matrix(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist) - 1)])))

cent_dist2<-cent_dist_yrs%>%
  separate(trt_year, into=c("treatment", "calendar_year2"), sep = '_')%>%
  mutate(calendar_year = as.numeric(cent_dist2$calendar_year2)-1,
         calendar_year2 = as.numeric(calendar_year2))%>%
  filter(calendar_year!=2001)


#comparing with function
pplots_allyears<-dat%>%
  filter(project_name=="pplots", treatment=="N1P0"|treatment=="N2P3"|treatment=="N2P0")

mult_change_allyears<-multivariate_change(pplots_allyears, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id", treatment.var = "treatment")


###putting together
alldat<-mult_change_allyears%>%
  left_join(cent_dist2)

ggplot(data=alldat, aes(x = composition_change, y = composition_change_all))+
         geom_point()+
  xlab("calculated for each treatment sepeartely")+
  ylab("calculated fro all treatments together")+
  geom_abline(slope = 1)

#package output
package<-ggplot(data=mult_change_allyears, aes(x=calendar_year2, y=composition_change, group=treatment, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Compositional Change")+
  xlab("Year-Pair")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  ggtitle("Treatments seperate")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())

#all together output
other<-ggplot(data=cent_dist2, aes(x=calendar_year2, y=composition_change_all, group=treatment, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Compositional Change")+
  xlab("Year-Pair")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  ggtitle("all data together")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())

grid.arrange(package, other, ncol=2)
