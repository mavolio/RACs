library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

###CLEANING CODYN DATASET
#restrict to species that are present in an experiment
splist<-codyndat%>%
  group_by(site_code, project_name, community_type, species)%>%
  summarize(present=sum(abundance))%>%
  filter(present!=0)%>%
  select(-present)

#merge back and will drop species that do not exist in a dataset
codyndat_clean<-merge(codyndat, splist, by=c("site_code","project_name","community_type","species"))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))

##richness
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

codyndat_diversity <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            E_simp=E_simp(abundance))

#Whittacker beta diversity (gamma / average alpha)
gammadiv<-codyndat_clean%>%
  filter(abundance!=0)%>%
  group_by(site_project_comm, experiment_year, species)%>%
  summarize(splist=mean(abundance))%>%
  ungroup()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(gamma=length(species))

alphadiv<-codyndat_diversity%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(alpha=mean(S))

beta_div<-merge(gammadiv, alphadiv, by=c("site_project_comm","experiment_year"))%>%
  mutate(wbeta=gamma/alpha)


#turnover
turnover<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="total")%>%
  separate(id, c("site_project_comm", "plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(total=mean(total))

hist(codyndat_diversity$S)
hist(codyndat_diversity$E_simp)
hist(codyndat_diversity$E_q)
hist(beta_div$wbeta)
hist(turnover$total)
