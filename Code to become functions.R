library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)

###get corre dataset to work on and focus on familiar examples
corre<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Dec2016.csv")%>%
  filter(project_name=="pplots"|project_name=="herbdiv"|project_name=="WENNDEx")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type,sep="_"))%>%
  select(-X)%>%
  mutate(id=paste(site_project_comm, plot_id, treatment, sep="::"))

#function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

###
### Functions that work on a plot level, no comparison
###

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

###Calculating the metrics. No need to do anything different here.
corre_diversity <- group_by(corre, site_project_comm, calendar_year, treatment, plot_id) %>% 
  summarize(S=S(relcov),
            E_q=E_q(relcov))%>%
  tbl_df()%>%
  group_by(site_project_comm, treatment, calendar_year)%>%
  summarize(S=mean(S),
            E_Q=mean(E_q))

##comparing plots through time (RACS). This is through time, do not need a treatment, can do on any dataset but a treatment column is also possible.

##use codyn package to calculate appearances and dissapearances. To do this need to incorporate year and treatment into the replicate var. There should be a better way to do this.

corre_id<-corre%>%
  mutate(id=paste(site_project_comm, treatment, plot_id, sep="::"))

loss<-turnover(df=corre_id, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="id", metric="disappearance")

gain<-turnover(df=corre_id, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="id", metric="appearance")

gains_loss<-merge(loss, gain, by=c("calendar_year","id"))%>%
  separate(id, c("site_project_comm","treatment", "plot_id"), sep="::")%>%
  group_by(site_project_comm, calendar_year, treatment)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

ggplot(data=subset(gains_loss,site_project_comm=="KNZ_pplots_0"), aes(x=calendar_year, y=gain, group=treatment))+
  geom_point(size=3, aes(color=treatment))+
  scale_color_manual(values=c("green","blue","blue","blue","red","purple","purple","purple"))+
  geom_line()

ggplot(data=subset(gains_loss,site_project_comm=="KNZ_pplots_0"), aes(x=calendar_year, y=loss, group=treatment))+
  geom_point(size=3, aes(color=treatment))+
  scale_color_manual(values=c("green","blue","blue","blue","red","purple","purple","purple"))+
  geom_line()

##reordering through time
reordering=data.frame(id=c(), calendar_year=c(), MRSc=c())

explist<-unique(corre$site_project_comm)

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corre%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  wide<-subset%>%
    gather(genus_species, relcov,11:ncol(subset))
  
  ##add ranks dropping zeros
  rank_pres<-wide%>%
    filter(relcov!=0)%>%
    tbl_df()%>%
    group_by(site_project_comm, calendar_year, treatment, plot_id)%>%
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-wide%>%
    filter(relcov==0)
  ##get species richness for each year
  rich<-group_by(wide, site_project_comm, calendar_year, plot_id)%>%
    summarize(S=S(relcov))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year","plot_id"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  #get uniuque id's
  spc_id<-unique(subset$id)
  
  for (i in 1:length(spc_id)){
    subset2<-rank%>%
      filter(id==spc_id[i])
    id<-spc_id[i]
    
    #now get all timestep within an experiment
    timestep<-sort(unique(subset2$calendar_year))    
    
    for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
      subset_t1<-subset2%>%
        filter(calendar_year==timestep[i])
      
      subset_t2<-subset2%>%
        filter(calendar_year==timestep[i+1])
      
      subset_t12<-merge(subset_t1, subset_t2, by=c("genus_species","id"), all=T)%>%
        filter(relcov.x!=0|relcov.y!=0)
      
      MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
      
      metrics<-data.frame(id=id, calendar_year=timestep[i+1], MRSc=MRSc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering=rbind(metrics, reordering)  
    }
  }
}

reorder_raw<-reordering%>%
  separate(id, c("site_project_comm","plot_id", "treatment"), sep="::")

summary(aov(MRSc ~ treatment*calendar_year + Error(plot_id/calendar_year), data=subset(reorder_raw, site_project_comm=="KNZ_pplots_0")))

reorder_means<-reorder_raw%>%
  group_by(site_project_comm, calendar_year, treatment)%>%
  summarise(MRSc=mean(MRSc))

ggplot(data=subset(reorder_means,site_project_comm=="KNZ_pplots_0"), aes(x=calendar_year, y=MRSc, group=treatment))+
  geom_point(size=3, aes(color=treatment))+
  scale_color_manual(values=c("green","blue","blue","blue","red","purple","purple","purple"))+
  geom_line()

ggplot(data=subset(reorder_means,site_project_comm=="SEV_WENNDEx_0"), aes(x=calendar_year, y=MRSc, group=treatment))+
  geom_point(size=3, aes(color=treatment))+
  scale_color_manual(values=c("green","red","blue","purple","black","red","blue","purple"))+
  geom_line()

ggplot(data=subset(reorder_means,site_project_comm=="NIN_herbdiv_0"), aes(x=calendar_year, y=MRSc, group=treatment))+
  geom_point(size=3, aes(color=treatment))+
  scale_color_manual(values=c("purple","green","purple","black","purple","black","purple","black", "purple","black"))+
  geom_line()

####comparing control versus treatment plots

###need a contorl and treatment column
corre2<-corre%>%
  mutate(trt=ifelse(site_code=="KNZ"&treatment=="N1P0", "C", ifelse(site_code=="NIN"&treatment=="1NF","C", ifelse(site_code=="SEV"&treatment=="C","C","T"))))

reordering_ct=data.frame(site_project_comm=c(), treatment=c(), calendar_year=c(), MRSc_diff=c(), spdiffc=c())

explist<-unique(corre2$site_project_comm)

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corre2%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  spc<-explist[i]
  
  ##make wide and get averages of each species by treatment
  wide<-subset%>%
    gather(genus_species, relcov,12:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, trt,genus_species)%>%
    summarize(relcov=mean(relcov))
  
  ##add ranks dropping zeros
  rank_pres<-wide%>%
    filter(relcov!=0)%>%
    tbl_df()%>%
    group_by(site_project_comm, calendar_year, treatment, trt)%>%
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-wide%>%
    filter(relcov==0)
  ##get species richness for each year
  rich<-group_by(wide, site_project_comm, calendar_year, treatment)%>%
    summarize(S=S(relcov))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year", "treatment"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  timestep<-sort(unique(rank$calendar_year)) 
  for(i in 1:(length(timestep))){
    
    time<-rank%>%
    filter(calendar_year==timestep[i])
  
    time_id<-timestep[i]
    
  #fitler out control plots
  control<-time%>%
    filter(trt=="C")
  
  treat_list<-unique(subset(time, trt!="C")$treatment)
  
  for (i in 1:length(treat_list)){
    treat<-time%>%
      filter(treatment==treat_list[i])
    
    treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("site_project_comm", "calendar_year","genus_species"), all=T)%>%
        filter(relcov.x!=0|relcov.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      spdiff<-subset_ct%>%
        filter(relcov.x==0|relcov.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      metrics<-data.frame(site_project_comm=spc, treatment=treat_id, calendar_year=time_id, MRSc_diff=MRSc_diff, spdiffc=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering_ct=rbind(metrics, reordering_ct)  
    }
  }
}

reorder_ct_raw<-reordering_ct

ggplot(subset(reorder_ct_raw, site_project_comm=="KNZ_pplots_0"), aes(x=calendar_year, y=MRSc_diff))+
  geom_point(aes(color=treatment))+
  geom_line(aes(group=treatment))+
  scale_color_manual(values=c("blue","blue","blue", "red","purple","purple","purple"))

ggplot(subset(reorder_ct_raw, site_project_comm=="KNZ_pplots_0"), aes(x=calendar_year, y=spdiffc))+
  geom_point(aes(color=treatment))+
  geom_line(aes(group=treatment))+
  scale_color_manual(values=c("blue","blue","blue", "red","purple","purple","purple"))

