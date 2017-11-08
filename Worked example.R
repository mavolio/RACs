##worked example with pplots

setwd("~/Dropbox/converge_diverge/datasets/Longform")

library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)

theme_set(theme_bw(20))

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
  

ggplot(data=subset(average, treatment=="N1P0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=1)+
  scale_color_manual(values=c("black","gray"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
    ylab("Cumulative Relative Abundance")+
    xlab("Relative Rank")
ggplot(data=subset(average, treatment=="N2P3"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=1)+
  scale_color_manual(values=c("black","gray"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")

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

##richness and evenness changes
diversity <- group_by(pplots, treatment, calendar_year, plot_id) %>% 
  summarize(S=S(relcov),
            E_q=E_q(relcov))%>%
  arrange(treatment, plot_id, calendar_year)%>%
  group_by(plot_id)%>%
  mutate(S_diff=c(NA,diff(S)),
         E_diff=c(NA,diff(E_q)))%>%
  na.omit

##gains and losses
loss<-turnover(df=pplots, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="disappearance")
gain<-turnover(df=pplots, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="appearance")

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

reordering=data.frame(plot_id=c(), calendar_year=c(), MRSc=c())#expeiment year is year of timestep2

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
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(plot_id=id, calendar_year=timestep[i+1], MRSc=MRSc)#spc_id
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

merge1<-merge(diversity, gain, by=c("plot_id", "calendar_year"))
merge2<-merge(merge1, loss, by=c("plot_id", "calendar_year"))
merge3<-merge(merge2, reordering, by=c("plot_id", "calendar_year"))
allmetrics<-merge(merge3, d_output, by=c("plot_id", "calendar_year"))%>%
  group_by(treatment)%>%
  summarize(sdiff=mean(S_diff),
            ediff=mean(E_diff),
            gain=mean(appearance),
            loss=mean(disappearance),
            reorder=mean(MRSc),
            cc=mean(Dstar),
            n=length(plot_id),
            sdsdiff=sd(S_diff),
            sdediff=sd(E_diff),
            sdgain=sd(appearance),
            sdloss=sd(disappearance),
            sdreorder=sd(MRSc),
            sdcc=sd(Dstar))%>%
  mutate(seS=sdsdiff/sqrt(n),
         seE=sdsdiff/sqrt(n),
         se=sdsdiff/sqrt(n),
         seS=sdsdiff/sqrt(n),
         seS=sdsdiff/sqrt(n),
         )
            


# theme_set(theme_bw(12))
# ###read in datasets where we already have this
# 
# vpplots<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\CORRE_RAC_Metrics_Oct2017_allyears_2.csv")%>%
#   filter(site_project_comm=="KNZ_pplots_0")%>%
#   filter(treatment=="N1P0"|treatment=="N2P3")%>%
#   gather(metric, value, S_diff:dispersion_diff)
# 
# ggplot(data=vpplots, aes(x=calendar_year, y=value, color=treatment))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~metric, ncol=3, scale="free")
# 
# hpplots<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\CORRE_ContTreat_Compare_OCT2017.csv")%>%
#   filter(site_project_comm=="KNZ_pplots_0")%>%
#   filter(treatment=="N2P3")%>%
#   select(-plot_mani)%>%
#   gather(metric, value, PCSdiff:disp_diff)
# 
# ggplot(data=hpplots, aes(x=calendar_year, y=value))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~metric, ncol=3, scale="free")
