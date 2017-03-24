library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)

###this shows that curves are more similar to one another than we would expect if we randomly shuffle the values. this if for even curves that look very different to me. if using raw abundance data do not need to relativize, if using relcov data need to re-relativize the data in the loop. This is why there are code differences.

###get corre dataset to work on and focus on familiar examples
corre<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Dec2016.csv")%>%
  filter(project_name=="pplots"|project_name=="herbdiv"|project_name=="WENNDEx")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type,sep="_"))%>%
  select(-X)%>%
  mutate(id=paste(site_project_comm, plot_id, treatment, sep="::"))


##creating code to compare two curves following the general methods outlined in the ROAP paper.

#1. get two curves to compare.
pplots<-corre%>%
  filter(project_name=="pplots" & treatment=="N1P0"|treatment=="N2P3")%>%
  filter(calendar_year==2014)%>%
  filter(plot_id==25|plot_id==13)%>%
  group_by(plot_id)%>%
  mutate(rank=rank(-relcov,ties.method="average"))%>%
  arrange(rank)

#2. just plot these to compare - at quick glance, I think these are very different!
ggplot(data=pplots, aes(x=rank, y = relcov, group=plot_id))+
  geom_point(aes(color=genus_species))+
  geom_line(aes(linetype=treatment))


# #3. figure out the actual difference between the two curve
real_rank<-pplots%>%
  group_by(plot_id)%>%
  mutate(rank=rank(-relcov, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(-relcov)%>%
  mutate(cumabund=cumsum(relcov))%>%
  arrange(rank)%>%
  ungroup()

ggplot(data=real_rank, aes(x=relrank, y = cumabund, group=plot_id))+
  geom_point()+
  geom_step()

result_actual <- real_rank %>%
  do({
    y <- unique(.$plot_id)###assumption this is a length 2 list
    df1 <- filter(., plot_id==y[[1]])
    df2 <- filter(., plot_id==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

#4. set up bootstrapping 

#get blank dataset to fill in
blank<-pplots%>%
  select(plot_id)%>%
  mutate(relcov=NA)

cc_output=data.frame(Darea=c())#

for (i in 1:1000){
  blank$relcov<-sample(pplots$relcov, 22, replace = F)#need to re-relative to 1
  
  blank_totcov<-blank%>%
    group_by(plot_id)%>%
    summarize(totcov=sum(relcov))
  
  blank_relcov<-merge(blank_totcov, blank, by="plot_id")%>%
    mutate(relcov2=relcov/totcov)#%>%
    #group_by(plot_id)%>%
    #summarize(totcov=sum(relcov2))
  
  ranks<-blank_relcov%>%
    filter(relcov2!=0)%>%
    group_by(plot_id)%>%
    mutate(rank=rank(-relcov2, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-relcov2)%>%
    mutate(cumabund=cumsum(relcov2))%>%
    arrange(rank)%>%
    ungroup()
  
  ggplot(data=ranks, aes(x=relrank, y = cumabund, group=plot_id))+
    geom_point()+
    geom_step()
  
  #this does not have to be looped.
    result <- ranks %>%
        do({
        y <- unique(.$plot_id)###assumption this is a length 2 list
        df1 <- filter(., plot_id==y[[1]])
        df2 <- filter(., plot_id==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(Darea=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(Darea=result$Darea)
    
    cc_output<-rbind(cc_output, d_output1)
}

hist(cc_output$Darea)

#calculating p-value

cc_output_stat<-cc_output%>%
  mutate(count=ifelse(Darea>=0.1113217, 0, 1))

sum(cc_output_stat$count)
########
#######
########
######okay doing this on real pplot data

#1. get two curves to compare
pplots<-read.csv("~/Dropbox/pplots/Compiled Datasets/public data/KNZ_PPLOTS_Data_2002_2014_spp.csv")%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  filter(calendar_year==2014)%>%
  filter(plot_id==25|plot_id==27)%>%
  filter(abundance!=0)%>%
  group_by(plot_id)%>%
  mutate(rank=rank(-abundance,ties.method="average"))%>%
  arrange(rank)

##creating code to compare two curves following the general methods outlined in the ROAP paper.

#2. just plot these to compare - at quick glance, I think these are very different!
ggplot(data=pplots, aes(x=rank, y = abundance, group=plot_id))+
  geom_point(aes(color=genus_species))+
  geom_line(aes(linetype=treatment))


#3. figure out the actual difference between the two curve
real_rank<-pplots%>%
  group_by(plot_id)%>%
  mutate(rank=rank(-abundance, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(-abundance)%>%
  mutate(cumabund=cumsum(abundance))%>%
  arrange(rank)%>%
  ungroup()

ggplot(data=real_rank, aes(x=relrank, y = cumabund, group=plot_id))+
  geom_point()+
  geom_step()

result_actual <- real_rank %>%
  do({
    y <- unique(.$plot_id)###assumption this is a length 2 list
    df1 <- filter(., plot_id==y[[1]])
    df2 <- filter(., plot_id==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

#4. set up bootstrapping 

#get blank dataset to fill in
blank<-pplots%>%
  select(plot_id)%>%
  mutate(relcov=NA)

cc_output=data.frame(Darea=c())#

for (i in 1:1000){
    blank$relcov<-sample(pplots$abundance, 27, replace = F)#need to re-relative to 1

  ranks<-blank%>%
    filter(relcov!=0)%>%
    group_by(plot_id)%>%
    mutate(rank=rank(-relcov, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-relcov)%>%
    mutate(cumabund=cumsum(relcov))%>%
    arrange(rank)%>%
    ungroup()
  
  ggplot(data=ranks, aes(x=relrank, y = cumabund, group=plot_id))+
    geom_point()+
    geom_step()
  
  #this does not have to be looped.
  result <- ranks %>%
    do({
      y <- unique(.$plot_id)###assumption this is a length 2 list
      df1 <- filter(., plot_id==y[[1]])
      df2 <- filter(., plot_id==y[[2]])
      sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
      sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
      r <- sort(unique(c(0, df1$relrank, df2$relrank)))
      h <- abs(sf1(r) - sf2(r))
      w <- c(diff(r), 0)
      data.frame(Darea=sum(w*h))#do has to output a dataframe
    })
  
  d_output1=data.frame(Darea=result$Darea)
  
  cc_output<-rbind(cc_output, d_output1)
}

ggplot(data=cc_output, aes(x=Darea))+
  geom_histogram()+
  geom_vline(xintercept = 8.62)


