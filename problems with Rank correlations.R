###get codyndat_rank in the measuremnts and analysis file

###trouble shooting
sgs<-codyndat_rank%>%
  filter(site_project_comm=="SGS_UNUN_0")

reordering1=data.frame(id=c(), experiment_year=c(), MRS=c(), RC=c())#expeiment year is year of timestep2

spc_id<-unique(sgs$id)

for (i in 1:length(spc_id)){
  subset<-sgs%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$experiment_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(experiment_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(experiment_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRS<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    RC<-cor(subset_t12$rank.x, subset_t12$rank.y, method="kendall")
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRS=MRS, RC=RC)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering1=rbind(metrics, reordering1)  
  }
}




###further trouble shooting
#Here is an example where there is a rank shift but an r of 0? Not sure why.
subset_t1<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_7_2"&experiment_year==1998)

subset_t2<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_7_2"&experiment_year==1999)

subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
  filter(abundance.x!=0|abundance.y!=0)

MRS<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)

RC<-cor(subset_t12$rank.x, subset_t12$rank.y, method="kendall")


#Here is an example where there is a rank shift but an r of NA not sure why
##give error of standard deviation is zero because one of the inputs have all the same rank
subset_t1<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_5A_15"&experiment_year==2007)

subset_t2<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_5A_15"&experiment_year==2008)

subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
  filter(abundance.x!=0|abundance.y!=0)

MRS<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)

RC<-cor(subset_t12$rank.x, subset_t12$rank.y, method="kendall")


#Here is an example where there is no rank shift and an r of NA not sure why
##here there is no error, just only 1 species so can't correlate with 1 datapoint
subset_t1<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_5A_4"&experiment_year==1997)

subset_t2<-sgs%>%
  filter(id=="SGS_UNUN_0::SGS_UNUN_5A_4"&experiment_year==1998)

subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
  filter(abundance.x!=0|abundance.y!=0)

MRS<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)

RC<-cor(subset_t12$rank.x, subset_t12$rank.y, method="kendall")


