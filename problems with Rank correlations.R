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


#trouble shooting E_Q
test1 <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(S=S(abundance),
            E_q=E_q(abundance),
            Gini=Gini(abundance),
            E_simp=E_simp(abundance))

##when there is only 1 species get NA
##when a community is perfectly even, also get NA, that is a problem EX: SGS_UNUN_5A_15!!

x<-c(10, 10, 10)
z<-4


#first if statment says return a 1 if an even community.
  
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
E_q(x)

test<-codyndat%>%
  filter(plot_id=="PUPL5", abundance!=0, experiment_year==1993)


y<-test$abundance
E_q(y)

