
########comparing control versus treatment plots

#' @param sim foo dataset with columns for time, replicate, species, abundance, a treatment column for grouping replicates, and a column speficfying if the treatment is a control or a treatment
add_ranks_treatment_control_sppools <- function(df) {
  
  ###add zeros and average up species pool for control and treatment plots
  wide<-df%>%
    spread(species, abundance, fill=0)
  
  ##make long and get averages of each species by treatment
  long<-wide%>%
    gather(species, abundance, 5:ncol(wide))%>%
    group_by(time, treatment, species, C_T)%>%
    summarize(abundance=mean(abundance))
  
  ##add ranks dropping zeros
  rank_pres<-long%>%
    filter(abundance!=0)%>%
    tbl_df()%>%
    group_by(time, treatment, C_T)%>%
    mutate(rank=rank(-abundance, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-long%>%
    filter(abundance==0)
  ##get species richness for each year
  rich<-group_by(long, time, treatment, C_T)%>%
    summarize(S=S(abundance))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("time", "treatment","C_T"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  ct_rank<-rbind(rank_pres, zero_rank)
  return(ct_rank)
}

calculate_SERSp <- function(ct_rank){
  SERSp=data.frame(treatment=c(), time=c(), Sd=c(), Ed=c(), Rd=c(), spd=c())      
  timestep<-sort(unique(ct_rank$time)) 
  
  for(i in 1:(length(timestep))){
    
    time<-ct_rank%>%
      filter(time==timestep[i])
    
    time_id<-timestep[i]
    
    #need to do all comparisions NOT SURE HOW TO PROCEED FROM HERE TO SUBSET ALL POSSIBLE COMBINATIONS
    comparison<-combn(unique(time$treatment),2)
    
    
    #filter out control plots
    control<-time%>%
      filter(C_T=="Control")
    
    treat_list<-unique(subset(time, C_T=="Treatment")$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("time","species"), all=T)%>%
        filter(abundance.x!=0|abundance.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      
      spdiff<-subset_ct%>%
        filter(abundance.x==0|abundance.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      ##eveness richness
      s_c <- S(subset_ct$abundance.x)
      e_c <- E_q(subset_ct$abundance.x)
      s_t <- S(subset_ct$abundance.y)
      e_t <- E_q(subset_ct$abundance.y)
      
      sdiff<-abs(s_c-s_t)/nrow(subset_ct)
      ediff<-abs(e_c-e_t)/nrow(subset_ct)
      
      metrics<-data.frame(treatment=treat_id, time=time_id, Sd=sdiff, Ed=ediff, Rd=MRSc_diff, spd=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      SERSp=rbind(metrics, SERSp)  
    }
  }
  return(SERSp)
}
