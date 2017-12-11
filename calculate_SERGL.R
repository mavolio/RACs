
##calculating SERGL
# @param sim foo dataset with columns for time, plot, species, and abundance
calculate_SERGL <- function(df, replicate.var, species.var, abundance.var, time.var) {
  rankdf <- add_ranks(df, replicate.var, species.var, abundance.var, time.var)
  SERGL=data.frame(replicate=c(), time=c(), S=c(), E=c(), R=c(), G=c(), L=c())#expeiment year is year of timestep2
  
  replist<-unique(rankdf[[replicate.var]])
  
  for (i in 1:length(replist)){
    #THIS BREAKS AT THIS STEP - see if single bracket fixes this.
    subber <- subset(rankdf, rankdf[replicate.var]==replist[i])
    
    replicate<-replist[i]
    
    #now get all timestep within an experiment
    timestep<-sort(unique(subber[[time.var]]))    
    
    #NEED TO REMOVE THE LOOP
    for(i in 1:(length(timestep)-1)) {
      subset_t1<-subset(subset, subset[[time.var]]==timestep[i])
      
      subset_t2<-subset(subset, subset[[time.var]]==timestep[i+1])
      
      #HOW TO MERGE BY COLUMN NAMES?
      subset_t12<-merge(subset_t1, subset_t2, by=c("species","replicate"), all=T)
      subset_t12<-subset(subset_t12, abundance.x!=0|abundance.y!=0)
      
      #reordering
      MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
      #ricness and evenness differences
      s_t1 <- S(subset_t12$abundance.x)
      e_t1 <- E_q(as.numeric(subset_t12$abundance.x))
      s_t2 <- S(subset_t12$abundance.y)
      e_t2 <- E_q(as.numeric(subset_t12$abundance.y))
      
      sdiff<-abs(s_t1-s_t2)/nrow(subset_t12)
      ediff<-abs(e_t1-e_t2)/nrow(subset_t12)
      
      #gains and losses
      subset_t12$gain<-ifelse(subset_t12$abundance.x==0, 1, 0)
      subset_t12$loss<-ifelse(subset_t12$abundance.y==0, 1, 0)
      
      gain<-sum(subset_t12$gain)/nrow(subset_t12)
      loss<-sum(subset_t12$loss)/nrow(subset_t12)
      
      metrics<-data.frame(replicate=replicate, time=timestep[i+1], S=sdiff, E=ediff, R=MRSc, G=gain, L=loss)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      SERGL=rbind(metrics, SERGL)  
    }
  }
  return(SERGL)
}

SERGL_func <- function(df) {
  rank <- add_ranks(df)
  result <- calculate_SERGL(rank)
  return(result)
}