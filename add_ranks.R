
##FUNCTION TO ADD RANKS
#' @param for a dataset with columns for time, replicate, species, and abundance, (and optionally an id column for grouping and second column for defining the groupings)
#' 

add_ranks <- function(df, replicate.var, species.var, abundance.var, time.var) {
  ##SHOULD THE FIRST STEP TO BE ONLY SELECT THE RELEVANT COLUMNS?
  df<-subset(df, select=c("time", "abundance","species","replicate"))
  
  ##add ranks for present species
  rank_pres<-subset(df, df[[abundance.var]]!=0)
  rank_pres$rep_time<-paste(rank_pres[[replicate.var]], rank_pres[[time.var]], sep="_")
  rank_pres$rank<-ave(rank_pres[[abundance.var]], rank_pres$rep_time, FUN=function(x) rank(-x, ties.method = "average"))
  rank_pres<-subset(rank_pres, select=-rep_time)
  
  #adding zeros
  #NEED TO REMOVE THE LOOP HERE
  replist<-unique(df[[replicate.var]])
  
  allsp<-data.frame()
  for (i in 1:length(replist)){
    subset <- subset(df[[replicate.var]]==replist[i])
    
    replicate<-replist[i]
    
    #FOR RESHAPE HOW TO DEAL WITH COLUMN NAMES? ###
    subset2<-subset(subset, select=c("time","species","abundance"))
    wide<- reshape(subset2, idvar="time", timevar="species", direction="wide")
    wide[is.na(wide)] <- 0
    
    long<-reshape(wide, idvar="time", ids="time", time=names(wide), timevar="abundance", direction="long")
    long$replicate<-replist[i]
    colnames(long)[3]<-"abundance"
    
    allsp<-rbind(long, allsp)  
  }
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-subset(allsp, all.sp[[abundance.var]]==0)
  
  ##get species richness for each year
  SpR<-aggregate([[abundance.var]]~[[replicate.var]]+[[time.var]], FUN=S, data=allsp)
  colnames(SpR)[4]<-"S"
  
  ##merge together make zero abundances rank S+1
  #HOW TO DEAL WITH COLUMN NAMES IN MERGE?? - toes it work to not have in quotes?
  zero_rank<-merge(zeros, SpR, by=c(time.var,replicate.var))
  zero_rank$rank<-zero_rank$S+1
  zero_rank<-subset(zero_rank, select=-S)
  
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  return(rank)
}
