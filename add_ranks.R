## a function to fill zeros
fill_zeros <- function(df, time.var, species.var, abundance.var){
  df2<-subset(df, select=c(time.var,species.var,abundance.var))
  wide<- reshape(df2, idvar=time.var, timevar=species.var, direction="wide")
  wide[is.na(wide)] <- 0
  
  long<-reshape(wide, idvar=time.var, ids=time.var, time=names(wide), timevar=abundance.var, direction="long")
  colnames(long)[3]<-abundance.var
  return(long)
}


##FUNCTION TO ADD RANKS
#' @param for a dataset with columns for time, replicate, species, and abundance, (and optionally an id column for grouping and second column for defining the groupings)
#' 

add_ranks <- function(df, replicate.var, species.var, abundance.var, time.var) {
  ##SHOULD THE FIRST STEP TO BE ONLY SELECT THE RELEVANT COLUMNS?
  df<-subset(df, select=c(time.var, abundance.var, species.var, replicate.var))
  
  ##add ranks for present species
  rank_pres<-subset(df, df[[abundance.var]]!=0)
  rank_pres$rep_time<-paste(rank_pres[[replicate.var]], rank_pres[[time.var]], sep="_")
  rank_pres$rank<-ave(rank_pres[[abundance.var]], rank_pres$rep_time, FUN=function(x) rank(-x, ties.method = "average"))
  rank_pres<-subset(rank_pres, select=-rep_time)
  
  #adding zeros
  
  # sort and apply fill_zeros to all replicates
  df <- df[order(df[[replicate.var]]),]
  X <- split(df, df[replicate.var])
  out <- lapply(X, FUN=fill_zeros, time.var, species.var, abundance.var)
  ID <- unique(names(out))
  out <- mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                out, ID, SIMPLIFY = FALSE)
  allsp <- do.call("rbind", out)

 
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros <- subset(allsp, allsp[[abundance.var]]==0)
  
  ##get species richness for each year
  ## Note to Meghan: This uses a function that I thought was only for community_structure
  ## Might make a separate file of shared diversity functions
  myformula <- as.formula(paste(abundance.var, "~", replicate.var, "+", time.var))
  SpR<-aggregate(myformula, FUN=S, data=allsp)
  colnames(SpR)[3]<-"S"
  
  ##merge together make zero abundances rank S+1
  #HOW TO DEAL WITH COLUMN NAMES IN MERGE?? - toes it work to not have in quotes?
  zero_rank<-merge(zeros, SpR, by=c(time.var,replicate.var))
  zero_rank$rank<-zero_rank$S+1
  zero_rank<-subset(zero_rank, select=-S)
  
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  return(rank)
}
