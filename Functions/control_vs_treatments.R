
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




##### DO THIS AT THE TOP OF THE PRIMARY FUNCTION #####
## create a dataframe of all unique treatment combinations
namesvector = unique(df[[treatment.var]])

myperms <- as.data.frame(cbind(cola=as.character(), colb = as.character()))
for (i in 1:length(namesvector)) {
  cola <- as.character(namesvector[i])
  for (j in i:length(namesvector)) {
    if(i != j){
      colb <- as.character(namesvector[j])
      suba <- cbind(cola, colb)
      myperms <- rbind(suba, myperms)
    }
  }
}
names(myperms) <- c(treatment.var, paste(treatment.var, "2", sep = ""))

## A more parsimonious utils way to do this,
## if we go this way need to shift to a df 
# myperms <- combn(unique(time$treatment),2)

## Average values within each treatment, species and year
## Then take ranks
## Notes for Meghan:
## 1) will we ever want to make time.var optional? People might often just be comparing treatments...
## 2) double check that 0s have been filled in before the summing
myformula = as.formula(paste(abundance.var, "~", treatment.var, "+", species.var, "+", time.var))
sumdf <- aggregate(myformula, data = df, "mean")

## need to add treatment into this!
rankdf <- add_ranks(sumdf, time.var, species.var, abundance.var, treatment.var)

rankdf2 <- rankdf
rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
rankdf2[[treatment.var]] <- NULL

rankdfall <- merge(rankdf, myperms, all.y = T)

#rankdf2all <- merge(rankdf2, myperms, all.y = T)

ranktog <- merge(rankdfall, rankdf2, by=c(time.var, species.var, paste(treatment.var, "2", sep="")))
ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="_")
ranktog <- subset(ranktog, ranktog[[paste(abundance.var, ".x", sep = "")]]!=0 | ranktog[[paste(abundance.var, ".y", sep = "")]]!=0)
ranktog <- subset(ranktog, !is.na(ranktog[[paste(abundance.var, ".x", sep = "")]]) & !is.na(ranktog[[paste(abundance.var, ".y", sep = "")]]))


X <- split(ranktog, ranktog$splitvariable)

out <- lapply(X, FUN=aggfunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
ID <- unique(names(out))
out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
              out, ID, SIMPLIFY = FALSE)
output <- do.call("rbind", out)  

outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
names(outnames) = c(time.var, treatment.var, paste(treatment.var, "2", sep=""))
output$splitvariable <- NULL
output <- cbind(outnames, output)



subrank <- subset(ranktog,  treatment == "N1P3" & treatment2 == "N2P2")
subrank <- subset(subrank, subrank[[paste(abundance.var, ".x", sep = "")]]!=0 | subrank[[paste(abundance.var, ".y", sep = "")]]!=0)
subrank <- subset(subrank, !is.na(subrank[[paste(abundance.var, ".x", sep = "")]]) & !is.na(subrank[[paste(abundance.var, ".y", sep = "")]]))


#### ## SPLIT BY YEAR AND DO EACH YEAR ####


SERSp=data.frame(treatment=c(), time=c(), Sd=c(), Ed=c(), Rd=c(), spd=c())      

for (i in 1:length(myperms)){
  
df1 <- subset(rankdf, rankdf[[treatment.var]]== as.character(myperms[[i,1]]))
df2 <- subset(rankdf, rankdf[[treatment.var]]== as.character(myperms[[i,2]]))

df12 <- merge(df1, df2, by=c(time.var,species.var), all=T)
df12 <- subset(df12, df12[[paste(abundance.var, ".x", sep = "")]]!=0 | df12[[paste(abundance.var, ".y", sep = "")]]!=0)
df12 <- subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))

out <- aggfunc(df12, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = ""))



# out <- lapply(X, FUN=aggfunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
# ID <- unique(names(out))
# out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
#               out, ID, SIMPLIFY = FALSE)
# output <- do.call("rbind", out) 




MRSc_diff <- mean(abs(df12$rank.x-df12$rank.y))/nrow(df12)

## Meghan: I think there's a problem up the chain because there aren't any 0s
spdiff <- subset(df12, df12[[paste(abundance.var, ".x", sep = "")]] == 0 | df12[[paste(abundance.var, ".y", sep = "")]] == 0)

spdiffc <- nrow(spdiff)/nrow(df12)

##eveness richness
s_c <- S(subset_ct$abundance.x)
e_c <- EQ(subset_ct$abundance.x)
s_t <- S(subset_ct$abundance.y)
e_t <- EQ(subset_ct$abundance.y)

sdiff<-abs(s_c-s_t)/nrow(df12)
ediff<-abs(e_c-e_t)/nrow(df12)

metrics<-data.frame(treatment=treat_id, time=time_id, Sd=sdiff, Ed=ediff, Rd=MRSc_diff, spd=spdiffc)#spc_id

##calculate differences for these year comparison and rbind to what I want.

SERSp=rbind(metrics, SERSp)  

}



  