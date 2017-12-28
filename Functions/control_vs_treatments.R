## Notes for Meghan:
## 1) will we ever want to make time.var optional? People might often just be comparing treatments... Same for replicate.var, don't have the optional ifelse in yet
## 2) double check that 0s have been filled in before the summing

##calculating RAC changes between treatments
#' @title Rank Abundance Curve Changes Between Treatments
#' @description 
#' @param df A data frame containing time, species and abundance, treatment and replicate columns 
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param treatment.var The name of the treatment column
#' @param replicate.var The name of the optional replicate column 
#' 
#' 
########comparing control versus treatment plots
SERGL_treatment <- function(df, time.var, species.var, abundance.var, treatment.var, replicate.var=NULL){

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

## A more parsimonious way to do this is in "utils"
## if we go this way we'd need to shift the output to a df to be consistent with below
# myperms <- (combn(unique(time$treatment,2))

## Average values within each treatment, species and year
myformula = as.formula(paste(abundance.var, "~", treatment.var, "+", species.var, "+", time.var))
sumdf <- aggregate(myformula, data = df, "mean")

## Take ranks
rankdf <- add_ranks_time(sumdf, time.var, species.var, abundance.var, treatment.var)

## Create a second rankdf with a renamed treatment.var column
rankdf2 <- rankdf
rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
rankdf2[[treatment.var]] <- NULL

## Merge rankdf with all possible permutations of treatment combinations
rankdfall <- merge(rankdf, myperms, all.y = T)

## Merge the data together (for all possible permutations of treatments)
ranktog <- merge(rankdfall, rankdf2, by=c(time.var, species.var, paste(treatment.var, "2", sep="")))

## Create a variable to split on (time as well as unique treatment combos)
ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="_")

## Remove 0s (Meghan - I copied this from SERGL, is this what we should do?)
ranktog <- subset(ranktog, ranktog[[paste(abundance.var, ".x", sep = "")]]!=0 | ranktog[[paste(abundance.var, ".y", sep = "")]]!=0)

## Remove instances where abundance.var is NA, probably not necessary since we are doing over space and not time but in as a check
ranktog <- subset(ranktog, !is.na(ranktog[[paste(abundance.var, ".x", sep = "")]]) & !is.na(ranktog[[paste(abundance.var, ".y", sep = "")]]))

## Split the dataframe
X <- split(ranktog, ranktog$splitvariable)

## Apply the basic RAC function
out <- lapply(X, FUN=aggfunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
ID <- unique(names(out))
out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
              out, ID, SIMPLIFY = FALSE)
output <- do.call("rbind", out)  

## Add in the identifying column names
outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
names(outnames) = c(time.var, treatment.var, paste(treatment.var, "2", sep=""))
output$splitvariable <- NULL
output <- cbind(outnames, output)

return(output)
}
