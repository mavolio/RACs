#####CALCULATING DIVERSITY METRICS
##still need to make a funciton to encapsulate all these
#' @title Community structure
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 

##make for diversity
##rename columns

community_structure <- function(df, replicate.var, abundance.var, time.var, diversity="shannons") {
  if(is.null(replicate.var)){
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    
    if(diversity=="shannons"){
      comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),diversity=shannons(x)))
      names(comstruct)[3]<-"shannons"
      
    }
    else{
      comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),diversity=simpsons(x)))
      names(comstruct)[3]<-"simpsons"

    }
  } 
  else {
    
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
    
    if(evenness=="shannons"){
      comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),diversity=shannons(x)))
      names(comstruct)[3]<-"shannons"
    } 
    else{
      comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),diversity=simpsons(x)))
      names(comstruct)[3]<-"simpsons"
    }
  }
  return(comstruct)
}

#Problems
#only outputs richness
# does not work with replicate missing

##test
test<- community_structure(df, replicate.var = "replicate", abundance.var = "abundance", time.var = "time")
test2<-community_structure(df2, abundance.var = "abundance", time.var = "time")

#### PRIVATE FUNCTIONS ####


#1) function to calculate Simpson's Divsersity from Smith and Wilson 1996
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
simpsons<-function(x, S=length(x[x!=0 & !is.na(x)]), N=sum(x[x!=0&!is.na(x)]), ps=x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D<-sum(p2)
  (1/D)/S
}

#2) function to calculate Shannon's Divsersity 
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
Shannons<-function(x, S=length(x[x!=0 & !is.na(x)]), N=sum(x[x!=0&!is.na(x)]), ps=x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D<-sum(p2)
  (1/D)/S
}



