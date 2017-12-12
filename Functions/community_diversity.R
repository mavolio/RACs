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

community_diversity <- function(df, replicate.var, abundance.var, time.var, diversity="shannons") {
  if(is.null(replicate.var)){
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    
    if(diversity=="shannons"){
      comdiv <- aggregate(myformula, data=df, FUN=function(x)diversity=shannons(x))
      names(comdiv)[3]<-"shannons"
      
    }
    else{
      comdiv <- aggregate(myformula, data=df, FUN=function(x)diversity=simpsons(x))
      names(comdiv)[3]<-"simpsons"

    }
  } 
  else {
    
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
    
    if(diversity=="shannons"){
      comdiv <- aggregate(myformula, data=df, FUN=function(x)diversity=shannons(x))
      names(comdiv)[3]<-"shannons"
    } 
    else{
      comdiv <- aggregate(myformula, data=df, FUN=function(x)diversity=simpsons(x))
      names(comdiv)[3]<-"simpsons"
    }
  }
  return(comdiv)
}

#Problems
# does not work with replicate missing

##test
test<- community_diversity(df, replicate.var = "replicate", abundance.var = "abundance", time.var = "time")
test3<- community_diversity(df, replicate.var = "replicate", abundance.var = "abundance", time.var = "time", diversity = "simpson")
test2<-community_structure(df2, abundance.var = "abundance", time.var = "time")

#### PRIVATE FUNCTIONS ####


#1) function to calculate Simpson's Divsersity from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
simpsons<-function(x, N=sum(x[x!=0&!is.na(x)]), ps=x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D<-sum(p2)
  1/D
}

#2) function to calculate Shannon's Divsersity 
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
shannons<-function(x, N=sum(x[x!=0&!is.na(x)]), ps=x[x!=0&!is.na(x)]/N ){
  -sum(ps*log(ps))
}




