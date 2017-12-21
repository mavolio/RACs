#####CALCULATING DIVERSITY METRICS
#' @title Community structure
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 

community_diversity <- function(df,  time.var, abundance.var, replicate.var=NULL,  diversity = "Shannon") {
  if(is.null(replicate.var)){
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    
    if(diversity == "Shannon"){
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
      names(comdiv)[2] <- "Shannon"
      
    }
    else{
      comdiv <- aggregate(myformula, data=df, FUN = function(x) diversity = Simpson(x))
      names(comdiv)[2] <- "Simpson"

    }
  } 
  else {
    
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
    
    if(diversity == "Shannon"){
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
      names(comdiv)[3] <- "Shannon"
    } 
    else{
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Simpson(x))
      names(comdiv)[3] <- "Simpson"
    }
  }
  return(comdiv)
}


#### PRIVATE FUNCTIONS ####


#1) function to calculate Simpson's Divsersity from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
Simpson <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D <- sum(p2)
  1/D
}

#2) function to calculate Shannon's Divsersity 
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
Shannon <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N ){
  -sum(ps*log(ps))
}




