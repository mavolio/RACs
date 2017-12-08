
#####CALCULATING DIVERSITY METRICS
##still need to make a funciton to encapsulate all these
#' @title Community structure
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 


community_structure <- function(df, replicate.var, abundance.var, time.var, evenness="E_q") {
  if(is.null(replicate.var)){
    
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),evenness=E_q(x)))
  
    } else {
      
  myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
  comstruct <- aggregate(myformula, data=df, FUN=function(x)c(SpR=S(x),evenness=E_q(x)))
    }
  return(comstruct)
}




#1) function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}


S<-function(df, species.var){
  
  x1<-x[x!=0]
  length(x1)
}


# 2) function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
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


#function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
E_simp<-function(x, S=length(x[x!=0]), N=sum(x[x!=0]), ps=x[x!=0]/N, p2=ps*ps ){
  D<-sum(p2)
  (1/D)/S
}



