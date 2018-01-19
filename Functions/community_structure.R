
#####CALCULATING DIVERSITY METRICS
#' @title Community structure
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 


community_structure <- function(df,  time.var=NULL, 
                                abundance.var, replicate.var = NULL, 
                                evenness = "EQ") {
  
  if(is.null(replicate.var)){
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    
    if(evenness == "EQ"){
    comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = EQ(x))))
    names(comstruct)[2] <- "richness"
    names(comstruct)[3] <- "evenness_EQ"
    
    } else {
      
    comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = SimpEven(x))))
    names(comstruct)[2] <- "richness"
    names(comstruct)[3] <- "evenness_Simpson"
    
    } 
  } else {
    
    if(is.null(time.var)){
      myformula <- as.formula(paste(abundance.var, "~", replicate.var))
      
      if(evenness == "EQ"){
        comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = EQ(x))))
        names(comstruct)[2] <- "richness"
        names(comstruct)[3] <- "evenness_EQ"
        
      } else {
        comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = SimpEven(x))))
        names(comstruct)[2] <- "richness"
        names(comstruct)[3] <- "evenness_Simpson"
      }
      
    } else {
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
  
    if(evenness == "EQ"){
    comstruct <- do.call(data.frame, aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = EQ(x))))
    names(comstruct)[3] <- "richness"
    names(comstruct)[4] <- "evenness_EQ"
    
  } else {
    comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = SimpEven(x))))
    names(comstruct)[3] <- "richness"
    names(comstruct)[4] <- "evenness_Simpson"
  }
    }
  }
  
    return(comstruct)
}

if(is.null(time.var)) {
  myformula <- as.formula(paste(abundance.var, "~", replicate.var))
  
  if(evenness == "EQ"){
    comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = EQ(x))))
    names(comstruct)[2] <- "richness"
    names(comstruct)[3] <- "evenness_EQ"
    
  } else {
    comstruct <- do.call(data.frame,aggregate(myformula, data = df, FUN = function(x) c(SpR = S(x), evenness = SimpEven(x))))
    names(comstruct)[2] <- "richness"
    names(comstruct)[3] <- "evenness_Simpson"
  }
} 

#### PRIVATE FUNCTIONS ####


#'  A function to calculate richness
#' @param x the vector of abundances of each species
S <- function(x){
  x1 <- x[x!=0 & !is.na(x)]
  stopifnot(x1 == as.numeric(x1))
  length(x1)
  }

#'  A function to calculate EQ evenness from Smith and Wilson 1996
#' @param x the vector of abundances of each species
#' if all abundances are equal it returns a 1
EQ <- function(x){
  x1 <- x[x!=0 & !is.na(x)]
  if (length(x1) == 1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r <- rank(x1, ties.method = "average")
  r_scale <- r/max(r)
  x_log <- log(x1)
  fit <- lm(r_scale~x_log)
  b <- fit$coefficients[[2]]
  2/pi*atan(b)
}


#'  A function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @param S the number of species in the sample
#' @param x the vector of abundances of each species
#' @param N the total abundance
#' @param p the vector of relative abundances of each species
SimpEven <- function(x, S = length(x[x!=0 & !is.na(x)]), N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N, p2 = ps*ps ){
  D <- sum(p2)
  (1/D)/S
}



