##calculating RAC changes
#' @title Rank Abundance Curve Changes
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' 
#' 
#' TO DO: Add time.var_pair, add for null replicates?

RAC_changes <- function(df, time.var, species.var, abundance.var, replicate.var=NULL) {
  if(is.null(replicate.var)){
  
  rankdf <- add_ranks(df, time.var, species.var, abundance.var, replicate.var)
  
  # current year rankdf
  df2 <- rankdf
  
  # previous year rank df
  df1 <- rankdf
  df1[[time.var]] <- df1[[time.var]] + 1
  
  # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
  df12 <- merge(df1, df2,  by=c(species.var, time.var), all=T)
  df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]]!=0|df12[[paste(abundance.var, ".y", sep = "")]]!=0)
  df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
  
  # sort and apply turnover to all replicates
  df12 <- df12[order(df12[[time.var]]),]
  X <- split(df12, df12[[time.var]])
  

  out <- lapply(X, FUN=aggfunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
  ID <- unique(names(out))
  out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                out, ID, SIMPLIFY = FALSE)
  output <- do.call("rbind", out)  

  outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
  names(outnames) = time.var
  
  output$splitvariable <- NULL
  output <- cbind(outnames, output)
  }
  else{
    
    rankdf <- add_ranks(df,  time.var, species.var, abundance.var, replicate.var)
    
    # current year rankdf
    df2 <- rankdf
    
    # previous year rank df
    df1 <- rankdf
    df1[[time.var]] <- df1[[time.var]] + 1
    
    # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
    df12 <- merge(df1, df2,  by=c(species.var,replicate.var, time.var), all=T)
    df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]]!=0|df12[[paste(abundance.var, ".y", sep = "")]]!=0)
    df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
    df12$splitvariable <- paste(df12[[replicate.var]], df12[[time.var]], sep="_") 
    
    # sort and apply turnover to all replicates
    df12 <- df12[order(df12$splitvariable),]
    X <- split(df12, df12$splitvariable)
    
    
    out <- lapply(X, FUN=aggfunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    output <- do.call("rbind", out)  
    
    outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
    names(outnames) = c(replicate.var, time.var)
    
    output$splitvariable <- NULL
    output <- cbind(outnames, output)
    
  }
  return(output)
}

test2<-RAC_changes(pdata2, time.var="time",species.var = "species", abundance.var = "abundance")#not working
test1<-RAC_changes(pdata, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")

### PRIVATE FUNCTIONS ###


## function for the richness and evenness differences, gains and losses, and rankshifts returning a dataframe with those and the MRSc output
#rename this
aggfunc <- function(df, rank.var1, rank.var2, abundance.var1, abundance.var2){
  #ricness and evenness differences
  s_t1 <- S(df[[abundance.var1]])
  e_t1 <- EQ(as.numeric(df[[abundance.var1]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
  
  sdiff <- abs(s_t1-s_t2)/nrow(df)
  ediff <- abs(e_t1-e_t2)/nrow(df)
  
  #gains and losses
  df$gain <- ifelse(df[[abundance.var1]]==0, 1, 0)
  df$loss <- ifelse(df[[abundance.var2]]==0, 1, 0)
  
  gain <- sum(df$gain)/nrow(df)
  loss <- sum(df$loss)/nrow(df)
  
  mrsc <- mean(abs(df[[rank.var1]]-df[[rank.var2]])/nrow(df))
  
  metrics <- data.frame(Richness_change=sdiff, Evenness_change=ediff, Rank_change=mrsc, Gains=gain, Losses=loss)
  return(metrics)
}

S<-function(x){
  x1 <- x[x!=0 & !is.na(x)]
  stopifnot(x1==as.numeric(x1))
  length(x1)
}

# 2) function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
EQ<-function(x){
  x1<-x[x!=0 & !is.na(x)]
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
