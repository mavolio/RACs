##calculating abundance change of each species
#' @title Species Abundance Changes
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' 

abundance_change <- function(df, time.var, species.var, abundance.var, replicate.var=NULL) {
  if(is.null(replicate.var)){
    
    rankdf <- add_ranks(df, time.var, species.var, abundance.var, replicate.var=NULL)
    
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
    
    
    out <- lapply(X, FUN=abundchangefunc, paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    output <- do.call("rbind", out)  

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
    
    
    out <- lapply(X, FUN=abundchangefunc, paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    output <- do.call("rbind", out)  
    
    outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
    names(outnames) = c(replicate.var, time.var)
    outnames <- outnames[1]
    
    output$splitvariable <- NULL
    output <- cbind(outnames, output)
    
  }
  return(output)
}



### PRIVATE FUNCTIONS ###

abundchangefunc <- function(df, abundance.var1, abundance.var2){

  df$abund_change <- df[[abundance.var1]]-df[[abundance.var2]]
  df$time_pair<-paste(unique(df[[time.var]])-1, unique(df[[time.var]]), sep="-")
  
  df <- subset(df, select = c("time_pair", species.var, "abund_change"))
  
  colnames(df)[1]<-paste(time.var, "pair", sep="_")
  return(df)
}
