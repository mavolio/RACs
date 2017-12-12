calculate_abundance_change <- function(df, replicate.var, species.var, abundance.var, time.var) {
  rankdf <- add_ranks(df, replicate.var, species.var, abundance.var, time.var)
  
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
  
  
  out <- lapply(X, FUN=abundchangefunc, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
  ID <- unique(names(out))
  out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                out, ID, SIMPLIFY = FALSE)
  output <- do.call("rbind", out)  
  
  outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'_',fixed=TRUE)))
  names(outnames) = c(replicate.var, time.var)
  
  output$splitvariable <- NULL
  output <- cbind(outnames, output)
  
  return(output)
}


### PRIVATE FUNCTIONS ###

## function for the richness and evenness differences, gains and losses, and returning a dataframe with those and the MRSc output
abundchangefunc <- function(df, abundance.var1, abundance.var2){
  #ricness and evenness differences
  delta_abundance<-df[[abundance.var.1]]-df[[abundance.var2]]
  
  return(delta_abundance)
}
