#'@title Multivariate changes in composition and dispersion
#' @description 
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column


multivariate_change <- function(df, time.var, species.var, abundance.var, replicate.var, treatment.var = NULL){

    if(is.null(treatment.var)){
  
      mult_change(df, time.var, species.var, abundance.var, replicate.var)
      
      mult_com_change <- distances
    }
  
  else {
    
  # calculate change for each treatment
  df <- df[order(df[[treatment.var]]),]
  X <- split(df, df[treatment.var])
  out <- lapply(X, FUN = mult_change, time.var, species.var, abundance.var, replicate.var)
  ID <- unique(names(out))
  out <- mapply(function(x, y) "[<-"(x, treatment.var, value = y) ,
                out, ID, SIMPLIFY = FALSE)
  mult_com_change <- do.call("rbind", out)
  
  
  }
  
  return(mult_com_change)
}

#######Private Functions

mult_change <- function(df, time.var, species.var, abundance.var, replicate.var){
  require(vegan)

  #get years
timestep <- sort(unique(df[[time.var]]))

#transpose data
df2<-subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="_")
species<-codyn:::transpose_community(df2, 'id', species.var, abundance.var)
species$id <- row.names(species)
speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="_"))
colnames(speciesid)[1] <- time.var
colnames(speciesid)[2] <- replicate.var
species2 <- cbind(speciesid, species)
species3 <- subset(species2, select = -id)

#calculate bray-curtis dissimilarities
bc <- vegdist(species3[,3:ncol(species3)], method="bray")

#calculate distances of each plot to year centroid (i.e., dispersion)
disp <- betadisper(bc, species3[[time.var]], type="centroid")

#getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))

##extracting only the comparisions we want year x to year x=1.
###year x+1
cent_dist_yrs <- data.frame(time = timestep[2:length(timestep)],
                         composition_change = diag(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist)-1)]))

#collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
disp2 <- data.frame(time=species2[[time.var]],
                 dist = disp$distances)

myformula <- as.formula(paste("dist", "~", time.var))
disp2.2<-aggregate(myformula, mean, data=disp2)

##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
disp_yrs <- data.frame(time = timestep[2:length(timestep)],
                    dispersion_change = diff(disp2.2$dist))

#merge together change in mean and dispersion data
distances <- merge(cent_dist_yrs, disp_yrs, by=time.var)

distances$time_pair<-paste(distances[[time.var]]-1, distances[[time.var]], sep="_")

distances<-subset(distances, select = c("time_pair", "composition_change", "dispersion_change"))
colnames(distances)[1]<-paste(time.var, "pair", sep="_")

return(distances)
}