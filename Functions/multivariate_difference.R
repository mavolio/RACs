#'@title Multivariate differences in composition and dispersion
#' @description 
#' @param df A data frame containing an optional time column, species, abundance and replicate, and treatment columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the name of the treatment column
#' 
multivariate_difference <- function(df, time.var=NULL, species.var, abundance.var, replicate.var, treatment.var){
  
  if(is.null(time.var)){
    
    
    
  }
  
  else{
    
  }
  
  return()
}

###private functions

mult_diff <- function(df, species.var, abundance.var, replicate.var, treatment.var){
  require(vegan)
  
  #get treatments
  labels <- sort(unique(df[[treatment.var]]))
  
  #transpose data
  df2<-subset(df, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="_")
  species<-codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="_"))
  colnames(speciesid)[1] <- treatment.var
  colnames(speciesid)[2] <- replicate.var
  species2 <- cbind(speciesid, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissimilarities
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
  
  #getting distances between treatments; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist <- as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  #extracting all treatment differences
  cent_dist2 <- as.data.frame(cbind(rownames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,1]],
        colnames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,2]],
        cent_dist[lower.tri(cent_dist, diag=T)]))
  cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
  
  colnames(cent_dist3)[1] <- paste(treatment.var, 2, sep="")
  colnames(cent_dist3)[2] <- treatment.var
  colnames(cent_dist3)[3] <- "comm_diff"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2 <- data.frame(treatment=species3[[treatment.var]],
                      dist = disp$distances)
  
  myformula <- as.formula(paste("dist", "~", treatment.var))
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