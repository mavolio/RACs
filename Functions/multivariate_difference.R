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
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist <- as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###year x+1
  cent_dist_yrs <- data.frame(treatment.var = timestep[2:length(timestep)],
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
  
  
  year<-unique(df$time)

#makes an empty dataframe
Mult_Comp_Disp_Diff=data.frame() 
##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(year)) {
  
  #subsets out each dataset
  subset<-df%>%
    filter(time==year[i])%>%
    select(treatment, time, species, abundance, replicate, C_T)
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(C_T, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
    spread(species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(time=year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$C_T=="Control"],]))
  
  #renaming column
  colnames(cent_C_T)[3]<-"comp_diff"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2=data.frame(time=year[i],
                   treatment=species$treatment,
                   C_T=species$C_T,
                   replicate=species$replicate,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(time, treatment, C_T)%>%
    summarize(dispersion=mean(dist))
  
  control<-disp2$dispersion[disp2$C_T=="Control"]
  
  ##subtract control from treatments
  disp_treat=disp2%>%
    mutate(disp_diff=dispersion-control)%>%
    select(-dispersion)
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_C_T, disp_treat, by=c("time","treatment"))
  
  #pasting dispersions into the dataframe made for this analysis
  Mult_Comp_Disp_Diff=rbind(Mult_Comp_Disp_Diff, distances)  
}

return(Mult_Comp_Disp_Diff)
}