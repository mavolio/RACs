multivariate_change <- function(df){
  #make a new dataframe with just the label;
  treatlist<-unique(df$treatment)
  
  #makes an empty dataframe
  Mult_Comp_Disp_Change=data.frame(replicate=c(), time=c(), compositional_change=c(), dispersion_change=c()) 
  
  #Calculating bc mean change and dispersion
  for(i in 1:length(treatlist)) {
    
    #subsets out each dataset
    subset=subset(df, treatment==treatlist[i])  
    #get years
    timestep<-sort(unique(subset$time))
    
    #transpose data HOW TO DO THIS WITHOUT ANOTHER PACAGE
    species=subset%>%
      spread(species, abundance, fill=0)
    
    #calculate bray-curtis dissimilarities
    bc=vegdist(species[,5:ncol(species)], method="bray") #this is 4 here b/c there is a treatment column, but if obs data only would be ,3:ncol()
    
    #calculate distances of each plot to year centroid (i.e., dispersion)
    disp=betadisper(bc, species$time, type="centroid")
    
    #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
    cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
    
    ##extracting only the comparisions we want year x to year x=1.
    ###(experiment_year is year x+1
    cent_dist_yrs=data.frame(time=timestep[2:length(timestep)],
                             mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
    
    #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
    disp2=data.frame(time=species$time,
                     dist=disp$distances)
    
    disp2.2<-aggregate(dist~time, mean, data=disp2)
    
    
    ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
    disp_yrs=data.frame(time=timestep[2:length(timestep)],
                        dispersion_diff=diff(disp2.2$dist))
    
    #merge together change in mean and dispersion data
    distances<-merge(cent_dist_yrs, disp_yrs, by=c("time"))
    
    #pasting dispersions into the dataframe made for this analysis
    Mult_Comp_Disp_Change=rbind(distances, Mult_Comp_Disp_Change)  
  }
  return(Mult_Comp_Disp_Change)
}