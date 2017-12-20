#how to require vegan?
#make if else if a treatment column is there or not

#'@title Multivariate changes in composition and dispersion
#' @description 
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column


#do i need to remove the loop?
#problem with dropping columns
multivariate_change <- function(df, time.var, species.var, abundance.var, replicate.var, treatment.var=NULL){
  require(vegan)
  
    if(is.null(treatment.var)){
  
    #Calculating bc mean change and dispersion
        #get years
        timestep<-sort(unique(df[[time.var]]))
        
        #transpose data
        df2<-subset(df, select=c(time.var,species.var,abundance.var,replicate.var))
        df2$id<-paste(df2[[time.var]], df2[[replicate.var]], sep="_")
        species<-codyn:::transpose_community(df2, 'id',species.var, abundance.var)
        species$id<-row.names(species)
        speciesid<-do.call(rbind.data.frame, strsplit(species$id, split="_"))
        colnames(speciesid)[1]<-time.var
        colnames(speciesid)[2]<-replicate.var
        species2<-cbind(speciesid, species)
        
        #calculate bray-curtis dissimilarities
        bc=vegdist(species2[,4:ncol(species2)], method="bray")
        
        #calculate distances of each plot to year centroid (i.e., dispersion)
        disp=betadisper(bc, species2$time, type="centroid")
        
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
        Mult_Comp_Disp_Change<-distances
      
    }
  else{
    
    
  #make a new dataframe with just the label;
  treatlist<-unique(df[[treatment.var]])
  
  #makes an empty dataframe
  # Mult_Comp_Disp_Change=data.frame(replicate=c(), treatment=c(),time=c(), compositional_change=c(), dispersion_change=c()) 
  # 
  

  
  # # calculate change for each treatment
  # df <- df[order(df[[replicate.var]]),]
  # X <- split(df, df[replicate.var])
  # out <- lapply(X, FUN=fill_zeros, time.var, species.var, abundance.var)
  # ID <- unique(names(out))
  # out <- mapply(function(x, y) "[<-"(x, treatment.var, value = y) ,
  #               out, ID, SIMPLIFY = FALSE)
  # allsp <- do.call("rbind", out)
  
  #Calculating bc mean change and dispersion PUT ALL THIS IN A FUNCTION
  
  df2<-subset(df, select=c(time.var,replicate.var, species.var,abundance.var, treatment.var))
  df2$id<-paste(df2[[time.var]], df2[[replicate.var]], df2[[treatment.var]], sep="_")
  df2<-subset(df2, select=-c(time.var,replicate.var))#why doesn't this work?
  species<- reshape(df2, idvar="id", timevar=species.var, direction="wide")
  species[is.na(wide)] <- 0
  speciesid<-do.call(rbind.data.frame, strsplit(species$id, split="_"))
  colnames(speciesid)[1]<-"time"
  colnames(speciesid)[2]<-"replicate"
  colnames(speciesid)[3]<-"treatment"
  species2<-cbind(speciesid, species)
  
  for(i in 1:length(treatlist)) {
    
    #subsets out each dataset
    subset=subset(species2, treatment==treatlist[i])  
    #get years
    timestep<-sort(unique(subset$time))
    
    #calculate bray-curtis dissimilarities
    bc=vegdist(species[,5:ncol(species)], method="bray") 
    
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
                        dispersion_diff=diff(disp2.2$dist),
                        treatment=treatlist[i])
    
    #merge together change in mean and dispersion data
    distances<-merge(cent_dist_yrs, disp_yrs, by=c("time"))
   
    return(distances) 
    #pasting dispersions into the dataframe made for this analysis
    Mult_Comp_Disp_Change=rbind(distances, Mult_Comp_Disp_Change)  
  }
  }
  return(Mult_Comp_Disp_Change)
}

test2<-multivariate_change(df, time.var="time", replicate.var = "replicate", treatment.var = "treatment", species.var = "species", abundance.var = "abundance")
test1<-multivariate_change(df3, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")
