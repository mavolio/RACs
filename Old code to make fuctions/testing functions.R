library(tidyverse)
library(codyn)
library(vegan)


setwd("~/Dropbox/SESYNC/SESYNC_RACs/For R Package")
setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")
pdata<-read.csv('pplots_example.csv')

#not working


#two steps working
cd2<-curve_diff_func(df)
sergl<-SERGL_func(df)
sersp<-SERSp_func(df)
spchange<-sp_abund_change_func(df)
spdiff<-sp_abund_diff_func(df)


#working
mult_change<-multivariate_change_func(df)
mult_diff<-multivariate_diff_func(df)
cc<-curve_change_func(df)


####testing multivariate with old code
#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps.
##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.
###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

df<-pplots%>%
  filter(treatment=="N1P0")

  #get years
  experiment_years<-sort(unique(df$year))
  
  #transpose data
  species=df%>%
    spread(species, relative_cover, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$year, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(experiment_year=experiment_years[2:length(experiment_years)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(experiment_year=species$year,
                   plot_id=species$plot,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(experiment_year)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(experiment_year=experiment_years[2:length(experiment_years)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by="experiment_year")
  
####differences
  time.var <- 'year'
  species.var <- 'species'
  replicate.var <- 'plot'
  treatment.var <- 'treatment'
  abundance.var <- 'relative_cover'
  
  df1<-pplots%>%
    filter(year==2002)
  
  time.var <- 'calendar_year'
  species.var <- 'genus_species'
  replicate.var <- 'plot_id'
  treatment.var <- 'treatment'
  abundance.var <- 'relcov'
  
  df<-subset(corredat, site_project_comm=="KNZ_pplots_0"&calendar_year==2003&treatment %in% c("N1P0","N2P0","N2P3"))
  df2<-subset(corredat, site_project_comm=="KNZ_pplots_0"&calendar_year==2003)
  
  
  
 
  ############
  #removing relevent info from bc matrix
  #calculate bray-curtis dissimilarities
  df2<-subset(df1, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  
  bc <- as.data.frame(as.matrix(vegdist(species, method="bray")))
  
  #extracting lower diagonal
  bc2 <- as.data.frame(cbind(rownames(bc)[which(lower.tri(bc), arr.ind=T)[,1]],
                             colnames(bc)[which(lower.tri(bc), arr.ind=T)[,2]],
                             bc[lower.tri(bc)]))
  c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
  c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
  
  bc3 <- cbind(bc2, c1, c2)
  bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
  colnames(bc3)[4] <- paste(treatment.var, 2, sep="")
  colnames(bc3)[6] <- treatment.var
  bc3$compare <- ifelse(bc3[[treatment.var]] == bc3[[paste(treatment.var, 2, sep="")]], 1, 2)
  
  #within treatment differences
  bc_within <- subset(bc3, compare == 1)
  myformula <- as.formula(paste("bc_dissim", "~", treatment.var))
  bc_within_ave <- aggregate(myformula, mean, data=bc_within)
  colnames(bc_within_ave)[2] <- "BC_dissim_within"
  
  #between treatment differences
  bc_between <- subset(bc3, compare == 2)
  myformula2 <- as.formula(paste("bc_dissim", "~", treatment.var, "+", paste(treatment.var, 2, sep = "")))
  bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
  colnames(bc_between_ave)[3] <- "BC_disism_between_diff"
  
  #mege into get bc_within differences for each treatment
  bc_dis1 <- merge(bc_between_ave, bc_within_ave, by = treatment.var)
  bc_dis <- merge(bc_dis1, bc_within_ave, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  
  #calculate absolute difference
  bc_dis$BC_within_dissim_diff <- bc_dis$BC_dissim_within.y - bc_dis$BC_dissim_within.x
  
  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL
  
  ############
  ########### CHANGE
  ###########
  ############
  data("pplots")
  time.var <- 'year'
  species.var <- 'species'
  replicate.var <- 'plot'
  treatment.var <- 'treatment'
  abundance.var <- 'relative_cover'
  
  df1<-pplots%>%
    filter(treatment == "N1P0")
  
  #removing relevent info from bc matrix
  #calculate bray-curtis dissimilarities
  df2<-subset(df1, select = c(time.var, species.var, abundance.var, replicate.var))
  df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  
  bc <- as.data.frame(as.matrix(vegdist(species, method="bray")))
  
  #extracting lower diagonal
  bc2 <- as.data.frame(cbind(rownames(bc)[which(lower.tri(bc), arr.ind=T)[,1]],
                             colnames(bc)[which(lower.tri(bc), arr.ind=T)[,2]],
                             bc[lower.tri(bc)]))
  c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
  c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
  
  bc3 <- cbind(bc2, c1, c2)
  bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
  colnames(bc3)[4] <- paste(time.var, 2, sep="")
  colnames(bc3)[6] <- time.var
  bc3$compare <- ifelse(bc3[[time.var]] == bc3[[paste(time.var, 2, sep="")]], 1, 2)
  
    #within treatment differences
  bc_within <- subset(bc3, compare == 1)
  myformula <- as.formula(paste("bc_dissim", "~", time.var))
  bc_within_ave <- aggregate(myformula, mean, data=bc_within)
  colnames(bc_within_ave)[2] <- "BC_dissim_within"
  
  #between treatment differences
  bc_between <- subset(bc3, compare == 2)
  myformula2 <- as.formula(paste("bc_dissim", "~", time.var, "+", paste(time.var, 2, sep = "")))
  bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
  colnames(bc_between_ave)[3] <- "BC_between_change"
  
  bc_between_ave$yr1 <- as.integer(as.factor(bc_between_ave[[time.var]]))
  bc_between_ave$yr2 <- as.integer(as.factor(bc_between_ave[[paste(time.var, 2, sep = "")]]))
  bc_between_ave$diff <- bc_between_ave$yr2 - bc_between_ave$yr1
  bc_between_ave2 <- subset(bc_between_ave, diff==1)
  bc_between_ave2$yr1 <- NULL
  bc_between_ave2$yr2 <- NULL
  bc_between_ave2$diff <- NULL
  
  #mege into get bc_within differences for each treatment
  bc_dis1 <- merge(bc_between_ave2, bc_within_ave, by = time.var)
  bc_dis <- merge(bc_dis1, bc_within_ave, by.x = paste(time.var, 2, sep = ""), by.y = time.var)
  
  #calculate absolute difference
  bc_dis$BC_within_change <- bc_dis$BC_dissim_within.y - bc_dis$BC_dissim_within.x
  
  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL
  
  
#######testing curve_change/diff
  #change
  df<-subset(df, replicate ==25 & time %in% c(2002, 2003))%>%
    filter(abundance!=0)%>%
    group_by(time, replicate)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  
  df <- df[order(df$time, df$cumabund),]
  
  timestep2 <- unique(df$time)#assumes this is a length of 2
  
  df1 <- df[df$time == timestep2[1],]
  df2 <- df[df$time == timestep2[2],]
  sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
  sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
  r <- sort(unique(c(0, df1$relrank, df2$relrank)))
  h <- abs(sf1(r) - sf2(r))
  w <- c(diff(r), 0)
  CC=sum(w*h)
  
  
  #difference just two reps
  df<-subset(pdata, replicate %in% c(6,25) & time == 2002)%>%
    filter(abundance!=0)%>%
    group_by(time, replicate)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  
  df <- df[order(df$time, df$cumabund),]
  
  replicate2 <- unique(df$replicate)#assumes this is a length of 2
  
  df1 <- df[df$replicate == replicate2[1],]
  df2 <- df[df$replicate == replicate2[2],]
  sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
  sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
  r <- sort(unique(c(0, df1$relrank, df2$relrank)))
  h <- abs(sf1(r) - sf2(r))
  w <- c(diff(r), 0)
  CC=sum(w*h)
  
  #difference just in a block
  df<-subset(pdata, replicate %in% c(29,25) & time == 2002)%>%
    filter(abundance!=0)%>%
    group_by(time, replicate)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  
  df <- df[order(df$time, df$cumabund),]
  
  replicate2 <- unique(df$replicate)#assumes this is a length of 2
  
  df1 <- df[df$replicate == replicate2[1],]
  df2 <- df[df$replicate == replicate2[2],]
  sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
  sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
  r <- sort(unique(c(0, df1$relrank, df2$relrank)))
  h <- abs(sf1(r) - sf2(r))
  w <- c(diff(r), 0)
  CC=sum(w*h)

  ######  difference pooling
  
  
  pdata2<-pdata%>%
    spread(species, abundance, fill = 0)%>%
    gather( species, abundance, 5:102)%>%
    group_by(time, treatment, species)%>%
    summarize(abundance = mean(abundance))
  
  df<-subset(pdata2, treatment %in% c("N2P0", "N1P0") & time == 2002)%>%
    filter(abundance!=0)%>%
    group_by(treatment)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  trt <- unique(df$treatment)#assumes this is a length of 2
  
  df1 <- df[df$treatment == trt[1],]
  df2 <- df[df$treatment == trt[2],]
  sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
  sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
  r <- sort(unique(c(0, df1$relrank, df2$relrank)))
  h <- abs(sf1(r) - sf2(r))
  w <- c(diff(r), 0)
  CC=sum(w*h)
  
  
  ####why is it not working with culardoc?
  
  
  time.var <- 'calendar_year'
  species.var <- 'genus_species'
  replicate.var <- 'plot_id'
  treatment.var <- 'treatment'
  abundance.var <- 'relcov'
  
  
  df1<-subset(corredat, site_project_comm == "SEV_Nfert_0" & calendar_year == 2004)
  df<-df1[order(df1[[treatment.var]], df[[replicate.var]]),]
  df2<-subset(df1, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  
  bc <- as.data.frame(as.matrix(vegdist(species, method="bray")))
  
  #extracting lower diagonal
  bc2 <- as.data.frame(cbind(rownames(bc)[which(lower.tri(bc), arr.ind=T)[,1]],
                             colnames(bc)[which(lower.tri(bc), arr.ind=T)[,2]],
                             bc[lower.tri(bc)]))
  c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
  c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
  
  bc3 <- cbind(bc2, c1, c2)
  bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
  colnames(bc3)[4] <- paste(treatment.var, 2, sep="")
  colnames(bc3)[6] <- treatment.var
  bc3$compare <- ifelse(bc3[[treatment.var]] == bc3[[paste(treatment.var, 2, sep="")]], 1, 2)
  
  #within treatment differences
  bc_within <- subset(bc3, compare == 1)
  myformula <- as.formula(paste("bc_dissim", "~", treatment.var))
  bc_within_ave <- aggregate(myformula, mean, data=bc_within)
  colnames(bc_within_ave)[2] <- "BC_dissim_within"
  
  #between treatment differences
  bc_between <- subset(bc3, compare == 2)
  myformula2 <- as.formula(paste("bc_dissim", "~", treatment.var, "+", paste(treatment.var, 2, sep = "")))
  bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
  colnames(bc_between_ave)[3] <- "BC_disism_between_diff"
  
  #mege into get bc_within differences for each treatment
  bc_dis1 <- merge(bc_between_ave, bc_within_ave, by = treatment.var)
  bc_dis <- merge(bc_dis1, bc_within_ave, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  
  #calculate absolute difference
  bc_dis$BC_within_dissim_diff <- bc_dis$BC_dissim_within.y - bc_dis$BC_dissim_within.x
  
  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL
  
  
  
  #################trying to get this into 1 function
  
  df2<-subset(df1, select = c(time.var, species.var, abundance.var, replicate.var))
  df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  species.1 <- species
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- "time_forfixxyz"
  colnames(speciesid)[2] <- replicate.var
  speciesid[[time.var]]<-as.numeric(as.character(speciesid$time_forfixxyz))
  speciesid.2 <- subset(speciesid, select = -time_forfixxyz)
  species2 <- cbind(speciesid.2, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissmilarity
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
         
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[time.var]], type="centroid")
         
         
    #### doing compositional differences
    if(comp_metric == "ave_BC_dissim"){
      
    bc1 <- vegdist(species, method="bray")
    
    #extracting lower diagonal
    bc2 <- as.data.frame(cbind(rownames(bc1)[which(lower.tri(bc1), arr.ind=T)[,1]],
                               colnames(bc1)[which(lower.tri(bc1), arr.ind=T)[,2]],
                               bc1[lower.tri(bc1)]))
    c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
    c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
    
    bc3 <- cbind(bc2, c1, c2)
    bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
    colnames(bc3)[4] <- paste(time.var, 2, sep="")
    colnames(bc3)[6] <- time.var
    bc3$compare <- ifelse(bc3[[time.var]] == bc3[[paste(time.var, 2, sep="")]], 1, 2)
    
    #between time differences
    bc_between <- subset(bc3, compare == 2)
    myformula2 <- as.formula(paste("bc_dissim", "~", time.var, "+", paste(time.var, 2, sep = "")))
    bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
    colnames(bc_between_ave)[3] <- "BC_dissim_change"
    
    #select only consecutive years
    bc_between_ave$yr1 <- as.integer(as.factor(bc_between_ave[[time.var]]))
    bc_between_ave$yr2 <- as.integer(as.factor(bc_between_ave[[paste(time.var, 2, sep = "")]]))
    bc_between_ave$diff <- bc_between_ave$yr2 - bc_between_ave$yr1
    bc_between_ave2 <- subset(bc_between_ave, diff==1)
    bc_between_ave2$yr1 <- NULL
    bc_between_ave2$yr2 <- NULL
    bc_between_ave2$diff <- NULL
    
    comp <- bc_between_ave2
    
  } else {
    timestep <- sort(unique(df1[[time.var]]))
    
    #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
    cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))
    
    ##extracting only the comparisions, year x to year x+1.
    cent_dist_yrs <- data.frame(
      time1 = timestep[1:length(timestep)-1],
      time2 = timestep[2:length(timestep)],
      centroid_distance_change = diag(
        as.matrix(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist)-1)])))
    
    comp <- cent_dist_yrs
    colnames(comp)[1] <- time.var
    colnames(comp)[2] <- paste(time.var, 2, sep="")
  }
  
  ##doing dispersion
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2 <- data.frame(time=species2[[time.var]],
                      dist = disp$distances)
  colnames(disp2)[1] <-time.var
  
  myformula <- as.formula(paste("dist", "~", time.var))
  disp2.2<-aggregate(myformula, mean, data=disp2)

  ##merge together  
  bc_dis1 <- merge(comp, disp2.2, by = time.var)
  bc_dis <- merge(bc_dis1, disp2.2, by.x = paste(time.var, 2, sep = ""), by.y = time.var)
  
  #calculate absolute difference
  bc_dis$disp_change <- bc_dis$dist.y - bc_dis$dist.x
  
  bc_dis$dist.y <- NULL
  bc_dis$dist.x <- NULL
  
  #########getting difference into 1 functino
  mult_diff <- function(df, time.var, species.var, abundance.var, replicate.var, comp_metric) {
    
    df1<-df[order(df[[treatment.var]], df[[replicate.var]]),]
    df2<-subset(df1, select = c(treatment.var, species.var, abundance.var, replicate.var))
    df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
    species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
    species1<-species
    species$id <- row.names(species)
    speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
    colnames(speciesid)[1] <- treatment.var
    colnames(speciesid)[2] <- replicate.var
    species2 <- cbind(speciesid, species)
    species3 <- subset(species2, select = -id)
    
    #calculate bray-curtis dissmilarity
    bc <- vegdist(species3[,3:ncol(species3)], method="bray")
    
    #calculate distances of each plot to year centroid (i.e., dispersion)
    disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
    
    
    #### doing compositional differences
    if(comp_metric == "ave_BC_dissim"){
      bc1<-as.data.frame(as.matrix(vegdist(species1, method = "bray")))
      #extracting lower diagonal
      bc2 <- as.data.frame(cbind(rownames(bc1)[which(lower.tri(bc1), arr.ind=T)[,1]],
                                 colnames(bc1)[which(lower.tri(bc1), arr.ind=T)[,2]],
                                 bc1[lower.tri(bc1)]))
      c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
      c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
      
      bc3 <- cbind(bc2, c1, c2)
      bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
      colnames(bc3)[4] <- paste(treatment.var, 2, sep="")
      colnames(bc3)[6] <- treatment.var
      bc3$compare <- ifelse(bc3[[treatment.var]] == bc3[[paste(treatment.var, 2, sep="")]], 1, 2)
      
      #between time differences
      bc_between <- subset(bc3, compare == 2)
      myformula2 <- as.formula(paste("bc_dissim", "~", treatment.var, "+", paste(treatment.var, 2, sep = "")))
      bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
      colnames(bc_between_ave)[3] <- "BC_dissim_diff"
      
      comp <- bc_between_ave
      
    } else {
      #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
      cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))
      
      cent_dist2 <- as.data.frame(cbind(rownames(cent_dist)[which(lower.tri(cent_dist, diag=T), 
                                                                  arr.ind=T)[,1]],
                                        colnames(cent_dist)[which(lower.tri(cent_dist, diag=T), 
                                                                  arr.ind=T)[,2]],
                                        cent_dist[lower.tri(cent_dist, diag=T)]))
      cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
      cent_dist3[3]<-as.numeric(as.character(cent_dist3[[3]]))
      
      colnames(cent_dist3)[1] <- paste(treatment.var, 2, sep="")
      colnames(cent_dist3)[2] <- treatment.var
      colnames(cent_dist3)[3] <- "centroid_distance_diff"
      
      comp <- cent_dist3
    }
    
    ##doing dispersion
    disp2 <- data.frame(treatment=species3[[treatment.var]],
                        dist = disp$distances)
    
    myformula <- as.formula(paste("dist", "~", treatment.var))
    disp2.2 <- aggregate(myformula, mean, data=disp2)
    
    #mege into get dispersion for each treatment
    cent_dist_disp <- merge(comp, disp2.2, by = treatment.var)
    cent_dist_disp2 <- merge(cent_dist_disp, disp2.2, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
    
    cent_dist_disp2$dispersion_diff <- cent_dist_disp2$dist.y - cent_dist_disp2$dist.x
    
    cent_dist_disp2$dist.x <- NULL
    cent_dist_disp2$dist.y <- NULL
    
    return(cent_dist_disp2)
  }
  