library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

setwd("~/Dropbox/SESYNC/SESYNC_RACs/For R Package")
setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")
df<-read.csv('pplots_example.csv')

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
    filter(year==2003)
  
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
  bc_dis$BC_within_dissim_diff <- bc_dis$BC_dissim_within.x - bc_dis$BC_dissim_within.y

  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL
  
  ############
  ########### CHANGE
  ###########
  ############
  
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
  
  #mege into get bc_within differences for each treatment
  bc_dis1 <- merge(bc_between_ave, bc_within_ave, by = time.var)
  bc_dis <- merge(bc_dis1, bc_within_ave, by.x = paste(time.var, 2, sep = ""), by.y = time.var)
  
  #calculate absolute difference
  bc_dis$BC_within_change <- bc_dis$BC_dissim_within.x - bc_dis$BC_dissim_within.y
  
  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL
  
  

  