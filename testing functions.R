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
  
  
  
  df2<-subset(df, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- treatment.var
  colnames(speciesid)[2] <- replicate.var
  species2 <- cbind(speciesid, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissimilarities
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
  
  #getting distances between treatments with euclidean distances
  cent_dist <- as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  #extracting all treatment differences
  cent_dist2 <- as.data.frame(cbind(rownames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,1]],
                                    colnames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,2]],
                                    cent_dist[lower.tri(cent_dist, diag=T)]))
  cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
  cent_dist3[3]<-as.numeric(as.character(cent_dist3[[3]]))
  
  colnames(cent_dist3)[1] <- paste(treatment.var, 2, sep="")
  colnames(cent_dist3)[2] <- treatment.var
  colnames(cent_dist3)[3] <- "composition_diff"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2 <- data.frame(treatment=species3[[treatment.var]],
                      dist = disp$distances)
  
  myformula <- as.formula(paste("dist", "~", treatment.var))
  disp2.2 <- aggregate(myformula, mean, data=disp2)
  
  #mege into get dispersion for each treatment
  cent_dist_disp <- merge(cent_dist3, disp2.2, by = treatment.var)
  cent_dist_disp2 <- merge(cent_dist_disp, disp2.2, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  
  #calculate absolute difference
  cent_dist_disp2[[treatment.var]] <- as.character(cent_dist_disp2[[treatment.var]])
  cent_dist_disp2[[paste(treatment.var, 2, sep = "")]] <- as.character(cent_dist_disp2[[paste(treatment.var, 2, sep = "")]])
  
  cent_dist_disp2$abs_dispersion_diff <- abs(cent_dist_disp2$dist.x - cent_dist_disp2$dist.y)
  cent_dist_disp2$trt_greater_disp <- as.character(ifelse(cent_dist_disp2$dist.x > cent_dist_disp2$dist.y, cent_dist_disp2[[treatment.var]], cent_dist_disp2[[paste(treatment.var, 2, sep = "")]]))
  
  cent_dist_disp2$dist.x <- NULL
  cent_dist_disp2$dist.y <- NULL
  
#######################
  ########################
  ##########################
  #####code from converge-diverge
  
  df<-subset(corredat, site_project_comm=="KNZ_pplots_0"&calendar_year==2003&treatment %in% c("N1P0","N2P0","N2P3"))%>%
    mutate(plot_mani = ifelse(treatment=="N1P0",0,1))
  df2<-subset(corredat, site_project_comm=="KNZ_pplots_0"&calendar_year==2003)%>%
    mutate(plot_mani = ifelse(treatment=="N1P0",0,1))
  

    #need this to keep track of plot mani
  labels=df2%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=df2%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,11:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  #merging back with labels to get back plot_mani
  centroid=merge(cent_C_T, labels, by="treatment")
  
  #collecting and labeling distances to centroid from betadisper
  trt_disp=data.frame(data.frame(exp_year=exp_year$exp_year[i], 
                                 plot_id=species$plot_id,
                                 treatment=species$treatment,
                                 dist=disp$distances))%>%
    tbl_df()%>%
    group_by(exp_year, treatment)%>%
    summarize(dispersion=mean(dist))
  
  
  ###testing with function
  pplots<-subset(corredat, site_project_comm == 'KNZ_pplots_0'&calendar_year==2003)
  test<- multivariate_difference(pplots, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment')
  