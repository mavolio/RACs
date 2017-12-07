library(tidyverse)
library(codyn)
library(vegan)



#work
setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")
df<-read.csv('pplots_example.csv')


#####CALCULATING DIVERSITY METRICS
##still need to make a funciton to encapsulate all these

com_struct <- function(df, replicate.var, abundance.var, time.var, evenness="E_q") {

  ###NO CLUE WHY THIS IS NOT GIVING ME MULTIPLE COLUMNS
  comstruct<-aggregate(abundance~time+replicate, data=df, FUN=function(x)c(SpR=S(x),evenness=E_q(x)))
  
  #works here
  head(aggregate(abundance~time+replicate, data=df, FUN=function(x)c(SpR=S(x),evenness=E_q(x))))
    
  #comstruct<-aggregate(df[[abundance.var]]~df[[replicate.var]]+df[[time.var]], FUN=function(x)c(SpR=S(x),evenness=E_q(x)), data=df)
  
  return(comstruct)
}
  
#1) function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

# 2) function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
E_q<-function(x){
  x1<-x[x!=0]
  if (length(x1)==1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r<-rank(x1, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x1)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]
  2/pi*atan(b)
}


#function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
E_simp<-function(x, S=length(x[x!=0]), N=sum(x[x!=0]), ps=x[x!=0]/N, p2=ps*ps ){
  D<-sum(p2)
  (1/D)/S
}


##FUNCTION TO ADD RANKS
#' @param for a dataset with columns for time, replicate, species, and abundance, (and optionally an id column for grouping and second column for defining the groupings)
#' 

add_ranks <- function(df, replicate.var, species.var, abundance.var, time.var) {
  ##SHOULD THE FIRST STEP TO BE ONLY SELECT THE RELEVANT COLUMNS?
  df<-subset(df, select=c("time", "abundance","species","replicate"))
  
  ##add ranks for present species
  rank_pres<-subset(df, df[[abundance.var]]!=0)
  rank_pres$rep_time<-paste(rank_pres[[replicate.var]], rank_pres[[time.var]], sep="_")
  rank_pres$rank<-ave(rank_pres[[abundance.var]], rank_pres$rep_time, FUN=function(x) rank(-x, ties.method = "average"))
  rank_pres<-subset(rank_pres, select=-rep_time)
  
  #adding zeros
  #NEED TO REMOVE THE LOOP HERE
  replist<-unique(df[[replicate.var]])
  
  allsp<-data.frame()
  for (i in 1:length(replist)){
    subset <- subset(df[[replicate.var]]==replist[i])
    
    replicate<-replist[i]
    
    #FOR RESHAPE HOW TO DEAL WITH COLUMN NAMES?
    subset2<-subset(subset, select=c("time","species","abundance"))
    wide<- reshape(subset2, idvar="time", timevar="species", direction="wide")
    wide[is.na(wide)] <- 0
      
    long<-reshape(wide, idvar="time", ids="time", time=names(wide), timevar="abundance", direction="long")
    long$replicate<-replist[i]
    colnames(long)[3]<-"abundance"
    
    allsp<-rbind(long, allsp)  
  }
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-subset(allsp, all.sp[[abundance.var]]==0)
  
  ##get species richness for each year
  SpR<-aggregate([[abundance.var]]~[[replicate.var]]+[[time.var]], FUN=S, data=allsp)
  colnames(SpR)[4]<-"S"

  ##merge together make zero abundances rank S+1
  #HOW TO DEAL WITH COLUMN NAMES IN MERGE??
  zero_rank<-merge(zeros, SpR, by=c("time","replicate","C_T"))
  zero_rank$rank<-zero_rank$S+1
  zero_rank<-subset(zero_rank, select=-S)
  
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  return(rank)
}

##calculating SERGL
# @param sim foo dataset with columns for time, plot, species, and abundance
calculate_SERGL <- function(rank, replicate.var, species.var, abundance.var, time.var) {
  
  SERGL=data.frame(replicate=c(), time=c(), S=c(), E=c(), R=c(), G=c(), L=c())#expeiment year is year of timestep2
  
  replist<-unique(rank[[replicate.var]])
  
  for (i in 1:length(replist)){
    #THIS BREAKS AT THIS STEP
    subset<-subset(rank, rank[[replicate.var]]==replist[i])
    
    replicate<-replist[i]
    
    #now get all timestep within an experiment
    timestep<-sort(unique(subset[[time.var]]))    
    
    #NEED TO REMOVE THE LOOP
    for(i in 1:(length(timestep)-1)) {
      subset_t1<-subset(subset, subset[[time.var]]==timestep[i])
      
      subset_t2<-subset(subset, subset[[time.var]]==timestep[i+1])
      
      #HOW TO MERGE BY COLUMN NAMES?
      subset_t12<-merge(subset_t1, subset_t2, by=c("species","replicate"), all=T)
      subset_t12<-subset(subset_t12, abundance.x!=0|abundance.y!=0)
      
      #reordering
      MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
      #ricness and evenness differences
      s_t1 <- S(subset_t12$abundance.x)
      e_t1 <- E_q(as.numeric(subset_t12$abundance.x))
      s_t2 <- S(subset_t12$abundance.y)
      e_t2 <- E_q(as.numeric(subset_t12$abundance.y))
      
      sdiff<-abs(s_t1-s_t2)/nrow(subset_t12)
      ediff<-abs(e_t1-e_t2)/nrow(subset_t12)
      
      #gains and losses
      subset_t12$gain<-ifelse(subset_t12$abundance.x==0, 1, 0)
      subset_t12$loss<-ifelse(subset_t12$abundance.y==0, 1, 0)
      
      gain<-sum(subset_t12$gain)/nrow(subset_t12)
      loss<-sum(subset_t12$loss)/nrow(subset_t12)
      
      metrics<-data.frame(replicate=replicate, time=timestep[i+1], S=sdiff, E=ediff, R=MRSc, G=gain, L=loss)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      SERGL=rbind(metrics, SERGL)  
    }
  }
  return(SERGL)
}

SERGL_func <- function(df) {
  rank <- add_ranks(df)
  result <- calculate_SERGL(rank)
  return(result)
}

#we can ignore this, once calcualte_SERGL is working this is basically the same thing and we can copy the code
# calc_sp_abund_changes<- function(rank){
#   sp_abund_changes<-data.frame()
#   
#   replist<-unique(rank$replicate)
#   
#        for (i in 1:length(replist)){
#         subset<-subset(rank, replicate==replist[i])
#         
#         replicate<-replist[i]
#         
#         #now get all timestep within an experiment
#         timestep<-sort(unique(subset$time))    
#         
#         for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
#           subset_t1<-subset(subset, time==timestep[i])
#           
#           subset_t2<-subset(subset, time==timestep[i+1])
#           
#           subset_t12<-merge(subset_t1, subset_t2, by=c("species","replicate"), all=T)
#           subset_t12<-subset(subset_t12, abundance.x!=0|abundance.y!=0)
#        ##top species changes
#       subset_t12$delta_abundance<-subset_t12$abundance.y-subset_t12$abundance.x
#       
#       abundchange<-subset_t12[,c(1,2, 8,13)]
#       colnames(abundchange)[3]<-"time"
#       
#       sp_abund_changes=rbind(sp_abund_changes, abundchange)
#     }
#   }
#   return(sp_abund_changes)
# }

# sp_abund_change_func <- function(df) {
#   rank <- add_ranks(df)
#   result <- calc_sp_abund_changes(rank)
#   return(result)
# }

# TODO define functions throughout whole file
# TODO rename data frame argument to df
# TODO add arguments to name the time, site, species, and abundance columns
# TODO refactor to allow users to use arbitrary column names (remove dplyr?)
# TODO refactor to optimize speed of function
# TODO add new function to codyn
# TODO document the functions fully
# TODO add tests to codyn for the new functions

#4/5)Bray-Curtis Mean Change & Bray-Curtis Dissimilarity 
##THIS NEEDS TO BE FIXED SO YOU HAVE 2 OPTIONS, IF OBS ONLY GROUP BY YEAR, IF TRT GROUP BY YEAR+TIME

multivariate_change_func <- function(df){
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


#############################################################
#############################################################
#Giving up here

#6) curve comparision
####Looking at the shape of the curve - cc
###this compares the areas of difference between two curves that are consequtive time steps for a plot.
curve_change_func <- function(df){
  curvechange=data.frame(replicate=c(), time=c(), CC=c())#expeiment year is year of timestep2
  
  ranks<-df%>%
    filter(abundance!=0)%>%
    group_by(time, replicate)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  timestep<-sort(unique(ranks$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(time==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(replicate)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(time==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(replicate)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="replicate")
    #dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
    ##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="replicate")
    
    #dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(time, replicate)%>%
      mutate(numplots=length(replicate))%>%
      ungroup()%>%
      group_by(replicate)%>%
      mutate(min=min(numplots))%>%
      select(replicate, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="replicate")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(time, replicate)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
      group_by(replicate) %>%
      do({
        y <- unique(.$time)###assumption this is a length 2 list
        df1 <- filter(., time==y[[1]])
        df2 <- filter(., time==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(CC=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(replicate=result$replicate, time=timestep[i+1], CC=result$CC)#expeiment year is year of timestep2
    
    curvechange<-rbind(curvechange, d_output1)
  }
  return(curvechange)
}

########comparing control versus treatment plots

#' @param sim foo dataset with columns for time, replicate, species, abundance, a treatment column for grouping replicates, and a column speficfying if the treatment is a control or a treatment
add_ranks_treatment_control_sppools <- function(df) {
  
  ###add zeros and average up species pool for control and treatment plots
  wide<-df%>%
    spread(species, abundance, fill=0)
  
  ##make long and get averages of each species by treatment
  long<-wide%>%
    gather(species, abundance, 5:ncol(wide))%>%
    group_by(time, treatment, species, C_T)%>%
    summarize(abundance=mean(abundance))
  
  ##add ranks dropping zeros
  rank_pres<-long%>%
    filter(abundance!=0)%>%
    tbl_df()%>%
    group_by(time, treatment, C_T)%>%
    mutate(rank=rank(-abundance, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-long%>%
    filter(abundance==0)
  ##get species richness for each year
  rich<-group_by(long, time, treatment, C_T)%>%
    summarize(S=S(abundance))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("time", "treatment","C_T"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  ct_rank<-rbind(rank_pres, zero_rank)
  return(ct_rank)
}

calculate_SERSp <- function(ct_rank){
  SERSp=data.frame(treatment=c(), time=c(), Sd=c(), Ed=c(), Rd=c(), spd=c())      
  timestep<-sort(unique(ct_rank$time)) 
  
  for(i in 1:(length(timestep))){
    
    time<-ct_rank%>%
      filter(time==timestep[i])
    
    time_id<-timestep[i]
    
    #need to do all comparisions NOT SURE HOW TO PROCEED FROM HERE TO SUBSET ALL POSSIBLE COMBINATIONS
    comparison<-combn(unique(time$treatment),2)
    
    
    #filter out control plots
    control<-time%>%
      filter(C_T=="Control")
    
    treat_list<-unique(subset(time, C_T=="Treatment")$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("time","species"), all=T)%>%
        filter(abundance.x!=0|abundance.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      
      spdiff<-subset_ct%>%
        filter(abundance.x==0|abundance.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      ##eveness richness
      s_c <- S(subset_ct$abundance.x)
      e_c <- E_q(subset_ct$abundance.x)
      s_t <- S(subset_ct$abundance.y)
      e_t <- E_q(subset_ct$abundance.y)
      
      sdiff<-abs(s_c-s_t)/nrow(subset_ct)
      ediff<-abs(e_c-e_t)/nrow(subset_ct)
      
      metrics<-data.frame(treatment=treat_id, time=time_id, Sd=sdiff, Ed=ediff, Rd=MRSc_diff, spd=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      SERSp=rbind(metrics, SERSp)  
    }
  }
  return(SERSp)
}

SERSp_func <- function(df) {
  ct_rank <- add_ranks_treatment_control_sppools(df)
  result <- calculate_SERSp(ct_rank)
  return(result)
}

calc_species_abund_differences<-function(ct_rank){
  spabunddiff<-data.frame()
  
  timestep<-sort(unique(ct_rank$time)) 
  
  for (i in 1:length(timestep)){
    time<-ct_rank%>%
      filter(time==timestep[i])
    
    time_id<-timestep[i]
    
    #fitler out control plots
    control<-time%>%
      filter(C_T=="Control")
    
    treat_list<-unique(subset(time, C_T=="Treatment")$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("time","species"), all=T)%>%
        filter(abundance.x!=0|abundance.y!=0)
      
      ##top species changes
      subset_ct$delta_abundance<-subset_ct$abundance.y-subset_ct$abundance.x
      
      spabund_diff<-subset_ct[, c(1:2,7,11)]
      colnames(spabund_diff)[3]<-"treatment"
      
      spabunddiff=rbind(spabunddiff, spabund_diff)
    }
  }
  return(spabunddiff)
}

sp_abund_diff_func <- function(df) {
  ct_rank <- add_ranks_treatment_control_sppools(df)
  result <- calc_species_abund_differences(ct_rank)
  return(result)
}

multivariate_diff_func <- function(df){
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

#curve cahnge
calculate_curve_diff <- function(rank){
  curvechange_diff=data.frame(treatment=c(), time=c(), CC=c())#expeiment year is year of timestep2
  
  ranks<-rank%>%
    filter(abundance!=0)%>%
    group_by(time, treatment, C_T)%>%
    mutate(maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  control<-ranks%>%
    filter(C_T=="Control")
  
  treat_list<-unique(subset(ranks, C_T=="Treatment")$treatment)
  
  for(i in 1:length(treat_list)) {
    subset_trt<-ranks%>%
      filter(treatment==treat_list[i])
    
    #dataset of two treatments    
    subset_t12<-rbind(control, subset_trt)
    
    result <- subset_t12 %>%
      group_by(time) %>%
      do({
        y <- unique(.$treatment)###assumption this is a length 2 list
        df1 <- filter(., treatment==y[[1]])
        df2 <- filter(., treatment==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(CC=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(treatment=treat_list[i], time=result$time, CC=result$CC)#expeiment year is year of timestep2
    
    curvechange_diff<-rbind(curvechange_diff, d_output1)
  }
  
  return(curvechange_diff)
}

curve_diff_func <- function(df) {
  ct_rank <- add_ranks_treatment_control_sppools(df)
  curve_change_diff <- calculate_curve_diff(ct_rank)
  return(curve_change_diff)
}
