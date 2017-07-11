############### SESYNC Research Support: Community metrics ########## 
## Community metrics processing.
## 
## DATE CREATED: 07/05/2017
## DATE MODIFIED: 07/05/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Fisheries by Meghan Avolio
## ISSUE: 
## TO DO:
##
## COMMIT: learning to write a function first step
##
## Links to investigate:

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)
library(lubridate)
library(dplyr)

###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

### Add function S
S <- function(x){
  x1<-x[x!=0]
  length(x1)
}

### Other functions ####

#function_processing_data <- "community_metrics_reordering_functions_07052017.R" #PARAM 1
#script_path <- "/nfs/bparmentier-data/Data/projects/community_metrics/scripts" #path to script #PARAM 
#source(file.path(script_path,function_processing_data)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

in_dir <- "/nfs/bparmentier-data/Data/projects/community_metrics/data" #parmam 1
out_dir <- "/nfs/bparmentier-data/Data/projects/community_metrics/outputs" #param 2

num_cores <- 2 #param 3: number of cores used for parallel processing
create_out_dir_param=TRUE # param 4

out_suffix <-"reordering_07052017" #output suffix for the files and ouptut folder #param 5
in_filename <- "SimCom_June.csv" #param 6

col_year_id_name <- "" #param 7
col_plot_id_name <-  "" #param 8
  
############## START SCRIPT ############################

######### PART 0: Set up the output dir ################

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

#df_data <- read.table(file.path(in_dir),in_filename,sep=",")
df_data <- read.csv(file.path(in_dir,in_filename),header=T)

dim(df_data)

#### start building the function

#col_abundance_id <- "abundance"  
#col_time_id <- "time"
#col_plot_id <- "site"
#col_experiment_id <- "id"

debug(reorder_fun)
test <- reorder_fun(df_data=df_data,
                    col_abundance_id="abundance",
                    col_time_id="iteration",
                    col_plot_id="site",
                    col_experiment_id="id",
                    col_species_id="species")

  
reorder_fun <- function(df_data,col_abundance_id,col_time_id,
                        col_plot_id,
                        col_experiment_id,col_species_id){
  ## This function adds zero to ...
  #
  ### INPUTS:
  #1)
  #2)
  #3)
  ### OUTPUTS
  #object as a lits with items:
  #1)
  #2)
  
  ###### Begin script #######
  
  ##SIM dataset
  ##add in zeros
  
  ##add ranks 
  columns_group_by <- c(col_experiment_id, col_time_id, col_plot_id)
  sim_rank_pres<-df_data%>%
    filter(col_abundance_id!=0)%>%
    tbl_df()%>%
    #group_by(columns_group_by)%>%
    group_by_(col_experiment_id, col_time_id, col_plot_id)%>%
    mutate(rank=rank(-df_data[[col_abundance_id]], ties.method = "average"))%>%
    tbl_df()
  
  #adding zeros
  sim_addzero <- df_data %>%
    group_by(col_experiment_id) %>%
    nest() %>%
    mutate(spread_df = purrr::map(data, ~spread(., key=col_species_id, value=col_abundance_id, fill=0) %>%
                                    gather(key=col_species_id, value=col_abundance_id, -col_plot_id, -col_time_id))) %>%
    unnest(spread_df)
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  sim_zeros<-sim_addzero%>%
    filter(col_abundance_id==0)
  ##get species richness for each year
  sim_S<-group_by(df_data, col_experiment_id, col_time_id, col_plot_id)%>%
    summarize(S=S(col_abundance_id))
  ##merge together make zero abundances rank S+1
  sim_zero_rank<-merge(sim_zeros, sim_S, by=c(col_experiment_id,col_time_id,col_plot_id))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  sim_rank<-rbind(sim_rank_pres, sim_zero_rank)
  
  #prepare return dim_df
  #reorder_obj<- list(dim_df,col_names)
  #names(reorder_obj)<- c("dim_data_frame","col_names")
  
  return(sim_rank)
}



#######


reordering=data.frame(id=c(), calendar_year=c(), MRSc=c())



explist<-unique(corre$site_project_comm)



for (i in 1:length(explist)){
  
  ##get zero abundances to be filled in for all species.
  
  ##this works the first time only
  
  subset<-corre%>%
    
    filter(site_project_comm==explist[i])%>%
    
    spread(genus_species, relcov, fill=0)
  
  wide<-subset%>%
    
    gather(genus_species, relcov,11:ncol(subset))
  
  
  
  ##add ranks dropping zeros
  
  rank_pres<-wide%>%
    
    filter(relcov!=0)%>%
    
    tbl_df()%>%
    
    group_by(site_project_comm, calendar_year, treatment, plot_id)%>%
    
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    
    tbl_df()
  
  
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  
  ##pull out zeros
  
  zeros<-wide%>%
    
    filter(relcov==0)
  
  ##get species richness for each year
  
  rich<-group_by(wide, site_project_comm, calendar_year, plot_id)%>%
    
    summarize(S=S(relcov))
  
  ##merge together make zero abundances rank S+1
  
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year","plot_id"))%>%
    
    mutate(rank=S+1)%>%
    
    select(-S)%>%
    
    tbl_df()
  
  ##combine all
  
  rank<-rbind(rank_pres, zero_rank)
  
  
  
  #get uniuque id's
  
  spc_id<-unique(subset$id)
  
  
  
  for (i in 1:length(spc_id)){
    
    subset2<-rank%>%
      
      filter(id==spc_id[i])
    
    id<-spc_id[i]
    
    
    
    #now get all timestep within an experiment
    
    timestep<-sort(unique(subset2$calendar_year))    
    
    
    
    for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds
      
      subset_t1<-subset2%>%
        
        filter(calendar_year==timestep[i])
      
      
      
      subset_t2<-subset2%>%
        
        filter(calendar_year==timestep[i+1])
      
      
      
      subset_t12<-merge(subset_t1, subset_t2, by=c("genus_species","id"), all=T)%>%
        
        filter(relcov.x!=0|relcov.y!=0)
      
      
      
      MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
      
      
      
      metrics<-data.frame(id=id, calendar_year=timestep[i+1], MRSc=MRSc)#spc_id
      
      ##calculate differences for these year comparison and rbind to what I want.
      
      
      
      reordering=rbind(metrics, reordering)  
      
    }
    
  }
  
}