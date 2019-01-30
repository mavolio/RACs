library(tidyverse)
library(vegan)
library(devtools)

install_github("NCEAS/codyn",  ref = github_pull(83))
install_github("NCEAS/codyn", ref = "anderson")
library(codyn)


#Files from home

corredat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0")

corredat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-corredat1%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-corredat1%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-corredat1%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-corredat1%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat<-rbind(corredat1, azi, jrn, knz, sak)


####files from work
corredat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0")

corredat1<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)
corredat<-rbind(corredat1, azi, jrn, knz, sak)


#problems
#gvn face - only 2 years of data so will only have one point for the dataset.


plotinfo<-corredat%>%
  select(site_project_comm, calendar_year, plot_id, treatment, treatment_year)%>%
  unique()


#####CALCULATING DIVERSITY METRICs

spc<-unique(corredat$site_project_comm)
div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  div_eq<-rbind(div_eq, out)
}

##testing diversity

diversity<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-community_diversity(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id', metric = "Shannon")
  out$site_project_comm<-spc[i]
  
  diversity<-rbind(diversity, out)
}

#####CALCULATING RAC changes
spc<-unique(corredat$site_project_comm)
delta_rac<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  delta_rac<-rbind(delta_rac, out)
}

#####CALCULATING abundance changes
spc<-unique(corredat$site_project_comm)
delta_abund<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-abundance_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  delta_abund<-rbind(delta_abund, out)
}


##calculating multivariate changes

spc<-unique(corredat$site_project_comm)
delta_mult<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-multivariate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment='treatment')
  out$site_project_comm<-spc[i]
  
  delta_mult<-rbind(delta_mult, out)
}

###calculating curve change
spc<-unique(corredat$site_project_comm)
delta_curve<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-curve_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  delta_curve<-rbind(delta_curve, out)
}

# ##getting the average for each treatment in a year
# 
# corre_diversity<-merge(plotinfo, diversity, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(S_diff=mean(S_diff), even_diff=mean(E_diff, na.rm=T))
# 
# corre_gainloss<-merge(plotinfo, gain_loss, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(gain=mean(appearance), loss=mean(disappearance))
# 
# corre_reordering<-merge(plotinfo, reordering, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(MRSc=mean(MRSc))
# 
# ####MERGING TO A SINGE DATASET and exporting
# 
# merge1<-merge(corre_diversity, corre_gainloss, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
# merge2<-merge(merge1, corre_reordering, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
# all_metrics<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","treatment"), all=T)
# 
# write.csv(all_metrics, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allyears_2.csv")
# 
# 
# write.csv(all_metrics, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allyears_2.csv")


# ###Getting b-C distnace of each plot to itself comparing t1 to t2.
# 
# corredat$expplot<-paste(corredat$site_project_comm, corredat$plot_id, sep="::")
# 
# exp_plot_list<-unique(corredat$expplot)
# 
# 
# #makes an empty dataframe
# bray_curtis_dissim=data.frame(site_project_comm_plot=c(), calendar_year=c(), bc_dissim=c()) 
# 
# ##calculating bray-curtis mean change and disperison differecnes
# for(i in 1:length(exp_plot_list)) {
#   
#   #subsets out each dataset
#   subset=corredat%>%
#     filter(expplot==exp_plot_list[i])%>%
#     select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id)
#   
#   #get years
#   experiment_years<-sort(unique(subset$calendar_year))
#   
#   #transpose data
#   species=subset%>%
#     spread(genus_species, relcov, fill=0)
#   
#   #calculate bray-curtis dissimilarities
#   bc=as.matrix(vegdist(species[,5:ncol(species)], method="bray"))
#   
#   ###experiment_year is year x+1
#   bc_dis=data.frame(site_project_comm_plot=exp_plot_list[i],
#                     calendar_year=experiment_years[2:length(experiment_years)],
#                     bc_dissim=diag(bc[2:nrow(bc),1:(ncol(bc)-1)]))
#   
#   #pasting dispersions into the dataframe made for this analysis
#   bray_curtis_dissim=rbind(bc_dis, bray_curtis_dissim)  
# }
# 
# corre_braycurtis<-bray_curtis_dissim%>%
#   separate(site_project_comm_plot, into=c("site_project_comm","plot_id"), sep="::")
# 
# ###merging to a single dataset and adding treatment information
# merge1<-merge(gain_loss, diversity, by=c("site_project_comm","calendar_year","plot_id"))
# merge2<-merge(merge1, reordering,by=c("site_project_comm","calendar_year","plot_id")) 
# merge3<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","plot_id"))
# corre_all<-merge(plotinfo, merge3, by=c("site_project_comm","calendar_year","plot_id"))
# 
# write.csv(corre_all, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allReplicates.csv")
# 
# write.csv(corre_all, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allReplicates_2.csv")


#######Doing difference

#problems
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE)
#zero-length inputs cannot be mixed with those of non-zero length 

# #####CALCULATING RAC differences with block
# 
# blocked<-corredat%>%
#   filter(block!=0)%>%
#   filter(site_project_comm!="ARC_MNT_0"&site_project_comm!="BAY_LIND_0"&site_project_comm!="dcgs_gap_0"&site_project_comm!="JRN_study278_0"&site_project_comm!="KLU_KGFert_0"&site_project_comm!="KLU_BFFert_0"&site_project_comm!=""&site_project_comm!="LATNJA_CLIP_Heath"&site_project_comm!="LATNJA_CLIP_Meadow"&site_project_comm!="NWT_bowman_DryBowman"&site_project_comm!="NWT_bowman_WetBowman"&site_project_comm!="NWT_snow_0"&site_project_comm!="TRA_Lovegrass_0")
# 
# spc<-unique(blocked$site_project_comm)
# diff_rac_block<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-blocked%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
#   out$site_project_comm<-spc[i]
#   
#   diff_rac_block<-rbind(diff_rac_block, out)
# }
# 
# #####CALCULATING RAC differences without blocks pooling up to treatment
# 
# trt_control<-corredat%>%
#   filter(block==0)
# 
# spc<-unique(trt_control$site_project_comm)
# diff_rac_ct<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-trt_control%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
#   out$site_project_comm<-spc[i]
#   
#   diff_rac_ct<-rbind(diff_rac_ct, out)
# }
# 
# #####CALCULATING abundance differences with block
# spc<-unique(blocked$site_project_comm)
# diff_abund_block<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-blocked%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
#   out$site_project_comm<-spc[i]
#   
#   diff_abund_block<-rbind(diff_abund_block, out)
# }
# #####CALCULATING abundance differences without blocks pooling up to treatment
# 
# trt_control<-corredat%>%
#   filter(block==0)
# 
# spc<-unique(trt_control$site_project_comm)
# diff_abund_ct<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-trt_control%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
#   out$site_project_comm<-spc[i]
#   
#   diff_abund_ct<-rbind(diff_abund_ct, out)
# }
# #####CALCULATING RAC differences without blocks pooling up to treatment for all datasets
# spc<-unique(corredat$site_project_comm)
# diff_rac_all<-data.frame()
# 
# #####CALCULATING RAC differences with block
# 
# blocked<-corredat%>%
#   filter(block!=0)%>%
#   filter(site_project_comm!="ARC_MNT_0"&site_project_comm!="BAY_LIND_0"&site_project_comm!="dcgs_gap_0"&site_project_comm!="JRN_study278_0"&site_project_comm!="KLU_KGFert_0"&site_project_comm!="KLU_BFFert_0"&site_project_comm!=""&site_project_comm!="LATNJA_CLIP_Heath"&site_project_comm!="LATNJA_CLIP_Meadow"&site_project_comm!="NWT_bowman_DryBowman"&site_project_comm!="NWT_bowman_WetBowman"&site_project_comm!="NWT_snow_0"&site_project_comm!="TRA_Lovegrass_0")
# 
# 
# spc<-unique(blocked$site_project_comm)
# diff_rac_block<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-blocked%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
#   out$site_project_comm<-spc[i]
#   
#   diff_rac_block<-rbind(diff_rac_block, out)
# }
# 
# #####CALCULATING RAC differences without blocks pooling up to treatment
# 
# trt_control<-corredat%>%
#   filter(block==0)
# 
# spc<-unique(trt_control$site_project_comm)
# diff_rac_ct<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-trt_control%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)
#   out$site_project_comm<-spc[i]
#   
#   diff_rac_ct<-rbind(diff_rac_ct, out)
# }
# 
# #####CALCULATING abundance differences with block
# spc<-unique(blocked$site_project_comm)
# diff_abund_block<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-blocked%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
#   out$site_project_comm<-spc[i]
#   
#   diff_abund_block<-rbind(diff_abund_block, out)
# }
#####CALCULATING RAC differences without blocks pooling up to treatment
spc<-unique(corredat$site_project_comm)
diff_rac_all<-data.frame()
for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)
  out$site_project_comm<-spc[i]
  
  diff_rac_all<-rbind(diff_rac_all, out)
}
# 
# pplots<-subset(corredat, site_project_comm == 'KNZ_pplots_0')
# test<- RAC_difference(pplots, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)


##test with reference treatment
test<-subset(corredat, site_project_comm == "SCL_TER_0")
testing<-RAC_difference(test, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE, reference.treatment = "OO")

#####CALCULATING abundance differences without blocks pooling up to treatment for all datasets
spc<-unique(corredat$site_project_comm)
diff_abund_all<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)
  out$site_project_comm<-spc[i]
  
  diff_abund_all<-rbind(diff_abund_all, out)
}

##calculating multivariate differences

spc<-unique(corredat$site_project_comm)
diff_mult<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment='treatment')
  out$site_project_comm<-spc[i]
  
  diff_mult<-rbind(diff_mult, out)
}

write.csv(diff_mult, "~/Dropbox/SESYNC/SESYNC_RACs/R Files/anderson_corre_multdiff_new.csv", row.names = F )"

# ##calculating dissimilarity differences
# 
# spc<-unique(corredat$site_project_comm)
# diff_dissim<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-corredat%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-dissimilarity_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment='treatment')
#   out$site_project_comm<-spc[i]
#   
#   diff_dissim<-rbind(diff_dissim, out)
# }

write.csv(diff_dissim,"C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\Bray_Curtis_Ave_dissim_03162018.csv" )

write.csv(test,"C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\BCave_Vs_centdist_03162018.csv" )

# ###why are there differences?
# diff_dissim$calendar_year <- as.integer(diff_dissim$calendar_year)
# diff_mult$calendar_year <- as.integer(diff_mult$calendar_year)
# 
# test <- diff_dissim%>%
#   full_join(diff_mult)
# 
# 
# plot(test$centroid_distance_diff, test$BC_between_diff)
# cor.test(test$centroid_distance_diff, test$BC_between_diff)
# 
# test2<-unique(diff_dissim$site_project_comm)
# 
# num_dissim<-diff_dissim%>%
#   select(site_project_comm, calendar_year, treatment, treatment2, BC_between_diff)%>%
#   unique()%>%
#   mutate(present=1)
# 
# num_mult<-diff_mult%>%
#   select(site_project_comm, calendar_year, treatment, treatment2, abs_dispersion_diff)%>%
#   unique()%>%
#   mutate(present2=1)%>%
#   full_join(num_dissim)


#####CALCULATING curve differences with pooling
diff_curve_block<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-curve_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)
  out$site_project_comm<-spc[i]
  
  diff_curve_block<-rbind(diff_curve_block, out)
}

#####CALCULATING curve differences without blocks pooling up to treatment
spc<-unique(corredat$site_project_comm)
diff_curve_ct<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-curve_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = TRUE)
  out$site_project_comm<-spc[i]
  
  diff_curve_ct<-rbind(diff_curve_ct, out)
}


# merge1<-merge(div_diff, reordering_ct, by=c("site_project_comm","calendar_year","treatment"))
# all_Cont_Treat_Compare<-merge(merge1, corre_braycurtis_control_treat,by=c("site_project_comm","calendar_year","treatment"))
# 
# write.csv(all_Cont_Treat_Compare, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_ContTreat_Compare_OCT2017.csv")
# 
