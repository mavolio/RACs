setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")
setwd("~/Dropbox/SESYNC/SESYNC_RACs/For R Package")
#all data
pdata<-read.csv('pplots_example.csv')

#no replicate
pdata2<-subset(pdata, replicate==1)
pdata2<-pdata2[,-3]

#no treatment
pdata3<-subset(pdata, treatment=="N1P0")
pdata3<-pdata3[-2]

replicate.var <- 'replicate'
treatment.var <- 'treatment'
species.var <- 'species'
time.var <- 'time'
abundace.var <- 'abundance'

#RAC_changes
#works with time and rep.
test1<-RAC_changes(pdata, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")
#works with no rep
test2<-RAC_changes(pdata2, time.var="time",species.var = "species", abundance.var = "abundance")

#community_diversity
#works with time and rep.
test1.shan<-community_diversity(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance")
test1.simp<-community_diversity(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance", diversity = "Simpson")
#works with no rep
test2.simp<-community_diversity(pdata2, time.var="time", abundance.var = "abundance", diversity = "Simpson")
test2.shan<-community_diversity(pdata2, time.var="time", abundance.var = "abundance")

#community_structure
test1.eq<-community_structure(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance")
test1.esimp<-community_structure(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance", evenness = "SimpEven")
#works with no rep
test2.esimp<-community_structure(pdata2, time.var="time", abundance.var = "abundance", evenness = "SimpEven")
test2.eq<-community_structure(pdata2, time.var="time", abundance.var = "abundance")

#multivariate change
#with treatment
test2<-multivariate_change(pdata, time.var="time", replicate.var = "replicate", treatment.var = "treatment", species.var = "species", abundance.var = "abundance")
#without treatment
test1<-multivariate_change(pdata3, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")

