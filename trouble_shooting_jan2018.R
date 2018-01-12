df<-subset(codyndat_clean, site_project_comm=="AND.control.0")
df<-subset(codyndat_clean, site_project_comm=='AND.postlog.0')

test<-data.frame()
rep<-unique(df$plot_id)

for (i in 1:length(rep)){
  
  subset<-subset(df, plot_id==rep[i])
  
  long<-fill_zeros(subset, 'experiment_year', 'species', 'abundance')
  
  long$plot_id<-rep[i]
  test<-rbind(test, long)
  
}
test<-RAC_change(df, 'experiment_year','species','abundance','plot_id')  
test<-add_ranks_time(df, 'experiment_year', 'species', 'abundance','plot_id')
test<-fill_zeros(subset(df, plot_id==5), 'experiment_year','species','abundance')
test<-abundance_change(df,'experiment_year','species','abundance','plot_id')


##looking at where curve_change is different
cc<-codyn_curvechange%>%
  separate(experiment_year_pair, into=c('year1','experiment_year'), sep = "-", remove=F)

test<-merge(d_output, cc, by=c('experiment_year','site_project_comm','plot_id','curve_change'), all=T)

test_cc<-curve_change(df, 'experiment_year', 'species', 'abundance','plot_id')

df_test<-subset(subset_t12_3, plot_id==11)

  timestep2t <- unique(df_test$experiment_year)#assumes this is a length of 2

df1t <- subset(df_test, experiment_year == timestep2[1])
df2t <- subset(df_test, experiment_year == timestep2[2])
sf1t <- stepfun(df1$relrank, c(0, df1$cumabund))
sf2t <- stepfun(df2$relrank, c(0, df2$cumabund))
rt <- sort(unique(c(0, df1t$relrank, df2t$relrank)))
ht <- abs(sf1t(rt) - sf2t(rt))
wt <- c(diff(rt), 0)
CCt=sum(wt*ht)

result <- subset_t12_3 %>%
    df1 <- filter(., experiment_year==y[[1]])
    df2 <- filter(., experiment_year==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
