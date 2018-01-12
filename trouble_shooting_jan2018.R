df<-subset(codyndat_clean, site_project_comm=="AND.control.0")

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
