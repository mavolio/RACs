###trying to learn functions for dataframes
J_longform<-function(df, time.var="year", species.var="species", abundance.var="abundance"){
  com.use<-codyn:::transpose_community(df, time.var, species.var, abundance.var)
  div.out<-diversity(com.use)
  rich.out<-specnumber(com.use)
  J<-div.out/log(rich.out)
  return(J)
}

df <- data.frame(
  plot=c(1, 1, 1, 2, 2, 2, 2),
  sp=c('a','b','c','a','b','c','d'),
  n=c(5,4,3,15,15,15,1)
)

plotsum<-function(df, plot_col, n_col) {
  group_by_(df, plot_col) %>%
    summarise_(sumabund=quote(sum(n_col)))
}

plotsum(df, plot_col='plot', n_col='n')
