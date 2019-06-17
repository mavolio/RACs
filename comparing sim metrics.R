library(tidyverse)

setwd("C:\\Users\\megha\\Dropbox\\Manuscripts\\RACs\\Submit Ecosphere\\Revision")

change<-read.csv("sim_change_metrics.csv")
diff<-read.csv("sim_diff_metrics.csv")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 1),
                   symbols = c("*", " "))
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=5, col="red")
}        


pairs(change[,c(7:8,11:14)], font.labels=0.5, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)


pairs(diff[,c(8:12)], font.labels=0.5, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)