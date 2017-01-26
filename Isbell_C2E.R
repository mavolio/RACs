library(vegan)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(reshape2)

set.seed(21)

#make rank-abundance relationship
x <- round(rexp(n=100,rate=0.1)) #100 species, rate affects total abundance
x <- x+1 #this keeps species with abundances that were rounded down to zero
plot(log(1:100),x[order(x,decreasing=TRUE)])
plot(1:100,x[order(x,decreasing=TRUE)])

df <- data.frame(sp=c(letters[1:25],LETTERS[1:25],paste(letters[1:25],LETTERS[1:25],sep=""),paste(LETTERS[1:25],letters[1:25],sep="")),x=x[order(x,decreasing=TRUE)])

sp2 <- sample(df$sp,50) #sample half the 'metacommunity' for the control
ctl <- df[df$sp%in%sp2,]
names(ctl)[2] <- "ctl"
df <- merge(df,ctl,all.x=T)
df <- df[order(df$x,decreasing=T),]

#first consider a range of species losses and gains
for(i in 1:45) {
	gns <- sample(df$sp[is.na(df$ctl)],i) #sample new species to gain: h
	df$gn <- df$x
	df$gn[!df$sp%in%c(as.character(gns),as.character(df$sp[!is.na(df$ctl)]))] <- NA
	names(df)[ncol(df)] <- paste("gn",i,sep=".")
	
	lss <- sample(df$sp[!is.na(df$ctl)],i) #sample new species to lose: f
	df$ls <- df$ctl
	df$ls[df$sp%in%c(as.character(lss))] <- NA
	names(df)[ncol(df)] <- paste("ls",i,sep=".")
}

df2 <- as.data.frame(t(as.matrix(df[,2:ncol(df)])))
colnames(df2) <- paste("sp",1:100,sep="")
df3 <- apply(df2,MARGIN = 1,FUN=as.logical)
df2$rich <- colSums(df3,na.rm=T)
df2$comm <- colnames(df[,2:ncol(df)])
df2$N <- rowSums(df2[names(df2)%in%paste("sp",1:100,sep="")],na.rm=TRUE)

#next consider high species evenness, lower dominance
df4 <- df2
#for(i in 1:100) {df4[!is.na(df4[,i]),i] <- round(df4$N/df4$rich)[!is.na(df4[,i])]}
for(i in 1:100) {df4[!is.na(df4[,i]),i] <- 10}
df4$even <- "max"
df2$even <- "none"
df5 <- rbind(df2,df4)
df44 <- df2
for(i in 1:30) {df44[,i] <- df44[,i]-round(((df44[,i]-10)/2))}
for(i in 33:100) {df44[,i] <- df44[,i]+round(((10-df44[,i])/2))}
df44$even <- "medium"
df5 <- rbind(df5,df44)

df4a <- df2
for(i in 1:30) {df4a[,i] <- df4a[,i]-round(((df4a[,i]-10)/4))}
for(i in 33:100) {df4a[,i] <- df4a[,i]+round(((10-df4a[,i])/4))}
df4a$even <- "low"
df5 <- rbind(df5,df4a)

df4b <- df2
for(i in 1:30) {df4b[,i] <- df4b[,i]-round(((df4b[,i]-10)/(4/3)))}
for(i in 33:100) {df4b[,i] <- df4b[,i]+round(((10-df4b[,i])/(4/3)))}
df4b$even <- "high"
df5 <- rbind(df5,df4b)

#next consider reordering within levels of richness and evenness
df7 <- df5 #shuffle ranks
for(i in 3:nrow(df7)) {df7[i,1:100][!is.na(df7[i,1:100])] <- sample(df7[i,1:100][!is.na(df7[i,1:100])],1,replace=F)}
#rbind(df5[3,],df7[3,])
df7$rank.shift <- "random"

df8 <- df5 #reverse ranks
for(i in 3:nrow(df8)) {df8[i,1:100][!is.na(df8[i,1:100])] <- sort(df8[i,1:100][!is.na(df8[i,1:100])])}
#rbind(df5[3,],df8[3,])
df8$rank.shift <- "high"

df9 <- df5 #minor tweak to ranks: shuffle abundances for top three species
for(i in 3:nrow(df9)) {df9[i,1:100][!is.na(df9[i,1:100])][1:3] <- df9[i,1:100][!is.na(df9[i,1:100])][sample(1:3,3)]}
#rbind(df5[3,],df9[3,])
df9$rank.shift <- "low"

df5$rank.shift <- "none"

#consolidate fake communities
d10 <- rbind(df5,df7)
d10 <- rbind(d10,df8)
d10 <- rbind(d10,df9)

d10$chEven <- factor(d10$even,levels=c("none","low","medium","high","max"))
d10$rank.shift <- factor(d10$rank.shift,levels=c("high","random","low","none"))

d10[is.na(d10)] <- 0
d10$bc <- as.matrix(vegdist(d10[,1:100]))[,2]
d10$S <- specnumber(d10[,1:100])
d10$InvSimpEven <- diversity(d10[,1:100],"inv")/d10$S


cg <- c(brewer.pal(7,"RdBu")[c(1,2,6,7)])
  
ggplot(d10,aes(x=S,y=bc,group=interaction(rank.shift,even),color=rank.shift,shape=chEven))+ 
  geom_point() + scale_colour_manual(values=cg, name="Rank Shift") + scale_shape_manual(values=c(19,6,5,2,1), name="Evenness Change") + scale_x_continuous("Richness") + scale_y_continuous("Bray-Curtis Dissimilarity") +
  theme_bw() + theme(panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_rect(colour="black"),plot.title=element_text(hjust=0.5)) +labs(title="Simulated Data")

pdf("CommChangeDiss.pdf",width=10, height=6)
ggplot(d10,aes(x=S,y=bc,group=interaction(rank.shift,chEven),color=rank.shift,shape=chEven))+ 
  geom_point() + facet_grid(rank.shift~chEven) + scale_colour_manual(values=cg, name="Rank Shift") + scale_shape_manual(values=c(19,6,5,2,1), name="Evenness Increase") + scale_x_continuous("Richness") + scale_y_continuous("Bray-Curtis Dissimilarity") +
  theme_bw() + theme(panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_rect(colour="black"),plot.title=element_text(hjust=0.5)) +labs(title="Simulated Data")
dev.off()


d10$InvSimpEven <- round(d10$InvSimpEven,1)
d.m <- melt(data=d10[d10$rank.shift=="none",],id.vars=c("S","InvSimpEven"),measure.vars="bc")
d.c <- dcast(d.m,S~InvSimpEven,fun.aggregate=mean,na.rm=T)
rownames(d.c) <- d.c[,1]
d.c <- as.matrix(d.c[,3:ncol(d.c)])
plot_ly(z= ~d.c, type="surface")

