library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

setwd("~/Dropbox/SESYNC/SESYNC_RACs/For R Package")

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


