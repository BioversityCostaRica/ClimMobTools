####################################################
###### MODELLING TRICOT DATA 
#Kauê de Sousa
#<k.desousa(at)cgiar.org>
#Bioversity International - Costa Rica
#First run 22Oct2017
#Updated 07DecNov2017
####################################################
############################################################################
#install.packages(c("parallel","ggplot2","psychotools", "grid", "partykit" ,"psychotree",
#                   "devtools" ,"igraph","RSpectra","PlackettLuce","car","vegan","lbfgs"))
library(tidyverse)
library(car)
library(vegan)
library(psychotools)
library(grid)
library(partykit)
library(psychotree)
library(devtools)
library(igraph)
library(RSpectra)
library(PlackettLuce)
library(lbfgs)
library(parallel)
library(ggplot2)

#set work directory
#mypath <- "/home/tucurrique2/allan/Documents/"#"C:/Users/Suelen/Dropbox/"
#wd <- "./bioversity_k_j/"
mypath <- "C:/Users/KAUE/Dropbox (Bioversity CR)/"
wd <- "./vanEtten_etal2018_replication_data"
setwd(paste0(mypath, wd))

#load extra functions
source("./script/ClimMobTools.R")
ls()
#read data set
mydata <- readr::read_csv("./input/climmob_ethiopia_19Dec2017.csv")
R <- readr::read_csv("./input/climmob_ethiopia_rankigs.csv")
dim(mydata)
# #create output directory
# output <- "./output/dorum_wheat_ETH/"
# dir.create(output, showWarnings = F)

#check consistency of paired evaluations
#remove wrong paired evaluations best == worst 
n=nrow(mydata)

mydata$Soil <- factor(mydata$Soil, levels = sort(unique(mydata$Soil)))
#mydata$season <- factor(mydata$season, levels = sort(unique(mydata$season)))
mydata$planting_date <- as.integer(mydata$planting_date)

#predict() is no applicable to an object of class "PlackettLuce"
#generate two pseudo-atrributes to run PLtree with no effect to results in null model
mydata$p1 <- 1
mydata$p2 <- 0
#add new location attributes for design model
mydata$xy <- mydata$lon + mydata$lat
mydata$yx <- mydata$lon - mydata$lat
names(mydata)

#verify if there is multicolinearity
#calculate diversity of attributes
names(mydata)
div <- as.data.frame(mydata[,c(20:104)])
div[is.na(div)] <- 0
div[1:ncol(div)] <- lapply(div[1:ncol(div)], as.numeric)
summary(div)
div <- as.matrix(vegan::diversity(div, index = "shannon", MARGIN = 2))
div[,1] <- exp(div[,1])
div <- as.matrix(div[div[,1] > max(div[,1]) * .5, ])

# mod <- lm(as.formula(paste0("GDD ~ " , paste(cbind(dimnames(div)[[1]]),collapse = " + "))),
#           data = mydata)
# keep <- as.vector(names(subset(mod$coefficients, !is.na(mod$coefficients))))[-1]
# 
# div <- div[keep,]
# 
# mod <- lm(as.formula(paste0("GDD ~ " , 
#                             paste(cbind(names(div)),collapse = " + "))),
#           data = mydata)
# vif <- vif(mod)

#define list of attributes variables ####
attrib <- list(null = c("p1","p2"),
               desi = c("lat","lon","xy","yx"),
               clim = c("Soil", dimnames(div)[[1]]))

mydata <- merge(mydata, R, by="Farmer_No")

items <- names(R)[-1]

R <- as.matrix(mydata[,items])
colnames(R) <- items

G <- grouped_rankings(as.rankings(R), seq_len(n))




############
n = nrow(mydata)
mydata <- mydata[,as.vector(unlist(attrib))]
#Make sure to put the response variable in the dataset
mydata$G <- G

tree <- pltree(as.formula(paste0("G ~ " , paste(cbind(attrib[[3]]),collapse = " + "))),
               data = mydata, minsize = minsize, alpha = alpha)
plot(tree)



#Create cluster to do parallel computation to speed things up
cluster <- makeCluster(2)
registerDoParallel(cluster)
getDoParWorkers()

# k-fold cross-validation
kfolds <- list()
for(i in c(1:3)){
  kfolds[[i]] <-  crossvalidation_PLTE(G ~ ., 
                                       d = mydata[c("G", attrib[[i]])], 
                                       k = 10, 
                                       folds = NULL, 
                                       minsize = 30, 
                                       alpha = 0.05, 
                                       npseudo = 0.5, 
                                       ntrees = 25, 
                                       cluster = cluster)
}

kfolds[[1]]$average_tau
kfolds[[2]]$average_tau
kfolds[[3]]$average_tau

#blocked cross-validation
#Make the folds first and calculate k
folds <- as.integer(mydata$season)
k <- max(folds)

blocked <- list()
for(i in c(1:3)){
  blocked[[i]] <- crossvalidation_PLTE(G ~ ., 
                                       d = mydata[c("G", attrib[[i]])], 
                                       k = k, 
                                       folds = folds, 
                                       minsize = 30, 
                                       alpha = 0.05, 
                                       npseudo = 0.5, 
                                       ntrees = 25, 
                                       cluster = cluster)
  
}

blocked[[1]]$average_tau
blocked[[2]]$average_tau
blocked[[3]]$average_tau

save(blocked, kfolds, mydata, G, file = "./output/ethiopia_wheat.RData")


stopCluster(cluster)



stopCluster(cluster)












