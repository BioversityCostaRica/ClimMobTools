####################################################
###### MODELLING CROWDSOURCING DATA 
#Kaue de Sousa
#<k.desousa(at)cgiar.org>
#Bioversity International - Costa Rica
#First run 22Oct2017
#Updated 15Mar2017
####################################################

library(tidyverse)
library(lubridate)
library(car)
library(caret)
library(RCurl)
library(PlackettLuce)
library(psychotree) #
library(partykit) #
library(e1071)
library(densratio)
library(foreach)
library(doParallel)
library(abind)

#add ClimMob tools for data learning
tools <- RCurl::getURL("https://raw.githubusercontent.com/kauedesousa/ClimMobTools/pilot/ClimMobLearning.R", ssl.verifypeer = FALSE)
eval(parse(text = tools))
#cores for parallel
n_cpu <- 3
#working directory
#wd <- "C:/Users/Jacob Van Etten/Dropbox/vanEtten_etal2018_replication_data/"
wd <- "C:/Users/KAUE/Dropbox (Bioversity CR)/vanEtten_etal2018_replication_data/"
#wd<-"/home/tucurrique2/ALLAN/vanEtten_etal2018_replication_data/"
setwd(wd)

#create folder for outputs
output <- "./output/"
dir.create(output, showWarnings = FALSE)

output <- "./output/commonbean_NIC/"
dir.create(output, showWarnings = FALSE)


#read data
mydata <- read_csv("./input/NIC/climmob_nicaragua.csv")

#select rows for analysis
keep <- !is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c)

mydata <- mydata[keep,]

#soil and season as factor
mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], as.factor)

#add new location time and spatial covariates 
mydata$xy <- mydata$lon + mydata$lat
mydata$yx <- mydata$lon - mydata$lat

#the planting day (decimal day of the year)
mydata$planting_day <- as.integer(format(mydata$planting_date,'%j'))

#CHECK correlation between overall performance and yield
overall <- mydata[,c(2:4,13:15)]


yield <- mydata[,c(2:4,20:22)]
names(yield) <- names(overall)


#remove NA's in yield
keep <- !is.na(yield$best) & !is.na(yield$worst)
#keep those rows with no NA's
yield <- yield[keep,]
overall <- overall[keep,]

#get rankings for yield
yield <- grouped_rankings(as.PL(yield, local = F, additional.rank = F ), seq_len(nrow(yield)))
yield <- yield[1:length(yield),, as.grouped_rankings = FALSE]

#get rankings for overall
overall <- grouped_rankings(as.PL(overall, local = F, additional.rank = F ), seq_len(nrow(overall)))
overall <- overall[1:length(overall),, as.grouped_rankings = FALSE]

#calculate Kendall's correlation
#export this output
capture.output(get_tau(yield, overall)[1] * -1 , file = paste0(output,"Kendall_tau_yield_vs_overall.txt"))


#This analysis will evaluate the performance of Plackett-Luce models comparing whith follwing models
# PLT-null  without covariates null model 
# PLT-design  with covariates related to experiental design (season, planting date, location) 
# PLT-climate  with climatic covariates
# PLT-climdesign with covariates from design and climate
approach <- c("PL-null","PL-design","PL-climate","PLT-climdesign")

### PCA to identify clusters in covariates
covar <- cbind(mydata [, c("season","soil","xy","yx","lon","lat","planting_day","planting_date")],
               mydata [, c(28:209)])

#identify covariates with near 0 variance
varout <- caret::nearZeroVar(covar)
#remove those whith near 0 variance
covar <- covar[,-varout]


#apply BoxCox transformation
trans <- caret::preProcess(covar[-c(1)],
                           method = c("BoxCox", "center",
                                      "scale", "pca"))
#get principal components predictions from preprocessed object "trans"
mypca <- predict(trans, newdata = covar[-c(1)])

# #plot PCA
# ggplot(mypca) +
#    geom_point(aes(x = PC1, y = PC2, colour = covar$season, shape = covar$season ))

### Run PL model using covariates shifts as weights 
#add all elements for modelling, the principal components, rankings and covariates 
mydata <- cbind(mypca[,c("PC1","PC2")], mydata[, c(2:4,13:19) ], covar)
n <- nrow(mydata) 


#generate PL grouped rankings
G <- grouped_rankings(as.PL(mydata, local = T, additional.rank = T ), rep(seq_len(n), 4))
#add to mydata 
mydata <- cbind(G, mydata)

## Perform a covariate shifts analysis
#define folds
folds <- as.integer(as.factor(as.character(mydata$season)))
k <- max(folds) #number of folds

#get weights using principal components for the training data
myweights <- matrix(NA, nrow  = n, ncol = k, dimnames = list(seq_len(n), NULL))
#run over folds
for(i in 1:k){
  fold_i <- i
  train <- mydata[folds != fold_i,]
  test <- mydata[folds == fold_i,]
  
  #estimate density ratio by Kullback-Leibler Importance Estimation Procedure (KLIEP)
  i.dens <- densratio::KLIEP(train[,c("PC1","PC2")], test[,c("PC1","PC2")], verbose = FALSE)
  #compute density ratio
  i.weights <- i.dens$compute_density_ratio(train[,c("PC1","PC2")])
  
  # plot(density(i.weights_lambda))
  # plot(density(i.weights))
  
  #add weights to list
  #myweights[[i]] <- i.weights
  myweights[row.names(train) ,i] <- i.weights
}

rm(test,train)


#calculate the weights represented by observations per season 
wseason <- as.vector(summary(mydata$season))
wseason <- sqrt(wseason/n)
wseason <- wseason/sum(wseason)

# ## Perform a forward variable selection on PL-climate covariates
# to be able to predict a PLtree object from a null model I will create a NULL variable and add to mydata 
mydata$P1 <- 1 

#run over climate covariates and select those who contribute to improve PLT-null 
climcov <- names(mydata)[22:170]
#remove covariates related to the firts quarter of crop timespan
climcov <- climcov[ !grepl("_1st", climcov) ]
climcov <- climcov[ !grepl("DT", climcov) ]
climcov <- climcov[ !grepl("GDD", climcov) ]
lambda <- c(0, 0.035) #seq(0.030, 0.040, 0.001) #to smooth weigths


#climcov <- c("maxNT_2nd","minNT_2nd","P1","R10mm_4th","SDII_2nd_anomaly")

#combinations lambda x covariates 
CLC <- rbind(rep(climcov, each = length(lambda)), rep(lambda, times = length(climcov)))

#model parameters
minsize <- 150 #pltree 
bonferroni <- TRUE #pltree
alpha <- 0.01 #pltree
npseudo <- 0.5 #pltree
weights <- myweights

#AIC from null model (used as baseline)
par_n <-  0
#vector to keep best covariates
var_keep <- NULL
#if TRUE the routine will keep running, this parameters is updated at the end of each "while" routine
best <- TRUE
#keep number of runs
runs <- FALSE
#vector with best parameters (loglik, vars and lambda) in each run
best_parameters <- NULL

rm(covar,mypca,trans)

#Create cluster to do parallel computation to speed things up
cluster <- makeCluster(n_cpu)
registerDoParallel(cluster)
getDoParWorkers()

#models with weights
while(best){
  
  cat("Starting run ", sum(runs)+1, ". Time: ", date(), "\n")
  
  nCLC <- ncol(CLC)
  
  #get predictions for test from nodes and put in matrix (foreach)
  models <- foreach(i = 1:nCLC,
                       .combine = acomb,
                       .packages = c("PlackettLuce","psychotree"),
                       .export = ls(globalenv())) %dopar% (
                         f1(formula = as.formula(paste0("G ~ ", paste(unique(c(var_keep, CLC[1,i])), collapse = " + "))), 
                            d = mydata, 
                            k = k, folds = folds, minsize = minsize, 
                            alpha = alpha, bonferroni = bonferroni,
                            weights = weights, npseudo = npseudo, lambda = as.numeric(CLC[2,i]))
                         )

  #calculate Akaike weights along covariates per season
  #get weights in a transposed matrix 
  AW <- as.matrix(models[,1:k])
  AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))
  AW <- apply(AW, 2, function(X) akaike.weights(X)[[3]] )
  
  #then take the weighted average per season Stouffer's method with weights
  for(i in seq_len(k)) AW[,i] <- AW[,i] * wseason[i]
  
  #take the Stouffer's weighted mean 
  meanAW <- rowSums(AW)
  
  #take a simple mean of inernodes generated by each model
  nodes <- models[,c((k+1):(k+k))]
  nodes <- apply(nodes, 2, as.integer )
  nodes <- rowMeans(nodes)
  
  #take Kendall tau
  tau <- models[,(ncol(models)-1)]
  
  #take call from each model
  call <- models[,ncol(models)]
  
  #take maximum parameter from Akaike weights
  par_max <- max(meanAW)
  
  #take the position of par_max in the CLC matrix
  index_par_max <- which.max(meanAW)
  
  #Is par_max best (higher) than par_n?
  best <- par_max >= par_n #if not, stop
  
  #take the name of best variable
  best_var <- CLC[1, index_par_max]
  
  #if best, save the outputs 
  if(best){
    
    runs <- c(runs, best)
    
    #take the lambda used to smooth weights and obtain par_max
    best_lambda <- CLC[2, index_par_max]
    
    #keep parameters from this run in a .csv file 
    out <- as_tibble(cbind(call, CLC[2,], models[,c(1:k)], AW, meanAW = meanAW, inernodes = nodes, tau = tau))
    colnames(out) <- c("call", "lambda",paste0(rep("Deviance",k), 1:k), paste0(rep("AkaikeWeight",k), 1:k), "meanAW", "inernodes", "Kendall_tau")
    
    #remove COMBINATIONS of this covariable in CLC
    CLC <- CLC[ , - which(best_var == CLC[1,]) ]
    
    #remove null var from the first run, no longer necessary 
    CLC <- CLC[, !grepl("P1", CLC[1,])]
    
    #keep this model for the next run
    var_keep <- c(var_keep, best_var)
    
    #add best model to the next round 
    CLC <- cbind(CLC, matrix(rbind(paste(var_keep, collapse = ' + '), best_lambda)))
    
    #change the base for par_n (minimun accepted value)
    par_n <- par_max
    
    #keep the best parameters form this run
    best_parameters <- rbind(best_parameters, cbind(par_max, best_lambda, covars = toString(var_keep)))
    
    write.csv(best_parameters, paste0(output, "best_parameters.csv" )  )
    write.csv(out, paste0(output, "model_parameters_run", sum(runs),".csv" )  )
    cat(" #################################################################### \n Run ", sum(runs), "\n Lambda: ",  best_lambda,"\n Covars:", toString(var_keep), 
        "\n #################################################################### \n\n" )
  }
  
}

#return best model
best_model <- var_keep
best_lambda <- as.numeric(best_parameters[nrow(best_parameters), 2])

cat("\n\n END Forward Selection. Time:", date() )

#do forward without weights 
CLC <- as.matrix(rbind(climcov, rep(0, times = length(climcov))), dimnames = list(1:2, NULL))
par_n <-  0
var_keep <- NULL
best <- TRUE
runs <- FALSE
lambda <- 0

while(best){
  
  cat("Starting run without weigths ", sum(runs)+1, ". Time: ", date(), "\n")
  
  nCLC <- ncol(CLC)
  
  #get predictions for test from nodes and put in matrix (foreach)
  models <- foreach(i = 1:nCLC,
                    .combine = acomb,
                    .packages = c("PlackettLuce","psychotree"),
                    .export = ls(globalenv())) %dopar% (
                      f1(formula = as.formula(paste0("G ~ ", paste(unique(c(var_keep, CLC[1,i])), collapse = " + "))), 
                         d = mydata, 
                         k = k, folds = folds, minsize = minsize, 
                         alpha = alpha, bonferroni = bonferroni,
                         weights = weights, npseudo = npseudo, lambda = lambda)
                    )
  
  #calculate Akaike weights along covariates per season
  #get weights in a transposed matrix 
  AW <- as.matrix(models[,1:k])
  AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))
  AW <- apply(AW, 2, function(X) akaike.weights(X)[[3]] )
  
  #then take the weighted average per season Stouffer's method with weights
  for(i in seq_len(k)) AW[,i] <- AW[,i] * wseason[i]
  
  #take the Stouffer's weighted mean 
  meanAW <- rowSums(AW)
  
  #take a simple mean of inernodes generated by each model
  nodes <- models[,c((k+1):(k+k))]
  nodes <- apply(nodes, 2, as.integer )
  nodes <- rowMeans(nodes)
  
  #take Kendall tau
  tau <- models[,(ncol(models)-1)]
  
  #take call from each model
  call <- models[,ncol(models)]
  
  #take maximum parameter from Akaike weights
  par_max <- max(meanAW)
  
  #take the position of par_max in the CLC matrix
  index_par_max <- which.max(meanAW)
  
  #Is par_max best (higher) than par_n?
  best <- par_max >= par_n #if not, stop
  
  #take the name of best variable
  best_var <- CLC[1, index_par_max]
  
  #if best, save the outputs 
  if(best){
    
    runs <- c(runs, best)
    
    #keep parameters from this run in a .csv file 
    out <- as_tibble(cbind(call, CLC[2,], models[,c(1:k)], AW, meanAW = meanAW, inernodes = nodes, tau = tau))
    colnames(out) <- c("call", "lambda",paste0(rep("Deviance",k), 1:k), paste0(rep("AkaikeWeight",k), 1:k), "meanAW", "inernodes", "Kendall_tau")
    
    #remove COMBINATIONS of this covariable in CLC
    CLC <- CLC[ , - which(best_var == CLC[1,]) ]
    
    #remove null var from the first run, no longer necessary 
    CLC <- CLC[, !grepl("P1", CLC[1,])]
    
    #keep this model for the next run
    var_keep <- c(var_keep, best_var)
    
    #add best model to the next round 
    CLC <- cbind(CLC, matrix(rbind(paste(var_keep, collapse = ' + '), best_lambda)))
    
    #change the base for par_n (minimun accepted value)
    par_n <- par_max
    
    write.csv(out, paste0(output, "model_parameters_run", sum(runs),"_NO_WEIGHTS.csv" )  )
    cat(" #################################################################### \n Run ", sum(runs), "\n Lambda: ",  0,"\n Covars:", toString(var_keep), 
        "\n #################################################################### \n\n" )
  }
  
}

stopCluster(cluster)

#get outputs from each model approach in blocked cross-validation
blockedfolds <-  matrix(NA, nrow = length(approach), ncol = 4, dimnames = list(1:length(approach),c("Model","Call", "Ktau","Deviance")))
blockedfolds[,1] <- approach

#add PLT-null already calculated above
PLTnull <- unlist(crossvalidation_PLTE(as.formula(c("G ~ P1")), d = mydata, k = k, folds = folds, minsize = minsize, 
                                    alpha = alpha, bonferroni = bonferroni))

#cross-validation using climatic covariates
weights <- myweights

blockedfolds[3,2:4] <- unlist(crossvalidation_PLTE(as.formula(paste0(c("G ~ "), paste(best_model, collapse = " + "))),
                                                   d = mydata, k = k, folds = folds, minsize = minsize,
                                                   alpha = alpha, bonferroni = bonferroni, weights = weights, lambda = best_lambda)[1:3])

#PLT-designclim
blockedfolds[4,2:4] <- unlist(crossvalidation_PLTE(as.formula(paste0(c("G ~ soil + xy + yx + lon + lat +"), paste(best_model, collapse = " + "))),
                                                   d = mydata, k = k, folds = folds, minsize = minsize, 
                                                   alpha = alpha, bonferroni = bonferroni, weights = weights, lambda = best_lambda)[1:3])

print(blockedfolds)
write.csv(blockedfolds, paste0(output, "performance_PLT_blocked_crossvalidation.csv"))


#Do the same in a 10-fold cross validation 
kfolds <-  matrix(NA, nrow = length(approach), ncol = 4, dimnames = list(1:length(approach),c("Model","Call", "Ktau","Deviance")))
kfolds[,1] <- approach

#number of folds
k <- 10
#take sample
folds <- sample(rep(1:k, times = ceiling(n/k), length.out=n), replace=FALSE)

#PLT-null 
kfolds[1,2:4] <- unlist(crossvalidation_PLTE(as.formula(c("G ~ P1")),
                                             d = mydata, k = k, folds = folds, minsize = minsize, 
                                             alpha = alpha, bonferroni = bonferroni)[1:3])

#PLT-design
kfolds[2,2:4] <- unlist(crossvalidation_PLTE(as.formula(c("G ~ season + soil + xy + yx + lon + lat + planting_day + planting_date")),
                                      d = mydata, k = k, folds = folds, minsize = minsize, 
                                      alpha = alpha, bonferroni = bonferroni)[1:3])


#PLT-climate
#calculate new weights
myweights <- matrix(NA, nrow  = n, ncol = k, dimnames = list(seq_len(n), NULL))
#run over folds
for(i in 1:k){
  fold_i <- i
  train <- mydata[folds != fold_i,]
  test <- mydata[folds == fold_i,]
  
  #estimate density ratio by Kullback-Leibler Importance Estimation Procedure (KLIEP)
  i.dens <- densratio::KLIEP(train[,c("PC1","PC2")], test[,c("PC1","PC2")], verbose = FALSE)
  #compute density ratio
  i.weights <- i.dens$compute_density_ratio(train[,c("PC1","PC2")])
  #add to matrix
  myweights[row.names(train) ,i] <- i.weights
}

weights <- myweights

kfolds[3,2:4] <- unlist(crossvalidation_PLTE(as.formula(paste0(c("G ~ "), paste( best_model, collapse = " + "))),
                                             d = mydata, k = k, folds = folds, minsize = minsize, 
                                             alpha = alpha, bonferroni = bonferroni, weights = weights, lambda =  best_lambda)[1:3])

#PLT-designclim
kfolds[4,2:4] <- unlist(crossvalidation_PLTE(as.formula(paste0(c("G ~ soil + xy + yx + lon + lat +"), paste( best_model, collapse = " + "))),
                                             d = mydata, k = k, folds = folds, minsize = minsize, 
                                             alpha = alpha, bonferroni = bonferroni, weights = weights, lambda =  best_lambda)[1:3])


print(kfolds)
write.csv(kfolds, paste0(output, "performance_PLT_10fold_crossvalidation.csv"))

#calculate wheights for the entire dataset
myweights <- matrix(NA, nrow  = n, ncol = k, dimnames = list(seq_len(n), NULL))
#run over folds
for(i in 1:k){
   fold_i <- i
   #estimate density ratio by Kullback-Leibler Importance Estimation Procedure (KLIEP)
   i.dens <- densratio::KLIEP(mydata[folds != fold_i , c("PC1","PC2")], mydata[folds == fold_i , c("PC1","PC2")], verbose = TRUE)
   #compute density ratio
   i.weights <- i.dens$compute_density_ratio(mydata[,c("PC1","PC2")])

   myweights[,i] <- i.weights
}

myweights <- (rowMeans(myweights, na.rm=TRUE)) ^ best_lambda

tree <- pltree(as.formula(paste0(c("G ~ "), paste(best_model, collapse = " + "))),
                data = mydata, weights = myweights, alpha = alpha, minsize = minsize, bonferroni = bonferroni)


capture.output(tree, summary(tree), file = paste0(output, "pltree_clim.txt") )


#extract the terminal nodes
terminal_nodes <- partykit::nodeids(tree, terminal=TRUE)

for (j in 1:length(terminal_nodes)){
  qvcoeff <- qvcalc(tree[[terminal_nodes[j]]]$node$info$object)$qvframe
  nobs <- tree[[terminal_nodes[j]]]$node$info$nobs
  qvcoeff$player <- factor(rownames(qvcoeff), levels = rownames(qvcoeff))
  #plot estimates
  plot <- ggplot(qvcoeff, aes(x=estimate, y=player)) +
    geom_point(pch=21, size=3,fill="black",colour="black") +
    #geom_line(aes(group = 1)) +
    geom_errorbarh(aes(xmin=estimate-(quasiSE),
                       xmax=estimate+(quasiSE)),
                   colour="black",height = .2) +
    theme_bw() +
    labs(x = NULL, y = NULL, title=paste("Node ", terminal_nodes[j], " (n= ", n, ")", sep="")) +
    theme(plot.title = element_text(size=17),
          axis.text = element_text(size=12, colour="black"),
          legend.text= element_text(size=12, colour="black"),
          axis.text.x = element_text(size=17, angle = 0, hjust=1, vjust=1, face="plain"),
          axis.text.y = element_text(size=17, angle = 0, hjust=1, vjust=1, face="plain"),
          plot.background = element_blank())
  #export plot
  ggsave(paste0(output, "pltree_qvcalc_node_", terminal_nodes[j] , ".png"), plot = plot, dpi = 300,
         width = 20, height = 20, units = "cm")
}

