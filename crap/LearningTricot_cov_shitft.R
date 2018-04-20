####################################################
###### MODELLING CROWDSOURCING DATA 
#Kaue de Sousa <k.desousa(at)cgiar.org> & Jacob van Etten <j.vanetten(at)cgiar.org>
#Bioversity International - Costa Rica
#First run 22Oct2017
#Updated 15Mar2017
####################################################

library(tidyverse)
library(car)
library(caret)
library(RCurl)
library(PlackettLuce)
library(psychotree)
library(partykit)
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

#read data
NIC <- read_csv("./input/NIC/climmob_nicaragua.csv")
IND <- read_csv("./input/IND/climmob_india.csv")
#ETH <- ""

#bind all datasets
df <- as_tibble(rbind(NIC, IND)) #,ETH))
#put crop as first column
df <- df[,c(207,206, 1:205 )]
rm(IND, NIC)
#crop that will be analysed here 
crop <- c("commonbean","wheat") #,"durumwheat"))

#This analysis will evaluate the performance of Plackett-Luce models comparing whith following models
# PLT-null  without covariates null model 
# PLT-design  with covariates related to experiental design (season, planting date, location) 
# PLT-season with season as covariate
# PLT-climate  with climatic covariates
# PLT-clim+spatial with covariates from design and climate
approach <- c("PLT-null","PLT-design","PLT-season","PLT-climate","PLT-clim+spatial")

#run over crops
for(m in 1:length(crop)){
  
  cat("################################## \n Starting analysis:", toupper(crop[m]), " \n Time:", date())
  
  mydata <- df[df$crop == crop[m] , ]
  
  #create folder for outputs
  output <- "./output/"
  dir.create(output, showWarnings = FALSE)
  
  #an output for the crop 
  output <- paste0(output,unique(mydata$crop),"_",unique(mydata$country),"/")
  dir.create(output, showWarnings = FALSE)
  
  #remove possible storaged factors levels in the subseted data
  ##soil and season as factor
  mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], as.character)
  mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], as.factor)
  
  #select rows for analysis
  #deal with commonbean dataset
  if(crop[m] == "commonbean"){
    #remove NA's in rankings for overall_performance vs local
    keep <- !is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c)
    mydata <- mydata[keep,]
  }
  
  #deal with wheat dataset
  if(crop[m] == "wheat"){
    #remove wrong evaluations best == worst in overall_performance
    mydata <- mydata[mydata$best!=mydata$worst ,]
    
    #remove items (varieties) tested in only 1 season
    rm_item <- as.data.frame(table(mydata$variety_a, mydata$season))
    rm_item[rm_item > 0] <- 1
    rm_item <- aggregate(rm_item[,"Freq"], by=list(rm_item[,"Var1"]), sum)
    rm_item <- as.vector(c(as.character(rm_item[rm_item$x<2, 1 ]), "HUW468"))
    for(i in 1:length(rm_item)){
      mydata$variety_a <- ifelse(mydata$variety_a == rm_item[i], NA, mydata$variety_a)
      mydata$variety_b <- ifelse(mydata$variety_b == rm_item[i], NA, mydata$variety_b)
      mydata$variety_c <- ifelse(mydata$variety_c == rm_item[i], NA, mydata$variety_c)
    }
    #remove rankings with more than 1 NA per observer
    keep <- NULL
    for(i in 1:nrow(mydata)){
      if(sum(is.na(mydata[i,c("variety_a","variety_b","variety_c")])) > 1) keep[i] <- FALSE else keep[i] <- TRUE
    }
    mydata <- mydata[keep, ]
  }
  
  dim(mydata)
  
  #add new location time and spatial covariates 
  mydata$xy <- mydata$lon + mydata$lat
  mydata$yx <- mydata$lon - mydata$lat
  
  #the planting day (decimal day of the year)
  mydata$planting_day <- as.integer(format(mydata$planting_date,'%j'))
  
  #CHECK correlation between overall performance and yield
  #overall performance
  overall <- mydata[,c("variety_a","variety_b","variety_c","characteristic","best","worst")]
  #yield data 
  yield <- mydata[,c("variety_a","variety_b","variety_c","yield","best_yield","worst_yield")]
  #same names of overall in yield 
  names(yield) <- names(overall)
  
  #remove NA's in yield
  keep <- !is.na(yield$best) & !is.na(yield$worst) &  yield$best != yield$worst
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
  
  #PCA to identify clusters in covariates
  covar <- cbind(mydata [, c("season","planting_day","planting_date","soil","lon","lat","xy","yx")],
                 mydata [, c(26:207)])

  #identify covariates with near 0 variance
  varout <- caret::nearZeroVar(covar)
  #remove those whith near 0 variance
  covar <- covar[,-varout]
  
  #Here a covariate shift approach is applied.
  #The output of this analysis will be used as weights in PLT models 
  #First perform a principal component analysis.
  #Since anomalies implies in, some cases, negative values the normalisation of data is not possible with log() 
  #To overcome this issue we apply a BoxCox transformation
  trans <- caret::preProcess(covar[-c(1)],
                             method = c("BoxCox", "center",
                                        "scale", "pca"))
  #get principal components predictions from preprocessed object "trans"
  mypca <- predict(trans, newdata = covar[-c(1)])
  
  # #plot PCA
  # ggplot(mypca) +
  #    geom_point(aes(x = PC1, y = PC2, colour = covar$season, shape = covar$season ))
  
  #Add the first two principal component into the main dataset together with farmers rankings and covariates
  mydata <- cbind(mypca[,c("PC1","PC2")], 
                  mydata[, c(4:6, 12:18)],
                  covar)
  
  #keep the number of rows in this dataset
  n <- nrow(mydata) 
  
  # Perform a covariate shifts analysis
  #define folds
  folds <- as.integer(as.factor(as.character(mydata$season)))
  k <- max(folds) #number of folds
  
  
  cat("Calculating covariate shift in seasons \n")
  
  
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
  
  
  #Generate a Plackett-Luce grouped rankings
  #for common beans we add the comparison with the local var in partial rankings 
  #in this approach each comparison with local is a individual rank
  if(crop[m] == "commonbean"){
    G <- grouped_rankings(as.PL(mydata, local = TRUE, additional.rank = TRUE ), rep(seq_len(n), 4))
  }
  #information about local item is not available for the other datasets
  if(crop[m] == "wheat"){
    G <- grouped_rankings(as.PL(mydata, local = FALSE), seq_len(n))
  }
  
  #Add these rankings to the main dataset 
  mydata <- cbind(G, mydata)
  
  #Calculate the weights represented by seasons
  wseason <- as.vector(summary(mydata$season))
  wseason <- sqrt(wseason/n)
  wseason <- wseason/sum(wseason)
  
  # Perform a forward variable selection on PL-climate covariates
  # To be able to predict a PLtree object from a null model we create a NULL variable and add to the main dataset 
  mydata$P1 <- 1 
  
  #run over climate+spatial covariates and select those who contribute to improve PLT-null 
  climcov <- names(mydata)[17:ncol(mydata)]
  #remove covariates related to the firts quarter of crop timespan
  climcov <- climcov[ !grepl("_1st", climcov) & !grepl("DT", climcov) & !grepl("GDD", climcov) & !grepl("Rx1", climcov) ]
  #define values of lambda to be tested in this exercise
  if(crop[m] == "commonbean"){
    lambda <- seq(0.030, 0.040, 0.001)  
  }
  if(crop[m] == "wheat"){
    lambda <- seq(0.015, 0.035, 0.001)  
  }
  
  #climcov <- c("maxNT_2nd","minNT_2nd","P1","R10mm_4th","SDII_2nd_anomaly")#use this for tests
  
  #combinations lambda x covariates 
  CLC <- rbind(rep(climcov, each = length(lambda)), rep(lambda, times = length(climcov)))
  
  #Define parameters for the PL model 
  if(crop[m] == "commonbean"){
    minsize <- 200
  }
  if(crop[m] == "wheat"){
    minsize <- 2000
  }
  
  bonferroni <- TRUE
  alpha <- 0.01
  npseudo <- 0.5
  weights <- myweights
  
  
  #Define initial parameters for forward selection
  ##AIC from null model (used as baseline)
  par_n <-  0
  ##vector to keep best covariates
  var_keep <- NULL
  ##if TRUE the routine will keep running, this parameters is updated at the end of each "while" routine
  best <- TRUE
  ##keep number of runs
  runs <- FALSE
  ##vector with best parameters (loglik, vars and lambda) in each run
  best_parameters <- NULL
  
  #remove unused objects and reduce size of globalenv() in parallel export
  rm(covar,mypca,trans, i.weights, i.dens, G, keep, yield, overall)
  
  #Create cluster to do parallel computation to speed things up
  cluster <- makeCluster(n_cpu)
  registerDoParallel(cluster)
  getDoParWorkers()
  
  #Perform forward selection using weights and lambdas 
  while(best){
    
    cat("Starting Forward Selection using weigths. Run ", sum(runs)+1, "\n Time: ", date(), "\n")
    
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
  
  cat("End Forward Selection using weigths. \n Time: ", date(), "\n")
  
  #Perform forward without weights (lambda = 0)
  CLC <- as.matrix(rbind(climcov, rep(0, times = length(climcov))), dimnames = list(1:2, NULL))
  par_n <-  0
  var_keep <- NULL
  best <- TRUE
  runs <- FALSE
  lambda <- 0
  
  while(best){
    
    cat("Starting Forward Selection without weigths. Run ", sum(runs)+1, "\n Time: ", date(), "\n")
    
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
  
  cat("End Forward Selection without weigths. \n Time: ", date(), "\n")
  
  #Stop connection with cores
  stopCluster(cluster)
  
  #Define list of covariables to use in each approach
  attrib <- list(null = ("P1"),
                 desi = c("lat","lon","xy","yx","planting_date","planting_day","season","soil"),
                 seas = c("season"),
                 clim = best_model,
                 clsp = c(best_model, c("lat","lon","xy","yx","soil") ))
  
  
  #Get outputs from each model approach in blocked cross-validation
  blockedfolds <-  matrix(NA, nrow = length(approach), ncol = (6+(k*3)), 
                          dimnames = list(1:length(approach),
                                          c("Approach","Model","mean.Ktau","mean.AIC","mean.Deviance",
                                            paste0(rep("Ktau",times = k), 1:k),
                                            paste0(rep("AIC",times = k), 1:k),
                                            paste0(rep("Deviance",times = k), 1:k), "mean.AW")))
  
  #Add approach to the first column
  blockedfolds[,1] <- approach
  
  cat("Starting blocked cross-validation", date() ,"\n")
  #run blocked cross-validation over approaches 
  for(a in 1:length(approach)){
    
    formula_a <- as.formula(paste0("G ~ ", paste(attrib[[a]], collapse = " + ") ))
    
    #apply weights and lambda on PLT-clim and PLT-clim+spatial
    if(a >= 5) weights <- myweights else weights <- NULL
    if(a >= 5) lambda <- best_lambda else lambda <- 0
    
    blockedfolds[ a , c(2:(ncol(blockedfolds)-1))] <- unlist(crossvalidation_PLTE(formula_a, 
                                                              d = mydata, k = k, folds = folds, minsize = minsize, 
                                                              alpha = alpha, bonferroni = bonferroni, weights = weights, lambda = lambda,
                                                              verbose = FALSE)[c(1:7)])
  }
  
  #calculate Akaike weights along approaches over Deviance
  AW <- as.matrix(blockedfolds[,(ncol(blockedfolds)-(k)):(ncol(blockedfolds)-1)])
  AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))
  AW <- apply(AW, 2, function(X) akaike.weights(X)[[3]] )
  
  #then take the weighted average per season Stouffer's method with weights
  for(i in seq_len(k)) AW[,i] <- AW[,i] * wseason[i]
  
  #take the Stouffer's weighted mean 
  blockedfolds[,ncol(blockedfolds)] <- rowSums(AW)
  
  #write this matrix in a csv file 
  write.csv(blockedfolds, paste0(output, "performance_PLT_blocked_crossvalidation.csv"), row.names = FALSE)
  
  
  cat("Starting k-fold cross-validation",date(),"\n")
  
  #Do the same in a 10-fold cross validation 
  #Define number of folds
  k <- 10
  #Define samples
  folds <- sample(rep(1:k, times = ceiling(n/k), length.out=n), replace=FALSE)
  
  cat("Calculating new weights for blocked cross-validation \n")
  #Calculate new weights
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
  
  #matrix to keep the outputs 
  kfolds <-  matrix(NA, nrow = length(approach), ncol = (6+(k*3)), 
                    dimnames = list(1:length(approach),
                                    c("Approach","Model","mean.Ktau","mean.AIC","mean.Deviance",
                                      paste0(rep("Ktau",times = k), 1:k),
                                      paste0(rep("AIC",times = k), 1:k),
                                      paste0(rep("Deviance",times = k), 1:k), "mean.AW")))
  kfolds[,1] <- approach
  
  #run k-fold over approaches 
  for(a in 1:length(approach)){
    
    formula_a <- as.formula(paste0("G ~ ", paste(attrib[[a]], collapse = " + ") ))
    
    #apply weights and lambda on PLT-clim and PLT-clim+spatial
    if(a >= 5) weights <- myweights else weights <- NULL
    if(a >= 5) lambda <- best_lambda else lambda <- 0
    
    kfolds[ a , c(2:(ncol(kfolds)-1))] <- unlist(crossvalidation_PLTE(formula_a, 
                                                              d = mydata, k = k, folds = folds, minsize = minsize, 
                                                              alpha = alpha, bonferroni = bonferroni, weights = weights, lambda = lambda,
                                                              verbose = FALSE)[c(1:7)])
  }
  
  #calculate Akaike weights along approaches over Deviance
  AW <- as.matrix(kfolds[,(ncol(kfolds)-(k)):(ncol(kfolds)-1)])
  AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))
  AW <- apply(AW, 2, function(X) akaike.weights(X)[[3]] )
  
  #Calculate the weights represented by folds
  wfold <- as.vector(summary(as.factor(folds)))
  wfold <- sqrt(wfold/n)
  wfold <- wfold/sum(wfold)
  
  #then take the weighted average per season Stouffer's method with weights
  for(i in seq_len(k)) AW[,i] <- AW[,i] * wfold[i]
  
  #take the Stouffer's weighted mean 
  kfolds[,ncol(kfolds)] <- rowSums(AW)
  
  #write csv file 
  write.csv(kfolds, paste0(output, "performance_PLT_10fold_crossvalidation.csv"), row.names = FALSE)
  
  
  cat("Fit pltree with best climate model \n")
  
  #Calculate wheights for the entire dataset and fit a pltree
  myweights <- matrix(NA, nrow  = n, ncol = k, dimnames = list(seq_len(n), NULL))
  #run over folds
  for(i in 1:k){
    fold_i <- i
    #estimate density ratio by Kullback-Leibler Importance Estimation Procedure (KLIEP)
    i.dens <- densratio::KLIEP(mydata[folds != fold_i , c("PC1","PC2")], mydata[folds == fold_i , c("PC1","PC2")], verbose = FALSE)
    #compute density ratio
    i.weights <- i.dens$compute_density_ratio(mydata[,c("PC1","PC2")])
    
    myweights[,i] <- i.weights
  }
  
  myweights <- (rowMeans(myweights, na.rm=TRUE)) ^ best_lambda
  
  tree <- pltree(as.formula(paste0(c("G ~ "), paste(attrib$clim, collapse = " + "))),
                 data = mydata, weights = myweights, alpha = alpha, minsize = round(n*0.1, -1), bonferroni = bonferroni)
  
  #export tree as chart 
  png(paste(output, "PL_tree_plot.png", sep=""),
      units="cm", width=60,height=25, res=200)
  plot(tree, abbreviate = 2, worth = TRUE, ref = 1)
  dev.off()
  
  #export summary of fitted pltree
  capture.output(tree, summary(tree), file = paste0(output, "pltree_clim.txt") )
  
  #extract the terminal nodes
  terminal_nodes <- partykit::nodeids(tree, terminal = TRUE)
  
  lims <- NULL
  #identify limits in error bars
  for (j in 1:length(terminal_nodes)){
    lims <- rbind(lims, qvcalc(tree[[terminal_nodes[j]]]$node$info$object)$qvframe)
  }
  
  xmax <- round(max(lims[,"estimate"] + lims[,"quasiSE"]) + 0.5, digits = 1)
  xmin <- round(min(lims[,"estimate"] - lims[,"quasiSE"]) - 0.5, digits = 1)
  
  #generate plots with error bars
  for (j in 1:length(terminal_nodes)){
    qvcoeff <- qvcalc(tree[[terminal_nodes[j]]]$node$info$object)$qvframe
    nobs <- as.integer(tree[[terminal_nodes[j]]]$node$info$nobs)
    qvcoeff$player <- factor(rownames(qvcoeff), levels = rownames(qvcoeff))
    #plot estimates
    plot <- ggplot(qvcoeff, aes(x=estimate, y=player)) +
      geom_point(pch=21, size=3,fill="black",colour="black") +
      geom_errorbarh(aes(xmin=estimate-(quasiSE),
                         xmax=estimate+(quasiSE)),
                     colour="black",height = 0) +
      xlim(min = xmin, max = xmax) +
      theme_bw() +
      labs(x = NULL, y = NULL, title=paste("Node ", terminal_nodes[j], " (n= ", nobs, ")", sep="")) +
      theme(plot.title = element_text(size=17),
            axis.text = element_text(size=12, colour="black"),
            legend.text= element_text(size=12, colour="black"),
            axis.text.x = element_text(size=17, angle = 0, hjust=1, vjust=1, face="plain"),
            axis.text.y = element_text(size=17, angle = 0, hjust=1, vjust=1, face="plain"),
            plot.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    #export plot
    ggsave(paste0(output, "pltree_qvcalc_node_", terminal_nodes[j] , ".png"), plot = plot, dpi = 300,
           width = 20, height = 20, units = "cm")
  }
  
  cat("################################## \n End analysis:", toupper(crop[m]) , " \n Time:", date())

}
