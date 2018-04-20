####################################################
###### MODELLING CROWDSOURCING DATA 
#Kaue de Sousa <k.desousa(at)cgiar.org> & Jacob van Etten <j.vanetten(at)cgiar.org>
#Bioversity International - Costa Rica
#First run 22Oct2017
#Updated 02Apr2017
####################################################

library(Matrix) ###
library(reshape2) ###
library(tidyverse)
library(readr)
library(car)
library(caret)
library(RCurl)
library(PlackettLuce)
library(psychotree)
library(partykit)
library(e1071)
library(foreach)
library(doParallel)
library(abind)

#add ClimMob tools for data learning
tools <- RCurl::getURL("https://raw.githubusercontent.com/kauedesousa/ClimMobTools/pilot/ClimMobLearning.R", ssl.verifypeer = FALSE)
eval(parse(text = tools))
#cores for parallel
#detectCores()
n_cpu <- 3
#working directory
#wd <- "C:/Users/Jacob Van Etten/Dropbox/vanEtten_etal2018_replication_data/"
wd <- "C:/Users/KAUE/Dropbox (Bioversity CR)/vanEtten_etal2018_replication_data/"
#wd<-"/home/tucurrique2/ALLAN/vanEtten_etal2018_replication_data/"
setwd(wd)

# #read data from India, Ethiopia and Nicaragua
# eth <- read_csv("./input/ETH/climmob_ethiopia.csv", na = c("NA",""," "))
# ind <- read_csv("./input/IND/climmob_india.csv", na = c("NA",""," "))
# nic <- read_csv("./input/NIC/climmob_nicaragua.csv", na = c("NA",""," "))
# #merge data 
# df <- rbind(nic, ind, eth)
# #remove rice data 
# df <- df[df$crop!="rice",]
# #export file
# write.csv(df, "./input/tricot_data.csv", row.names = FALSE)

#read data
df <- read_csv("./input/tricot_data.csv", na = c("NA",""," "))

#Crop that will be analysed 
summary(as.factor(df$crop))
crop <- sort(unique(df$crop))

#This analysis will evaluate the performance of Plackett-Luce models under following models:
# PLT-null  null model (without explanatory variables)
# PLT-design  with explanatory variables related to experiental design (season, planting date, location, soil) 
# PLT-season with season as explanatory variable
# PLT-climate  with climatic explanatory variables
# PLT-clim+loc with explanatory variables from climate and location
approach <- c("PLT-null","PLT-design","PLT-climate","PLT-clim+loc")

names(df)

#run over crops
for(m in 1:length(crop)){
  
  cat("################################## \n Starting analysis:", toupper(crop[m]), " \n Time:", date(), "\n")
  
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
  
  #Add new spatial explanatory variables 
  mydata$xy <- mydata$lon + mydata$lat
  mydata$yx <- mydata$lon - mydata$lat
  
  #The planting day (decimal day of the year)
  mydata$planting_day <- as.integer(format(mydata$planting_date,'%j'))
  
  #CHECK correlation between overall performance and yield
  #only available for common wheat (IND) and common beans (NIC)
  if(crop[m] != "durumwheat"){
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
    
    #Calculate Kendall's correlation
    #export this output
    capture.output(get_tau(yield, overall)[1] * -1 , file = paste0(output,"Kendall_tau_yield_vs_overall.txt"))
  }
  
  #Take explanatory variables in a separate data frame and check variance
  covar <- cbind(mydata [, c("season","planting_day","planting_date","soil","lon","lat","xy","yx")],
                 mydata [, c(24:205)])

  #Identify explanatory variables with near 0 variance
  varout <- caret::nearZeroVar(covar)
  #Remove those whith near 0 variance
  covar <- covar[,-varout]
  
  #Keep only the rankings in mydata
  if(crop[m] != "durumwheat"){
    mydata <- cbind(mydata[, c("variety_a","variety_b","variety_c",
                               "characteristic","best","worst",
                               "overall_vs_local","var_a","var_b","var_c")])
  }
  
  if(crop[m] == "durumwheat"){
    mydata <- cbind(mydata[, c("package","variety_a","variety_b","variety_c","variety_d",
                               "rank_variety_a","rank_variety_b","rank_variety_c","rank_variety_d")])
  }
  
  #Take the number of rows in this dataset
  n <- nrow(mydata) 
  
  #Generate a Plackett-Luce grouped rankings
  #for common beans we add the comparison with the local var in partial rankings 
  #in this approach each comparison with local is a individual rank
  if(crop[m] == "commonbean"){
    G <- grouped_rankings(as.PL(mydata, local = TRUE, additional.rank = TRUE ), rep(seq_len(n), 4))
  }
  #information about local item is not available for the other datasets
  #a simple grouped rankings is created
  if(crop[m] == "wheat"){
    G <- grouped_rankings(as.PL(mydata, local = FALSE), seq_len(n))
  }
  
  #ethiopia
  if(crop[m] == "durumwheat"){
    
    #convert data in long format 
    R <- cbind(melt(mydata[,c(1,2:5)], id.vars = c("package")),
               melt(mydata[,c(1,6:9)], id.vars = c("package")))
    #keep item names and rankings
    R <- R[,c(1,3,6)]
    #check the frequency of observations per accession 
    rmitem <- as.matrix(summary(as.factor(R[,2])))
    #remove those tested in less than 10 plots 
    rmitem <- rmitem[rmitem[,1] < 10, ] 
    rmitem <- names(rmitem)
    
    #arrange observations per package number
    R <- arrange(R, package)
    #rankings as integer
    R[,3] <- as.integer(R[,3])
    #take item names
    items <- sort(unique(R[,2]))
    #match those to remove 
    rmitem <- match(rmitem, items)
    
    #convert in a sparsed matrix
    R <- sparseMatrix(i = rep(1:n, each = 4),
                      j = match(R[,2], items),
                      x = R[,3], 
                      dimnames = list(as.integer(unique(R[,1])), items ))
    #Plackett-Luce rankings 
    R <- as.rankings(R)
    #exclude selected items from rankings
    R <- R[, -rmitem]
    #PL rankings into grouped rankings
    G <- grouped_rankings(R, 1:n)
  }
  
  #Merge these rankings with covariables
  mydata <- cbind(G, covar)
  
  #Define folds based on the season where this crop was evaluated
  folds <- as.integer(as.factor(as.character(mydata$season)))
  #number of folds
  k <- max(folds)
  
  #Calculate the weights represented by seasons
  wseason <- as.vector(summary(mydata$season))
  wseason <- sqrt(wseason/n)
  wseason <- wseason/sum(wseason)
  
  ## Perform a forward variable selection on climate+spatial explanatory variables and select those who contribute to improve predictions between seasons
  #to be able to predict a PLtree object from a null model we create a NULL variable and add to the main dataset 
  mydata$P1 <- 1 
  
  #select explanatory variables
  expvar <- names(mydata)[5:ncol(mydata)]
  #remove explanatory variables related to the firts quarter of crop 
  expvar <- expvar[ !grepl("_1st", expvar) & !grepl("GDD", expvar) ]
  #expvar <- c("P1","maxNT_2nd","minNT_2nd","Rx1day_2nd","Rx5day_2nd","R10mm_4th","lat")#use this for tests
  
  #Define parameters for the PL model 
  if(crop[m] == "commonbean"){
    minsize <- 200
  }
  if(crop[m] == "wheat"){
    minsize <- 2000
  }
  if(crop[m] == "durumwheat"){
    minsize <- 250
  }
  bonferroni <- TRUE
  alpha <- 0.01
  npseudo <- 0.5
  weights <- NULL
  lambda <- 0
  
  #Define initial parameters for forward selection
  ##AIC from null model (used as baseline)
  par_n <-  0
  ##vector to keep best explanatory variables
  var_keep <- NULL
  ##if TRUE the routine will keep running, this parameters is updated at the end of each "while" routine
  best <- TRUE
  ##keep number of runs
  runs <- FALSE
  ##vector with best parameters (loglik, vars and lambda) in each run
  best_parameters <- NULL
  
  #remove unused objects and reduce size of globalenv() in parallel export
  rm(covar,G, keep, yield, overall)
  
  #Create cluster to do parallel computation to speed things up
  cluster <- makeCluster(n_cpu)
  registerDoParallel(cluster)
  getDoParWorkers()
  
  ## Perform forward selection 
  while(best){
    
    cat("Starting Forward Selection. Run ", sum(runs)+1, "\n Time: ", date(), "\n")
    
    fs <- length(expvar)
    
    #get predictions for test from nodes and put in matrix (foreach)
    models <- foreach(i = 1:fs,
                      .combine = acomb,
                      .packages = c("PlackettLuce","psychotree"),
                      .export = ls(globalenv())) %dopar% (
                        f1(formula = as.formula(paste0("G ~ ", paste(unique(c(var_keep, expvar[i])), collapse = " + "))), 
                           d = mydata, 
                           k = k, folds = folds, minsize = minsize, 
                           alpha = alpha, bonferroni = bonferroni,
                           weights = weights, npseudo = npseudo, lambda = lambda)
                      )
    
    #calculate Akaike weights along explanatory variables per season
    AW <- as.matrix(models[,1:k])
    AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))
    AW <- apply(AW, 2, function(X) akaike.weights(X)[[3]] )
    
    #then take the weighted average per season Stouffer's method with weights
    for(i in seq_len(k)) AW[,i] <- AW[,i] * wseason[i]
    
    #take the Stouffer's weighted mean 
    meanAW <- rowSums(AW)
    
    #Take a simple mean of inernodes generated by each model
    nodes <- models[,c((k+1):(k+k))]
    nodes <- apply(nodes, 2, as.integer )
    nodes <- rowMeans(nodes)
    
    #take Kendall tau
    tau <- models[,(ncol(models)-1)]
    
    #take call from each model
    call <- models[,ncol(models)]
    
    #take maximum parameter from Akaike weights
    par_max <- max(meanAW)
    
    #take the position of par_max in expvar vector
    index_par_max <- which.max(meanAW)
    
    #Is par_max best (higher) than par_n?
    best <- par_max > par_n #if not, stop
    
    #if best, save the outputs 
    if(best){
      
      #take the name of best variable
      best_var <- expvar[index_par_max]
      
      #sum runs 
      runs <- c(runs, best)
      
      #ensemble outputs from this run and export to a .csv file 
      out <- as_tibble(cbind(call, models[,c(1:k)], AW, meanAW = meanAW, inernodes = nodes, tau = tau))
      colnames(out) <- c("call", paste0(rep("Deviance",k), 1:k), paste0(rep("AkaikeWeight",k), 1:k), "meanAW", "inernodes", "Kendall_tau")
      
      #remove best_var from this run 
      expvar <- expvar[!grepl(best_var, expvar)]
      
      #remove null var from the first run, no longer necessary 
      expvar <- expvar[!grepl("P1", expvar)]
      
      #keep this model for the next run
      var_keep <- c(var_keep, best_var)
      
      #add best model to the next round 
      expvar <- c(expvar, paste(var_keep, collapse = ' + '))
      
      #change the base for par_n (minimun accepted value)
      par_n <- par_max
      
      #keep the best parameters form this run
      best_parameters <- rbind(best_parameters, cbind(par_max, model = toString(var_keep)))
      
      write.csv(out, paste0(output, "model_parameters_run", sum(runs),".csv" )  )
      cat(" ####### End run", sum(runs), "\n" )
    }
    
  }
  
  #return best model
  best_model <- var_keep
 
  cat("End Forward Selection. \n Time: ", date(), "\n")
  
  #Stop connection with cores
  stopCluster(cluster)
  
  #Define list of explanatory variables to use in each approach
  attrib <- list(null = ("P1"),
                 desi = c("lat","lon","xy","yx","planting_date","planting_day","season","soil"),
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
    
    cat("#### blocked cross-validation in", approach[a], "\n")
    
    formula_a <- as.formula(paste0("G ~ ", paste(attrib[[a]], collapse = " + ") ))
    
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
  
  #matrix to keep the outputs 
  kfolds <-  matrix(NA, nrow = length(approach), ncol = (6+(k*3)), 
                    dimnames = list(1:length(approach),
                                    c("Approach","Model","mean.Ktau","mean.AIC","mean.Deviance",
                                      paste0(rep("Ktau",times = k), 1:k),
                                      paste0(rep("AIC",times = k), 1:k),
                                      paste0(rep("Deviance",times = k), 1:k), "sum.Deviance")))
  kfolds[,1] <- approach
  
  #run k-fold over approaches 
  for(a in 1:length(approach)){
    cat("#### k-fold in", approach[a], "\n")
    
    formula_a <- as.formula(paste0("G ~ ", paste(attrib[[a]], collapse = " + ") ))
    
    kfolds[ a , c(2:(ncol(kfolds)-1))] <- unlist(crossvalidation_PLTE(formula_a, 
                                                              d = mydata, k = k, folds = folds, minsize = minsize, 
                                                              alpha = alpha, bonferroni = bonferroni, weights = weights, lambda = lambda,
                                                              verbose = FALSE)[c(1:7)])
  }
  
  #calculate Akaike weights along approaches over Deviance
  AW <- as.matrix(kfolds[,(ncol(kfolds)-(k)):(ncol(kfolds)-1)])
  AW <- apply(AW, 2, function(x) round(as.numeric(x), digits = 10))

  #take the sum of Deviances
  kfolds[,ncol(kfolds)] <- rowSums(AW)
  
  #write csv file 
  write.csv(kfolds, paste0(output, "performance_PLT_10fold_crossvalidation.csv"), row.names = FALSE)
  
  
  #Fit pltree with best climate model
  cat("Fit pltree with best climate model \n")
  
  tree <- pltree(as.formula(paste0(c("G ~ "), paste(attrib$clim, collapse = " + "))),
                 data = mydata, weights = weights, alpha = alpha, minsize = round(n*0.15, -2), bonferroni = bonferroni)
  
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
  
  lims[,3] <- ifelse(is.nan(lims[,3]), mean(lims[,3], na.rm=TRUE), lims[,3])
  
  xmax <- round(max(lims[,"estimate"] + lims[,"quasiSE"]) + 0.05, digits = 1)
  xmin <- round(min(lims[,"estimate"] - lims[,"quasiSE"]) - 0.05, digits = 1)

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
    
    h <- if(crop[m]=="durumwheat") 45 else 20
    w <- 20
    
    ggsave(paste0(output, "pltree_qvcalc_node_", terminal_nodes[j] , ".png"), plot = plot, dpi = 300,
           width = w, height = h, units = "cm")
  }
  
  cat("################################## \n End analysis:", toupper(crop[m]) , " \n Time:", date())

}
