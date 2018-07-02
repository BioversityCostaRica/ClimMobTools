####################################################
###### MODELLING CROWDSOURCING DATA 
#Kaue de Sousa <k.desousa (at) cgiar.org> & Jacob van Etten <j.vanetten (at) cgiar.org>
#Bioversity International - Costa Rica
#First run 22Oct2017
#Updated 02May2017
####################################################

library(Matrix)
library(reshape2)
library(tidyverse)
library(magrittr)
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

#Define number of cores for parallelisation
n_cpu <- 3
#Define working directory
wd <- ""
setwd(wd)

#read data
df <- read_csv("./input/tricot_data.csv", na = c("NA",""))

#Crop that will be analysed 
crop <-  sort(unique(df$crop))

#This analysis will evaluate the performance of Plackett-Luce models under following models:
# PLT-null  null model (without explanatory variables)
# PLT-design  with explanatory variables related to experiental design (season, planting date, location, soil) 
# PLT-season with season as explanatory variable
# PLT-climate  with climatic explanatory variables
# PLT-clim+loc with explanatory variables from climate and location
approach <- c("PLT-null","PLT-design","PLT-climate","PLT-clim+loc")

#Run over crops
for(m in seq_along(crop)){

  cat("################################## \n Starting analysis:", toupper(crop[m]), " \n Time:", date(), "\n")

  #Subset data for the m crop
  mydata <- df[df$crop == crop[m] , ]

  #Create folder for outputs
  output <- "./output/"
  dir.create(output, showWarnings = FALSE)

  #An output folder for the crop
  output <- paste0(output,unique(mydata$crop),"_",unique(mydata$country),"/")
  dir.create(output, showWarnings = FALSE)

  #Reclassify factors levels to fit in the subseted data
  mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], as.character)
  mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], as.factor)

  #Check consistency of data and remove inconsistent rows
  if(crop[m] == "commonbean"){
    #remove NA's in rankings for overall_performance vs local
    keep <- !is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c)
    mydata <- mydata[keep,]
  }

  if(crop[m] == "wheat"){
    #for wheat
    #remove items (varieties) tested in only 1 season
    #here we remove items, not the entire row, rows with more than 1 NA will be removed next
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
    
    #remove wrong evaluations best == worst in overall_performance
    keep <- ifelse(mydata$best==mydata$worst, FALSE, keep)
    
    mydata <- mydata[keep, ]
  }

  #Check correlation between overall performance and yield
  #this information is only available for wheat (IND) and common beans (NIC)
  if(crop[m] != "durumwheat"){
    #rankings from overall performance
    overall <- mydata[,c("variety_a","variety_b","variety_c","characteristic","best","worst")]
    #rankings from yield
    yield <- mydata[,c("variety_a","variety_b","variety_c","yield","best_yield","worst_yield")]
    #make sure that both datasets has the same names
    names(yield) <- names(overall)

    #remove NA's in yield
    ykeep <- !is.na(yield$best) & !is.na(yield$worst) &  yield$best != yield$worst
    #keep those rows with no NA's
    yield <- yield[ykeep,]
    overall <- overall[ykeep,]

    #get rankings for yield
    yield <- grouped_rankings(as.PL(yield, local = F, additional.rank = F ), seq_len(nrow(yield)))
    yield <- yield[1:length(yield),, as.grouped_rankings = FALSE]

    #get rankings for overall performance
    overall <- grouped_rankings(as.PL(overall, local = F, additional.rank = F ), seq_len(nrow(overall)))
    overall <- overall[1:length(overall),, as.grouped_rankings = FALSE]

    #Calculate Kendall's correlation
    #export this output
    capture.output(get_tau(yield, overall)[1] * -1 , file = paste0(output,"Kendall_tau_yield_vs_overall.txt"))
  }

  #Take explanatory variables in a separate dataframe
  covar <- cbind(mydata [, c("season","lat","lon","xy","yx","soil", "planting_day","planting_date")],
                 mydata [, c(30:(ncol(mydata)-4))])

  #Identify explanatory variables with near zero variance
  varout <- caret::nearZeroVar(covar)
  cat("Removing these variables with near zero variance: \n",sort(names(covar)[varout]),"\n")
  #Remove variables with near zero variance
  covar <- covar[,-varout]
  
  #Take rankings from mydata 
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
  ## for common beans we add the comparison with the local var as partial rankings
  ## in this approach the three variaties from tricot are compared together in the same rank
  ## then each variety is compared with local as an individual rank (additional.rank)
  if(crop[m] == "commonbean"){
    G <- grouped_rankings(as.PL(mydata, local = TRUE, additional.rank = TRUE ), rep(seq_len(n), 4))
  }
  
  ## information about local variety is not available for wheat
  ## then a simple grouped rankings is created
  if(crop[m] == "wheat"){
    G <- grouped_rankings(as.PL(mydata, local = FALSE), seq_len(n))
  }

  ## durum data is another format
  ## reorganise the dataset to convert into grouped rankings
  if(crop[m] == "durumwheat"){

    #convert data in long format
    R <- bind_cols(melt(mydata[,c(1,2:5)], id.vars = c("package")),
                   melt(mydata[,c(1,6:9)], id.vars = c("package")))
    
    #keep item names and rankings
    R %<>% 
      dplyr::select(c(1,3,6)) %>%
      set_colnames(c("package","var","rank")) %>%
      arrange(., package)
    
    #check the number of observations per variety
    #remove those tested in less than 20 farms 
    rmitem <- R %>% 
      group_by(var) %>%  
      count(var) %>%
      filter(n < 20) %>%
      dplyr::select(var) %>%
      as.matrix() %>%
      as.vector()

    #take variety names
    items <- sort(unique(R[,"var"]))
    
    #match those to remove
    rmitem <- match(rmitem, items)

    #convert dataset from long format into a sparsed matrix
    R <- sparseMatrix(i = rep(1:n, each = 4),
                      j = match(R[,"var"], items),
                      x = R[,"rank"],
                      dimnames = list(as.integer(unique(R[,"package"])), items ))
    #sparsed matrix into a ordinary matrix from R package base
    R <- as.matrix(R)
    #remove varieties selected previously
    R <- R[,-rmitem]
    #remove rows with less tham 2 items
    R <- R[-which(rowSums(R > 0) < 2), ]

    #row names here are the packages codes, get this to merge rankings and explanatory variables
    keep <- row.names(R)
 
    #add package codes to covar
    covar <- cbind(package = mydata[,"package"],covar)
    
    #subset covar
    covar <- covar[covar$package %in% keep , ]
    #reorder by packages
    covar %<>% arrange(., package)

    #define new n
    n <- nrow(R)

    #rename rows the rankings
    row.names(R) <- 1:n
    row.names(covar) <- row.names(R)

    pack <- mydata[,"package"]
    
    #convert to a grouped rankings
    G <- grouped_rankings(as.rankings(R), seq_len(n))

    #remove packages codes
    covar <- covar[,-1]
    
  }

  cat("This analysis will use",n,"observations.\n")
  
  #Merge grouped rankings with explanatory variables
  mydata <- cbind(G, covar)

  #Define folds based on the season where this crop was evaluated
  folds <- as.integer(as.factor(as.character(mydata$season)))
  #number of folds
  k <- max(folds)

  #Calculate the weights of each season based on the square root of n observations per season divided by total n
  wseason <- as.vector(summary(mydata$season))
  wseason <- sqrt(wseason/n)
  wseason <- wseason/sum(wseason)

  #Perform a forward variable selection on explanatory variables and select those who better contribute to improve predictions between seasons
  ## also compare the performance of explanatory variables with a null model 
  ## to be able to predict a PLtree object from a null model we create a NULL variable (variable with no variance) and added to the main dataset
  mydata$P1 <- 1

  #Select explanatory variables
  expvar <- names(mydata)[5:ncol(mydata)]

  #Define PLT parameters 
  ## minimum size of each node 
  minsize <- round(n*0.3, -1)
  ## if bonferroni correction will be applied
  bonferroni <- TRUE
  ## the significance level for spliting the data into nodes
  alpha <- 0.01

  #Define initial parameters for forward selection
  ## baseline deviance
  par_n <-  0
  ## vector to keep best explanatory variables
  var_keep <- NULL
  ## if TRUE the routine will keep running, this parameters is updated at the end of each "while" routine
  best <- TRUE
  ## number of runs
  runs <- FALSE
  ## vector with best parameters (loglik, vars and lambda) in each run
  best_parameters <- NULL

  #Remove unused objects and reduce size of globalenv() in parallel export
  rm(covar, yield, ykeep, overall)

  #Create cluster to do parallelisation
  cluster <- makeCluster(n_cpu)
  registerDoParallel(cluster)

  #Perform forward selection
  #run untial the model reach its best performance
  while(best){

    cat("Starting Forward Selection. Run ", sum(runs)+1, "\n Time: ", date(), "\n")

    fs <- length(expvar)

    #get predictions from nodes and put in matrix (foreach)
    models <- foreach(i = 1:fs,
                      .combine = acomb,
                      .packages = c("PlackettLuce","psychotree"),
                      .export = ls(globalenv())) %dopar% (
                        f1(formula = as.formula(paste0("G ~ ", paste(unique(c(var_keep, expvar[i])), collapse = " + "))),
                           d = mydata,
                           k = k, folds = folds, minsize = minsize,
                           alpha = alpha, bonferroni = bonferroni)
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

    #take the model call from each model
    call <- models[,ncol(models)]

    #take maximum parameter from Akaike weights
    par_max <- max(meanAW)

    #take the position of par_max in expvar vector
    index_par_max <- which.max(meanAW)

    #Is par_max best (higher) than par_n?
    best <- par_max > par_n #if not, the forward selection will stop

    #if best, save the outputs
    if(best){

      #take the name of best variable
      best_var <- expvar[index_par_max]
      cat("### Best covariate found:", best_var, "\n")

      #sum runs
      runs <- c(runs, best)

      #take outputs from this run and export to a .csv file
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

  #Stop cluster connection
  stopCluster(cluster)

  #take the best model
  best_model <- var_keep

  #Save parameters used in this analysis
  write.table(rbind(n = n, minsize = minsize, bonferroni = bonferroni,
                    alpha = alpha, model = paste(best_model, collapse = " + ")),
              file = paste0(output, "PLT_parameters.txt" ))

  cat("End Forward Selection. \n Time: ", date(), "\n Best model will use:",  best_model ,"\n")

  
  #Define list of explanatory variables to use in each approach
  attrib <- list(null = c("P1"),
                 desi = c("lon","lat","xy","yx","planting_day","season","soil"),
                 clim = c(best_model),
                 clsp = c(best_model, c("lon","lat","xy","yx","soil") ))
  
  
  #Run the best model against the PLT-null, PLT-clim+loc and PLT-design
  cat("Starting blocked cross-validation", date() ,"\n")
  
  blockedfolds <-  matrix(NA, nrow = length(approach), ncol = (6+(k*3)),
                          dimnames = list(1:length(approach),
                                          c("Approach","Model","mean.Ktau","mean.AIC","mean.Deviance",
                                            paste0(rep("Ktau",times = k), 1:k),
                                            paste0(rep("AIC",times = k), 1:k),
                                            paste0(rep("Deviance",times = k), 1:k), "mean.AW")))

  #Add approach to the first column
  blockedfolds[,1] <- approach

  #run blocked cross-validation over approaches
  for(a in 1:length(approach)){

    cat("#### blocked cross-validation in", approach[a], "\n")

    formula_a <- as.formula(paste0("G ~ ", paste(attrib[[a]], collapse = " + ") ))

    blockedfolds[ a , c(2:(ncol(blockedfolds)-1))] <- unlist(crossvalidation_PLTE(formula_a,
                                                              d = mydata, k = k, folds = folds, minsize = minsize,
                                                              alpha = alpha, bonferroni = bonferroni,
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



  #Run average season using historical data
  cat("Starting average season with blocked cross-validation", date() ,"\n")

  #load climatology data
  load(paste0("./input/", unique(df$country[df$crop==crop[m]]),"/", crop[m],"_climatology.RData"))
  #define number of predictions to be generated n.years * n.estimated planting dates
  npreds <- length(climatology) * length(climatology[[1]])

  #list into array
  nr <- nrow(climatology[[1]][[1]])
  nc <- ncol(climatology[[1]][[1]])

  climatology %<>%
   unlist(.) %>%
   array(., dim = c(nr, nc, npreds),
         dimnames = list(1:nr, names(climatology[[1]][[1]]), 1:npreds ))

  #durumwheat have a particular order, filter data to fit this order
  if(crop[m]=="durumwheat"){

    tf <- as.data.frame(cbind(pack, climatology[[1]][[1]]))

    keep <- tf %<>%
      filter(., pack %in% keep ) %>%
      rownames_to_column() %>%
      arrange(., pack) %>%
      dplyr::select(1) %>%
      as.matrix() %>%
      as.vector()

    rm(tf)
    }


  #matrix to keep results
  avgseason <-  matrix(NA, nrow = npreds, ncol = (3+(k*3)),
                          dimnames = list(1:npreds,
                                          c("mean.Ktau","mean.AIC","mean.Deviance",
                                            paste0(rep("Ktau",times = k), 1:k),
                                            paste0(rep("AIC",times = k), 1:k),
                                            paste0(rep("Deviance",times = k), 1:k))))


  pb <- txtProgressBar(min = 1, max = npreds, style = 3)

  for(a in seq_len(npreds)){

    #temporary dataframe with G and explanatory variables
    a_df <- cbind(G, as.data.frame(climatology[keep, , a]))


    avgseason[ a , ] <- unlist(crossvalidation_PLTE(as.formula(paste0("G ~ ", paste(attrib$clim, collapse = " + ") )),
                                                    d = a_df, k = k, folds = folds, minsize = minsize,
                                                    alpha = alpha, bonferroni = bonferroni, 
                                                    verbose = FALSE, mean.method = "stouffer")[c(2:7)])

    setTxtProgressBar(pb, a)

  }

  close(pb)

  cat("End of average season, average Deviance of:",mean(avgseason[,3]),"\n")

  #export data
  write.csv(avgseason, paste0(output, "average_season.csv"), row.names = FALSE)

  #Run a 10-fold cross validation
  cat("Starting k-fold cross-validation",date(),"\n")
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
                                                              alpha = alpha, bonferroni = bonferroni,
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
                 data = mydata, alpha = alpha, minsize = minsize)

  #export tree as chart
  png(paste(output, "PL_tree_plot.png", sep=""),
      units="cm", width=60,height=25, res=200)
  plot(tree, abbreviate = 2, worth = TRUE, ref = 1)
  dev.off()

  #export summary of fitted pltree
  capture.output(tree, summary(tree), file = paste0(output, "pltree_clim.txt") )

  save(tree, mydata, file = paste0(output,"pltree_fit.RData"))
  
  #extract the terminal nodes
  nodes <- partykit::nodeids(tree, terminal = TRUE)

  lims <- NULL
  
  #identify limits in error bars
  for (j in seq_along(nodes)){

    nj <- tree[[ nodes[j] ]]$node$info$object

    lims <- nj %<>%
      itempar(., vcov = TRUE, alias = TRUE) %>%
      qvcalc::qvcalc.itempar(.) %>%
      extract2(2) %>%
      rbind(., lims )


    
    
  }

  xmax <- round(max(lims[,"estimate"] + lims[,"quasiSE"]) + 0.01, digits = 2)
  xmin <- round(min(lims[,"estimate"] - lims[,"quasiSE"]) - 0.01, digits = 2)

  probs <- NULL
  
  #Generate plots with error bars using qvcal
  for (j in seq_along(nodes)){

    #Get node information
    nj <- tree[[ nodes[j] ]]$node$info$object
    #Get item names
    items <- as.factor(c(names(tree[[nodes[j]]]$node$info$coefficients)))
    #Get number of observers in this node
    nobs <- as.integer(tree[[nodes[j]]]$node$info$nobs)

    #Calculate coefficients using qvcalc
    qvcoeff <- nj %<>%
      itempar(., vcov = TRUE, alias = TRUE) %>%
      qvcalc::qvcalc.itempar(.) %>%
      extract2(2) %>%
      as.matrix() %>%
      as_tibble() %>%
      mutate(items = items,
             bmin = estimate-(quasiSE),
             bmax = estimate+(quasiSE))

    qvcoeff$bmin <- ifelse(qvcoeff$bmin <= 0, 0.001, qvcoeff$bmin)

    #Plot coefficients
    plot <- ggplot(qvcoeff, aes(x=estimate, y=items)) +
      geom_point(pch=21, size=3,fill="black",colour="black") +
      geom_errorbarh(aes(xmin=bmin,
                         xmax=bmax),
                     colour="black",height = 0.17) +
      geom_vline(xintercept = 1/length(items), colour="#E5E7E9", size = 0.5) +
      theme_bw() +
      scale_x_continuous(limits=c(0, xmax), breaks = seq(0,xmax,by=0.02) ) +
      labs(x = NULL, y = NULL, title = paste0("Node ", nodes[j], " (n= ", nobs, ")")) +
      theme(plot.title = element_text(size=17),
            axis.text = element_text(size=12, colour="black"),
            legend.text= element_text(size=12, colour="black"),
            axis.text.x = element_text(size=17, angle = 0, hjust=0.5,vjust=1, face="plain"),
            axis.text.y = element_text(size=17, angle = 0, hjust=1, vjust=0.5, face="plain"),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

    #Export plot
    #define width and height
    if(crop[m]=="durumwheat") h <- 55
    if(crop[m]=="wheat") h <- 30
    if(crop[m]=="commonbean") h <- 20

    w <- 20

    ggsave(paste0(output, "pltree_qvcalc_node_", nodes[j] , ".png"), plot = plot, dpi = 400,
           width = w, height = h, units = "cm")
    
    qvcoeff$node <- nodes[j]
    
    probs <- rbind(probs, qvcoeff[,c(1:5,8)])
    
  }

  #Export probabilities of winning per node
  write_csv(probs, paste0(output, "win_probabilies.csv"))

  #Use historical data to get predictions using fitted PLT-clim model
  avgpred <- array(NA, dim = c(n, length(items), npreds), dimnames = list(1:n, as.character(items), 1:npreds) )

  for(i in seq_len(npreds)) {

    clim_i <- as_tibble(climatology[ keep , , i])

    avgpred[ , , i] <- predict(tree, newdata = clim_i)

  }

  save(avgpred, file = paste0(output, "average_predictions.RData"))

  cat("################################## \n End analysis:", toupper(crop[m]) , " \n Time:", date() , "\n ################################## \n")
  
  
}
