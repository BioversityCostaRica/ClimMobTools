#FUNCTIONS FOR DATA LEARNING AND PORTFOLIOS CONTRUCTIONS USING CROWDSOURCING DATA
#Jacob van Etten, Kauê de Sousa, Heather Turner 
#First run 31 Nov 2017
#Updated in 01 Feb 2018


# Tools for data handling ####
#Prepare the data for analysis with PlackettLuce.
# require(PlackettLuce)
# x = a tibble with tricot data
# local = add comparison with local item 
# additional.rank = comparison with local item will be done as addional rank comparing all items with local individually
# function return a PlackettLuce rank (see PlackettLuce::as.rankings)
as.PL <- function(x, local = FALSE, additional.rank = TRUE)
  {
  #get nrow in x
  n <- nrow(x)
  
  #First we consider the best and worst rankings. 
  #These give the item the observer though was best or worst, coded as A, B or C for the first, second or third variety assigned to the farmer respectively.
  #Convert these to numeric values, allowing us to impute the middle-ranked variety (a strict ranking is assumed here, so the sum of each row should be 6)
  x <- within(x, {
    best <- match(best, c("A", "B", "C")) 
    worst <- match(worst, c("A", "B", "C"))
    middle <- 6 - best - worst
  })
  
  #if there is any NA in item names add a pseudo-item which will be removed later
  if(sum(is.na(x [c("variety_a", "variety_b", "variety_c")])) > 0)  x[is.na(x)] <- "ZZZ"
  
  #This gives an ordering of the three items each observer was given.
  #The names of these items are stored in separate columns
  items <- as.matrix( x [c("variety_a", "variety_b", "variety_c")])
  
  #So we can convert the itemsIDs to the item names
  x <- within(x, {
    best <- items[cbind(seq_len(n), best)]
    worst <- items[cbind(seq_len(n), worst)]
    middle <- items[cbind(seq_len(n), middle)]
  })
  
  #get vector with item names
  itemnames = c(sort(unique(as.vector(unlist(x[c("variety_a", "variety_b", "variety_c")])))))
  #Convert the orderings of the items given to each observer to sub-rankings of the full set of varieties:
  R <- PlackettLuce::as.rankings(x[c("best","middle","worst")], input = "ordering", labels = itemnames)
  #Remove pseudo-item 
  if(any(grepl("ZZZ", itemnames))) R <- R[,-match("ZZZ", itemnames)]
  
  if(local){
    if(additional.rank){
      #Treat the paired comparisons as additional rankings.
      #First we can convert the orderings of the trial varieties to sub-rankings of the full set of items
      #including the local as an additional item, so that we can add the paired comparisons shortly
      #The comparisons with the local item are stored in another set of columns
      #add local to itemnames
      itemnames <- c("Local", itemnames)
      paired <- list()
      for (id in c("a", "b", "c")){
        ordering <- matrix("Local", nrow = n, ncol = 2)
        worse <- x[[paste0("var_", id)]] == "Worse"
        ## name of winner
        ordering[!worse, 1] <- x[[paste0("variety_", id)]][!worse]
        ## name of loser
        ordering[worse, 2] <- x[[paste0("variety_", id)]][worse]
        paired[[id]] <- ordering
      }
      #Again we convert these orderings to sub-rankings of the full set of items and combine them with the rankings of order three:
      paired <- lapply(paired, as.rankings, input = "ordering", labels = itemnames)
      R <- rbind(R, paired[["a"]], paired[["b"]], paired[["c"]])
      
    }
    if(!additional.rank){
      #incorporate the local item into the ranking to capture the full preference information from each observer. 
      #First we count the number of trial items that the observers ranks as better than local item to infer the rank of the local item
      nbetter <- rowSums(x[c("var_a","var_b","var_c")] == "Better")
      add_local <- nbetter + 1
      #Then we adjust the ranks of the trial items that the observer ranks as worse than the local and add the local to the rankings
      worse <- R >= add_local
      R[worse] <- R[worse] + 1
      R <- as.rankings(cbind(R, Local = add_local))
    }
    
  }
  return(R)
}

# Tools for data learning and modelling ####
#Jacob van Etten
#k-fold cross-validation with PL trees
library(PlackettLuce)
library(foreach)
library(doParallel)
library(abind)

#Calculate correlation coefficient
get_tau <- function(predictedcoef, observedrank, ...)
{
  
  nr1 <- nrow(predictedcoef)
  nr2 <- nrow(observedrank)
  if(nr1 != nr2){stop("predictedrank and observedrank do not have equal number of rows")}
  
  nc1 <- ncol(predictedcoef)
  nc2 <- ncol(observedrank)
  if(nc1 != nc2){stop("predictedrank and observedrank do not have equal number of columns")}
  
  nc <- nc1
  
  # The following apply produces a matrix with two rows
  # First row: Kendall tau value for that row
  # Second row: how many pairs entered the comparison as a weight
  # To make the predicted probability of winning match the ranks (1 = best), we take the negative probability
  # No need to convert the probability to ranks -- the cor() function does that internally
  tau_N <- apply(cbind(observedrank, predictedcoef), 1, 
                 function(x) {
                   
                   tau_cor <- cor(x[1:nc], -x[(nc+1):(2*nc)], method = "kendall", use = "pairwise.complete.obs")
                   n <- sum(x[1:nc] > 0, na.rm=TRUE)
                   weight <- n*(n-1)/2
                   return(c(tau_cor, weight))
                   
                 }
  )
  
  # Extract the values from the matrix
  tau <- tau_N[1,]
  N <- tau_N[2,]
  
  # N are very small per observer, so simple averaging will have very little bias
  # No need for z transformation at this stage
  # And z tranformation would always give -Inf or Inf with N=2
  tau_average <- sum(tau * N) / sum(N)
  
  # Effective N is the equivalent N needed if all were compared to all
  # N_comparisons = ((N_effective - 1) * N_effective) / 2
  # This is used for significance testing later
  N_effective <- .5 + sqrt(.25 + 2 * sum(N)) 
  
  return(c(tau_average, N_effective))
  
}

#Function that sends i-fold vector (assigning to inside and outside of fold) of length n to nodes and gets back predictions for the i-fold test set
#It then aggregates the predictions and calculates the tau
pltree_ensemble <- function(train, cluster, ntrees)
{
  
  #Export train, which is what changes (different fold)
  clusterExport(cluster, "train", envir = environment())
  
  #Create function to combine prediction matrices into 3-dimensional array
  #Idea obtained from here:
  #https://stackoverflow.com/questions/17570806/parallel-for-loop-with-an-array-as-output
  acomb <- function(...) abind(..., along = 3)
  
  #Function that creates a prediction array for the test fold
  #Array of 3 dimensions: farms (rows), varieties (columns) and trees (third dimension)
  f1 <- function(formula, dat, minsize, npseudo, alpha, newdata) {
    
    m <- do.call("pltree", list(formula = formula, 
                                data = dat,
                                minsize = minsize,
                                npseudo = npseudo,
                                alpha = alpha))
    
    pp <- predict(object=m, newdata = newdata)
    
    return(pp)
    
  }
  
  #get predictions for test from nodes and put in matrix (foreach)
  preds <- foreach(j = 1:ntrees, 
                   .combine=acomb,
                   .packages=c("PlackettLuce")) %dopar% (
                     
                     f1(formula = formula, 
                        dat = d[sample(which(train), size = sum(train), replace=TRUE),], 
                        minsize = minsize,
                        npseudo = npseudo,
                        alpha = alpha,
                        newdata = d[!train,])
                     
                   )
  
  #get average prediction for test fold 
  #averaging probability of winning across trees
  preds_matrix <- apply(preds, c(1,2), mean)
  
  return(preds_matrix)
  
}

crossvalidation_PLTE <- function(formula, d, k = 10, folds = NULL, minsize = 50, alpha = 0.05,
                                 npseudo=0, ntrees = 25, cluster, ...)
{
  
  #Create folds if needed, check length if given
  n <- nrow(d)
  if(is.null(folds)){
    
    folds <- sample(rep(1:k, times=ceiling(n/k), length.out=n), replace=FALSE)
    
  } else{
    
    if(length(folds) != n){stop("length of folds is not the same as nrow(data)")}
    
  }
  
  #Setting up things for the loop
  #rs as an empty vector, which will be filled in the loop with correlation coefficients
  taus <- rep(NA, times=k)
  Ns_effective <- rep(NA, times=k)
  
  #Prepare matrix with observations to later compare with predictions for each fold
  observed <- d$G[1:length(d$G),, as.grouped_rankings = TRUE]
  index_observed <- attr(observed, "index")
  observed <- observed[1:length(observed),, as.grouped_rankings = FALSE]
  observed[observed == 0] <- NA
  
  #These objects are sent to the worker nodes only once, and then used in every iteration
  #Only the train vector changes and is sent to the nodes for every fold (done inside the pltree_ensemble function)
  clusterExport(cluster, c("formula", "d", "minsize", "npseudo", "alpha"), envir = environment())
  
  for(i in 1:k){
    
    cat("Starting fold ", i, " of ", k, " folds. Time: ", date(), "\n")
    
    #The train and test part change in every iteration
    train <- folds != i
    
    #Use PL tree ensemble on the train data to create a prediction of the test part of the data
    #Output is a matrix with predictions with cases in the rownames
    #Some rows with lacking data may have dropped, so we spend some effort on reformatting
    #observed and predicted matrices to make them match
    preds <- pltree_ensemble(train = train, cluster = cluster, ntrees = ntrees)
    preds_index <- as.integer(rownames(preds))
    
    #"shrink" observed ranks
    obs_index <- index_observed[index_observed %in% preds_index] #inside i-th fold and inside preds
    observed_shrunk <- observed[index_observed %in% preds_index,]
    
    #"expand" predicted coef 
    match_index <- match(obs_index, preds_index)
    predicted_expanded <- preds[match_index,]
    
    #calculate correlation by comparing predicted test data with observed test data
    r <- get_tau(predicted_expanded, observed_shrunk) 
    
    taus[i] <- r[1]
    
    Ns_effective[i] <- r[2]
    
    cat("Fold ", i, " completed. Kendall's correlation: ", r[1], "\n")
    
  }
  
  # Tau is then averaged weighted by number of predicted cases
  # This is equivalent to calculating tau over the whole set (tallying discordant and concordant pairs)
  # Since each sample is small (2 or 3 varieties), this is without much bias, not warrantin a z transformation
  # before averaging.
  # In blocked cross-validation, folds can be of unequal size
  # In k-fold cross-validation, this will slightly correct if n/k is not a round number
  # It will not affect the tau value if n/k is a round number
  foldsize <- table(folds)
  average_tau <- sum(taus * as.vector(foldsize)) / sum(foldsize)
  
  # The Z score is calculated from the average tau
  # N is taken as total N_effective, taking into account that not all were compared to all
  # Calculation of Z score follows: Abdi, H., 2007. The Kendall rank correlation coefficient. Encyclopedia of Measurement and Statistics. Sage, Thousand Oaks, CA, pp.508-510.
  N_effective <- sum(Ns_effective)
  z_score <- average_tau / sqrt((4*N_effective + 10)/((9*N_effective)*(N_effective-1)))
  
  result <- list(folds=folds, taus=taus, foldsize = foldsize, average_tau = average_tau, z_score=z_score, N_effective=N_effective)
  
  return(result)
  
}


#Predict function that works for a PL tree ensemble
#Aggregating is done by taking the mean
#x is a list of PL trees (ensemble)
#newdata is a dataframe with the variables used to predict for new cases
#... passes arguments to the predict function
predict_ensemble <- function(x, newdata, aggregate_function = mean, ...){
  
  ntrees <- length(x)
  vars <- names(coef(x[[1]]))
  nvars <- length(vars)
  array_result <- array(NA, dim= c(nrow(newdata), nvars, ntrees))
  
  for(i in 1:ntrees){
    tree <- x[[i]]
    pr <- predict(tree, newdata = newdata, ...)
    array_result[,,i] <- pr
    
  }
  
  result <- t(apply(array_result, c(1,2), aggregate_function))
  rownames(result) <- vars
  colnames(result) <- rownames(newdata)
  return(result)
  
}


# Tools for portfolios construction ####


#Calculate minimum regret
minRegretn <- function(loss, n){
  
  regretn <- function(x, loss){
    
    n <- length(x)
    for(i in 1:n){loss[i,] <- loss[i,] * x[i]}
    loss <- apply(loss, 2, sum)
    result <- max(loss) + (sum(x) - 1)^2 * 1e10 #penalty if the contributions don´t add up to 1.
    return(result)
    
  }
  
  #optim can sometimes get stuck in a local optimum
  #so we give it random starting values ten times for each n
  #and take the best solution
  for(i in 1:(n*10)){
    
    v <- 1
    r1 <- optim(runif(n), 
                regretn, 
                method="L-BFGS-B", 
                lower=rep(0,3), 
                upper=rep(1,3), 
                loss=loss)
    if(r1$value < v){r2 <- r1; v <- r2$value}
    
  }
  result <- NULL
  result$loss <- r2$value
  result$contributions <- r2$par
  return(result)
  
}

create_portfolio <- function(ensemble, newdata, n, digits=5){
  
  #Predict variety probability of winning
  pred <- predict_ensemble(ensemble, newdata)
  
  #For the loss matrix, get the winner varieties first
  winners <- apply(pred, 2, max)
  loss <- t(1 - t(pred) / winners) #Varieties in rows, seasonal climates in columns
  
  vars <- names(coef(x$trees[[1]]))
  nvars <- length(vars)
  
  if(n == 1){
    
    #Prepare results dataframes
    result <- data.frame(vars = vars, loss = NA)
    
    #Create "portfolios" of one variety
    for(i in 1:nvars){result$loss[i] <- round(max(loss[i,]), digits=digits)}
    colnames(result) <- c("Var", "loss_portfolio")
    
    return(result)
    
  } else {
    
    cc <- combn(1:nvars, n)
    loss_portfolio <- NULL
    contributions <- matrix(NA, nrow=ncol(cc), ncol=n)
    colnames(contributions) <- paste("Contr_Var", 1:n, sep="_")
    
    for(i in 1:ncol(cc)){
      
      s <- minRegretn(loss=loss[cc[,i],], n=n)
      loss_portfolio[i] <- round(s$loss, digits=digits)
      contributions[i,] <- round(s$contributions, digits=digits)
      
    }
    
    varieties <- matrix(NA, nrow=ncol(cc), ncol=n)
    for(i in 1:n){varieties[,i] = vars[cc[i,]]}
    colnames(varieties) <- paste("Var", 1:n, sep="_")
    
    result <- data.frame(varieties, contributions, loss_portfolio)
    return(result)
    
  }
  
}





#crap code ####

# #Take samples for bagging
# #require(PlackettLuce)
# # G = a grouped rankings (see PlackettLuce::grouped_rankings)
# # data = a dataframe (tibble) with explanatory variables (attributes)
# # nboot = an integer indicating the number of samples for bagging procedure
# #function return a list with n (nboot) train data and test data
# bagging.samples <- function(G, data, nboot=100)
# {
#   n <- nrow(data)
#   samples <- list()
#   for(b in 1:nboot){
#     s <- sample(1:n, replace = TRUE)
#     
#     samples[[b]] <- list(train = data[s, ], test = data[-s, ], G = G[s], s = s)
#   }
#   return(samples)
# }
# 
# #PL model in a bagging procedure
# bagging.PL <- function(B, A, minsize = 50, alpha = 0.05, backup = FALSE)
# {
#   n <- nrow(B$train)
#   items <- as.character(dimnames(attr(B$G, "rankings"))[[2]])
#   Y <- B$G
#   tree <- PlackettLuce::pltree(as.formula(paste0("Y ~ ", paste(A , collapse = " + "))), 
#                                data = B$train[A], 
#                                minsize = minsize, alpha = alpha, 
#                                bonferroni = TRUE, method = "L-BFGS")
#   
#   coeff <- array(dim = c(n, length(items), 1), dimnames = list(NULL, items, NULL))
#   coeff[ - B$s , , ] <- predict(tree, newdata = B$test[A] )
#   
#   out <- list(coeff = coeff, tree = tree)
#   
#   if(backup) {save(out, file = "./backup_baggingPL.RData")}
#   
#   return(out)
# }
# 
# extract.coeff <- function(X) 
# {
#   s <- X[[1]][[1]]
#   coeff <- array(dim = c(dim(s)[1:2], length(X)) , dimnames = dimnames(s))
#   for(i in 1:length(X)){ coeff[ , , i] <- X[[i]][[1]] }
#   
#   coeff <- apply(coeff, c(1,2), function(x) { mean(x, na.rm = TRUE)})
#   coeff <- coeff * -1
#   coeff <- t(apply(coeff, 1, function(x){ rank(x, ties.method = "random") } ))
#   
#   return(coeff)
# }
# 
# ##Cross-validation with bagging
# # B <- list()
# # for(i in 1:k){
# #   data <- mydata
# #   data$set <- ifelse(data$folds == i, "test","train")
# #   samples <- bagging.samples(G, data[ data$folds != i , as.vector(unlist(attrib))], nboot = 3 )
# #   B[[i]] <- list(data = data, samples = samples)
# # }
# 
# crossvalidation.PL <- function(B, A, minsize = 50, alpha = 0.05)
# {
#   n <- nrow(B$data)
#   items <- as.character(dimnames(attr( B$samples[[1]]$G , "rankings"))[[2]])
#   test <- B$data[ B$data$set == "test", ]
#   
#   cvtrees  <- parLapply(cl, X = B$samples, fun = function(X) {  bagging.PL(X , A, minsize = minsize, alpha = alpha )$tree })
#   
#   cvpred   <- lapply(cvtrees, function(x) {predict(x, newdata = test)})
#   
#   cvpred   <- array(unlist(cvpred), dim = c(nrow(cvpred[[1]]), ncol(cvpred[[1]]), length(cvpred) ))
#   
#   coeff    <- array( dim = c(n, length(items), 1 ), dimnames = list(NULL, items, NULL))
#   
#   coeff[B$data$set == "test",  , ] <- apply(cvpred, c(1,2), function(x) {mean(x, na.rm=TRUE)} )
#   
#   return(coeff)
# }
# 
# extract.coeff2 <- function(X) 
# {
#   s <- X[[1]][[1]]
#   X <- array(unlist(X), dim = c(dim(X[[1]])[1:2], length(X) ), dimnames = dimnames(X[[1]]))
#   coeff <- apply(X, c(1,2), function(x) { mean(x, na.rm = TRUE)})
#   coeff <- coeff * -1
#   coeff <- t(apply(coeff, 1, function(x){ rank(x, ties.method = "random") } ))
#   return(coeff)
# }
# 
# 
# #Calculate Kendall tau
# K.tau <- function(model.rank = NULL, observed.rank = NULL, drop.local = FALSE, ...)
# {
#   if(drop.local) observed.rank <- observed.rank[,dimnames(observed.rank)[[2]]!="Local"]
#   
#   n <- nrow(observed.rank)
#   
#   tau <- rep(NA, times = n)
#   
#   for(b in 1:nrow(model.rank)){
#     o <- observed.rank[b, observed.rank[b,] > 0] 
#     tau[b] <- cor(as.vector(o), model.rank[b, names(o)], method = "kendall" )
#   }
#   
#   return(mean(tau, na.rm=TRUE)) 
# }
# 
# #extract splitting attributes from nodes 
# extract.from.nodes <- function(X, depth = 5)
# {
#   l <- length(X)
#   attr <- names(X[[1]]$data)
#   infonodes <- array(NA, dim = c(l, depth , 2 ))
#   for(n in 1:l){
#     ids <- partykit::nodeids(X[[n]], terminal = F)[-nodeids(X[[n]], terminal = T)]
#     if(length(ids)==0) next
#     for(i in 1:depth){
#       infonodes[n,i,1] <- attr[X[[n]][[ids[i]]]$node$split$varid] 
#       infonodes[n,i,2] <- if(is.null(X[[n]][[ids[i]]]$node$split$breaks)) NA else X[[n]][[ids[i]]]$node$split$breaks
#       if(length(ids) < depth & i==length(ids)) break
#     }
#   }
#   return(infonodes)
# }







# bagging.PL <- function(formula, G = NULL, data = NULL, nboot = 100, alpha = 0.05, minsize = 25,  ...)
# {
#   n <-  nrow(data)
#   e <- all.names(formula[[3]], functions = FALSE, max.names = -1L)
#   items <- as.character(dimnames(attr(G, "rankings"))[[2]])
#   
#   samples <- matrix(NA, ncol = nboot, nrow = n) 
#   
#   aggr.coeff <- array(NA, dim = c(n, length(items), nboot), dimnames = list(c(1:n),items,c(1:nboot)))
#   models <- list()
#   for(b in 1:nboot){
#     s <- sample(1:n,  replace = TRUE) 
#     
#     Y <- G[s]
#     
#     tree <- PlackettLuce::pltree(as.formula(paste0("Y ~ ", paste(e, collapse = " + "))), 
#                                  data = data[ s, e ], 
#                                  minsize = minsize, alpha = alpha, 
#                                  bonferroni = TRUE, method = "L-BFGS")
#     
#     aggr.coeff[ -s , , b] <- predict(tree, newdata = data[ -s, ])
#     models[[b]] <- tree
#     samples[,b] <- s
#   }
#   
#   aggr.rank <- apply(aggr.coeff, c(1,2), function(x) {mean(x, na.rm=TRUE)} )
#   aggr.rank <- aggr.rank * -1
#   aggr.rank <- t(apply(aggr.rank, 1, function(x){ rank(x, ties.method = "random") } ))
#   
#   outputs <- list(aggr.rank = aggr.rank, trees = models, coefficients = aggr.coeff, samples = samples, call = match.call())
#   
#   class(outputs) <- "PlackettLuce_Bagging"
#   return(outputs)
#   
# }


# #Bagging with k-fold cross validation
# kfold.PL <- function(formula, G = NULL, nboot = 100, data = NULL, alpha = 0.05,  minsize = NULL, backup = TRUE, ... )
# {
#   n <- nrow(data)
#   k <- max(data$folds)
#   e <- all.names(formula[[3]], functions = FALSE, max.names = -1L)
#   items <- as.character(dimnames(attr(G, "rankings"))[[2]])
#   
#   aggr.coeff <- array(NA, dim = c(n, length(items), k), dimnames = list(c(1:n),items,c(1:k)))
#   models <- list()
#   
#   for(b in 1:k){
#     train <- data[data["folds"] != b, ]
#     test  <- data[data["folds"] == b, ]
#     
#     kG  <- grouped_rankings(as.PL(train, local = T, additional.rank = T), rep(seq_len(nrow(train)), 4))
#     
#     tree <- bagging.PL(as.formula(paste0("kG ~ ", paste(c(e), collapse = " + "))),
#                        G = kG, data = train, nboot = nboot, alpha = 0.05, minsize = 25)
#     
#     models[[b]] <- tree$trees
#     
#     kpredict <- lapply(tree$trees, function(x) {predict(x, newdata = test)})
#     kpredict <- array(unlist(kpredict), dim = c(nrow(kpredict[[1]]), ncol(kpredict[[1]]), length(kpredict) ))
#     
#     aggr.coeff[ data["folds"] == b ,  ,  b ] <- apply(kpredict, c(1,2), function(x) {mean(x, na.rm=TRUE)} )
#     
#     if(backup) {save(aggr.coeff, b, file = "./backup_kfold.RData")}
#     
#   }
#   
#   aggr.rank <- apply(aggr.coeff, c(1,2), function(x) {mean(x, na.rm=TRUE)} )
#   aggr.rank <- aggr.rank * -1
#   aggr.rank <- t(apply(aggr.rank, 1, function(x){ rank(x, ties.method = "random") } ))
#   
#   
#   outputs <- list(aggr.rank = aggr.rank, trees = models, coefficients = aggr.coeff, call = match.call())
#   
#   class(outputs) <- "PlackettLuce_Bagging"
#   
#   return(outputs)
# }


# #Attribute bagging
# attr.bagging.PL <- function(formula, nboot = 100, G = NULL, data = NULL, alpha = 0.05,  minsize = NULL, ... )
# {
#   data = as.data.frame(data)
#   n = nrow(data)
#   e = all.names(formula[[3]], functions = FALSE, max.names = -1L)
#   le = length(e)
#   items = as.character(dimnames(attr(G, "rankings"))[[2]])
#   models = list()
#   samples <- matrix(NA, ncol = round(le*0.333), nrow = nboot)
#   PLcoeff <- array(NA, dim = c(n, length(items), nboot), dimnames = list(c(1:n),items,c(1:nboot)))
#   for(b in 1:nboot){
#     s <- sample(1:le, size = round(le*0.333), replace = FALSE)
#     newformula <- paste0("G ~ ", paste(e[s], collapse = " + "))
#     tree <- PlackettLuce::pltree(newformula, data = data[ , e[s] ], minsize = minsize, alpha = alpha)
#     models[[b]] <- tree
#     PLcoeff[,,b] <- predict(tree, data)
#     samples[b,] <- e[s]
#   }
#   aggr.rank <- apply(PLcoeff, c(1,2), function(x) {mean(x, na.rm=TRUE)} )
#   aggr.rank <- aggr.rank * -1
#   aggr.rank <- t(apply(aggr.rank, 1, function(x){ rank(x, ties.method = "random") } ))
#   
#   outputs <- list(aggr.rank = aggr.rank, trees = models, coefficients = PLcoeff, samples = samples, call = match.call())
#   class(outputs) <- "Attribute_bagging"
#   return(outputs)
# }
# 
# 
# #Calculate the percentual of times the model predicted the winning item 
# win.PL <- function(model.rank = NULL, observed.rank = NULL, drop.local = FALSE)
# {
#   if(drop.local) observed.rank <- observed.rank[,dimnames(observed.rank)[[2]]!="Local"]
#   
#   n = nrow(observed.rank)
#   
#   win <- rep(NA, times = n)
#   
#   for(k in 1:n){
#     o = observed.rank[k, observed.rank[k,] > 0] 
#     w <- o[model.rank[k,names(o)] == min(model.rank[k,names(o)])] == 1
#     win[k] <- (w/length(w))[1]
#   }
#   
#   return(mean(win, na.rm=TRUE)) 
# }