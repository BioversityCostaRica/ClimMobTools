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
  
  #alternative: averaging the log values (= geometric mean)
  # preds_matrix <- apply(preds, c(1,2), function(x) exp(mean(log(x))))
  
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
    
    cat("Starting fold ", i, "Time: ", date(), "\n")
    
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
    
    cat("Fold ", i, " out of ", k, " completed. Kendall's correlation: ", r[1], "\n")

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


#### Example, using bean data from PL package
example("beans", package = "PlackettLuce")
G <- grouped_rankings(R, rep(seq_len(nrow(beans)), 4))
format(head(G, 2), width = 50)
beans$year <- factor(beans$year)

#Make sure to put the response variable in the dataset
beans$G <- G

#Create cluster to do parallel computation to speed things up
cluster <- makeCluster(2)
registerDoParallel(cluster)
getDoParWorkers()

#add vars
beans$ev1 <- 0
beans$ev2 <- 1
s <- beans$season
levels(s) <- list(dry=c("Ap - 15", "Ap - 16"), wet=c("Po - 15", "Po - 16", "Pr - 16"))
beans$s <- s

# crossvalidation_PLTE(G ~  .,
#                      d = beans[c("G", "ev1", "ev2")],
#                      k = 5,
#                      minsize = 30,
#                      alpha = 0.05,
#                      npseudo = 0.5,
#                      ntrees = 25,
#                      cluster = cluster)

# # k-fold cross-validation
# crossvalidation_PLTE(G ~ ., 
#                           d = beans[c("G", "season", "year", "maxTN")], 
#                           k = 10, 
#                           folds = NULL, 
#                           minsize = 30, 
#                           alpha = 0.05, 
#                           npseudo = 0.5, 
#                           ntrees = 25, 
#                           cluster = cluster)

#blocked cross-validation
#Make the folds first and calculate k
folds <- as.integer(beans$season)
k <- length(unique(folds))

folds1 <- folds
folds1[folds1 == 1 | folds1 == 4] <- 1
folds1[folds1 == 2 | folds1 == 3 | folds1 == 5] <- 2
k1 <- length(unique(folds1))

climmodel_blocked <- crossvalidation_PLTE(G ~ ., 
                     d = beans[c("G", "maxTN")], 
                     k = k, 
                     folds = folds, 
                     minsize = 30, 
                     alpha = 0.01, 
                     npseudo = 0.5, 
                     ntrees = 50, 
                     cluster = cluster)

nullmodel_blocked <- crossvalidation_PLTE(G ~ ., 
                                  d = beans[c("G", "ev1", "ev2")], 
                                  k = k, 
                                  folds = folds1, 
                                  minsize = 30, 
                                  alpha = 0.01, 
                                  npseudo = 0.5, 
                                  ntrees = 25, 
                                  cluster = cluster)

round(cbind(nullmodel_blocked$taus, climmodel_blocked$taus) *100)

stopCluster(cluster)

# Comparing two models
# Following Grandvalet, Y. and Bengio, Y., 2006. Hypothesis testing for cross-validation. Montreal Universite de Montreal, Operationnelle DdIeR, 1285.
# And Santafe, G., Inza, I. and Lozano, J.A., 2015. Dealing with the evaluation of supervised classification algorithms. Artificial Intelligence Review, 44(4), pp.467-508.

