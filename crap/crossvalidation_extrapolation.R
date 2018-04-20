#Jacob van Etten
#k-fold cross-validation with PL trees
library(PlackettLuce)
library(BradleyTerry2)
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

# Calculate distance to data used to create the model
d_convexhull <- function(x, y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  xt <- t(x)
  
  h <- solve(xt %*% x)
  n <- nrow(y)
  r <- rep(NA, times=n)
  for(i in 1:n){
    yi <- y[i,]
    yit <- t(yi)
    
    r[i] <- yit %*% h %*% yi
  }
  
  return(r)
}


#### Example, using bean data from PL package
example("beans", package = "PlackettLuce")
G <- grouped_rankings(R, rep(seq_len(nrow(beans)), 4))
format(head(G, 2), width = 50)
beans$year <- factor(beans$year)

#Make sure to put the response variable in the dataset
beans$G <- G

#add vars
beans$ev1 <- 0
beans$ev2 <- 1
s <- beans$season

#Prepare matrix with observations to later compare with predictions for each fold
observed <- beans$G[1:length(beans$G),, as.grouped_rankings = TRUE]
index_observed <- attr(observed, "index")
observed <- observed[1:length(observed),, as.grouped_rankings = FALSE]
observed[observed == 0] <- NA

folds <- levels(beans$season)
tau <- rep(NA, times=length(folds))
N <- rep(NA, times=length(folds))

#Crossvalidation of observations within convex hull
for(i in 1:length(folds)){
  
  fold_i <- folds[i]
  train <- beans[beans$season != fold_i,]
  #w <- 1 / table(train$season)[as.character(train$season)] #Weights break predict...
  mod <- pltree(G ~ maxTN, data=train)
  test <- beans[beans$season == fold_i,]
  d_max <- max(d_convexhull(train["maxTN"], train["maxTN"]))
  d <- d_convexhull(train["maxTN"], test["maxTN"])
  #test <- test[d < d_max,] #Comment and uncomment this line to see the difference
  preds <- predict(mod, newdata=test)
  preds_index <- as.integer(rownames(preds))
  
  #"shrink" observed ranks
  obs_index <- index_observed[index_observed %in% preds_index] #inside i-th fold and inside preds
  observed_shrunk <- observed[index_observed %in% preds_index,]
  
  #"expand" predicted coef 
  match_index <- match(obs_index, preds_index)
  predicted_expanded <- preds[match_index,]
  
  tau_i <- get_tau(observed_shrunk, predicted_expanded)[1]
  N_i <- sum(d < d_max)
  
  cat("N: ", N_i, " tau:", tau_i, "\n")

  tau[i] <- tau_i 
  N[i] <- N_i
  
}

sum(tau * N) / sum(N)
