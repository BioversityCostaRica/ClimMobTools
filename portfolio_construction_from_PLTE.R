# Portfolio construction based on ensembles of PL trees
# Jacob van Etten, Nov 2017
library(PlackettLuce)

#load("C:/Users/Jacob Van Etten/Dropbox/@Other/ClimMob_CentralAmerica/output/red_beans_nicaragua_PL/sub_bagging_CentralAmerica.RData")
load("C:/Users/KAUE/Dropbox (Bioversity CR)/ClimMob_CentralAmerica/output/red_beans_nicaragua_PL/sub_bagging_CentralAmerica.RData")

#Predict function that works for any ensemble (list of models)
#Aggregating is done by taking the mean
#x is a list of PL trees (ensemble)
#newdata is a dataframe with the variables used to predict for new cases
#... passes arguments to the predict function
predict_ensemble <- function(x, newdata, aggregate_function = mean, ...){
  
  ntrees <- length(x$trees)
  vars <- names(coef(x$trees[[1]]))
  nvars <- length(vars)
  array_result <- array(NA, dim= c(nrow(newdata), nvars, ntrees))
  
  for(i in 1:ntrees){
    tree <- x$trees[[i]]
    pr <- predict(tree, newdata = newdata, ...)
    array_result[,,i] <- pr
    
  }
  
  result <- t(apply(array_result, c(1,2), aggregate_function))
  rownames(result) <- vars
  colnames(result) <- rownames(newdata)
  return(result)
  
}

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

#For testing, make a more manageable object
x <- NULL
x$trees <- sub_bagging_clim$trees[1:10] #10 trees, 11 varieties
rm(sub_bagging_clim)

#Create tiny dataset of ten rows
#These are now farms
#For robust portfolios, these should become different seasons for one farm or pixel
# newdata <- read.csv("C:/Users/Jacob Van Etten/Dropbox/@Other/ClimMob_CentralAmerica/processing/climmob_nicaragua28Nov2017.csv")
newdata <- read.csv("C:/Users/KAUE/Dropbox (Bioversity CR)/ClimMob_CentralAmerica//processing/climmob_nicaragua28Nov2017.csv")
newdata <- newdata[sample(1:nrow(newdata), 9),] #9 rows (= farms)

#Make the portfolios of 1, 2, and 3 varieties
pf1 <- create_portfolio(x, newdata, 1)  
pf2 <- create_portfolio(x, newdata, 2)
pf3 <- create_portfolio(x, newdata, 3)
