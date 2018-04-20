#Tools for environmental data gathering
#Kauê de Sousa 
#First run 31 Nov 2017
#Updated in 01 Feb 2018


# Tools for environmental data gathering ####

#Get time series environmental data from a defined timespan within sites
# X = a matrix with time series environmental data 
# pdata = a vector with sowing dates name 
# ts = an integer corresponding the timespan since the sowing date
# days.before = an integer corresponding how many days before sowing date will be included in the timespan
# function return a matrix with time series data from the chosen period
get.timespan <- function(X, pdate, ts = 50, days.before = 0)
  {
  n <- nrow(X)
  date <- as.data.frame(pdate)
  date <- match(as.character(date[,1]), as.character(dimnames(X)[[2]]))
  date <- date - days.before
  Y <- NULL
  for(i in 0:ts) Y <- cbind(Y, X[cbind(1:n, date+i)]) #need a better solution
  return(Y[,c(1:ts)])
}

#Get growing degree days using temperature data
#GDD is a heat index that can be used to predict when a crop will reach maturity
# X = an array with two dimensions containing time series temperature data
#     1st dimension contains the day temperature and 2nd dimension the night temperature
# gdd = an integer indicating the average degree days for the crop species
# gdd.base = an integer indicating the base for GDD calculation
# function return a vector with days the site took to reach the GDD for the crop
get.GDD <- function(X, gdd = NULL, gdd.base = 10)
  {
  Y <- (((X[,,1] + X[,,2]) / 2) - gdd.base)
  
  Y <- apply(Y, 1, function(x){
    for(d in 1:length(x)){
      i=d
      if(sum(x[1:d]) > gdd) break}
    return(i)
  })
  
  return(Y)
}

#Calculate time series temperature indices
# X = an array with two dimensions containing time series temperature data
#     1st dimension contains the day temperature and 2nd dimension the night temperature
#ts = an integer indicating the timespan for MLDS calculation
#index = character especifying the required indices, if null all available indices are calculated 
#return a data frame with especified temperature indices
# maxDT: maximun day temperature (°C)
# minDT: minimum day temperature (°C)
# maxNT: maximum night temperature (°C)
# minNT: minimum night temperature (°C)
# DTx25: days with temperature 25 >= T < 30 °C (days)
# DTx30: days with temperature 30 >= T < 35 °C
# DTX35: days with temperature T >= 35 °C
# NTb10: nights with temperature T < 10 °C
# NTx10: nights with temperature 10 >= T < 15 °C
# NTx15: nights with temperature 15 >= T < 20 °C
# NTx20: nights with temperature T >= 20 °C
# DTR: daily temperature range (max temperature - min temperature)
# NTR: night temperature range (max temperature - min temperature)

temp.index <- function(X, ts = NULL, index = NULL)
  {
  n <- dim(X)[1]
  if(is.null(index)) {
    ind <- as_tibble(matrix(nrow = n, ncol = 13,
                     dimnames = list(NULL, c("maxDT","minDT","maxNT","minNT","DTx25","DTx30","DTX35","NTb10","NTx10","NTx15","NTx20","DTR","NTR"))))
  }else{
      ind <- as_tibble(matrix(nrow = n, ncol = length(index),
                              dimnames = list(NULL, index)))
    }
  
  if(!is.na(match("maxDT", names(ind) ))) ind["maxDT"] <- apply(X[,,1], 1, function(X) max(X[1:ts], na.rm=TRUE )) #maximun day temperature (°C)
  if(!is.na(match("minDT", names(ind) ))) ind["minDT"] <- apply(X[,,1], 1, function(X) min(X[1:ts], na.rm=TRUE )) #minimum day temperature
  if(!is.na(match("maxNT", names(ind) ))) ind["maxNT"] <- apply(X[,,2], 1, function(X) max(X[1:ts], na.rm=TRUE )) #maximum night temperature
  if(!is.na(match("minNT", names(ind) ))) ind["minNT"] <- apply(X[,,2], 1, function(X) min(X[1:ts], na.rm=TRUE )) #minimum night temperature
  if(!is.na(match("DTx25", names(ind) ))) ind["DTx25"] <- apply(X[,,1], 1, function(X) sum(X[1:ts] >= 25 & X[1:ts] < 30, na.rm=TRUE)) #days with maximum temperature > 25 < 30 C
  if(!is.na(match("DTx30", names(ind) ))) ind["DTx30"] <- apply(X[,,1], 1, function(X) sum(X[1:ts] >= 30 & X[1:ts] < 35, na.rm=TRUE)) #days with maximum temperature > 30 < 35 C
  if(!is.na(match("DTx35", names(ind) ))) ind["DTx35"] <- apply(X[,,1], 1, function(X) sum(X[1:ts] >= 35, na.rm=TRUE)) #days with maximum temperature > 35 C
  if(!is.na(match("NTb10", names(ind) ))) ind["NTb10"] <- apply(X[,,2], 1, function(X) sum(X[1:ts] < 10, na.rm=TRUE)) #nights with temperature bellow 10 C
  if(!is.na(match("NTx10", names(ind) ))) ind["NTx10"] <- apply(X[,,2], 1, function(X) sum(X[1:ts] >= 10 & X[1:ts] < 15, na.rm=TRUE)) #nights with temperature between 10-15 C
  if(!is.na(match("NTx15", names(ind) ))) ind["NTx15"] <- apply(X[,,2], 1, function(X) sum(X[1:ts] >= 15 & X[1:ts] < 20, na.rm=TRUE)) #nights with temperature between 15-20 C
  if(!is.na(match("NTx20", names(ind) ))) ind["NTx20"] <- apply(X[,,2], 1, function(X) sum(X[1:ts] >= 20, na.rm=TRUE)) #nights with temperature above 20 C
  if(!is.na(match("DTR", names(ind) ))) ind["DTR"] <- apply(X[,,1], 1, function(X) max(X[1:ts], na.rm=TRUE) - min(X[1:ts], na.rm=TRUE)) #daily temperature range max temperature - min temperature
  if(!is.na(match("NTR", names(ind) ))) ind["NTR"] <- apply(X[,,2], 1, function(X) max(X[1:ts], na.rm=TRUE) - min(X[1:ts], na.rm=TRUE)) #night temperature range max temperature - min temperature
  
  ind[is.na(ind)] <- 0 
  
  return(ind)
}

#Calculate the MLDS rainfall index
#MLDS is a rainfall index indicating the maximum length of consecutive dry days (r < 1 mm) during a given timespan
#X = a vector matrix with time series rainfall data 
#ts = an integer indicating the timespan for MLDS calculation
#function return an integer with MLDS index
get.MLDS <- function(X, ts = NULL)
  {
   MLDS <- apply(X, 1, function(X){
     Y <- X[1:ts]
     Y <- ifelse(Y < 1, 0, Y)
     return(max(rle(Y)$lengths))
     }  )
  return(MLDS)
}

#Calculate the MLWS rainfall index
#MLWS is a rainfall index indicating the maximum length of consecutive wet days (r >= 1 mm) during a given timespan
#X = a vector matrix with time series rainfall data 
#ts = an integer indicating the timespan for MLWS calculation
#function return an integer with MLWS index
get.MLWS <- function(X, ts = NULL)
  {
  MLWS <- apply(precip, 1, function(X){
    Y <- X[1:ts]
    Y <- ifelse(Y >= 1, 0, runif(length(Y), min=2, max=10))
    return(max(rle(Y)$lengths))
    }  )
  return(MLWS)
}

#Calculate Rx5day rainfall index
#Rxday is a rainfall index indicating the maximum 5-day precipitation (mm) during a given timespan
#X = a vector matrix with time series rainfall data 
#ts = an integer indicating the timespan for Rx5day calculation
#function return an integer with Rx5day index

get.Rx5day <- function(X, ts = NULL)
  {
  Y <- apply(X, 1, function(X){
    r5day <- NULL
    for(i in 1:(ts-4)){
      r5day <- cbind(r5day, sum(X[i:(i+4)]))}
    return(max(r5day))
    })
  return(Y)
}


#Calculate time series rainfall indices
# X = an matrix with rainfall data
#ts = an integer indicating the timespan for MLDS calculation
#index = character especifying the required indices, if null all available indices are calculated 
#return a data frame with especified rainfall indices
# MLDS: maximum length of consecutive dry days (< 1 mm)
# MLWS: maximum length of consecutive wet days (>= 1 mm)
# R5mm: days with  rainfall 5 >= R < 10 mm
# R10mm: days with  rainfall 10 >= R < 15 mm
# R20mm: days with  rainfall R > 20 mm
# SDII: simple rainfall intensity index (mean of wet days / total rainfall)
# Rx1day: maximum 1-day rainfall (mm)
# Rx5day: maximum 5-day rainfall (mm)
# Rtotal: total rainfall (mm) in wet days (R >= 1)


rainfall.index <- function(X, ts = NULL, index = NULL)
  {
  n <- dim(X)[1]
  
  if(is.null(index)) {
    ind <- as_tibble(matrix(nrow = n, ncol = 9,
                            dimnames = list(NULL, c("MLDS","MLWS","R5mm","R10mm","R20mm","SDII","Rx1day","Rx5day","Rtotal"))))
  }else{
    ind <- as_tibble(matrix(nrow = n, ncol = length(index),
                            dimnames = list(NULL, index)))
  }
  
  if(!is.na(match("MLDS", names(ind) ))) ind["MLDS"]     <- get.MLDS(X, ts) #maximum length of consecutive dry days (< 1 mm)
  if(!is.na(match("MLWS", names(ind) ))) ind["MLWS"]     <- get.MLWS(X, ts) #maximum length of consecutive wet days
  if(!is.na(match("R5mm", names(ind) ))) ind["R5mm"]     <- apply(X, 1, function(X) sum(X[1:ts] >= 5 & X[1:ts] < 10) ) #days with  rainfall between 5-10 mm
  if(!is.na(match("R10mm", names(ind) ))) ind["R10mm"]   <- apply(X, 1, function(X) sum(X[1:ts] >= 10 & X[1:ts] < 15) ) #days with  rainfall between 10-15 mm
  if(!is.na(match("R20mm", names(ind) ))) ind["R20mm"]   <- apply(X, 1, function(X) sum(X[1:ts] >= 20) ) #days with  rainfall > 20 mm
  if(!is.na(match("SDII", names(ind) ))) ind["SDII"]     <- apply(X, 1, function(X) (sum(X[1:ts])/length(X[1:ts] > 0)) ) #simple rainfall intensity index (mean of rainy days / total rainfall)
  if(!is.na(match("Rx1day", names(ind) ))) ind["Rx1day"] <- apply(X, 1, function(X) max(X[1:ts])) #maximum 1-day rainfall
  if(!is.na(match("Rx5day", names(ind) ))) ind["Rx5day"] <- get.Rx5day(X, ts) #maximum 5-day rainfall
  if(!is.na(match("Rtotal", names(ind) ))) ind["Rtotal"] <- apply(X, 1, function(X) sum(X[1:ts], na.rm = T) ) #total rainfall (mm) in wet days (r >= 1)
  
  ind[is.na(ind)] <- 0
  
  return(ind)
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