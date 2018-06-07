#Tools for environmental data gathering
#Kaue de Sousa 
#First run 31 Nov 2017
#Updated in 11 May 2018


# Tools for environmental data gathering ####

#Get time series environmental data from a defined timespan within sites
# X = a matrix with time series environmental data 
# start.date = a vector with the starting date to capture the climate information (generally the planting date)
# ts = an integer corresponding the timespan since the sowing date
# days.before = an integer corresponding how many days before starting date will be included in the timespan
# function return a matrix with time series data from the chosen period (starting date to the lenght of timespan)

get.timespan <- function(X, start.date, ts = 50, days.before = 0)
  {
  
  if(length(ts) > 1) maxts <- max(ts) else maxts <- ts
  if(is.matrix(X)) n <- nrow(X) else n <- 1
  if(is.matrix(X)) days <- as.character(dimnames(X)[[2]]) else days <- names(X)
  
  date <- as.data.frame(start.date)
  date <- match(as.character(date[,1]), days)
  date <- date - days.before
  Y <- NULL
  if(is.matrix(X)) {for(i in 0:maxts) Y <- cbind(Y, X[cbind(1:n, date+i)]) } #need a better solution
  if(!is.matrix(X)) {for(i in 0:maxts) Y <- cbind(Y, X[cbind(date+i)]) } #need a better solution
  
  if(length(ts) > 1){
    
    for(i in seq_along(ts)){
      Y[i,ts[i]:ncol(Y)] <- NA
    }
    
  }
  
  return(Y[,c(1:maxts)])
}


#Get growing degree days using temperature data
#GDD is a heat index that can be used to predict when a crop will reach maturity
# X = an array with two dimensions containing time series temperature data
#     1st dimension contains the day temperature and 2nd dimension the night temperature
# gdd = an integer indicating the average degree days for the crop species
# gdd.base = an integer indicating the base for GDD calculation
# function return a vector with days the site took to reach the GDD for the crop
get.GDD <- function(X, gdd = NULL, base = 10)
  {
  Y <- (((X[,,1] + X[,,2]) / 2) - base)
  
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
# DTR: diurnal temperature range, mean difference between DT and NT (degree Celsius)
# SU: summer days, number of days with maximum temperature > 30 C
# TR: tropical nights, number of nights with maximum temperature > 25 C


temp.index <- function(X, ts = NULL, index = NULL)
  {
  n <- dim(X)[1]
  if(is.null(index)) {
    ind <- as_tibble(matrix(nrow = n, ncol = 7,
                     dimnames = list(NULL, 
                                     c("maxDT","minDT","maxNT","minNT","SU","TR","DTR"))))
  }else{
      ind <- as_tibble(matrix(nrow = n, ncol = length(index),
                              dimnames = list(NULL, index)))
    }
  
  if(!is.na(match("maxDT", names(ind) ))) ind["maxDT"] <- apply(X[,,1], 1, function(X) max(X[1:ts], na.rm=TRUE )) #maximun day temperature (degree Celsius)
  if(!is.na(match("minDT", names(ind) ))) ind["minDT"] <- apply(X[,,1], 1, function(X) min(X[1:ts], na.rm=TRUE )) #minimum day temperature (degree Celsius)
  if(!is.na(match("maxNT", names(ind) ))) ind["maxNT"] <- apply(X[,,2], 1, function(X) max(X[1:ts], na.rm=TRUE )) #maximum night temperature (degree Celsius)
  if(!is.na(match("minNT", names(ind) ))) ind["minNT"] <- apply(X[,,2], 1, function(X) min(X[1:ts], na.rm=TRUE )) #minimum night temperature (degree Celsius)
  if(!is.na(match("DTR", names(ind) ))) ind["DTR"] <- apply((X[,,1]-X[,,2]), 1, function(X) mean(X[1:ts], na.rm = TRUE) ) #diurnal temperature range, mean difference between DT and NT (degree Celsius)
  if(!is.na(match("SU", names(ind) ))) ind["SU"] <- apply(X[,,1], 1, function(X) sum(X[1:ts] > 30, na.rm=TRUE)) #summer days, number of days with maximum temperature > 30 C
  if(!is.na(match("TR", names(ind) ))) ind["TR"] <- apply(X[,,2], 1, function(X) sum(X[1:ts] > 25, na.rm=TRUE)) #tropical nights, number of nights with maximum temperature > 25 C

  return(ind)
}


#Calculate evapotranspiration (ETo)
#details in http://www.fao.org/docrep/S2022E/s2022e07.htm#3.1.4%20calculation%20example%20blaney%20criddle
#X = an array with two dimensions containing time series temperature data
#     1st dimension contains the day temperature and 2nd dimension the night temperature
#p = mean daily percentage of annual daytime hours for different latitudes
#ts = the time span to calculate the total evapotranspitation
#Kc = crop factor for watter requirement

#function return the ETo in mm/day

#FURTHER improvements:
#implement p using latitude and planting date 
get.ETo <- function(X = NULL, p = 0.27,  lat = NULL, Kc = 1){
  
  #calculate Tmean
  Tmean <- rowMeans(cbind(apply(X[,,1], 1, function(X) max(X, na.rm=TRUE )), apply(X[,,2], 1, function(X) min(X, na.rm=TRUE ))))
  #reference evapotranspiration 
  eto <- p * (0.46 * Tmean + 8) * Kc
  
  return(eto)
  
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
     return(max(rle(Y)$lengths, na.rm = TRUE))
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
  MLWS <- apply(X, 1, function(X){
    Y <- X[1:ts]
    Y <- ifelse(Y >= 1, 0, runif(length(Y), min=2, max=10))
    return(max(rle(Y)$lengths, na.rm = TRUE))
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
      r5day <- cbind(r5day, sum(X[i:(i+4)], na.rm = TRUE))}
    return(max(r5day, na.rm = TRUE))
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
  if(!is.na(match("R5mm", names(ind) ))) ind["R5mm"]     <- apply(X, 1, function(X) sum(X[1:ts] >= 5 & X[1:ts] < 10, na.rm=TRUE) ) #days with  rainfall between 5-10 mm
  if(!is.na(match("R10mm", names(ind) ))) ind["R10mm"]   <- apply(X, 1, function(X) sum(X[1:ts] >= 10 & X[1:ts] < 15, na.rm=TRUE) ) #days with  rainfall between 10-15 mm
  if(!is.na(match("R20mm", names(ind) ))) ind["R20mm"]   <- apply(X, 1, function(X) sum(X[1:ts] >= 20, na.rm = TRUE) ) #days with  rainfall > 20 mm
  if(!is.na(match("SDII", names(ind) ))) ind["SDII"]     <- apply(X, 1, function(X) (sum(X[1:ts], na.rm = TRUE)/length(X[1:ts] > 0)) ) #simple rainfall intensity index (mean of rainy days / total rainfall)
  if(!is.na(match("Rx1day", names(ind) ))) ind["Rx1day"] <- apply(X, 1, function(X) max(X[1:ts], na.rm = TRUE)) #maximum 1-day rainfall
  if(!is.na(match("Rx5day", names(ind) ))) ind["Rx5day"] <- get.Rx5day(X, ts) #maximum 5-day rainfall
  if(!is.na(match("Rtotal", names(ind) ))) ind["Rtotal"] <- apply(X, 1, function(X) sum(X[1:ts], na.rm = TRUE) ) #total rainfall (mm) in wet days (r >= 1)
  
  ind[is.na(ind)] <- 0
  
  return(ind)
}


