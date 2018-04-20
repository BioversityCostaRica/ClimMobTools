####################################################
###### PREDICT PLANTING DATES
#Kauê de Sousa
#<k.desousa(at)cgiar.org>
#Bioversity International - Costa Rica
#First run 30Nov2017
####################################################
library(readr)
library(lubridate)
library(dplyr)
library(raster)
library(survival)
library(ggplot2)
library(ggfortify)
library(RCurl)

#load ClimMob Tools from GitHub
#repository https://github.com/kauedesousa/ClimMobTools/
tools <- getURL("https://raw.githubusercontent.com/kauedesousa/ClimMobTools/pilot/ClimMobCovariates.R", ssl.verifypeer = FALSE)
eval(parse(text = tools))

#working directory
wd <- "C:/Users/KAUE/Dropbox (Bioversity CR)/vanEtten_etal2018_replication_data/"
setwd(wd)
#read data
mydata <- read_csv("./input/climmob_nicaragua31Jan2018.csv")
head(mydata)
#change planting date into as.Date
mydata$planting_date <- as.Date(mydata$planting_date, "%Y-%m-%d")
summary(mydata$planting_date)

#add the week
mydata$week <- as.integer(strftime(mydata$planting_date, format = "%V"))
#get day of the year 
mydata$day <- yday(mydata$planting_date)
#this data has status 1, meaning the sowing was observed
mydata$status <- 1
#add rainfall data
load("./input/chirps_nicaragua.RData")
#get precipitation for the timespan 
ts = 21
precip <- get.timespan(chirps, pdate = mydata$planting_date, ts = ts, days.before = ts)
precip <- rainfall.index(precip, ts = ts)

#load temperature
load("./input/modis_nicaragua.RData")
#temperature for the timespan
temp <- array(NA, c(nrow(mydata), ts, 2), dimnames = list(NULL,  NULL, c(1:2) ))
temp[,,1] <- get.timespan(modis_approx[,,1], pdate = mydata$planting_date, ts = ts)
temp[,,2] <- get.timespan(modis_approx[,,2], pdate = mydata$planting_date, ts = ts)
temp <- temp.index(temp, ts = ts, index = c("maxDT","minDT","maxNT","minNT"))

#remove NA's in observer's rankings
#precip <- precip[!is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c), ]
#temp <- temp[!is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c), ]
#mydata <- mydata[!is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c), ]

#select covariates for modelling exercise
mydata <- mydata[,c("status","day","week","lon","lat","planting_date","season","year","soil")]

#add rainfall covariates
mydata <- as_tibble(cbind(mydata, precip, temp))
rm(precip)
rm(temp)

#add spatial covariates
mydata$xy <- mydata$lon + mydata$lat
mydata$yx <- mydata$lon - mydata$lat

#filter data
#Rtotal == 0
train <- mydata[mydata$Rtotal > 0, ] 

#fit Cox Proportional Hazards Model
fit <- coxph(Surv(day, status) ~ season + maxDT + minDT + maxNT + minNT +
               SDII + MLDS + MLWS + R10mm + R5mm + Rx5day + Rx1day  + soil + xy + yx,  data = train  )
summary(fit)
cox_fit <- survfit(fit)
summary(cox_fit)
autoplot(cox_fit)
predict(fit, type="risk", se.fit = F)

#use this to predict planting dates backward
n = nrow(mydata)
years = c(2003:2016)
pdates <- as_tibble(matrix(NA, nrow = n, ncol = 15))
pdates[,15] <- mydata$planting_date
for(i in c(1:n)){
  print(i)
  for(j in 1:length(years)){
    #identify the season per observer
    if(as.character(mydata[i, "season"]) == "Pr") s <- seq(as.Date(paste0(years[j], "-06-01"), "%Y-%m-%d"), as.Date(paste0(years[j], "-07-30"), "%Y-%m-%d"), 1)
    if(as.character(mydata[i, "season"]) == "Po") s <- seq(as.Date(paste0(years[j], "-08-01"), "%Y-%m-%d"), as.Date(paste0(years[j], "-11-15"), "%Y-%m-%d"), 1)
    if(as.character(mydata[i, "season"]) == "Ap") s <- seq(as.Date(paste0(years[j], "-11-16"), "%Y-%m-%d"), as.Date(paste0(years[j]+1, "-01-31"), "%Y-%m-%d"), 1)
    #get precipitation and temperature for the timespan within season
    #X days (see ts above) previuous each estimated days within "s"
    precip <- get.timespan(chirps[rep(i,length(s)) ,], pdate = s, ts = ts, days.before = ts)
    precip <- rainfall.index(precip, ts = ts)
    temp <- array(NA, c(length(s), ts, 2), dimnames = list(NULL,  NULL, c(1:2) ))
    temp[,,1] <- get.timespan(modis_approx[rep(i,length(s)), ,1], pdate = s, ts = ts)
    temp[,,2] <- get.timespan(modis_approx[rep(i,length(s)), ,2], pdate = s, ts = ts)
    temp <- temp.index(temp, ts = ts, index = c("maxDT","minDT","maxNT","minNT"))
    #dataframe as input for predict
    test <- cbind(precip, temp, season = mydata[i,"season"], soil = mydata[i, "soil"], xy = mydata[i, "xy"], yx = mydata[i, "yx"])
    #predict and take the max "risk" 
    pdates[i, j ] <- s[which.max(predict(fit, newdata = test, type = "risk"))]
    
  }
}

pdates[1:14] <- lapply(pdates[1:14], function(X) as.Date(X, origin = "1970-01-01" ))

#get environmental data based on these predictions 
n = nrow(mydata)
gdd = 585 #average growing degree-days to panicle initiation
quant <- cbind( V1 = seq(15, 60, 15), V2 = c("_1st","_2nd","_3rd","_4th"))
climatology <- list()

for(i in c(1:14)){
  #temperature for the timespan of each site
  temp <- array(NA, c(n, 80, 2), dimnames = list(NULL,  NULL, c(1:2) ))
  temp[,,1] <- get.timespan(modis_approx[,,1], pdate = pdates[i], ts = 80)
  temp[,,2] <- get.timespan(modis_approx[,,2], pdate = pdates[i], ts = 80)
  precip <- get.timespan(chirps, pdate = pdates[i], ts = 80, days.before = 10)
  #add planting dates 
  climatology[[i]] <- as_tibble(cbind( pdates[i], GDD = get.GDD(temp, gdd, 10)))
  
  for(j in c(1:4)){
    x <- temp.index(temp, ts = as.integer(quant[j,1]))
    names(x) <- paste0(names(x),quant[j,2]) 
    y <- rainfall.index(precip, ts = as.integer(quant[j,1]) )
    names(y) <- paste0(names(y),quant[j,2])
    
    climatology[[i]] <- as_tibble(cbind(climatology[[i]], x, y))
    
  }
}

save(climatology, fit, pdates, file = "./processing/climatology_nicaragua.RData")



# #Other approach, returning the response (day) in predictions
# #fit Regression for a Parametric Survival Model
# fit <- survreg(Surv(day, status) ~ season + maxDT + minDT + maxNT + minNT + 
#                SDII + MLDS + MLWS + R10mm + R5mm + Rx5day + Rx1day  + soil + xy + yx,
#              data = mydata , scale = 1 , dist = "logistic")
# 
# summary(fit)
# 
# #predict as quantile
# pct <- 1:99/100
# pred_quant <- predict(fit, newdata = mydata[mydata$season=="Po" & mydata$year == 2015,], type="quantile", p = pct, se=T)
# 
# plot(colMeans(pred_quant$fit), pct, col = "black", type = "l",
#      xlab = "Day", ylab = "Proportion")# ,xlim=c(250,365))
# lines(colMeans(pred_quant$fit + 2*pred_quant$se.fit), pct, lty = 2)
# lines(colMeans(pred_quant$fit - 2*pred_quant$se.fit), pct, lty = 2)
# abline(h=0.5, col = "red")
# 
# as.Date(as.integer(mean(pred_quant$fit[,51])), origin = "2015-01-01")

# #predict as day (type = "response")
# fit_predict <- predict(fit, type="response")
# as.Date(fit_predict-1, origin = "2004-01-01") #day-of-year is zero based in this case


# #crap code 
# m <- as_tibble(rbind(cbind(date = c(seq(as.Date(paste0(i, "-06-01"), "%Y-%m-%d"), as.Date(paste0(i, "-07-30"), "%Y-%m-%d"), 1)),
#                            season = "Pr"),
#                      cbind(date = c(seq(as.Date(paste0(i, "-08-01"), "%Y-%m-%d"), as.Date(paste0(i, "-11-15"), "%Y-%m-%d"), 1)),
#                            season = "Po"),
#                      cbind(date = c(seq(as.Date(paste0(i, "-11-16"), "%Y-%m-%d"), as.Date(paste0(i+1, "-01-31"), "%Y-%m-%d"), 1)),
#                            season = "Ap")))
# m$date <- as.Date(as.integer(m$date), origin = "1970-01-01")

# #add year to season
# mydata$season <- as.factor(ifelse(grepl(2015, mydata$project), paste0(mydata$season,"-15"),
#                                   ifelse(grepl(2016, mydata$project), paste0(mydata$season,"-16"),
#                                          NA)))
# summary(mydata$season)

# mydata$time1 <- 30
# mydata$time2 <- 45
# head(mydata)
# #I'll generate data with status 0, in the 30 days previous the planting date in intervals of 10
# mpseudo <- NULL
# for(i in c(30,15)){
#   x <- mydata
#   x$planting_date <- x$planting_date - i
#   x$time1 <- if(i == 30) 0 else 15
#   x$time2 <- if(i == 30) 15 else 30
#   x$status <- 0
#   mpseudo <- as_tibble(rbind(mpseudo, x))
# }
# 
# mydata <- as_tibble(rbind(mydata, mpseudo))

# mydata$date <- as.integer(mydata$planting_date)
# #matrix for predicted dates
# aggr_predic <- matrix(NA, nrow = n, ncol = 100)
# for(b in 1:100){
#   s <- sample(1:n, replace = TRUE)
#   train <- mydata[ s, ]
#   test  <- mydata[-s, ]
#   X <- rpart(date ~ SDII + Rx1day + Rx5day + lon + lat + week + season, data = train)
#   aggr_predic[-s,b] <- round(predict(X, test))
# }
# 
# Y <- as.data.frame(rbind(cbind(date = mydata$date,
#                                Type = "Observed"),
#                          cbind(date = round(rowMeans(aggr_predic, na.rm=TRUE)),
#                                Type = "RPART")))
# 
# Y$date <- as.integer(as.character(Y$date))
# Y$date <- as.Date(Y$date, origin = "1970-01-01")
# Y$week <- format(Y$date, format = "%U")
# Y$year <- paste(as.integer(strftime(Y$date, "%Y")), Y$week, sep="-")
# 
# ggplot(data = Y) +
#   geom_line(aes(x = date, y = index , color = Type ))
# 

# #generate pseudo-absence data for planting dates
# names(mydata)
# absence <- data.frame(planting_date= as.Date(as.integer(runif(n*3, (min(mydata$planting_date)-60), (max(mydata$planting_date)+60))), origin = "1970-01-01"),
#                       lon=rep(mydata$lon,3), 
#                       lat=rep(mydata$lat, 3),
#                       season="No season")
# 
# #bind thw two dataframes
# mydata <- rbind(mydata[,names(absence)], absence)
# mydata$sowing <- ifelse(mydata$season=="No season",0,1)
# head(mydata)
# mydata$year <- as.integer(strftime(mydata$planting_date, "%Y"))
# summary(factor(mydata$year))
# min(mydata$planting_date)
# max(mydata$planting_date)
# plot(summary(as.factor(absence$planting_date)))
# #bind thw two dataframes
# mydata <- rbind(mydata[,names(absence)], absence)
# mydata$sowing <- ifelse(mydata$season=="No season",0,1)
# head(mydata)
# mydata$year <- as.integer(strftime(mydata$planting_date, "%Y"))
# summary(factor(mydata$year))
# min(mydata$planting_date)
# max(mydata$planting_date)


# # Logistic Regression Model
# logit = glm(season ~ date + year + sdii_before + rx1day_before  + mlds_before + mlws_before, 
#             data = train, family=binomial)
# summary(logit) # AIC = 1824.3
# logit$aic
# predict <- predict(logit, type="response")
# predict
# 
# # Predictions on the test set
# predictTest = predict(logit, type="response", newdata=test)
# 
# 
# # Confusion matrix with threshold of 0.5
# tableTest <- table(test$season, predictTest > 0.35)
# tableTest
# 
# # To check how accurate is your model, divide the
# # true positives + true negatives by your N
# (tableTest[1,1]+tableTest[2,2]) / sum(tableTest)
# 
# # The baseline model would be to always predict the most frequent
# # outcome, in this case 0.
# (tableTest[1,1]+tableTest[1,2]) /  sum(tableTest)
# 
# library(ROCR)
# ROCRpred <- ROCR::prediction(predictTest, test$season)
# ROCRperf <- ROCR::performance(ROCRpred, "tpr","fpr")
# plot(ROCRperf, colorize = T, text.adj = c(-0.2,1.7))
# 
# 
