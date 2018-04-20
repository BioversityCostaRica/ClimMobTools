#Jacob van Etten
#Make Mantel correlogram to determine if spatial blocks are needed
#And the size of these blocks if they are needed

library(PlackettLuce)
library(geosphere)
library(dplyr)
library(readr)
library(ncf)

#Prepare ClimMob data for analysis with PlackettLuce.
as.full <- function(x, local = FALSE, additional.rank = TRUE)
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
 
  result <- matrix(NA, nrow=length(itemnames), ncol=nrow(x))
  rownames(result) <- itemnames
  result[cbind(match(items[,1], itemnames), 1:nrow(items))] <- 1
  result[cbind(match(items[,2], itemnames), 1:nrow(items))] <- 2
  result[cbind(match(items[,3], itemnames), 1:nrow(items))] <- 3
  
  return(result)
  
}

cutlower <- c(0, 1000, 3000, 9000, 18000,54000)
cutupper <- c(cutlower[-1], 164000)

x <- read_csv("C:/Users/Jacob Van Etten/Dropbox/@Other/ClimMob_CentralAmerica/input/Data_shared_with_IoannisKosmodis/nicaragua_tricot.csv")
x <- x[!is.na(x$var_a),] 
x <- x[!is.na(x$var_b),] 
x <- x[!is.na(x$var_c),]
x <- dplyr::mutate(x, var_a = ifelse(var_a == "Peor", "Worse", "Better"), var_b = ifelse(var_b == "Peor", "Worse", "Better"), var_c = ifelse(var_c == "Peor", "Worse", "Better"))

dd <- as.full(x, local = TRUE)
ll <- as.matrix(x[c("lon", "lat")])
geodist <- distm(ll)
cd <- cor(dd, use="pairwise.complete.obs", method="kendall")
mm <- ncf::mantel.correlog(geodist, cd, resamp=999, increment=10000)
print(mm$p)
#Nicaragua:
#[1] 0.141 0.115 0.394 0.268 0.288 0.133 0.343 0.351 0.478 0.325
#[11] 0.282

