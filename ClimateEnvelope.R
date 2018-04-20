#Jacob van Etten, December 2017
#Script to plot climate envelopes
#Superimposing number of varieties within x% of optimum

# Portfolio construction based on ensembles of PL trees
# Jacob van Etten, Nov 2017
library(PlackettLuce)
library(partykit)
library(bamlss)
library(akima)
library(ggplot2)
library(plotly)
library(devtools)
install_github("ropensci/plotly")
library(sp)

load("C:/Users/Jacob Van Etten/Dropbox/@Other/ClimMob_CentralAmerica/output/red_beans_nicaragua_PL/sub_bagging_CentralAmerica.RData")
#load("C:/Users/KAUE/Dropbox (Bioversity CR)/ClimMob_CentralAmerica/output/red_beans_nicaragua_PL/sub_bagging_CentralAmerica.RData")

source("C:/Users/Jacob Van Etten/Dropbox/vanEtten_etal2018_replication_data/script/ClimMobTools.R")

#For testing, make a more manageable object
x <- sub_bagging_clim$trees
#rm(sub_bagging_clim)

# What are the two most important explanatory variables?
tvars <- table(unlist(lapply(x, function(x) {a <- as.vector(extract.from.nodes(x)[,,1]); a <- a[!is.na(a)]; return(a)})))
important_vars <- names(tvars)[order(tvars, decreasing=TRUE)][1:10]

# Get data and plot multi-dimensional scaling of variables
dd <- sub_bagging_clim$trees[[1]]$data
mds1 <- cmdscale(dist(t(dd[,colnames(dd)[colnames(dd) %in% important_vars]])))
plot(mds1)
text(mds1[,1], mds1[,2], rownames(mds1), cex = 0.9, xpd = TRUE)

#Two most distinctive important variables, one temperature, one rainfall
selected_vars <- c("maxNT_30", "Rx5day")

# Make prediction
pp <- predict_ensemble(x, dd)

#Combine the datasets of variety prob of winning and two climate variables
ddc <-  cbind(t(pp), dd[selected_vars])
ddc[,13] <- ddc[,13] * (7/120)

# Make interpolations of prob of winning for each of the varieties
s <- list(length=11)

for(i in 1:11) {
  
  ddc_xy <- ddc[,c(12,13,i)]
  coordinates(ddc_xy) <- ~ maxNT_30 + Rx5day
  s_i <- interp(ddc_xy, z=colnames(ddc)[i], nx=200, ny=200, linear=TRUE, duplicate="mean")
  s_i <- interp2xyz(s_i, data.frame=TRUE)
  s[[i]] <- s_i$z
  
}

s <- as.data.frame(s)
s$x <- s_i$x
s$y <- s_i$y
colnames(s) <- c(paste("Var", 1:11, sep="_"), "x", "y")

s_ma <- as.matrix(s[c("x", "y", "Var_1")])

p <- plot_ly(showscale = FALSE) %>%
  add_surface(z = ~s_ma, opacity = 0.98) 
p

ggplot(s) + 
  aes(x = x, y = y/(7/120), z = z, fill = z) + 
  geom_tile() + 
# geom_contour(color = "white", alpha = 0.5) + 
  scale_fill_distiller(palette="Spectral", na.value="white") + 
  theme_bw()

library(ggplot2)

geom_smooth(mapping = NULL, data = NULL, stat = "smooth",
            position = "identity", ..., method = "auto", formula = y ~ x,
            se = TRUE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE)
