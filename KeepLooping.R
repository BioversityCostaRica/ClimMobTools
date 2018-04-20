library(Matrix)

## really large {length(<dense equivalent>) is beyond R's limits}:
op <- options(warn = 2) # warnings (e.g. integer overflow!) become errors:
n <- 50000L

Lrg <- new("dgTMatrix", Dim = c(n,n))
l0 <- as(as(Lrg, "lMatrix"), "lgCMatrix")

for(i in 1:10)  dl0 <- try(as(l0, "denseMatrix"))
