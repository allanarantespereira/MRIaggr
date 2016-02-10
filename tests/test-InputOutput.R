
require(data.table)

path.code <- "C:/Users/hpl802/Documents/Projects/Creation_package/Package_MRIaggr/MRIaggr"
path.code <- "/home/brice/Bureau/Creation_package/Package_MRIaggr/MRIaggr"


excludeRfiles <- c("Generic_Functions.R","ClassMRIaggr_A_object.R","RcppExports.R")
source(file.path(path.code,"R","Generic_Functions.R"))
source(file.path(path.code,"R","ClassMRIaggr_A_Object.R"))
vecRfiles <- setdiff( list.files(file.path(path.code,"R")), excludeRfiles)
sapply(vecRfiles, function(x){source(file.path(path.code,"R",x))})

Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Filtering.cpp"))
Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Hemisphere.cpp"))
Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Potential.cpp"))
Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Potts.cpp"))
Rcpp:::sourceCpp(file.path(path.code,"src/Functions_W.cpp"))

#### array2dt ####

#### matrix2array ####

set.seed(10)
n <- 4
Y <- rnorm(n^2)

## conversion
res1 <- df2array(expand.grid(1:n + 0.5, 1:n + 0.5), contrast = Y)
res2 <- df2array(expand.grid(1:n + 0.5, 1:n + 0.5), contrast = Y, format = "matrix")
res3 <- df2array(expand.grid(2 * (1:n), 2 * (1:n)), contrast = Y)
res4 <- df2array(expand.grid(2 * (1:n), 2 * (1:n)), contrast = cbind(Y ,Y, Y), range.coords = c(10,10))

## display
par(mfrow = c(2,2), mar = rep(2,4), mgp = c(1.5,0.5,0))
fields::image.plot(unique(res1$coords[[1]]), unique(res1$coords[[2]]), res1$contrast[[1]],
                   xlab = "", ylab = "")
fields::image.plot(unique(res2$coords[[1]]), unique(res2$coords[[2]]), res2$contrast,
                   xlab = "", ylab = "")
fields::image.plot(res3$contrast[[1]])
fields::image.plot(res4$contrast[[2]])



#### df2array ####

#### dt2array ####
