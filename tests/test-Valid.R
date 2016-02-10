
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

#### validDimension ####

validDimension(value1 = 1:5, value2 = 1:5, method = "test")
validDimension(value1 = 1:5, value2 = 1:6, method = "test")
validDimension(value1 = 1:5, value2 = 1:5, type = "length" , method = "test")
validDimension(value1 = 1:5, value2 = 1:6, type = "length" , method = "test")

validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), type = "length" , method = "test")
validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), type = "nrow" , method = "test")
validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), method = "test")
validDimension(value1 = matrix(1:4,2,2), value2 = matrix(1,2,2), method = "test")