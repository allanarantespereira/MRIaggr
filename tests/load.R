require(data.table)
require(Rcpp)

#path.git <- "/home/brice/Documents/GitHub/MRIaggr"
path.git <- getwd()

excludeRfiles <- c("Generic_Functions.R","ClassMRIaggr_A_object.R","RcppExports.R")
source(file.path(path.git,"R","Generic_Functions.R"))
source(file.path(path.git,"R","ClassMRIaggr_A_Object.R"))
vecRfiles <- setdiff( list.files(file.path(path.git,"R")), excludeRfiles)
sapply(vecRfiles, function(x){source(file.path(path.git,"R",x))})

Rcpp:::sourceCpp(file.path(path.git,"src/Functions_Filtering.cpp"))
Rcpp:::sourceCpp(file.path(path.git,"src/Functions_Hemisphere.cpp"))
Rcpp:::sourceCpp(file.path(path.git,"src/Functions_Potential.cpp"))
Rcpp:::sourceCpp(file.path(path.git,"src/Functions_Potts.cpp"))
Rcpp:::sourceCpp(file.path(path.git,"src/Functions_W.cpp"))
