# require(data.table)
# 
# path.code <- "C:/Users/hpl802/Documents/Projects/Creation_package/Package_MRIaggr/MRIaggr"
# 
# excludeRfiles <- c("Generic_Functions.R","Class_MRIaggr_A_object.R","RcppExports.R")
# source(file.path(path.code,"R","Generic_Functions.R"))
# source(file.path(path.code,"R","Class_MRIaggr_A_object.R"))
# vecRfiles <- setdiff( list.files(file.path(path.code,"R")), excludeRfiles)
# sapply(vecRfiles, function(x){source(file.path(path.code,"R",x))})
# 
# Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Filtering.cpp"))
# Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Hemisphere.cpp"))
# Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Potential.cpp"))
# Rcpp:::sourceCpp(file.path(path.code,"src/Functions_Potts.cpp"))
# Rcpp:::sourceCpp(file.path(path.code,"src/Functions_W.cpp"))

#### data ####

path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")

nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0.nii"), format = "nifti")
nifti.Pat1_DWI_t0 <- readMRI(file.path(path.Pat1, "DWI_t0.nii"), format = "nifti")
nifti.Pat1_MASK_DWI_t0 <- readMRI(file.path(path.Pat1, "MASK_DWI_t0.nii"), format = "nifti")
nifti.Pat1_MASK_T2_FLAIR_t2 <- readMRI(file.path(path.Pat1, "MASK_T2_FLAIR_t2.nii"), 
                                       format = "nifti")


MRIaggr.Pa1 <- constMRIaggr(list(nifti.Pat1_TTP_t0, nifti.Pat1_DWI_t0,
                                 nifti.Pat1_MASK_DWI_t0, nifti.Pat1_MASK_T2_FLAIR_t2),
                            format = "MRIaggr",
                            ls.MergeParam = list(Lesion = c("MASK_t0","MASK_t2")),
                            identifier= "Pat1", default_value = "first",
                            param=c("TTP_t0","DWI_t0","MASK_t0","MASK_t2")
)

#### allocClinic ####

allocClinic(MRIaggr.Pa1) <- data.frame(age = 55, Sex = "male", stroke = TRUE, stringsAsFactors = FALSE)
allocClinic(MRIaggr.Pa1, overwrite = TRUE) <- data.frame(age = 56, stroke = TRUE, stringsAsFactors = FALSE)

#### allocContrast ####

allocContrast(MRIaggr.Pa1, param = "noise", overwrite = TRUE) <- rnorm(selectN(MRIaggr.Pa1))
allocContrast(MRIaggr.Pa1, param = c("noise","noise2"), overwrite = TRUE) <- cbind(rnorm(selectN(MRIaggr.Pa1)),
                                                                                   rnorm(selectN(MRIaggr.Pa1)))

region1 <- rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
region2 <- 5*rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
allocContrast(MRIaggr.Pa1, 
              param = c("noise3","noise5"),
              ls.MergeParam = list(nini = c("noise3","noise5")),
              overwrite = TRUE) <- cbind(region1,region2)
              

allocContrast(MRIaggr.Pa1, param = c("noise4","noise7"),
              ls.MergeParam = list(nini = c("noise4","noise7")),
              overwrite = TRUE) <- cbind(region1,
                                         5*rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
              )
unlist(selectContrast(MRIaggr.Pa1, param = "nini"))
selectContrast(MRIaggr.Pa1, param = "nini", subset = "nini")

selectContrast(MRIaggr.Pa1, param = "nini", subset = list(nini = "noise5"), rowNumber = TRUE, coords = TRUE)

#### allocDescStats ####

#### allocHemisphere ####
#### allocNormalization ####
#### allocTable ####
#### allocW ####

#### supprContrast ####
supprContrast(MRIaggr.Pa1) <- "nini"
# selectRegion(MRIaggr.Pa1)
supprContrast(MRIaggr.Pa1) <- "noise"

#### supprDescStats ####
