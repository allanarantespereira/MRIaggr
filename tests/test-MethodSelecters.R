
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


#### selectClinic ####
#### selectContrast ####

selectContrast(MRIaggr.Pa1)
selectContrast(MRIaggr.Pa1, coords = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100)
selectContrast(MRIaggr.Pa1, slice_i = 2, coords = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100, slice_i = 3, coords = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100, slice_i = 1, coords = TRUE)

selectContrast(MRIaggr.Pa1, rowNumber = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100, rowNumber = TRUE)
selectContrast(MRIaggr.Pa1, slice_i = 2, rowNumber = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100, slice_i = 3, rowNumber = TRUE)
selectContrast(MRIaggr.Pa1, subset = 5:100, slice_i = 1, rowNumber = TRUE)

selectContrast(MRIaggr.Pa1,param = "DWI_t0")
selectContrast(MRIaggr.Pa1, norm_mu = "default_value", param = "DWI_t0")

selectContrast(MRIaggr.Pa1,param = c("DWI_t0", "TTP_t0"))
selectContrast(MRIaggr.Pa1, norm_mu = "default_value", param = c("DWI_t0", "TTP_t0"))

selectContrast(MRIaggr.Pa1, subset = "Lesion")
selectContrast(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0","MASK_t2")))
selectContrast(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t2")))
selectContrast(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0","MASK_t2")), operator.withinR = "intersect")

which(unlist(lapply(MRIaggr.Pa1@contrast$Lesion,length)) > 0)

selectContrast(MRIaggr.Pa1, param = "TTP_t0", subset = list(TTP_t0 = c(">",5)))
selectContrast(MRIaggr.Pa1, param = c("TTP_t0","Lesion"), subset = list(TTP_t0 = c(">",5), 
                                                            Lesion = NULL), 
               operator.betweenR = "intersect")

list(TTP_t0 = c(">",5), 
     Lesion = NULL)

#### selectCoords ####

selectCoords(MRIaggr.Pa1)
selectCoords(MRIaggr.Pa1, subset = 5:100)
selectCoords(MRIaggr.Pa1, slice_i = 2)
selectCoords(MRIaggr.Pa1, subset = 5:100, slice_i = 3)
selectCoords(MRIaggr.Pa1, subset = 5:100, slice_i = 1)

selectCoords(MRIaggr.Pa1, rowNumber = TRUE)
selectCoords(MRIaggr.Pa1, subset = 5:100, rowNumber = TRUE)
selectCoords(MRIaggr.Pa1, slice_i = 2, rowNumber = TRUE)
selectCoords(MRIaggr.Pa1, subset = 5:100, slice_i = 3, rowNumber = TRUE)
selectCoords(MRIaggr.Pa1, subset = 5:100, slice_i = 1, rowNumber = TRUE)

selectCoords(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0","MASK_t2")))

#### selectDefault_value  ####

selectDefault_value(MRIaggr.Pa1)
selectDefault_value(MRIaggr.Pa1, param = "Lesion", format = "vector")

#### selectDescStats ####

selectDescStats(MRIaggr.Pa1)
selectDescStats(MRIaggr.Pa1, name = "rowNumber")

#### selectFieldDim ####

selectFieldDim(MRIaggr.Pa1)
selectFieldDim(MRIaggr.Pa1, coords = "i", format = "vector")

#### selectHemispheres ####

selectHemispheres(MRIaggr.Pa1)

#### selectHistory ####

selectHistory(MRIaggr.Pa1)

#### selectIdentifier ####

selectIdentifier(MRIaggr.Pa1)

#### selectMidplane ####

selectMidplane(MRIaggr.Pa1)

#### selectSubset ####

selectSubset(MRIaggr.Pa1, rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, slice_i = 2:3, rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, subset = 5:100, rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, subset = 5:100, slice_i = 5:7, rowNumber = TRUE)

selectSubset(MRIaggr.Pa1, subset = "Lesion", rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, subset = "MASK_t2", rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0","MASK_t2")), rowNumber = TRUE)
selectSubset(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0")), rowNumber = TRUE)

selectSubset(MRIaggr.Pa1, param = "TTP_t0", subset = list(TTP_t0 = c(">",5)), rowNumber = TRUE)

#### selectN ####

selectN(MRIaggr.Pa1)
selectN(MRIaggr.Pa1, slice_i = 2:3)
selectN(MRIaggr.Pa1, subset = 5:100)
selectN(MRIaggr.Pa1, subset = 5:100, slice_i = 5:7)

selectN(MRIaggr.Pa1, subset = "Lesion")
selectN(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0","MASK_t2")))
selectN(MRIaggr.Pa1, subset = list(Lesion = c("MASK_t0")))


#### selectNormalization ####
#### selectParameter ####

selectParameter(MRIaggr.Pa1)
selectParameter(MRIaggr.Pa1, type = "coords")
selectParameter(MRIaggr.Pa1, type = "ls_descStats")
selectParameter(MRIaggr.Pa1, type = "clinic")

#### selectRegion ####

selectRegion(MRIaggr.Pa1)
selectRegion(MRIaggr.Pa1, type = "names")
selectRegion(MRIaggr.Pa1, region = "Lesion", type = "names")

selectRegion(MRIaggr.Pa1, region = "Lesion")
selectRegion(MRIaggr.Pa1, region = "Lesion", region.value = "MASK_t0")

selectRegion(MRIaggr.Pa1, region = "TTP_t0", region.value = c(">",5))


#### selectTable ####
#### selectVoxelDim ####

selectVoxelDim(MRIaggr.Pa1)
selectVoxelDim(MRIaggr.Pa1, format = "vector")

#### selectW ####


