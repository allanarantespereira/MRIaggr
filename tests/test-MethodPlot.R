require(data.table)

path.code <- "C:/Users/hpl802/Documents/Projects/Creation_package/Package_MRIaggr/MRIaggr"
# path.code <- "/home/brice/Bureau/Creation_package/Package_MRIaggr/MRIaggr"


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

region1 <- rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
region2 <- 5*rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
allocContrast(MRIaggr.Pa1, 
              param = c("noise3","noise5"),
              ls.MergeParam = list(nini = c("noise3","noise5")),
              overwrite = TRUE) <- cbind(region1,region2)

# selectContrast(MRIaggr.Pa1, subset = list(Lesion = "MASK_t2"))[["Lesion"]]
#### multiplotMRI ####

multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("i","j","k"))
multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("i","j","k"), xlim = c(10,45))
multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("i","j","k"), xlim = c(10,45), asp = NULL)
multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("i","j","k"), xlim = c(10,45), ylim = c(10,45))

orthoplot(MRIaggr.Pa1, param = "TTP_t0")
orthoplot(MRIaggr.Pa1, param = "DWI_t0")

multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("i","j","k"), index1 = "Lesion")
multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_k = 2:3, slice_var = c("i","j","k"), index1 = "Lesion")
multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_k = 2:3, slice_var = c("i","j","k"), 
          index1 = list(subset = "Lesion", outline = TRUE))

multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_k = 2:3, slice_var = c("i","j","k"), 
          index1 = list(subset = list(Lesion = "MASK_t0"), outline = TRUE),
          index2 = list(subset = list(Lesion = "MASK_t2"), outline = TRUE)
)


df.contrast <- selectContrast(MRIaggr.Pa1, param = "TTP_t0", format = "vector", slice_k = 1)
df.coords <- selectCoords(MRIaggr.Pa1, slice_k = 1)

selectContrast(MRIaggr.Pa1, param = "TTP_t0", subset = c("nini","Lesion"), slice_i = 5)

# multiplot(df.coords, data.table(df.contrast), slice_var = c("i","j","k"), index1 = "Lesion")

multiplot(MRIaggr.Pa1, param = "TTP_t0", slice_var = c("j","k","i"),
          breaks = seq(0,38,1))

# multiplot(MRIaggr.Pa1, param = "TTP_t0", index

#### plotMRI ####

df.contrast <- selectContrast(MRIaggr.Pa1, param = "TTP_t0", format = "vector", slice_k = 1)
df.coords <- selectCoords(MRIaggr.Pa1, slice_k = 1)

plotMRI(data.table(df.contrast), df.coords, breaks = quantile(df.contrast), palette = heat.colors(4),
        col = NULL, asp = 1, xlim = NULL, ylim = NULL, pch = NULL, cex = NULL, axes = TRUE,
        col.NA = "blue", pch.NA = 21, xlab = NULL, ylab = NULL, main = "xxx", cex.main = 1)

