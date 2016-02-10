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

MRIaggr.Pa1

#### initCoords ####
#### initMask ####
#### initParameter ####
#### initSlice_var ####

#### initSubset ####

region1 <- rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)
region2 <- 5*rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.001)

allocContrast(MRIaggr.Pa1, 
              param = c("R1","R2","R3"),
              ls.MergeParam = list(nini = c("R1","R2","R3")),
              overwrite = TRUE) <- cbind(region1,region2,region1)
allocContrast(MRIaggr.Pa1, 
              param = c("R4","R5","R6"),
              ls.MergeParam = list(nini = c("R4","R5","R6")),
              overwrite = TRUE) <- cbind(region1,region2,region1)
# MRIaggr.Pa1@region$contrast

initSubset(MRIaggr.Pa1, subset = "Lesion")
initSubset(MRIaggr.Pa1, subset = "Lesion", operator.withinR = "none", operator.betweenR = "none")
initSubset(MRIaggr.Pa1, subset = "Lesion", operator.withinR = "intersect", operator.betweenR = "intersect")

initSubset(MRIaggr.Pa1, subset = "nini")
initSubset(MRIaggr.Pa1, subset = "nini", operator.withinR = "none", operator.betweenR = "none")
initSubset(MRIaggr.Pa1, subset = "nini", operator.withinR = "intersect", operator.betweenR = "intersect")

res <- initSubset(MRIaggr.Pa1, subset = list(nini = c("R1","R2","R3"),
                                             Lesion = c("MASK_t0","MASK_t2")), 
                  operator.withinR = "union", operator.betweenR = "union")
setdiff(unique(unlist(unlist(selectRegion(MRIaggr.Pa1)))), res )

res <- initSubset(MRIaggr.Pa1, subset = list(nini = c("R1","R2","R3"),
                                             Lesion = c("MASK_t0","MASK_t2")), 
                  operator.withinR = "intersect", operator.betweenR = "union")

res == intersect(selectRegion(MRIaggr.Pa1, region = "Lesion", region.value = "MASK_t0")[[1]],
                 selectRegion(MRIaggr.Pa1, region = "Lesion", region.value = "MASK_t2")[[1]])


