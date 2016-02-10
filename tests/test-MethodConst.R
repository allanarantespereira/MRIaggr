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

#### constCompressMRIaggr ####

MRIaggr.Pa1_compressed <- constCompressMRIaggr(MRIaggr.Pa1, compression.factor = 2,
                                               param = "DWI_t0")
MRIaggr.Pa1_compressed

MRIaggr.Pa1_compressed <- constCompressMRIaggr(MRIaggr.Pa1, compression.factor = 2,
                                               param = c("DWI_t0","TTP_t0"))
MRIaggr.Pa1_compressed

selectRegion(MRIaggr.Pa1_compressed)

#### constReduceMRIaggr ####
mask <- rbinom(selectN(MRIaggr.Pa1), size = 1, prob = 0.8)

MRIaggr.Pa1_reduced <- constReduceMRIaggr(MRIaggr.Pa1, mask = mask, numeric2logical = TRUE)
MRIaggr.Pa1_reduced

selectRegion(MRIaggr.Pa1_reduced)

selectContrast(MRIaggr.Pa1_reduced, param = "Lesion", subset = "Lesion")
selectContrast(MRIaggr.Pa1_reduced, param = "Lesion", subset = list(Lesion = "MASK_t0"))
selectContrast(MRIaggr.Pa1_reduced, param = "Lesion", subset = list(Lesion = "MASK_t2"))
