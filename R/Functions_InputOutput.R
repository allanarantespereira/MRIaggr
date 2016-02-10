#**********************************************************************
#**********************************************************************
#*************         Input Output function        *******************
#**********************************************************************
#**********************************************************************
#

# readMRI:         ok
# extractMRI:      ok
# constMRIaggr:    ok
# array2dt:        TO BE UPDATED
# dt2array:        TO BE UPDATED to enable one argument dt and optional arguments for compatibility with df vector ? split into 2 functions ??
# writeMRI:        TO BE UPDATED


#### 1- Read ####

####>>> readMRI ####
readMRI <- function (filename, format = NULL, na.value = 0, 
                     what = "numeric", size = "NA_integer_", dimensions = NULL, 
                     SPM = FALSE, reorient = FALSE, flipud = FALSE, recover.header = FALSE, 
                     checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose")){
  
  #### preparation
  validFormat <-  c("rawb.gz", "analyze", "nifti", "dicom")
  validExtension <-  c("nii", "nii.gz", "img.gz", "img", "dcm", "rawb.gz")
  
  fct_matching <- function(x){
    format <- switch(x,
                     "nii" = "nifti",
                     "nii.gz" = "nifti",
                     "img.gz" = "analyze",
                     "img" = "analyze",
                     "dcm" = "dicom",
                     "rawb.gz" = "rawb.gz",
                     stop("readMRI: invalid file extension \n",
                          "proposed extension: \"",extension,"\" \n",
                          "valid extensions : \"nii\" \"img.gz\" \"img\" \"dcm\" \"txt\" \n")
    )
  }
  
  Ofilename <- filename
  
  #### differentiate filename from path
  Ofilename_splitSlash <- strsplit(Ofilename, split = .Platform$file.sep, fixed = TRUE)[[1]]
  n.Ofilename_splitSlash <- length(Ofilename_splitSlash)
  if(n.Ofilename_splitSlash > 1){
    directory <- paste(Ofilename_splitSlash[-n.Ofilename_splitSlash], collapse = .Platform$file.sep)
  }else{
    directory <- "."
  }
  directory <- paste(Ofilename_splitSlash[-n.Ofilename_splitSlash], collapse = .Platform$file.sep) # name of the directory
  Ffilename <- Ofilename_splitSlash[n.Ofilename_splitSlash] # full filename
  
  #### extract file extension
  res_tempo <- getFileExtention(Ffilename)
  filename <- res_tempo$filename
  extension <- res_tempo$extension
  extension.compression <- res_tempo$extension.compression
  
  # check extension
  if(!is.na(extension)){
    extension.format <- fct_matching(extension)
  }else if(checkArguments == TRUE){
    
    candidates <- list.files(directory, include.dirs = FALSE)
    candidates <- candidates[grep(pattern = filename, x = candidates, fixed = TRUE)] # candidates with the right name
    candidates.extension <- unlist(sapply(candidates,function(x){
      res_tempo <- getFileExtention(x)
      if(res_tempo$extension %in% validExtension){
        return(res_tempo$extension)}
    }
    ))
    
    stop("readMRI: extension of argument \'filename\' is not specified \n",
         "propose filename: ",filename,"\n",
         if(length(candidates.extension)>0){paste("consider adding: \".",paste(candidates.extension, collapse = "\" or \"."),"\" \n", sep = "")},
         sep = ""
    )
  }
  
  #### recognise format 
  if(is.null(format)){
    format <- extension.format
  }
  
  #### check arguments
  if(checkArguments == TRUE){
    
    validCharacter(value = Ffilename, validLength = 1, refuse.NULL = TRUE, method = "readMRI")
    
    validCharacter(value = format, validLength = 1, validValues = validFormat,
                   refuse.NULL = FALSE, method = "readMRI")
    
    if(format != extension.format){
      stop("readMRI: mismatch between proposed \'format\' the extension of argument \'filename\' \n",
           "propose format: ",format,"\n",
           "filename extension: ",extension.format," (\".",extension,"\")\n")
    }
    
    if( file.exists(Ofilename) == FALSE){
      
      candidates <- list.files(directory, include.dirs = FALSE)
      candidates <- candidates[order(adist(x = Ffilename, y = candidates), decreasing = FALSE)]
      n.candidates <- length(candidates)
      
      stop("readMRI : file not found \n",
           "proposed file: \"",Ffilename,"\" \n",
           "in directory: \"",directory,"\" \n",
           "Best matching filenames in this directory (",min(5,n.candidates)," over ", n.candidates," files): \"",
           paste(candidates[1:min(5,n.candidates)], collapse = paste("\" \n",constBlanck(29 + nchar(n.candidates)),"\"", sep = "")),"\" \n",
           "current working directory: \"",getwd(),"\" \n") 
    }
  }
  
  #### header 
  if(format == "analyze"){
    
    header_file <- paste(filename, ".hdr", sep = "")
    headerGZ_file <- paste(filename, ".hdr.gz", sep = "")
    
    test.header_file <- file.exists(file.path(directory, header_file))
    test.headerGZ_file <- file.exists(file.path(directory, headerGZ_file))
    
    test.header <- (test.header_file == TRUE && extension.compression == FALSE) + (test.headerGZ_file == TRUE && extension.compression == TRUE)
    
    if(test.header == FALSE && recover.header == TRUE && (test.header_file || test.headerGZ_file) ){ # deal with (.img and .hdr.gz) or (.img.gz and .hdr)
      
      initPackage(package = "R.utils", argument = "recover.header = TRUE", method = "readMRI")
      
      if(dir.exists("MRIaggr_readMRI_tempoAnalyse")){
        stop("readMRI : please remove the \"MRIaggr_readMRI_tempoAnalyse\" directory of the current working directory \n")
      }else{
        if(verbose > 0){
          cat("* copy .img",if(extension.compression){".gz"}," and .hdr",if(test.headerGZ_file){".gz"}," files to a temporary directory \n", sep = "")
        }
        dir.create("MRIaggr_readMRI_tempoAnalyse")
        on.exit( unlink("MRIaggr_readMRI_tempoAnalyse", recursive = TRUE) )
        file.copy(Ofilename, 
                  file.path("MRIaggr_readMRI_tempoAnalyse",Ffilename))
        file.copy(file.path(directory, if(test.headerGZ_file){headerGZ_file}else{header_file}), 
                  file.path("MRIaggr_readMRI_tempoAnalyse",if(test.headerGZ_file){headerGZ_file}else{header_file}))
        
        if(extension.compression == FALSE){
          if(verbose > 0){cat("* uncompress the .hdr.gz file \n")}
          R.utils:::gunzip( filename = file.path("MRIaggr_readMRI_tempoAnalyse", headerGZ_file), ext = "gz")
        }else{
          if(verbose > 0){cat("* uncompress the .img.gz file \n")}
          R.utils:::gunzip( filename = file.path("MRIaggr_readMRI_tempoAnalyse", Ffilename), ext = "gz")
        }
      }
      
      directory <- "MRIaggr_readMRI_tempoAnalyse"
      extension <- "img"
      extension.compression <- FALSE
      test.header <- TRUE
    }
    
    if(test.header  == FALSE){
      stop("readMRI : header file not found for file: \"",filename,"\" \n",
           "for analyse format corresponding \"",if(extension.compression){".hdr.gz"}else{".hdr"},"\" file is need \n",
           "but was not found in the current directory: \"",directory,"\" \n",
           "current working directory: \"",getwd(),"\" \n",
           if(test.header_file || test.headerGZ_file){"consider setting recover.header to TRUE \n"})
    }
    
    
  }
  
  #### summary
  if(verbose > 0){
    cat("* orignal filename: \"",Ofilename, "\" \n",
        "* new filename    : \"",filename, "\" \n",
        "* extension       : \"",extension,"\" \n",
        if(format == "analyze"){paste("* header          : \"",if(extension.compression){".hdr.gz"}else{".hdr"},"\" \n", sep = "")},
        "* format          : \"",format,"\" \n", sep = "")
  }
  
  #### reading ####
  if (format == "rawb.gz") {
    
    validInteger(value = dimensions, validLength = 3, 
                 refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method = "readMRI")
    
    f <- gzfile(Ofilename, open = "rb")
    on.exit(close(f))
    data <- readBin(con = f, what = what, n = prod(dimensions), size = size)
    res <- array(data, dimensions)
    
  } else if (format == "analyze") {
    
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "readMRI")
    res <- oro.nifti::readANALYZE(file.path(directory,filename), SPM = SPM)
    
  } else if (format == "nifti") {
    
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "readMRI")
    res <- oro.nifti::readNIfTI(file.path(directory,filename), reorient = reorient)
    
  } else if (format == "dicom")  {
    # initPackage(package = "oro.dicom", argument = "format = \"oro.dicom\"", method = "readMRI")
    res <- oro.dicom::readDICOMFile(Ofilename, flipud = flipud)
  }
  
  if(!is.na(na.value)){
    
    if(format == "dicom"){
      test.na <- is.na(res$img)
      if(any(test.na)){res$img[test.na] <- na.value}        
    }else{
      test.na <- is.na(res)
      if(any(test.na)){res[test.na] <- na.value}
    }
    
  }
  
  return(res)
  
}

####>>> extractMRI ####

extractMRI <- function (MRI, fix = TRUE, as.vector = FALSE){ # internal
  
  class.MRI <- class(MRI)
  dim.MRI <- dim(MRI)
  
  #### extraction
  if(class.MRI == "anlz"){ # analyse object
    dim.array <- MRI@dim_[setdiff( which(MRI@dim_>0), 1)]
    unit <- MRI@vox_units
    dim.voxel <-  MRI@pixdim[MRI@pixdim != 0] ### WARNING may be wrong: >0???
    
    if(as.vector == TRUE){
      array <- as.vector(MRI@.Data)
    }else{
      array <- MRI@.Data
    }
    
  }else if(class.MRI %in% c("nifti", "niftiExtension", "niftiAuditTrail")){ # nifti object
    dim.array <- MRI@dim_[setdiff( which(MRI@dim_>0), 1)]
    unit <- oro.nifti::convert.units(oro.nifti::xyzt2space(MRI@xyzt_units))
    dim.voxel <- MRI@pixdim[MRI@pixdim > 0] 
    
    if(as.vector == TRUE){
      array <- as.vector(MRI@.Data)
    }else{
      array <- MRI@.Data
    }
    
  }else if(is.list(array) && length(array) == 2 && "img" %in% names(array)){ # dicom object
    dim.array <- rep(NA, length(dim.MRI))
    unit <- NA
    dim.voxel <-  rep(NA, length(dim.MRI))
    
    if(as.vector == TRUE){
      array <- as.vector(MRI$img)
    }else{
      array <- MRI$img
    }
    
  }else if( class.MRI == "matrix") { # matrix object
    
    dim.array <- dim(MRI)
    unit <- NA
    dim.voxel <-  rep(NA, length(dim.MRI))
    
    if(as.vector == TRUE){
      array <- as.vector(MRI)
    }else{
      array <- array(MRI, dim = dim.array)
    }
    
  }else if( class.MRI == "array") { # array object
    
    dim.array <- dim(MRI)
    unit <- NA
    dim.voxel <-  rep(NA, length(dim.MRI))
    
    if(as.vector == TRUE){
      array <- as.vector(MRI)
    }else{
      array <- MRI
    }
    
  }
  
  #### check
  if(fix == TRUE){ # fix empty response by using available information
    if(length(dim.array) == 0){
      dim.array <- dim(array)
    }
    if(length(dim.voxel) == 0){
      dim.voxel <- rep(NA,length(dim.array))
    }
    dim.voxel <- data.frame(rbind(dim.voxel), stringsAsFactors = FALSE)
    names(dim.voxel) <- letters[9:(8+length(dim.voxel))]
    
    dim.array <- data.frame(rbind(dim.array), stringsAsFactors = FALSE)
    names(dim.array) <- letters[9:(8+length(dim.array))]
  }
  
  
  #### export
  return(list(array = array,
              dim.array = dim.array,
              unit = unit,
              dim.voxel = dim.voxel))
}

#### 2- Const ####

####>>> constMRIaggr ####

constMRIaggr <- function(ls.MRI, param, drop1 =  TRUE, ls.MergeParam = NULL, refValue = 0, region.overall = "Stroke",
                         format = "MRIaggr", identifier, default_value = NULL, voxelDim = NULL, unit = NULL,
                         checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), rm.ls.MRI = FALSE){
  
  
  #### preparation
  name.ls.MRI <- as.character(substitute(ls.MRI))
  class.ls.MRI <- class(ls.MRI)
  
  ## convertion to a valid class
  if(class.ls.MRI %in% c("anlz", "nifti", "niftiExtension", "niftiAuditTrail") || (is.list(ls.MRI) && length(ls.MRI) == 2 && "img" %in% names(ls.MRI)) ){
    ls.MRI <- list(ls.MRI)
  }
  n.ls.MRI <- length(ls.MRI)
  
  if(is.null(param)){
    param <- names(ls.MRI)
  }
  
  dim.MRI <- dim(ls.MRI[[1]])
  if(drop1 == TRUE){
    dim.MRI <- dim.MRI[rev(cumprod(rev(dim.MRI) == 1)) == 0]
  }
  ndim.MRI <- length(dim.MRI)
  
  ## default values 
  if( identical(default_value, "first") ){
    
    default_value <- matrix(rep(1, ndim.MRI), nrow = 1)
    
  }else if( identical(default_value, "first") ){
    
    default_value <- matrix(dim.MRI, nrow = 1)
    
  }
  
  ## merge parameters
  res <- initMergeParam(ls.MergeParam = ls.MergeParam, param = param, refValue = refValue,
                        checkArguments = checkArguments, init = TRUE, method = "constMRIaggr")
  name.merge <-  res$name.merge
  name2merge <-  res$name2merge
  n.merge <-res$n.merge
  ls.indexMerge <- res$ls.indexMerge
  
  #### test ####
  if (checkArguments) {
    
    test.class <- unlist(lapply(ls.MRI, function(x){class(x) %in% c("anlz", "nifti", "niftiExtension", "niftiAuditTrail", "list", "array", "matrix")}))
    if(any(test.class == FALSE)){ # class
      stop("constMRIaggr: wrong specification of \'ls.MRI\' \n", 
           "\'ls.MRI\' must be a list of \"anlz\" or \"nifti\" or \"list\" or \"array\" \n", 
           "elements of \'ls.MRI\' that are not of type array : ", paste(which(test.class == FALSE), collapse = " "), " \n")
    }
    
    test1.dim <- lapply(ls.MRI, dim)
    if( length(unique(unlist(lapply(test1.dim, length)))) > 1 ){ # number of dimension i.e. 2D, 3D or 4D
      stop("constMRIaggr: wrong specification of \'ls.MRI\' \n", 
           "\'ls.MRI\' must have elements with the same number of dimensions \n", 
           "differing number of dimensions: ",paste(unique(unlist(lapply(test1.dim, length))), collapse = " "),"\n")
    }
    
    M_tempo <- matrix(unlist(test1.dim), nrow =  n.ls.MRI, byrow = TRUE)
    test2.dim <- lapply(1:ndim.MRI, function(x){ unique(M_tempo[,x])})
    if( any( unlist(lapply(test2.dim, length)) > 1) ){ # equal size in each dimension
      stop("constMRIaggr: wrong specification of \'ls.MRI\' \n", 
           "\'ls.MRI\' must have elements with the same length within each dimension \n", 
           "dimensions with differing size: \n",
           paste0(unlist(lapply(1:ndim.MRI, function(x){paste("dimension ",letters[8+x],": ", paste(test2.dim[[x]], collapse = " "), "\n", sep = " ")})))[unlist(lapply(test2.dim, length)) > 1]
      )
    }
    
    if(!is.null(voxelDim) && length(voxelDim) > ndim.MRI){ 
      stop("constMRIaggr: length of \'voxelDim\' is lower compared to the number of dimensions of the images \n", 
           "length(voxelDim): ",paste(voxelDim, collapse = " "),"\n", 
           "number of dimensions of the images: ",ndim.MRI,"\n")
    }
    
    if(!is.null(voxelDim) && any( c(voxelDim, rep(1, ndim.MRI - length(voxelDim)) - dim.MRI > 0) ) ){ 
      stop("constMRIaggr: specification of \'voxelDim\' does not match the dimension of the images\n", 
           "proposed voxelDim : ",paste(voxelDim, collapse = " "),"\n", 
           "augmented voxelDim: ",paste(voxelDim, rep(1, ndim.MRI - length(voxelDim)), sep = "", collapse = " "),"\n",
           "dimension of the images': ",paste(dim.MRI, collapse = " "),"\n")
    }
    
    validCharacter(value = format, validLength = 1, validValues = c("data.table", "data.frame", "matrix", "MRIaggr"), refuse.NULL = TRUE, method = "constMRIaggr")
    if(format == "MRIaggr"){
      validCharacter(value = identifier, validLength = 1, refuse.NULL = TRUE, method = "constMRIaggr")
    }
    
    validDimension(value1 = param, validDimension = n.ls.MRI, name1 = "param", name2 = "ls.MRI", type = "length", method  = "constMRIaggr")
    
    if(!is.null(default_value)){
      
      if(!is.matrix(default_value)){
        stop("constMRIaggr: wrong specification of \'default_value\' \n", 
             "\'default_value\' must be either \"first\" or \"last\" or a matrix \n", 
             "is(default_value): ",paste(is(default_value), collapse =  " "),"\n")
      }
      
      validDimension(value1 = default_value, validDimension = ndim.MRI, name1 = "default_value", name2 = "ls.MRI", type = "ncol", method  = "constMRIaggr")
      
      test.default <- sweep(default_value, MARGIN = 2, STATS = dim.MRI, FUN = "/")
      if( any(test.default > 1) ){
        stop("constMRIaggr: wrong specification of \'default_value\' \n", 
             "\'default_value\' indicates an observation outside the image \n", 
             "observation: ",paste( which( rowSums(test.default) > 0 ), collapse =  " "),"\n")
      }
      
    }
    
  }
  
  #### index of the default values
  if(!is.null(default_value)){
    default_value <- 1 + rowSums(sweep(default_value-1, MARGIN = 2, STATS = dim.MRI, FUN = "*"))
  }
  
  #### extractions of the values 
  ## coordinates
  if(verbose){cat("* Extraction of the coordinates: ",paste(letters[9:(8+ndim.MRI)], collapse = " "))}
  dt.MRI <- data.table::data.table(expand.grid(lapply(1:ndim.MRI,function(x){1:dim.MRI[x]})))
  data.table::setnames(dt.MRI, old = names(dt.MRI), new = letters[9:(8+ndim.MRI)])
  
  ## default values
  param.default <- c(setdiff(param, unlist(ls.MergeParam)), names(ls.MergeParam))
  defaultMRI <- setNames( data.frame( matrix(NA, nrow = 1, ncol =  length(param.default)), stringsAsFactors = FALSE ), param.default)
  
  ## index list
  if(!is.null(ls.MergeParam)){
    dt.MRI[,name.merge := .(list()), with = FALSE]
  }
  
  if(verbose){cat("\n")}
  
  ## Merging
  if(verbose){cat("* Merging : ")}
  
  for(iter_MRI in 1:n.ls.MRI){
    
    if(verbose){cat("(", iter_MRI, ") ",param[iter_MRI]," | ", sep = "")}
    
    if((iter_MRI>1) && (param[iter_MRI] %in% name2merge == FALSE) && is.null(default_value)){
      
      dt.MRI[, param[iter_MRI] :=  extractMRI(ls.MRI[[iter_MRI]], as.vector = TRUE)$array, with = FALSE] 
      
    }else{
      
      res_tempo <- extractMRI(ls.MRI[[iter_MRI]])
      
      if(iter_MRI == 1){
        if(is.null(unit)){unit <- res_tempo$unit}
        if(is.null(voxelDim)){dim.voxel <- res_tempo$dim.voxel[1, 1:ndim.MRI, drop = FALSE]}
        
      }
      
      if(param[iter_MRI] %in% name2merge) {
        
        index.Merge <- which(unlist(lapply(ls.MergeParam, function(x){param[iter_MRI] %in% x})))
        ls.indexMerge[[index.Merge]] <- c(ls.indexMerge[[index.Merge]], setNames(list(which(res_tempo$array!=refValue)),  param[iter_MRI]))
        
        #         dt.MRI[ls.indexMerge[[index.Merge]][[param[iter_MRI]]], name.merge[index.Merge] :=  list(lapply(.SD[[1]], FUN = function(x){c(x, param[iter_MRI])})),
        #                .SDcol = name.merge[index.Merge], with = TRUE] 
        
        dt.MRI[ls.indexMerge[[index.Merge]][[param[iter_MRI]]], name.merge[index.Merge] :=  list(lapply(1:.N, FUN = function(x){
          eval(parse(text = paste(
            "c(.SD[[1]][[x]],list(",param[iter_MRI],"= res_tempo$array[ ls.indexMerge[[index.Merge]][[param[iter_MRI]]][x] ] ))"
          )))
        })),
        .SDcol = name.merge[index.Merge], with = TRUE] 
        
        if(!is.null(default_value)){
          
          default_tempo <- res_tempo$array[default_value]
          mode <- names(which.max(table(default_tempo))[1, drop = FALSE])
          if(is.numeric(default_tempo)){ mode <- as.numeric( mode ) }
          
          if( param[iter_MRI] == ls.MergeParam[[index.Merge]][1]){
            defaultMRI[name.merge[index.Merge]] <- mode 
            
          }else{
            merge_tempo <- unique(c(defaultMRI[1,name.merge[index.Merge]], mode))
            if(length(merge_tempo) > 1){
              warning("constMRIaggr: different default values founded for ",name.merge[index.Merge]," \n", 
                      "merged variables: ",paste(names(ls.indexMerge[[index.Merge]]), collapse = " "),"\n", 
                      "default values: ",paste(merge_tempo, collapse =  " "),"\n")
            }
            defaultMRI[name.merge[index.Merge]] <- merge_tempo[1]
          }
        }
        
      }else{
        
        dt.MRI[, param[iter_MRI] :=  as.vector(res_tempo$array), with = FALSE] 
        
        if(!is.null(default_value)){
          default_tempo <- res_tempo$array[default_value]
          mode <- names(which.max(table(default_tempo))[1, drop = FALSE])
          if(is.numeric(default_tempo)){ mode <- as.numeric( mode ) }
          defaultMRI[param[iter_MRI]] <-  mode
        }
        
      }
      
    }
    
  }
  if(verbose){cat("\n")}
  
  
  
  # remaning parameters
  
  
  
  
  #### cleaning 
  if(rm.ls.MRI){
    rm(list = name.ls.MRI, envir = globalenv())
  }
  
  #### convertion
  if(format == "data.frame"){
    if(verbose){cat("* Convertion to data.frame \n")}
    dt.MRI <- as.data.frame(dt.MRI)
  }
  if(format == "matrix"){
    if(verbose){cat("* Convertion to matrix \n")}
    dt.MRI <- as.matrix(dt.MRI)
  }
  if(format == "MRIaggr" || format == "data.table"){
    if(verbose){cat("* Set coordinates as keys \n")}
    
    dt.MRI[, rowNumber := .I]
    data.table::setkeyv(dt.MRI, letters[9:(8+ndim.MRI)])
    rowNumber <- list(constMRIaggr = dt.MRI$rowNumber)
    dt.MRI[, rowNumber := NULL]
    
    # convert the index to take into account of the reorganisation of the rows
    if(!is.null(ls.indexMerge)){
      if(verbose){cat("* Reorganize region index \n")}
      
      ls.indexMerge <- lapply(ls.indexMerge, function(Region){
        lapply(Region, function(x){
          which(rowNumber$constMRIaggr %in% x)
        } )
      })
      
    }
    
  }
  
  dim.voxel <- data.frame(dim.voxel, unit = unit, stringsAsFactors = FALSE)
  rownames(dim.voxel) <- NULL
  
  dim.MRI <- setNames(data.frame( matrix(dim.MRI, nrow = 1, ncol = ndim.MRI), stringsAsFactors = FALSE),
                      letters[9:(8+ndim.MRI)])
  
  #### convertion to MRIaggr
  if(is.null(ls.indexMerge)){ ls.indexMerge <- list() }
  
  if(format == "MRIaggr"){ # ls.indexMerge = ls.indexMerge
    dt.MRI <- new(Class = "MRIaggr", 
                  identifier = identifier, 
                  contrast = dt.MRI, 
                  default_value = defaultMRI, 
                  fieldDim = dim.MRI, 
                  voxelDim = dim.voxel,
                  region = list(overall = region.overall,
                                contrast = ls.indexMerge),
                  ls_descStats = list(rowNumber = rowNumber) 
                  
    )
  }
  
  #### export
  if(format == "MRIaggr"){
    return(dt.MRI)
  }else{
    return(list(dt = dt.MRI,
                default_value = default_value, 
                fieldDim = dim.MRI, 
                voxelDim = dim.voxel,
                ls.indexMerge = ls.indexMerge))
  }
}  

#### 3- Convert ####
####>>> array2dt ####

methods::setMethod(f  = "array2dt", 
                   signature  = "array", 
                   definition = function(array, coords = NULL, format = "data.table",
                                         name_newparam = "res", names_coords = letters[9:(8 + ncol(coords))], na.rm = TRUE,
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### preparation 
                     p <- length(dim(array))  
                     n <- length(array)  
                     index_array <- arrayInd(1:n, .dim = dim(array))
                     if(na.rm){
                       index_array <- index_array[is.na(array) == FALSE,]
                     }
                     if(is.null(coords)){
                       coords <- index_array
                     }
                     
                     #### test
                     if( checkArguments ){
                       if(!is.null(coords) && any(dim(coords) != dim(index_array))){
                         stop("array2dt : incorrect dimension for \'coords\' : \n", 
                              "length(array[!is.na(array)]), length(dim(array))  = ", paste(dim(index_array), collapse = " "), "\n", 
                              "dim(coords) = ", paste(dim(coords), collapse = " "), "\n")
                       }
                       
                       validCharacter(format, validLength = 1, validValues = c("matrix", "data.frame", "data.table"), method  = "array2dt")
                       
                       validDimension(value1 = names_coords, validDimension = p, name1 = "names_coords", name2 = NULL, type = "length", method  = "array2dt")
                     }
                     
                     ### integration des donnees 
                     
                     
                     if(format == "data.table"){
                       data <- data.table(coords, array[index_array])
                       names(data) <- c(names_coords, name_newparam)
                       
                     }else if(format == "data.frame"){
                       
                       data <- data.frame(coords, array[index_array])
                       names(data) <- c(names_coords, name_newparam)
                       
                     }else if(format == "matrix"){
                       
                       data <- cbind(coords, array[index_array])
                       colnames(data) <- c(names_coords, name_newparam)
                       
                     }
                     
                     
                     ### export
                     return(data)
                   }
)


####>>> matrix2array ####
methods::setMethod(f  = "mat2array", 
                   signature  = "matrix", 
                   definition = function(coords, contrast, format = "any", default_value = NA, range.coords = NULL,
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### tests
                     if( checkArguments == TRUE ){
                       validClass(value = contrast, validClass = c("vector","matrix", "data.frame", "data.table"), superClasses = TRUE, method = "dt2array")
                       
                       validDimension(contrast, coords, name1 = "contrast", name2 = "coords", type = "NROW", method = "df2array")
                     }
                     
                     #### main
                     if(is.null(colnames(coords))){
                       colnames(coords) <- letters[9:(ncol(coords) + 1)]
                     }
                     
                     #### call and export
                     return(dt2array(as.data.table(coords), contrast = contrast, format = format, default_value = default_value, range.coords = range.coords,
                                     checkArguments = checkArguments))
                     
                   }
)

####>>> df2array ####
methods::setMethod(f  = "df2array", 
                   signature  = "data.frame", 
                   definition = function(coords, contrast, format = "any", default_value = NA, range.coords = NULL,
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### tests
                     if( checkArguments == TRUE ){
                       validClass(value = contrast, validClass = c("vector","matrix", "data.frame", "data.table"), superClasses = TRUE, method = "dt2array")
                       
                       validDimension(contrast, coords, name1 = "contrast", name2 = "coords", type = "NROW", method = "df2array")
                     }
                     
                     #### main
                     if(is.null(names(coords))){
                       names(coords) <- letters[9:(ncol(matrix) + 1)]
                     }
                     
                     #### call and export
                     return(dt2array(as.data.table(coords), contrast = contrast, format = format, default_value = default_value, range.coords = range.coords,
                                     checkArguments = checkArguments))
                     
                   }
)

####>>> dt2array ####
methods::setMethod(f  = "dt2array", 
                   signature  = "data.table", 
                   definition = function(coords, contrast, format = "any", default_value = NA, range.coords = NULL,
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### preparation 
                     if( is.null(contrast) ){
                       
                       names.coords <- key(coords)  
                       param <- setdiff(names.coords, names(coords))
                       contrast <- coords[param]
                       coords <- coords[names.coords]
                       
                     }else{
                       
                       names.coords <- names(coords)
                       if(!is.data.table(contrast)){
                         contrast <- as.data.table(contrast)
                       }
                       param <- names(contrast)
                       
                     }
                     n.coords <- length(names.coords)
                     n.param <- length(param)
                     n.data <- nrow(coords)
                     
                     ### test
                     if( checkArguments == TRUE ){
                       
                       # coords 
                       if(n.coords == 0){
                         stop("dt2array : wrong specification of \'coords\' \n", 
                              "key must be specified for \'coords\' in order to indicate the columns containing the coordinates \n")
                       }
                       
                       if(n.coords > 3){
                         stop("dt2array : wrong specification of \'coords\' \n", 
                              "\'coords\' can have up to 3 columns \n", 
                              "number of columns of \'coords\' : ", n.coords, "\n")
                       }
                       
                       if(any(duplicated(coords))){
                         stop("dt2array : wrong specification of \'coords\' \n", 
                              "there are duplicated values \n")
                       }
                       
                       # contrast
                       if(n.param == 0){
                         stop("dt2array : wrong specification of \'coords\' \n", 
                              "\'coords\' contains no parameter or the name of the parameter is missing in \'contrast\' \n")
                       }
                       
                       validClass(value = contrast, validClass = c("vector","matrix", "data.frame", "data.table"), superClasses = TRUE, method = "dt2array")
                       validDimension(contrast, coords, name1 = "contrast", name2 = "coords", type = "NROW", method = "dt2array")
                       
                       # format
                       validCharacter(value = format, validLength = 1, validValues = c("any", "matrix", "data.frame", "data.table", "list"), 
                                      refuse.NULL = TRUE, method = "dt2array")
                       
                       if(n.coords == 3 && length(unique(coords[[3]])) > 1 && format %in% c("matrix", "data.frame","data.table")){
                         stop("dt2array : wrong specification of \'format\' \n", 
                              "\'coords\' has several levels for the third coordinate \n", 
                              "the result cannot be converted in a 2D format \n",
                              "requested format: ",format,"\n")
                       }
                       
                       # range.coords
                       if(!is.null(range.coords) && n.coords != length(range.coords)){
                         stop("dt2array : wrong specification of \'range.coords\' \n", 
                              "if not NULL, \'range.coords\' must be a vector of length ", n.coords, " \n", 
                              "is(range.coords) : ", paste(is(range.coords), collapse = " "), " \n", 
                              "length(range.coords) : ", length(range.coords), "\n")
                       }  
                       
                     }
                     
                     #### definition des coordonnees des points dans le repere de la matrice
                     scale <- rep(0, n.coords)
                     
                     coords0 <- copy(coords)
                     
                     if(is.null(range.coords)){ 
                       
                       scale <- as.numeric(coords0[,lapply(.SD,min)]-1)
                       
                       coords0[, names.coords := lapply(1:n.coords, function(x){
                         .SD[[x]]-scale[x]
                       }), with = FALSE]
                       
                     }
                     
                     if(any(coords0 %% 1 != 0)){  
                       coords0[, names.coords := lapply(.SD, function(x){as.numeric(as.factor(x))}), with = FALSE]
                     }
                     
                     if(is.null(range.coords)){
                       range.coords <- as.numeric(coords0[,lapply(.SD,max)])
                     }else{
                       range.coords_tempo <- as.numeric(coords0[,lapply(.SD,max)])
                       if(any(range.coords < range.coords_tempo)){
                         stop("dt2array : wrong specification of \'range.coords\' \n", 
                              "\'range.coords\' must be at least ", paste(range.coords_tempo, collapse = " "), " \n", 
                              "requested \'range.coords\' : ", paste(range.coords, collapse = " "), "\n")    
                       }
                     }
                     
                     #### index correspondant aux points dans la matrice
                     Mindex <- coords0[[1]]
                     for(iter_dim in 2:n.coords){
                       Mindex <- Mindex + cumprod(range.coords)[iter_dim-1]*(coords0[[iter_dim]] - 1)
                     }
                     
                     ### integration des donnees
                     dataM <- list()
                     
                     for(iter_param in 1:n.param){
                       
                       dataM[[iter_param]] <- array(default_value, dim = range.coords)
                       
                       dataM[[iter_param]][Mindex] <- contrast[[iter_param]]
                       
                       if(format %in% c("matrix","data.frame","data.table")){
                         if(n.coords == 3){
                           dataM[[iter_param]] <- dataM[[iter_param]][,,1]
                         }
                         dataM[[iter_param]] <- do.call(paste0("as.",format),
                                                        list(dataM[[iter_param]])
                         )
                       }
                       
                     }
                     names(dataM) <- names(contrast)
                     
                     ### mise en forme
                     if(format %in% c("matrix", "data.frame", "data.table")  && n.param == 1){
                       dataM <- dataM[[1]]
                     }
                     
                     unique_coords <- lapply(1:ncol(coords0), function(x){scale[x]+seq(min(coords0[[x]]), max(coords0[[x]]), by = 1)})
                     names(unique_coords) <- names(coords)
                     
                     #### export
                     return(list(contrast = dataM, 
                                 coords = coords, 
                                 unique_coords = unique_coords))
                     
                     
                     
                   }
)

####>>> dt2array ####
methods::setMethod(f  = "dt2MRIaggr", 
                   signature  = "data.table", 
                   definition = function(coords, contrast = NULL, region.overall = "Stroke",
                                         identifier, default_value = NULL,
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     
                     #### preparation 
                     if( !is.null(contrast) ){
                       
                       names.coords <- names(coords)
                       if(!is.data.table(contrast)){
                         contrast <- as.data.table(contrast)
                       }
                       param <- names(contrast)
                       
                     }else{
                       names.coords <- key(coords)
                       param <- setdiff(names(coords),key(coords))
                     }
                     n.param <- length(param)
                     n.coords <- length(names.coords)
                     
                     ### test
                     if( checkArguments == TRUE ){
                       
                       # contrast
                       if(n.param == 0){
                         stop("dt2MRIaggr : wrong specification of \'coords\' \n", 
                              "\'coords\' contains no parameter or the name of the parameter is missing in \'contrast\' \n")
                       }
                       
                       # contrast
                       if(length(names.coords) == 0){
                         stop("dt2MRIaggr : wrong specification of \'coords\' \n", 
                              "key must be set in \'coords\' in order to identify the columns corresponding to the coordinates \n")
                       }
                       
                       if(!is.null(contrast)){
                         validClass(value = contrast, validClass = c("vector","matrix", "data.frame", "data.table"), superClasses = TRUE, method = "dt2MRIaggr")
                         validDimension(contrast, coords, name1 = "contrast", name2 = "coords", type = "NROW", method = "dt2MRIaggr")
                       }
                     }
                     
                     #### main
                     if( !is.null(contrast) ){
                       coords <- cbind(coords,contrast)
                       setkeyv(coords, names.coords)
                     }
                     
                     #### convertion to MRIaggr
                     dt.MRI <- new(Class = "MRIaggr", 
                                   identifier = identifier, 
                                   contrast = coords, 
                                   default_value = default_value, 
                                   region = list(overall = region.overall)
                     )
                     
                     #### export 
                     return(dt.MRI)
                     
                   }
)


#### 3- Write #####

####>>> writeMRI ####

writeMRI <- function (data, filename, format, gzipped = TRUE, verbose = optionsMRIaggr("verbose"), size = "NA_integer_"){
  
  validCharacter(value = format, validLength = 1, validValues = c("rawb.gz", "analyze", "nifti", "dicom"), 
                 refuse.NULL = TRUE, method = "writeMRI")
  
  objClass <- class(data)
  if (objClass != "array" || length(dim(data)) != 3){
    stop("writeMRI : incorrect \'data\' : \n", 
         "data has to be a 3 dimensional array  \n", 
         "proposed data (dimension) : ", paste(is(data), collapse = " "), " (", paste(dim(data), collapse = " "), ") \n")
  }
  
  if (format == "rawb.gz") {
    
    f <- gzfile(filename, open = "wb")
    writeBin(con = f, object = as.vector(data), size = size)
    on.exit(close(f))
    
  } else if (format == "analyze") {
    
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "writeMRI")
    oro.nifti::writeANALYZE(as(data, "anlz"), filename = filename, gzipped = gzipped, verbose = verbose)
    
  } else if (format == "nifti") {
    
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "writeMRI")
    oro.nifti::writeNIfTI(as(data, "nifti"), filename = filename, gzipped = gzipped, verbose = verbose)
    
  } else if (format == "dicom")  {
    
    # initPackage(package = "oro.dicom", argument = "format = \"oro.dicom\"", method = "writeMRI")
    cat("writeMRI for dicom files is not implemented \n")
    invisible(return(FALSE))    
    
  }
  
}

#### 4- initialisation functions ####

initMergeParam <- function(ls.MergeParam, param, refValue,
                           checkArguments, init, method){
  
  if(!is.null(ls.MergeParam)){
    name.merge <- names(ls.MergeParam)
    name2merge <- unlist(ls.MergeParam)
    
  }else{
    name.merge <- NULL
    name2merge <- NULL
    
  }
  
  #### tests 
  if(checkArguments && !is.null(ls.MergeParam)){
    
    if(!is.list(ls.MergeParam)){
      stop(method,": wrong specification of \'ls.MergeParam\' \n", 
           "\'ls.MergeParam\' must be a list \n", 
           "is(ls.MergeParam): ",paste(ls.MergeParam, collapse = " "),"\n")
    }
    
    validNames(value = ls.MergeParam, name = "ls.MergeParam", validLength = length(ls.MergeParam), method = method)
    
    if(any(name.merge %in% param)){
      stop(method,": wrong specification of \'ls.MergeParam\' \n", 
           "\'ls.MergeParam\' must not contain the same names as \'param\' \n", 
           "identical names: ",paste(intersect(name.merge,param), collapse = " "),"\n")
    }
    
    lapply(ls.MergeParam, function(x){
      validCharacter(value = x, validValues = param, name = "ls.MergeParam", validLength = NULL, refuse.NULL = TRUE, method = "constMRIaggr")
    })
    
    if(any(duplicated(name2merge))){
      stop(method,": wrong specification of \'ls.MergeParam\' \n", 
           "\'ls.MergeParam\' must not contain duplicated names \n", 
           "duplicated names: ",paste( unique(name2merge[duplicated(name2merge) == TRUE]), collapse = " "),"\n")
    }
    
    if(length(refValue) != 1 || (!is.character(refValue) && !is.numeric(refValue))){
      stop(method,": wrong specification of \'refValue\' \n", 
           "refValue must have length 1 and be a numeric or a character \n", 
           "lenght(refValue): ",length(refValue)," \n",
           "is(refValue): ",paste(is(refValue), collapse = " ")," \n",)
    }
    
  }
  
  #### init
  if(init){
    
    if(!is.null(ls.MergeParam)){
      
      n.merge <- length(name.merge)
      ls.indexMerge <- eval(parse(text = paste0("list(", paste(name.merge," = list()", collapse = ", \n "), ")")))
      
    }else{
      
      n.merge <- 0
      ls.indexMerge <- NULL
    }
    
    
    #### export ####
    return(list(name.merge = name.merge,
                name2merge = name2merge, 
                n.merge = n.merge,
                ls.indexMerge = ls.indexMerge))
  }
  
  
}
