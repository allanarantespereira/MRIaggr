#**********************************************************************
#**********************************************************************
#*************         0A Class MRIaggr             *******************
#**********************************************************************
#**********************************************************************
#

##### A) Object definition #############################################################

####>>> setClass ####
methods::setClass(
  
  Class = "MRIaggr", 
  
  representation(
    clinic = "data.frame", 
    contrast = "data.table",    # data.frame integrant toutes les cartos d interet de dim (L * nb, l)
    default_value = "data.frame", 
    fieldDim = "data.frame", 
    history = "list", 
    identifier = "character",   # id patient de la carto
    ls_descStats = "list",
    normalization = "list", 
    region = "list",
    table = "list",
    voxelDim = "data.frame", 
    W = "list"
    
  ), 
  
  ####>>> validity ####
  
  validity = function(object){
    
    if (optionsMRIaggr("checkArguments")) {
      
      #### @clinic
      
      #### @contrast
      contrast.keys <- key(object@contrast)
      n.keys <- length(contrast.keys)
      param <- setdiff(names(object@contrast), contrast.keys)
      
      if(!identical(contrast.keys, letters[9:10]) &&!identical(contrast.keys, letters[9:11]) && !identical(contrast.keys, letters[9:12])){
        stop("validity[MRIaggr] : wrong keys for \'@contrast\' \n", 
             "Keys of \'@contrast\' must be exactly: \"",paste(letters[9:10], collapse = "\" \""),"\" \n",
             "                                   or: \"",paste(letters[9:11], collapse = "\" \""),"\" \n", 
             "                                   or: \"",paste(letters[9:12], collapse = "\" \""),"\" \n", 
             "Proposed keys: \"",paste(contrast.keys, collapse = "\" \""),"\" \n")
      }
      
      if(identical(contrast.keys, letters[9:12])){ # Should be removed in further versions
        stop("validity[MRIaggr] : 4 dimensional datasets are not yet supported by MRIaggr \n",
             "if the size of the 4th dimension is small, you can still use several parameters to integrate it into the MRIaggr object \n")
      }
      
      if(any(duplicated(param))){
        stop("validity[MRIaggr] : \'@contrast\' contains duplicated parameters \n",
             "duplicated parameters: ",paste(unique(param[duplicated(param)]), collapse = " ")," \n",
             "columns: ",paste( which(names(object@contrast) %in% unique(param[duplicated(param)]) ), collapse = " ")," \n")
      }
      
      if("contrast" %in% names(object@ls_descStats)){
        test.index <- TRUE
        names.ContrastMerged <- names(object@ls_descStats$contrast)
        names.ContrastMerged_regions <- unlist(lapply(object@ls_descStats$contrast,names))
        
        if(any(duplicated(names.ContrastMerged_regions))){
          duplicated_tempo <- unique(names.ContrastMerged_regions[duplicated(names.ContrastMerged_regions)])
          stop("validity[MRIaggr] : \'@contrast\' contains duplicated values for merged regions \n",
               "duplicated values: ",paste(duplicated_tempo, collapse = " ")," \n",
               "corresponding columns names: ",paste( names.ContrastMerged[unlist(lapply(object@ls_descStats$contrast, function(x){any(names(x) %in% duplicated_tempo)}))], collapse = " ")," \n")
        }
        
      }else{
        test.index <- TRUE
      }
      
      if("rowNumber" %in% names(object@contrast)){
        stop("validity[MRIaggr] : \'@contrast\' must not contain a parameter named \"rowNumber\" \n")
      }
      
      ### @default_value
      param_default <- names(object@default_value)
      
      if(any(param %in% param_default == FALSE) || any(param_default%in% param == FALSE)){
        stop("validity[MRIaggr] : names of \'@default_value\' do not match those of @contrast \n", 
             "missing parameter in \'@default_value\'  : ", paste(param[param %in% param_default == FALSE], collapse = " "), "\n", 
             "missing parameter in  \'@contrast\' : ",paste(param_default[param_default %in% param == FALSE], collapse = " "), "\n")
      }
      
      #### @fieldDim
      validNames(value = object@fieldDim, name = "@fieldDim", validLength = n.keys, validValues = contrast.keys, method = "validity[MRIaggr]")
      
      #### @history
      
      #### @identifier
      
      #### @ls_descStats
      
      #### @normalization
      # !!!!
      
      #### @region
      validNames(value = object@region, name = "@region", validLength = NULL, 
                 requiredValues = c("overall","contrast"), method = "validity[MRIaggr]")
      validCharacter(value = object@region$overall, name = "@region$overall", validLength = 1, method = "validity[MRIaggr]")
      validClass(value = object@region$contrast, name = "@region$contrast", validClass = "list", method = "validity[MRIaggr]")
      validNames(value = object@region$contrast, name = "@region$region$contrast", refuse.NULL = FALSE, forbiddenValues = "overall", method = "validity[MRIaggr]")
      
      if(object@region$overall %in% c("Brain","Stroke")){
        
        validNames(value = object@region, name = "@region", validLength = NULL, 
                   requiredValues = c("hemispheres","midplane"), method = "validity[MRIaggr]")
        
        ## @hemispheres
        hemi_tempo <- object@region$hemispheres
        validClass(value = hemi_tempo, name = "@region$hemispheres", validClass = "data.frame", method = "validity[MRIaggr]")
        validNames(value = hemi_tempo, name = "@region$hemispheres", validLength = 2, 
                   validValues = c("left","right"), method = "validity[MRIaggr]")
        
        if(length(object@region$contrast)>0){
          names.Merged <- unlist(names(object@region$contrast))
          names.Merged_regions <- unlist(lapply(object@region$contrast, names))
        }else{
          names.ContrastMerged <- NULL
          names.ContrastMerged_regions <- NULL
        }
        
        if(length(hemi_tempo$left)>1 || !is.na(hemi_tempo$left)){
          
          if(!is.list(hemi_tempo) || length(hemi_tempo$left) != 1){
            stop("validity[MRIaggr] : wrong specfication of \'@region$hemispheres$left\' \n", 
                 "\'@region$hemispheres$left\' must contain one list of characters or be NA \n", 
                 "length(@region$hemispheres$left): ", length(hemi_tempo$left), "\n",
                 "is(@region$hemispheres$left): ", is(is(hemi_tempo$left), collapse = " "), "\n")
          }
          
          if( any(hemi_tempo$left[[1]] %in% c(names.Merged_regions,names.Merged) == FALSE) ){
            stop("validity[MRIaggr] : wrong specfication of \'@region$hemispheres$left\' \n", 
                 "valid values for \'@region$hemispheres$left[[1]]\': ",paste(c(names.Merged_regions,names.Merged), collapse = " ")," \n", 
                 "proposed invalid values: ",paste(hemi_tempo$left[[1]][object@hemispheres$left[[1]] %in% c(names.Merged_regions,names.Merged) == FALSE], collapse = " "), "\n")
          }
          
        }
        
        if(length(hemi_tempo$right)>1 || !is.na(hemi_tempo$right)){
          
          if(!is.list(hemi_tempo) || length(hemi_tempo$right) != 1){
            stop("validity[MRIaggr] : wrong specfication of \'@region$hemispheres$right\' \n", 
                 "\'@region$hemispheres$right\' must contain one list of characters or be NA \n", 
                 "length(@region$hemispheres$right): ", length(hemi_tempo$right), "\n",
                 "is(@region$hemispheres$right): ", is(is(hemi_tempo$right), collapse = " "), "\n")
          }
          
          if( any(hemi_tempo$right[[1]] %in% c(names.Merged_regions,names.Merged) == FALSE) ){
            stop("validity[MRIaggr] : wrong specfication of \'@region$hemispheres$right\' \n", 
                 "valid values for \'@region$hemispheres$right[[1]]\': ",paste(c(names.Merged_regions,names.Merged), collapse = " ")," \n", 
                 "proposed invalid values: ",paste(hemi_tempo$right[[1]][object@hemispheres$right[[1]] %in% c(names.Merged_regions,names.Merged) == FALSE], collapse = " "), "\n")
          }
          
        }
        
        #### @midplane
        midplane_tempo <- object@region$midplane
        validClass(value = midplane_tempo, name = "@region$midplane", validClass = "data.frame", method = "validity[MRIaggr]")
        validNames(value = midplane_tempo, name = "@region$midplane", validLength = n.keys, 
                   validValues = contrast.keys, method = "validity[MRIaggr]")
        
      }
      
      #### @table
      if(object@region$overall == "Stroke"){

        validNames(value = object@table, name = "@table", validLength = NULL, 
                   requiredValues = c("lesion","hypoperfusion","reperfusion"), method = "validity[MRIaggr]")
        
        validClass(object@table$lesion, name = "@table$lesion", validClass = "data.frame", method = "validity[MRIaggr]")
        validClass(object@table$hypoperfusion, name = "@table$hypoperfusion", validClass = "data.frame", method = "validity[MRIaggr]")
        validClass(object@table$reperfusion, name = "@table$reperfusion", validClass = "data.frame", method = "validity[MRIaggr]")
        
      }
      
      #### @voxelDim
      validNames(value = object@voxelDim, name = "@voxelDim", validLength = n.keys + 1, validValues = c(contrast.keys,"unit"), method = "validity[Carto3D]")
      
      #### @W
      validNames(value = object@W, name = "@W", validValues = c("Wmatrix", "Wblocks", "upper"), method = "validity[MRIaggr]")
      validClass(value = object@W$Wmatrix, name = "@W$Wmatrix", validClass = c("matrix", "dgCMatrix"), superClasses = TRUE, method = "validity[MRIaggr]")
      validLogical(value = object@W$upper, name = "@W$upper", validLength = 1, refuse.NULL = FALSE, refuse.NA = FALSE, method = "validity[MRIaggr]")
      
      if(is.null(object@W$upper) == FALSE && is.na(object@W$upper) == FALSE){
        validDimension(value1 = object@W$Wmatrix, validDimension = rep( nrow(object@contrast), 2), name1 = "@W$Wmatrix", name2 = NULL, type = c("nrow", "ncol"), method = "validity[MRIaggr]")
      }
      
      
      
      
    }
    
    return(TRUE)
    
  }   
)

####>>> initialize ####

methods::setMethod(
  f = "initialize", 
  signature = "MRIaggr", 
  definition = function(.Object, region.overall = NA, clinic, contrast, default_value, fieldDim,
                        history, identifier, ls_descStats, normalization, 
                        region, table, voxelDim, W){
    
    #### @clinic
    if(!missing(clinic)){
      
      if(!is.null(fieldDim)){
        .Object@clinic <- clinic
      }else{
        .Object@clinic <- data.frame()
      }
      
    }
    
    #### @contrast
    if(missing(contrast)){
      stop("initialize[MRIaggr] : \'@contrast\' must be specified \n")
    }
    .Object@contrast <- data.table::copy(contrast)
    contrast.keys <- key(.Object@contrast)
    n.keys <- length(contrast.keys)
    
    #### @default_value
    param <- setdiff(names(.Object@contrast), contrast.keys)
    
    if(!missing(default_value)){
      
      if(!is.null(default_value)){
        .Object@default_value <- default_value
      }else if(!is.null(contrast.keys)){
        .Object@default_value <- setNames(data.frame(  matrix(NA, nrow = 1L, ncol = length(param)), stringsAsFactors = FALSE), 
                                          param)
      }
      
    }
    
    #### @fieldDim
   if(!missing(fieldDim) && !is.null(fieldDim)){
        .Object@fieldDim <- fieldDim
    }else if(!is.null(contrast.keys)){
      .Object@fieldDim <- as.data.frame(.Object@contrast[, lapply(.SD,max),.SDcols = contrast.keys])
    }
    
    #### @history
    if(!missing(history)){
      .Object@history <- history
    }else{
      .Object@history <- list()
    }
    
    #### @identifier
    if(missing(identifier)){
      stop("initialize[MRIaggr] : \'@identifier\' must be specified \n")
    }
    
    .Object@identifier <- identifier
    
    #### @ls_descStats
    if(!missing(ls_descStats)){
      .Object@ls_descStats <- ls_descStats
    }else{
      .Object@ls_descStats <- list(contrast = NULL)
    }
    
    #### @normalization
    if(!missing(normalization)){
      .Object@normalization <- normalization
    }
    
    #### region
    if(!missing(region)){
      .Object@region <- region
    }else{
      .Object@region <- list()
    }
    
    if("overall" %in% names(.Object@region) ==  FALSE){
      Object@region <- c(.Object@region,
                         list(overall = region.overall)
      )
    }
    
    if("contrast" %in% names(.Object@region) ==  FALSE){
      .Object@region <- c(.Object@region,
                          list(contrast = list())
      )
    }
    
    if(.Object@region$overall %in% c("Brain","Stroke")){ ##> BRAIN MRIaggr
      
      if("hemispheres" %in% names(.Object@region) ==  FALSE){
        .Object@region <- c(.Object@region,
                            list(hemispheres = data.frame(left = NA, right = NA, stringsAsFactors = FALSE))
        )
      }
      
      if("midplane" %in% names(.Object@region) ==  FALSE){
        .Object@region <- c(.Object@region,
                            list(midplane = setNames(data.frame( matrix(NA, nrow = 1, ncol = n.keys), stringsAsFactors = FALSE), 
                                                     contrast.keys))
        )
      }
    }
    
    #### @table
    if(!missing(table)){
      .Object@table <- table
    }
    
    if(.Object@region$overall == "Stroke"){ ##> STROKE MRIaggr

      if("lesion" %in% names(.Object@table) ==  FALSE){
        .Object@table <- c(.Object@table,
                           list(lesion = data.frame())
        )
      }
      
      if("hypoperfusion" %in% names(.Object@table) ==  FALSE){
        .Object@table <- c(.Object@table,
                           list(hypoperfusion = data.frame())
        )
      }
      
      if("reperfusion" %in% names(.Object@table) ==  FALSE){
        .Object@table <- c(.Object@table,
                           list(reperfusion = data.frame())
        )
      }
    }
    
    #### @voxelDim
    if(!missing(voxelDim) && !is.null(voxelDim)){
      
      .Object@voxelDim <- voxelDim
      
    }else if(!is.null(contrast.keys)){
      .Object@voxelDim <- setNames(data.frame( matrix(NA, nrow = 1L, ncol = n.keys + 1L), stringsAsFactors = FALSE), 
                                   c(contrast.keys, "unit"))
    }
    
    #### @W
    if(!missing(W)){
      .Object@W <- W
    }else{
      .Object@W <- list(Wmatrix = matrix(), Wblocks = NULL, upper = NA)
    }
    
    #### export
    validObject(.Object)
    return(.Object)
  }
)




