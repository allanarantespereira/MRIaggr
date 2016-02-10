#**********************************************************************
#**********************************************************************
#*************         0B Class MRIaggr             *******************
#**********************************************************************
#**********************************************************************

##### B) Selecters ############################################
# selectClinic:         ok
# selectContrast:       ok (normalisation to check)
# selectCoords:         ok
# selectDefault_value:  ok
# selectDescStats:      ok
# selectFieldDim:       ok
# selectHemispheres:    ok
# selectHistory:        ok
# selectIdentifier:     ok
# selectMidplane:       ok
# selectSubset:         ok
# selectN:              ok
# selectNormalization:  to be updated
# selectParameter:      ok
# selectRegion:         ok
# selectTable:          oK
# selectVoxelDim:       ok
# selectW:              to be checked

####>>> selectClinic ####

methods::setMethod(f  = "selectClinic", 
                   signature  = "MRIaggr", 
                   definition = function(object, param = TRUE, 
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### tests
                     if (checkArguments) {
                       
                       validCharacter(param, validValues = NULL, validLength = NULL, refuse.NULL = TRUE, refuse.duplicates = TRUE, method = "selectClinic")
                       
                       if(!identical(param, TRUE) && any(param %in% names(object@clinic) == FALSE)){
                         
                         if(ncol(object@clinic) == 0){
                           stop("selectClinic[MRIaggr] : object@clinic is empty \n", 
                                "cannot extract a specif parameter \n",
                                "consider using allocClinic to fill the clinic slot of the object \n")
                         }
                         
                         dist_string <- stats::adist(param[param %in% names(object@clinic) == FALSE], names(object@clinic))
                         propositions <- as.vector(unlist(apply(dist_string, 1, function(y){names(object@clinic)[which(rank(y) < 5)]})))
                         
                         stop("selectClinic[MRIaggr] : wrong specification of \'param\' \n", 
                              "parameter(s) not found: ",paste(param[param %in% names(object@clinic) == FALSE], collapse = " "),"\n",
                              "best matching parameters  : ", paste(propositions, collapse = " "), " ... \n", 
                              "or use param = TRUE to select all parameters \n")
                       }
                       
                     }
                     
                     #### export
                     if(identical(param, TRUE)){
                       return(object@clinic)
                     }else{
                       return(object@clinic[,param, drop = FALSE])
                     }
                     
                   }
)


####>>> selectContrast ####

methods::setMethod(f  = "selectContrast", 
                   signature  = "MRIaggr", 
                   definition = function(object, param = TRUE, coords = FALSE, rowNumber = FALSE,
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, subset = NULL, na.rm = FALSE,
                                         format = "data.table", checkArguments = optionsMRIaggr("checkArguments"), ...){
                     
                     #### generic options
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c("norm_mu", "norm_sigma", "hemisphere", "operator.withinR", "operator.betweenR"), 
                                    refuse.NULL = FALSE, method = "multiplot[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     norm_mu <- optionsMRIaggr.eg$norm_mu
                     norm_sigma <- optionsMRIaggr.eg$norm_sigma
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     operator.withinR <- optionsMRIaggr.eg$operator.withinR
                     operator.betweenR <- optionsMRIaggr.eg$operator.betweenR
                     
                     #### initialisation
                     # coords
                     coords <- initCoords(object, coords = coords, checkArguments = checkArguments, init = TRUE,
                                          method = "selectCoords[MRIaggr]")
                     
                     # param
                     param <- initParameter(object = object, param = param, checkArguments = checkArguments, init = TRUE, accept.coords = FALSE, 
                                            method = "selectContrast") 
                     n.param <- length(param)
                     paramAll <- c(if(rowNumber){"rowNumber"}, coords, param)
                     coordsAll <- selectParameter(object, type = "coords")
                     
                     
                     test.hemiNorm_mu <- norm_mu %in% c("contralateral", "contralateral_1slice", "contralateral_3slices") 
                     test.hemiNorm_sigma <- norm_sigma %in% c("contralateral", "contralateral_1slice", "contralateral_3slices")
                     if(test.hemiNorm_mu || test.hemiNorm_sigma){
                       initParameter(object = object, param = "hemishere", checkArguments = checkArguments, init = FALSE, accept.coords = FALSE, 
                                     method = "selectContrast") 
                       param <- c("hemishere", param)
                     }
                     
                     # norm
                     test.default <- norm_mu == "default_value"
                     test.normalization <- (norm_mu %in% c(FALSE, "default_value") == FALSE) || (norm_sigma != FALSE)
                     
                     if(test.default){
                       default_value <- selectDefault_value(object, param = param)  
                     }
                     
                     if(test.normalization){
                       ls.normalization <- selectNormalization(object)
                     }
                     
                     #### tests
                     
                     if(checkArguments){
                       
                       # norm
                       validCharacter(value = norm_mu, validLength = 1, 
                                      validValues = c(FALSE, "global", "contralateral", "global_1slice", "contralateral_1slice", "global_3slices", "contralateral_3slices", "default_value"),
                                      refuse.NULL = FALSE, method = "selectContrast[MRIaggr]")
                       
                       validCharacter(value = norm_sigma, validLength = 1, 
                                      validValues = c(FALSE, "global", "contralateral", "global_1slice", "contralateral_1slice", "global_3slices", "contralateral_3slices"),
                                      refuse.NULL = FALSE, method = "selectContrast[MRIaggr]")
                       
                       if(test.default && length(default_value) == 0 ){
                         stop("selectContrast[MRIaggr] : cannot perform the normalization \n", 
                              "slot @default_value is empty \n")
                       }
                       
                       if(test.default && any(is.na(default_value)) ){
                         stop("selectContrast[MRIaggr] : cannot perform the normalization \n", 
                              "slot @default_value contains NA for parameter(s): ",paste(names(is.na(default_value)), collapse = " ")," \n")
                       }
                       
                       if(test.normalization && length(ls.normalization) == 0){
                         stop("selectContrast[MRIaggr] : cannot perform the normalization \n", 
                              "slot @default_value is empty \n",
                              "use the calcNormalization  function with argument \'update.object\' set to TRUE \n",
                              "to compute and allocate normalization values to the 'object' \n")
                       }
                       
                       # format
                       validCharacter(value = format, validLength = 1, validValues = c("data.table", "data.frame", "matrix", "vector"), method = "selectContrast[MRIaggr]")
                       
                       if( (format == "vector") && (length(paramAll) > 1) ){
                         stop("selectContrast[MRIaggr] : wrong specification of \'format\' \n", 
                              "vector format is not available when several parameters/ coords are requested \n", 
                              "requested parameters : \"", paste(paramAll, collpase = "\" \""), "\" \n")
                       }
                       
                     }
                     
                     #### subset
                     # select
                     subsetContrast <-  selectSubset(object, param = c(coordsAll, param), rowNumber = rowNumber,
                                                     slice_i = slice_i, slice_j = slice_j, slice_k = slice_k,
                                                     hemisphere = hemisphere, subset = subset,
                                                     operator.withinR = operator.withinR, operator.betweenR = operator.betweenR,
                                                     checkArguments = checkArguments)$subsetContrast
                     slice <- unique(subsetContrast$k)
                     n.slice <- length(slice)
                     
                     # remove NA values
                     if(na.rm){
                       subsetContrast <- na.omit(subsetContrast)
                     }
                     
                     # localise hemispheres
                     if(test.hemiNorm_mu || test.hemiNorm_sigma){
                       
                       index_hemiL <- as.vector(subsetContrast[.I(hemisphere == "left")])
                       index_hemiR <- as.vector(subsetContrast[.I(hemisphere == "right")])
                       
                       if(subsetContrast[hemisphere == "undefined", .N]){
                         warning("selectContrast[MRIaggr] : the hemisphere is undefined for some observations : by default they are assigned to the left hemisphere \n")
                         index_hemiL <- c(index_hemiL, as.vector(subsetContrast[.I(hemisphere == "undefined")]) )
                       }
                       
                     }
                     
                     #### normalisation
                     if(norm_mu != FALSE || norm_sigma != FALSE){
                       
                       ## tests 
                       if(checkArguments && any( subsetContrast[, lapply(.SD, is.numeric), .SDcols = param] == FALSE)){
                         stop("selectContrast[MRIaggr] : cannot perform the normalization \n", 
                              "some of the required contrast parameters are not numeric: ",paste(param[subsetContrast[, lapply(.SD, is.numeric), .SDcols = param] == FALSE], collapse = " ")," \n")
                       }
                       
                       
                       ## mean default
                       if( test.default ){
                         
                         subsetContrast[ , param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                          unlist(default_value), FUN = "-"), 
                                         with = FALSE]
                         
                       }
                       # mean normalisation
                       if(norm_mu == "global"){
                         
                         norm_tempo <- selectNormalization(object, type = "global", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "both")
                         
                         subsetContrast[ , param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                          as.matrix(norm_tempo), FUN = "-"), 
                                         with = FALSE]
                         
                       }else if(norm_mu == "contralateral"){
                         
                         if(length(index_hemiL) > 0){
                           norm_tempo <- selectNormalization(object, type = "global", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "right")
                           
                           subsetContrast[index_hemiL, param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                                      as.matrix(norm_tempo), FUN = "-"), 
                                          with = FALSE]
                           
                         }
                         
                         if(length(index_hemiR) > 0){
                           norm_tempo <- selectNormalization(object, type = "global", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "left")
                           
                           subsetContrast[index_hemiR, param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                                      as.matrix(norm_tempo), FUN = "-"), 
                                          with = FALSE]
                           
                         }
                         
                       }else if(norm_mu %in% c("global_1slice","global_3slices") ) {
                         
                         for(iter_slice in slice){
                           
                           index_slice <- as.vector(subsetContrast[.I(k == iter_slice)])
                           
                           norm_tempo <- switch(norm_mu,
                                                "global_1slice" = selectNormalization(object, type = "1slice", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "both", num = iter_slice),
                                                "global_3slices" = selectNormalization(object, type = "3slices", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "both", num = iter_slice)
                           )
                           
                           subsetContrast[index_slice , param := sweep(subsetContrast[index_slice, param, with = FALSE], 2, 
                                                                       as.matrix(norm_tempo), FUN = "-"), 
                                          with = FALSE]
                           
                         }
                         
                       }else if(norm_mu %in% c("contralateral_1slice","contralateral_3slices")){
                         
                         for(iter_slice in slice){
                           
                           index_slice <- as.vector(subsetContrast[.I(k == iter_slice)])
                           index_sliceL <- intersect(index_hemiL, index_slice)
                           index_sliceR <- intersect(index_hemiL, index_slice)
                           
                           if(length(index_sliceL) > 0){
                             
                             norm_tempo <- switch(norm_mu,
                                                  "contralateral_1slice" = selectNormalization(object, type = "1slice", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "right", num = iter_slice),
                                                  "contralateral_3slices" = selectNormalization(object, type = "3slices", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "right", num = iter_slice)
                             )
                             
                             subsetContrast[index_sliceL , param := sweep(subsetContrast[index_sliceL, param, with = FALSE], 2, 
                                                                          as.matrix(norm_tempo), FUN = "-"), 
                                            with = FALSE]
                           }
                           if(length(index_sliceR) > 0){
                             
                             norm_tempo <- switch(norm_mu,
                                                  "contralateral_1slice" = selectNormalization(object, type = "1slice", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "left", num = iter_slice),
                                                  "contralateral_3slices" = selectNormalization(object, type = "3slices", mu = TRUE, sigma = FALSE, param = paramMRI, hemisphere = "left", num = iter_slice)
                             )
                             
                             subsetContrast[index_sliceR , param := sweep(subsetContrast[index_sliceR, param, with = FALSE], 2, 
                                                                          as.matrix(norm_tempo), FUN = "-"), 
                                            with = FALSE]
                           }
                           
                         }
                         
                       }
                       
                       ## var normalisation
                       if(norm_sigma == "global"){
                         
                         norm_tempo <- selectNormalization(object, type = "global", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "both")
                         
                         subsetContrast[ , param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                          as.matrix(norm_tempo), FUN = "/"), 
                                         with = FALSE]
                         
                       }else if(norm_sigma == "contralateral"){
                         
                         if(length(index_hemiL) > 0){
                           norm_tempo <- selectNormalization(object, type = "global", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "right")
                           
                           subsetContrast[index_hemiL, param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                                      as.matrix(norm_tempo), FUN = "/"), 
                                          with = FALSE]
                           
                         }
                         
                         if(length(index_hemiR) > 0){
                           norm_tempo <- selectNormalization(object, type = "global", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "left")
                           
                           subsetContrast[index_hemiR, param := sweep(subsetContrast[, param, with = FALSE], 2, 
                                                                      as.matrix(norm_tempo), FUN = "/"), 
                                          with = FALSE]
                           
                         }
                         
                       }else if(norm_sigma %in% c("global_1slice","global_3slices") ) {
                         
                         for(iter_slice in slice){
                           
                           index_slice <- as.vector(subsetContrast[.I(k == iter_slice)])
                           
                           norm_tempo <- switch(norm_sigma,
                                                "global_1slice" = selectNormalization(object, type = "1slice", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "both", num = iter_slice),
                                                "global_3slices" = selectNormalization(object, type = "3slices", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "both", num = iter_slice)
                           )
                           
                           subsetContrast[index_slice , param := sweep(subsetContrast[index_slice, param, with = FALSE], 2, 
                                                                       as.matrix(norm_tempo), FUN = "/"), 
                                          with = FALSE]
                           
                         }
                         
                       }else if(norm_sigma %in% c("contralateral_1slice","contralateral_3slices") ){
                         
                         for(iter_slice in slice){
                           
                           index_slice <- as.vector(subsetContrast[.I(k == iter_slice)])
                           index_sliceL <- intersect(index_hemiL, index_slice)
                           index_sliceR <- intersect(index_hemiL, index_slice)
                           
                           if(length(index_sliceL) > 0){
                             
                             norm_tempo <- switch(norm_sigma,
                                                  "contralateral_1slice" = selectNormalization(object, type = "1slice", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "right", num = iter_slice),
                                                  "contralateral_3slices" = selectNormalization(object, type = "3slices", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "right", num = iter_slice)
                             )
                             
                             subsetContrast[index_sliceL , param := sweep(subsetContrast[index_sliceL, param, with = FALSE], 2, 
                                                                          as.matrix(norm_tempo), FUN = "/"), 
                                            with = FALSE]
                           }
                           if(length(index_sliceR) > 0){
                             
                             norm_tempo <- switch(norm_sigma,
                                                  "contralateral_1slice" = selectNormalization(object, type = "1slice", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "left", num = iter_slice),
                                                  "contralateral_3slices" = selectNormalization(object, type = "3slices", mu = FALSE, sigma = TRUE, param = paramMRI, hemisphere = "left", num = iter_slice)
                             )
                             
                             subsetContrast[index_sliceR , param := sweep(subsetContrast[index_sliceR, param, with = FALSE], 2, 
                                                                          as.matrix(norm_tempo), FUN = "/"), 
                                            with = FALSE]
                           }
                           
                         }
                         
                       }
                       
                     }
                     
                     #### remove columns that are not in param all (e.g. extra coordinates)
              
                     subsetContrast <- subsetContrast[, paramAll, with = FALSE] # [TO BE OPTIMIZED]
#                                                                index.rmCol <- setdiff(coordsAll, coords)
#                                                                if(length(index.rmCol) > 0){
#                                                                  subsetContrast[, index.rmCol := NULL, with = FALSE]
#                                                                }
                     
                     
                     #### format
                     if(format == "vector"){
                       subsetContrast <- subsetContrast[[1]]
                     }else if(format == "data.frame"){
                       subsetContrast <- as.data.frame(subsetContrast)
                     }else if(format == "matrix"){
                       subsetContrast <- as.matrix(subsetContrast)
                     }
                     
                     #### export 
                     return(subsetContrast)
                     
                   }
                   
)

####>>> selectCoords ####

methods::setMethod(f  = "selectCoords", 
                   signature  = "MRIaggr", 
                   definition = function(object, coords = c("i", "j", "k"), rowNumber = FALSE, spatial_res = c(1, 1, 1), 
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, subset = NULL,
                                         format = "data.table", checkArguments = optionsMRIaggr("checkArguments"), ...){
                     
                     
                     #### generic options
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c("hemisphere", "operator.withinR", "operator.betweenR"), 
                                    refuse.NULL = FALSE, method = "multiplot[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     operator.withinR <- optionsMRIaggr.eg$operator.withinR
                     operator.betweenR <- optionsMRIaggr.eg$operator.betweenR
                     
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     operator.withinR <- optionsMRIaggr.eg$operator.withinR
                     operator.betweenR <- optionsMRIaggr.eg$operator.betweenR
                     
                     coords <- initCoords(object, coords = coords, checkArguments = checkArguments, init = TRUE,
                                          method = "selectCoords[MRIaggr]")
                     
                     #### tests
                     if (checkArguments) {
                       
                       # spatial_res
                       validNumeric(value = spatial_res, validLength = 3, min = 0, method = "selectCoords[MRIaggr]")
                       
                       # format
                       validCharacter(value = format, validLength = 1, validValues = c("data.table", "data.frame", "matrix", "vector"), method = "selectCoords[MRIaggr]")
                       
                       if( (format == "vector") && (length(coords) > 1) ){
                         stop("selectCoords[MRIaggr] : wrong specification of \'coords\' \n", 
                              "vector coords is not available when several coordinates are requested \n", 
                              "requested coordinates : \"", paste(coords, collpase = "\" \""), "\"\n")}
                       
                     }
                     
                     #### extract subset
                     subsetCoords <-  selectSubset(object, param = coords, rowNumber = rowNumber,
                                                   slice_i = slice_i, slice_j = slice_j, slice_k = slice_k,
                                                   hemisphere = hemisphere, subset = subset, 
                                                   operator.withinR = operator.withinR, operator.betweenR = operator.betweenR,
                                                   checkArguments = checkArguments)$subsetContrast
                     
                     #### spatial resolution
                     if(spatial_res[1] != 1 && ("i" %in% coords) ){subsetCoords[,i := i * spatial_res[1]]}
                     if(spatial_res[2] != 1 && ("j" %in% coords) ){subsetCoords[,j := j * spatial_res[2]]}
                     if(spatial_res[3] != 1 && ("k" %in% coords) ){subsetCoords[,k := k * spatial_res[3]]}
                     
                     #### format
                     if(format == "vector"){
                       subsetCoords <- as.vector(subsetCoords)
                     }else if(format == "data.frame"){
                       subsetCoords <- as.data.frame(subsetCoords)
                     }else if(format == "matrix"){
                       subsetCoords <- as.matrix(subsetCoords)
                     }
                     
                     #### export 
                     return(subsetCoords)
                     
                   }
)


####>>> selectDefault_value ####

methods::setMethod(f  = "selectDefault_value", 
                   signature  = "MRIaggr", 
                   definition = function(object, param = TRUE, format = "data.frame", 
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     param <- initParameter(object = object, param = param, checkArguments = checkArguments, init = TRUE, 
                                            accept.coords = FALSE, method = "selectDefault_value")
                     
                     #### test
                     if ( checkArguments ) {
                       validCharacter(value = format, validLength = 1, validValues = c("data.frame", "vector"), method = "selectDefault_value[MRIaggr]")
                     }
                     
                     #### extraction 
                     if(format == "data.frame"){
                       default_value <- object@default_value[,param, drop = FALSE]
                     }else{
                       default_value <- unlist(object@default_value[1,param, drop = TRUE])
                     }
                     
                     #### export
                     return(default_value) }
)

####>>> selectDescStats ####

methods::setMethod(f  = "selectDescStats", 
                   signature  = "MRIaggr", 
                   definition = function(object, name = NULL, subset = NULL, 
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### test
                     if ( checkArguments ) {
                       validCharacter(value = name, validLength = 1, validValues = names(object@ls_descStats), refuse.NULL = FALSE, method = "selectDescStats[MRIaggr]")
                       
                       if(!is.null(subset) && !identical(name, "contrast")){
                         warning("selectDescStats[MRIaggr]: argument \'subset\' is disregarded if argument \'name\' is not set to \'contrast\' \n",
                                 "proposed name: ",name,"\n")
                       }
                       
                     }
                     
                     #### export
                     if(is.null(name)){
                       
                       return(object@ls_descStats) 
                       
                     }else if(is.null(subset)){       
                       
                       return(object@ls_descStats[[name]])  
                       
                     }else {
                       
                       return(initSubset(object, subset = subset, checkArguments = checkArguments, 
                                         arg_name = "subset", method = "selectDescStats"))
                       
                     }
                     
                   }
)

####>>> selectFieldDim ####

methods::setMethod(f  = "selectFieldDim", 
                   signature  = "MRIaggr", 
                   definition = function(object, coords = TRUE,
                                         format = "data.frame", checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### test
                     if(checkArguments){
                       validCharacter(format, validValues = c("data.frame", "vector"),validLength = 1, method = "selectVoxelDim[MRIaggr]")
                     }
                     
                     #### preparation
                     coords <- initCoords(object, coords = coords, checkArguments = checkArguments, init = TRUE,
                                          method = "selectVoxelDim[MRIaggr]")
                     
                     #### export
                     switch(format,
                            "data.frame" = return(object@fieldDim[,coords, drop = FALSE]),
                            "vector" = return(unlist(object@fieldDim[,coords, drop = TRUE]))
                     )
                   }
)

####>>> selectHemispheres ####

methods::setMethod(f  = "selectHemispheres", 
                   signature  = "MRIaggr", 
                   definition = function(object, hemisphere = "both", 
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### tests
                     if(checkArguments){
                       validCharacter(value = hemisphere, validLength = 1, validValues = c("both", "left", "right", "lesion", "contralateral"), method = "selectHemispheres[MRIaggr]")
                     }
                     
                     #### extraction
                     if(length(hemisphere) == 1L && hemisphere %in% c("both","left","right")){
                       
                       # export
                       switch(hemisphere,
                              "both" = return(object@region$hemispheres),
                              "left" = return(object@region$hemispheres$left),
                              "right" = return(object@region$hemispheres$right)
                       )
                       
                     } else {
                       
                       names.Merged <- selectRegion(object, type = "name")
                       
                       if(checkArguments && any(hemisphere %in% names.Merged == FALSE) ){
                         stop("selectHemispheres[MRIaggr] : wrong specfication of \'hemisphere\' \n", 
                              "if not \"both\" \"left\" \"right\" (only single value) \n",
                              "\'hemisphere\' must be in the following: ",paste(names.Merged, collapse = " ")," \n", 
                              "proposed invalid values: ",paste(hemisphere[hemisphere %in% names.Merged == FALSE], collapse = " "), "\n")
                       }
                       
                       test.left <- any(hemisphere %in% object@region$hemispheres$left[[1]])
                       test.right <- any(hemisphere %in% object@region$hemispheres$right[[1]])
                       
                       # export
                       if(test.left == TRUE && test.right == TRUE){
                         return("both")
                       }else if(test.left == TRUE){
                         return("left")
                       }else if(test.rigth == TRUE){
                         return("right")
                       }else{
                         return(NULL)
                       }
                       
                     }
                     
                   }
)

####>>> selectHistory ####

methods::setMethod(f  = "selectHistory", 
                   signature  = "MRIaggr", 
                   definition = function(object){
                     return(object@history)     
                   }
)

####>>> selectIdentifier ####

methods::setMethod(f  = "selectIdentifier", 
                   signature  = "MRIaggr", 
                   definition = function(object){
                     return(object@identifier)     
                   }
)

####>>> selectMidplane ####

methods::setMethod(f  = "selectMidplane", 
                   signature  = "MRIaggr", 
                   definition = function(object){
                     return(object@region$midplane) 
                   }
)

####>>> selectN ####

methods::setMethod(f  = "selectN", 
                   signature  = "MRIaggr", 
                   definition = function(object, slice_i = NULL, slice_j = NULL, slice_k = NULL,
                                         subset = NULL, checkArguments = optionsMRIaggr("checkArguments"), ...){
                     
                     #### generic options
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c("hemisphere", "operator.withinR", "operator.betweenR"), 
                                    refuse.NULL = FALSE, method = "multiplot[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     operator.withinR <- optionsMRIaggr.eg$operator.withinR
                     operator.betweenR <- optionsMRIaggr.eg$operator.betweenR
                     
                     #### define subset
                     subsetN <-  selectSubset(object, param = selectParameter(object, checkArguments =  FALSE)[1], rowNumber = FALSE, 
                                              slice_i = slice_i, slice_j = slice_j, slice_k = slice_k,
                                              hemisphere = hemisphere, subset = subset,
                                              operator.withinR = operator.withinR, operator.betweenR = operator.betweenR,
                                              checkArguments = checkArguments)$subsetContrast[,.N]
                     
                     #### export
                     return(subsetN)
                     
                   }
) 

####>>> selectNormalization ####

methods::setMethod(f  = "selectNormalization", 
                   signature  = "MRIaggr", 
                   definition = function(object, type = NULL, mu = TRUE, sigma = TRUE, hemisphere = "both", slice_k = NULL, param = NULL, 
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### test
                     if( checkArguments ){
                       validCharacter(value = type, validLength = 1, validValues = c("global", "slice", "3slices"), refuse.NULL = FALSE, method = "selectNormalization[MRIaggr]")  
                       
                       if( (mu == FALSE) && (sigma == FALSE) ){
                         stop("selectNormalization[MRIaggr] : \'mu\' and \'sigma\' cannot both be FALSE \n")
                       }
                       
                       validCharacter(value = hemisphere, validLength = 1, validValues = c("both", "left", "right"), method = "selectNormalization[MRIaggr]")
                     }
                     
                     
                     #### select all values
                     if(is.null(type)){
                       return(object@normalization)
                     }
                     
                     #### select specific values
                     
                     ## type global
                     if(type == "global"){
                       
                       res_norm <- object@normalization$norm_global
                       
                       ## test
                       validCharacter(value = param, validLength = NULL, validValues = names(res_norm), refuse.NULL = FALSE, method = "selectNormalization[MRIaggr]")
                       
                       if(is.null(hemisphere)){hemisphere <- c("both", "left", "right")}
                       
                       stat <- NULL
                       if(mu == TRUE){stat <- c(stat, "mu")}
                       if(sigma == TRUE){stat <- c(stat, "sigma")}
                       
                       nom_lignes <- as.vector(sapply(stat, function(object){paste(object, "_", hemisphere, sep = "")}))
                       res_norm <- res_norm[nom_lignes,, drop = FALSE]
                       
                       if(!is.null(param)){
                         res_norm <- res_norm[,param, drop = FALSE]
                       }
                       
                       return(res_norm)
                       
                     }else{
                       
                       if(mu == TRUE){stat <- "Mu"}
                       if(sigma == TRUE){stat <- "Sigma"}
                       nom <- paste("norm", stat, "_", type, "_", hemisphere, sep = "")
                       
                       ## test
                       if( checkArguments ){
                         if( (mu == TRUE) && (sigma == TRUE) ){
                           stop("selectNormalization[MRIaggr] :  arguments \'mu\' and \'sigma\' cannot be simultaneously TRUE \n")
                         }
                         
                         validCharacter(value = nom, name = "normalisation", validLength = 1, validValues = names(object@normalization), method = "selectNormalization[MRIaggr]")
                       }
                       
                       res_norm <- object@normalization[[nom]]
                       validCharacter(value = param, validLength = NULL, validValues = names(res_norm), method = "selectNormalization[MRIaggr]")
                       
                       if(!is.null(param)){
                         res_norm <- res_norm[,param, drop = FALSE]
                       }
                       
                       if(!is.null(slice_k)){					   
                         validInteger(value = slice_k, validLength = NULL, validValues = 1:nrow(res_norm), refuse.duplicates = TRUE,  method = "selectNormalization[MRIaggr]")
                         
                         res_norm <- res_norm[slice_k,, drop = FALSE]
                       }
                       
                       return(res_norm)
                     }
                     
                   }
)

####>>> selectParameter ####

methods::setMethod(f  = "selectParameter", 
                   signature  = "MRIaggr", 
                   definition = function(object, type = "contrast", 
                                         checkArguments = optionsMRIaggr("checkArguments")){                    
                     
                     #### tests
                     if(checkArguments){
                       validCharacter(value = type, validLength = 1, 
                                      validValues = c("clinic", "contrast", "contrastOnly", "region", "coords", "ls_descStats"), 
                                      method = "selectParameter[MRIaggr]")
                     }
                     
                     #### export
                     if(type == "contrast"){
                       
                       return( setdiff(names(object@contrast), key(object@contrast)) )  
                       
                     } else if(type == "contrastOnly"){
                       
                       return( setdiff(names(object@contrast), c(key(object@contrast), selectRegion(object, type = "names")) ) )  
                       
                     } else if(type == "region"){
                       
                       return( selectRegion(object, type = "names") )  
                       
                     } else if (type == "coords"){
                       
                       return(key(object@contrast))
                       
                     } else if (type == "clinic") {
                       
                       return(names(object@clinic))
                       
                     } else if(type == "ls_descStats"){
                       
                       return(names(object@ls_descStats))
                       
                     }
                   }
)

####>>> selectRegion ####

methods::setMethod(f  = "selectRegion", 
                   signature  = "MRIaggr", 
                   definition = function(object, region = NULL, region.value = NULL, type = "value",
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     if(!is.null(region)){
                       names.Contrast <- selectParameter(object, type = "contrastOnly")
                       names.Merged <- names(object@region$contrast) # must not  be written names.Merged <- selectRegion(object, type = "names") otherwise it is infine loop !!
                     }
                     
                     #### specific cases
                     if(identical(region,"overall")){
                       return(object@region$overall)
                     }
                     
                     #### test
                     if(checkArguments){
                       
                       if(!is.null(region)){
                         
                         if(any(region %in% c(names.Merged, names.Contrast) == FALSE) ){
                           stop("selectRegion[MRIaggr] : wrong specfication of \'region\' \n", 
                                "if not \"overall\" (only single value) or NULL \n",
                                "\'region\' must be a region stored in the object: ",paste(names.Merged, collapse = " ")," \n", 
                                "           or a contrast parameter: ",paste(names.Contrast, collapse = " ")," \n", 
                                "proposed invalid values: ",paste(region[region %in% c(names.Merged, names.Contrast) == FALSE], collapse = " "), "\n")
                         } 
                         
                         if(any(region %in% names.Contrast) && is.null(region.value)){
                           stop("selectRegion[MRIaggr] : wrong specfication of \'region\' \n", 
                                "if it contains a contrast parameter (here ",paste(region[region %in% names.Contrast], collapse = " ")," \n",
                                "argument \'region.value\' must be specified to indicate how to dichotomize the contrast parameter\n")
                         } 
                         
                         if(!is.null(region.value)){
                           
                           if(length(region) != 1 ){
                             stop("selectRegion[MRIaggr] : wrong specfication of \'region\' \n", 
                                  "if \'region.value\' argument is specifyied one and only one region/contrast parameter can be selected \n",
                                  "length(region): ",length(region)," \n")
                           }
                           
                           validNames(region, name = "region", validValues = c(names.Merged, names.Contrast), method = "selectRegion[MRIaggr]")
                           
                           if(region %in% names.Merged){
                             
                             validCharacter(region.value, name = paste0("region.value (region: ",region,")"),  validLength = NULL,
                                            validValues = names(object@region$contrast[[region]]), 
                                            method = "selectRegion[MRIaggr]")
                             
                           }else{ # if(names.subset[x] %in% names.Contrast)
                             
                             validCharacter(region.value, name = paste0("region.value (region: ",region,")"), validLength = 2, method = "selectRegion[MRIaggr]")
                             validCharacter(region.value[1], name = paste0("region.value[1] (region: ",region,")"), validLength = NULL, validValues = c("==",">=","<=",">","<"), method = "selectRegion[MRIaggr]")
                             if(!is.numeric(region.value[2]) && is.na(as.numeric(region.value[2])) ){
                               stop("selectRegion[MRIaggr]: wrong specification of \'region.value[2]\' (region: ",region,") \n ",
                                    "if not a numeric, it must be possible to convert it into a numeric using as.numeric \n",
                                    "propsed value: ",region.value[2],"\n")
                             }
                             
                           }
                         }
                         
                         if(!is.null(region.value) && type == "names"){
                           stop("selectRegion[MRIaggr] : wrong specfication of \'type\' \n", 
                                "there is no name to return when argument \'region.value\' is not NULL \n")
                         }
                         
                         validCharacter(type, name = "type", validLength = 1, validValues = c("value","names"), method = "selectRegion[MRIaggr]")
                         
                       }
                       
                     }
                     
                     ## main
                     if(is.null(region)){
                       
                       if(type == "names"){
                         return(names(object@region$contrast))
                       }else{
                         return(object@region$contrast)
                       }
                       
                     }else if(is.null(region.value)){
                       
                       if(type == "names"){
                         return(names(object@region$contrast[[region]]))
                       }else{
                         return(object@region$contrast[region])
                       }
                       
                     }else{
                       
                       if(region %in% names.Merged){
                         return(object@region$contrast[[region]][region.value])
                       }else{ # if(names.subset[x] %in% names.Contrast)
                         switch(region.value[1],
                                "==" = return(object@contrast[, as.numeric(na.omit(.I[.SD  == region.value[2] ])), .SDcols = region]),
                                ">=" = return(object@contrast[, as.numeric(na.omit(.I[.SD  >= region.value[2] ])), .SDcols = region]),
                                "<=" = return(object@contrast[, as.numeric(na.omit(.I[.SD  <= region.value[2] ])), .SDcols = region]),
                                ">" = return(object@contrast[, as.numeric(na.omit(.I[.SD  > region.value[2] ])), .SDcols = region]),
                                "<" = return(object@contrast[, as.numeric(na.omit(.I[.SD  < region.value[2] ])), .SDcols = region])
                         )
                       }
                       
                     }
                     
                   })

####>>> selectSubset ####

methods::setMethod(f  = "selectSubset", 
                   signature  = "MRIaggr", 
                   definition = function(object, rowNumber, param = NULL, 
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, hemisphere = "both", subset = NULL,
                                         operator.withinR = optionsMRIaggr("operator.withinR"), 
                                         operator.betweenR = optionsMRIaggr("operator.betweenR"),
                                         checkArguments = optionsMRIaggr("checkArguments")){
                     
                     test.slice <- (!is.null(slice_i) || !is.null(slice_j) || !is.null(slice_k))
                     
                     #### test
                     if(checkArguments){
                       validLogical(rowNumber, validLength = 1, method = "selectContrast[selectSubset]")    
                       
                       if( (missing(param) || is.null(param)) && rowNumber == FALSE){
                         stop("selectSubset[MRIaggr]: argument \'param\' must be specified or rowNumber set to TRUE \n")
                       }
                     }
                     
                     #### PART 1 - define subset ####
                     if(!is.null(subset)){
                       subset <- initSubset(object, subset = subset, 
                                            operator.withinR = operator.withinR, operator.betweenR = operator.betweenR,
                                            checkArguments = checkArguments, 
                                            arg_name = "subset", method = "selectSubset")
                     }
                     
                     if(hemisphere != "both"){
                       
                       ## search the hemisphere(s) containing the specifyied regions
                       if(hemisphere %in% c("left", "right") == FALSE){
                         hemisphere <- selectHemispheres(object, hemisphere = hemisphere)
                       }
                       
                       ## select the corresponding observations observations 
                       if(hemisphere == "right"){
                         if(is.null(subset)){
                           subset <- MRIaggr.Pa1@contrast[,.I[hemisphere == "right"]]
                         }else{
                           subset <- subset[MRIaggr.Pa1@contrast[subset,.I[hemisphere == "right"]]]  
                         }
                         
                       } else if(hemisphere == "left"){
                         if(is.null(subset)){
                           subset <- MRIaggr.Pa1@contrast[,.I[hemisphere == "left"]]
                         }else{
                           subset <- subset[MRIaggr.Pa1@contrast[subset,.I[hemisphere == "left"]]]  
                         }
                       }
                       
                     }
                     
                     test.subset <- !is.null(subset)
                     
                     #### PART 2 - define slice ####
                     if(test.slice == TRUE){
                      
                      n.coords <- length(key(object@contrast))
                       
                       slice_i <- initSlice_var(object = object, slice = slice_i, slice_var = "i", checkArguments = checkArguments, init = TRUE, method = "selectN[MRIaggr]")
                       slice_j <- initSlice_var(object = object, slice = slice_j, slice_var = "j", checkArguments = checkArguments, init = TRUE, method = "selectN[MRIaggr]")
                       
                       if(n.coords > 2){
                         slice_k <- initSlice_var(object = object, slice = slice_k, slice_var = "k", checkArguments = checkArguments, init = TRUE, method = "selectN[MRIaggr]")
                       }
                       
                       if(n.coords == 2){
                         CJ_slice <- CJ(i = slice_i, j = slice_j, sorted = FALSE)
                       }else{
                         CJ_slice <- CJ(i = slice_i, j = slice_j, k = slice_k, sorted = FALSE)
                       }
                      
                     }
                     
                     #### PART 3 - select subset ####
                    
                     if(test.slice == FALSE && test.subset == FALSE){ ### no subset
                       if(rowNumber == FALSE){
                         subsetContrast <- object@contrast[, .SD, .SDcols = param] 
                       }else if(rowNumber == TRUE){
                         subsetContrast <- object@contrast[, c(rowNumber = list(1:.N), .SD), .SDcols = param] 
                       }
                       
                     } else if(test.slice == TRUE && test.subset == FALSE){ ### slice
                       
                       subsetContrast <- object@contrast
                       
                     } else { #if(test.subset == TRUE){ ### subset
                       
                       if(test.slice == TRUE){
                         paramAll <- unique(c(param,selectParameter(object, type = "coords")))
                       }else{
                         paramAll <- param
                       }
                       
                       if(rowNumber == FALSE){
                         subsetContrast <- object@contrast[subset, .SD, .SDcols = paramAll, nomatch = 0]
                       }else{
                         subsetContrast <- object@contrast[, c(rowNumber = list(1:.N), .SD), .SDcols = paramAll][subset, nomatch = 0]
                       }
                       
                       if(test.slice == TRUE){
                         subsetContrast <- copy(subsetContrast)
                         setkeyv(subsetContrast, key(object@contrast))
                       }
                       
                     }
                     
                     ## slice
                     if(test.slice == TRUE){
                       
                       if(rowNumber == FALSE){
                         
                         subsetContrast <- subsetContrast[CJ_slice, .SD, .SDcols = param, nomatch = 0]
                         
                       }else{
                         
                         paramAll <- unique(c(param,selectParameter(object, type = "coords")))
                         subsetContrast <- subsetContrast[, c(rowNumber = list(1:.N), .SD), .SDcols = paramAll][CJ_slice, nomatch = 0]
                       }
                       
                     }
                     
                     #### export
                     return(list(subsetContrast = subsetContrast,
                                 slice_i = slice_i,
                                 slice_j = slice_j,
                                 slice_k = slice_k,
                                 subset = subset,
                                 test.slice = test.slice,
                                 test.subset = !is.null(subset)
                     ))
                   }
)


####>>> selectTable ####

methods::setMethod(f  = "selectTable", 
                   signature  = "MRIaggr", 
                   definition = function(object, type, size = FALSE, 
                                         checkArguments = optionsMRIaggr("checkArguments"))
                   {            
                     #### test
                     if(checkArguments){
                       validCharacter(value = type, validLength = 1, validValues = c("lesion", "reperfusion", "hypoperfusion"), method = "selectTable[MRIaggr]")
                     }
                     
                     #### extract
                     res <- switch(type, 
                                   "lesion" = object@table$lesion, 
                                   "reperfusion" = object@table$reperfusion, 
                                   "hypoperfusion" = object@table$hypoperfusion)
                     
                     #### conversion to volume
                     if(size == TRUE){
                       Volume <- prod(selectVoxelDim(object)[1:3])
                       if(type == "lesion"){
                         res <- res * Volume
                       }else{                
                         test.V <- sapply(names(res), function(x){substr(x, 1, 1)}) == "V"
                         res[,test.V] <- res[,test.V] * Volume                
                       }
                     }
                     
                     #### export
                     return(res) 
                   }
)

####>>> selectVoxelDim ####

methods::setMethod(f  = "selectVoxelDim", 
                   signature  = "MRIaggr", 
                   definition = function(object, coords = TRUE, unit = TRUE,
                                         format = "data.frame", checkArguments = optionsMRIaggr("checkArguments")){
                     
                     #### test
                     if(checkArguments){
                       validLogical(unit, validLength = 1, method = "selectVoxelDim[MRIaggr]")
                       
                       validCharacter(format, validValues = c("data.frame", "vector"),validLength = 1, method = "selectVoxelDim[MRIaggr]")
                     }
                     
                     #### preparation
                     coords <- initCoords(object, coords = coords, checkArguments = checkArguments, init = TRUE,
                                          method = "selectVoxelDim[MRIaggr]")
                     
                     if(unit == TRUE){
                       coords <- c(coords, "unit")
                     }
                     
                     #### export
                     switch(format,
                            "data.frame" = return(object@voxelDim[,coords, drop = FALSE]),
                            "vector" = return(unlist(object@voxelDim[,coords, drop = TRUE]))
                     )
                   }
)

####>>> selectW ####

methods::setMethod(f  = "selectW", 
                   signature  = "MRIaggr", 
                   definition = function(object, type = "Wmatrix",  upper = NULL, 
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, hemisphere = "both", subset = NULL,
                                         checkArguments = optionsMRIaggr("checkArguments")){ 
                     
                     #### tests
                     if(checkArguments){
                       # initPackage("spam", method = "selectW[MRIaggr]")
                       # initPackage("Matrix", method = "selectW[MRIaggr]")
                       validCharacter(type, validLength = 1, validValues = c("Wmatrix", "Wblocks", "upper"), method = "selectW")
                       
                       if(type == "Wmatrix"){
                         
                         if(!is.null(subset_W)){
                           validInteger(value = subset_W, validLength = NULL, min = 1, max = selectN(object), 
                                        refuse.duplicates = TRUE, method = "selectW[MRIaggr]")
                         }
                         
                         if(is.na(object@W$upper) ){
                           stop("selectW[MRIaggr] : no neighborhood matrix has been stored in the object \n", 
                                "use calcW to compute and store the neighborhood matrix in the object \n")
                         }
                         
                         if(is.logical(upper) && is.null(object@W$upper)){ # [TO BE FIXED]
                           stop("selectW[MRIaggr] : wrong specification of \'upper\' \n", 
                                "the neighborhood matrix was stored with all values \n", 
                                "cannot extract only upper or lower values, set \'upper\' to NULL \n")
                         }
                         
                         if(is.logical(upper) && is.logical(object@W$upper) && upper != object@W$upper){
                           stop("selectW[MRIaggr] : wrong specification of \'upper\' \n", 
                                "the neighborhood matrix has been stored with upper = ", object@W$upper, " \n", 
                                "cannot extract for upper = ", upper, " \n")
                         }
                       }
                     }
                     
                     #### extract
                     if(type == "Wblocks"){
                       
                       return(object@W$Wblocks)
                       
                     } else if(type == "upper"){
                       
                       return(object@W$upper)
                       
                     } else if (type == "Wmatrix"){
                       
                       test.subset <- !is.null(slice_i) || !is.null(slice_j) || !is.null(slice_k) || hemisphere != "both" || !is.null(subset)
                       
                       if(test.subset == TRUE){
                         subset_W <- as.vector(selectSubset(object, param = FALSE, rowNumber = TRUE,
                                                            slice_i = slice_i, slice_j = slice_j, slice_k = slice_k, hemisphere = hemisphere, subset = subset,
                                                            checkArguments = checkArguments)$subsetContrast)
                         subsetW <- object@W$Wmatrix[subset_W, subset_W]
                       }else{
                         subsetW <- object@W$Wmatrix
                       }
                       
                       ## upper or lower matrix to full matrix
                       if(is.null(upper) && is.logical(object@W$upper)){
                         subsetW <- subsetW + spam::t(subsetW)
                       }
                       
                       ## full matrix to upper matrix
                       if(identical(upper, TRUE) && is.null(object@W$upper)){
                         subsetW <- spam:::upper.tri(subsetW)
                       }
                       
                       ## full matrix to lower matrix
                       if(identical(upper, FALSE) && is.null(object@W$upper)){
                         subsetW <- spam:::lower.tri(subsetW)
                       }
                       
                       return(subsetW)
                     }
                     
                     
                     
                   }
)