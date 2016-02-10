#**********************************************************************
#**********************************************************************
#*************         0E Class MRIaggr            *******************
#**********************************************************************
#**********************************************************************
#

##### E) Init ############################################
# initCoords:       ok but not checked
# initIndex:        
# initParameter:   ok but not checked
# initSlice_var:   ok but not checked
# initSubset:      ok
#

####>>> initCoords ####
methods::setMethod(f  = "initCoords", 
                   signature  = "MRIaggr", 
                   definition = function(object, coords, checkArguments, init,
                                         arg_name = "coords", long_name = "coordinates", method)
                   { 
                     #### intialization
                     if(init == TRUE){
                       
                       if(identical(coords, TRUE)){
                         
                         coords <- selectParameter(object, type = "coords")
                         
                       }else if(is.numeric(coords)){
                         
                         if(checkArguments == TRUE){
                           validInteger(coords, validValues = 1:ncol(object@fieldDim), refuse.duplicates = TRUE, method = method)
                         }
                         
                         coords <- selectParameter(object, type = "coords")[coords]
                         
                       }else if(identical(coords, FALSE)){
                         
                         coords <- NULL
                         
                       }
                       
                     }
                     
                     #### test
                     if(checkArguments == TRUE){
                       
                       coords_data <- selectParameter(object, type = "coords")
                       
                       if(any(coords %in% coords_data == FALSE) || any(duplicated(coords)) ){
                         stop(method, "[MRIaggr] : wrong specification of \'", arg_name, "\' \n", 
                              "unknown ", long_name, " : ", paste(param[coords %in% coords_data == FALSE], collapse = " "), " \n", 
                              "duplicated ", long_name, " : ", paste(unique(coords[duplicated(coords)]), collapse = " "), "\n", 
                              "available ", long_name, " : ", paste(coords_data, collapse = " "), "\n")
                       }
                     }
                     
                     #### export
                     return(invisible(coords))
                   }
)

####>>> initIndex ####
methods::setMethod(f  = "initIndex", 
                   signature  = "MRIaggr", 
                   definition = function(object, index,
                                         slice_i, slice_j, slice_k, hemisphere,
                                         arg_name = "index", method)
                   { 
                     #### initialization ####
                     if(is.character(index)){
                       index <- list(subset = index)
                    }
                     
                     if( (is.list(index) && "subset" %in% names(index) ) ){
                      
                       index$coords <- selectCoords(object, 
                                                    slice_i = slice_i, slice_j = slice_j, slice_k = slice_k,
                                                    hemisphere = hemisphere, subset = index$subset,
                                                    format = "data.table")
                       index$subset <- NULL
                     
                      }else if(!is.list(index) || "coords" %in% names(index) == FALSE ){
                       stop(method, "[MRIaggr] : wrong specification of \'", arg_name, "\' ",
                            arg_name, " must be a character refering to regions stored in the object \n",
                            "or a list containing a subset element to be passed to selectCoords \n",
                            "or a list containing a coords element: a data frame containing the coordinates of the points to be displayed \n",
                            "names(",arg_name,"): ",paste(names(index), collapse =  " "),"\n")
                     }
                     
                     #### export ####
                     return(index)
                     
                   }
)


####>>> initParameter ####
methods::setMethod(f  = "initParameter", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, checkArguments, init, accept.coords,
                                         arg_name = "param", long_name = "parameters", method)
                   { 
                     #### intialization
                     if(init == TRUE){
                       
                       if(identical(param, TRUE)){
                         
                         param <- selectParameter(object, type = "contrast")
                         
                       }else if(is.numeric(param)){
                         
                         if(checkArguments == TRUE){
                           validInteger(coords, validValues = (ncol(object@fieldDim) + 1):ncol(object@contrast), refuse.duplicates = TRUE, method = method)
                         }
                         
                         param <- selectParameter(object, type = "contrast")[param]
                         
                       } else if(identical(param, FALSE)){
                         
                         param <- NULL
                         
                       }
                       
                     }
                     
                     
                     #### test
                     if(checkArguments == TRUE){
                       
                       param_data <- selectParameter(object, type = "contrast")
                       if(accept.coords == TRUE){
                         param_data <- c(param_data, selectParameter(object, type = "coords"))
                       }
                       
                       if(any(param %in% param_data == FALSE) || length(unique(param)) != length(param)){
                         stop(method, "[MRIaggr] : wrong specification of \'", arg_name, "\' \n", 
                              "unknown ", long_name, " : ", paste(param[param %in% param_data == FALSE], collapse = " "), " \n", 
                              "duplicated ", long_name, " : ", paste(unique(param[duplicated(param)]), collapse = " "), "\n", 
                              "available ", long_name, " : ", paste(param_data, collapse = " "), "\n")
                       }
                     }
                     
                     #### export
                     return(invisible(param))
                   }
)

####>>> initSlice_var ####
methods::setMethod(f  = "initSlice_var", 
                   signature  = "MRIaggr", 
                   definition = function(object, slice, slice_var, checkArguments, init, method){
                     
                     if(checkArguments == TRUE){	
                       validCharacter(value = slice_var, validLength = 1, validValues = names(object@fieldDim), method = paste(method,"[MRIaggr]",sep=""))
                     }
                     
                     if(init == TRUE && is.null(slice)){
                       slice <- seq(1, object@fieldDim[1, slice_var])
                     }
                     
                     if(checkArguments == TRUE){	
                       validInteger(value = slice, validLength = NULL, validValues = seq(1, object@fieldDim[1,slice_var], by = 1), 
                                    refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = TRUE, method = paste(method,"[MRIaggr]",sep=""))						 
                     }
                     
                     return(slice)
                   }
)

####>>> initSubset ####

methods::setMethod(f  = "initSubset", 
                   signature  = "MRIaggr", 
                   definition = function(object, subset, checkArguments = TRUE, operator.withinR = "union", operator.betweenR = "union",
                                         arg_name = "subset", method){
                     
                     n.subset <- length(subset)
                     names.Merged <- selectRegion(object, type = "names")
                     names.Contrast <- selectParameter(object, type = "contrastOnly")
                     if(is.character(subset) || is.list(subset) ){
                       subset.save <- subset
                       names.subset <- names(subset)
                     }
                     
                     #### tests 
                     if(checkArguments){
                      
                       if( is.list(subset) ){
                         # ok selectRegion will do the checks
                       }else if(is.character(subset)){
                         validCharacter(subset, name = "subset", validLength = NULL, validValues = names.Merged, method = method)
                       }else{
                         if(!is.null(subset)){
                           
                           test.integer <- is.integer(subset) == FALSE
                           test.min <- min(subset) < 1
                           test.max <- max(subset)>selectN(object)
                           test.duplicates <- any(duplicated(subset))
                           
                           if(test.integer || test.min|| test.max || test.duplicates ){
                             stop(method, ": wrong specification of \'",arg_name,"\' \n",
                                  "",subset," must be a vector of integer or characters (corresponding to regions stored in the object) \n",
                                  if(test.integer){paste0("is(",subset,"): ",paste( is(subset), collapse = " "),"\n")},
                                  if(test.min){paste0("duplicates: ",any(duplicated(subset)),"\n")},
                                  if(test.max){paste0("min(subset): ",min(subset),"\n")},
                                  if(test.duplicates){paste0("max(subset): ",max(subset),"\n")})
                           } 
                         }
                       }
                       
                       validCharacter(operator.withinR, name = "operator", validLength = 1, validValues = c("none","union","intersect"), method = method)
                       validCharacter(operator.betweenR, name = "operator", validLength = 1, validValues = c("none","union","intersect"), method = method)
                       if( sum(c(operator.withinR, operator.betweenR) == "none") == 1  ){
                         stop(method, ": wrong specification of \'operator.withinR\' or \'operator.betweenR\' \n",
                              "either both equals \"none\" or none \n",
                              "operator.withinR: ",operator.withinR,"\n",
                              "operator.betweenR: ",operator.betweenR,"\n")
                       }

                     }
                     
                     #### main
                     if(is.integer(subset)){
                       
                     } else if(operator.withinR == "none" && operator.betweenR == "none"){
                       
                       if(is.character(subset.save)){
                         
                       subset <- selectRegion(object, region = subset.save)
                       
                       }else{ # if(is.list(subset.save))
                         
                         subset <- lapply(1:n.subset, function(x){
                           if(is.null(subset.save[[x]])){
                             subset.save[[x]] <- selectRegion(object, region = names.subset[x], type = "names")
                             }
                           selectRegion(object, region = names.subset[x], region.value = subset.save[[x]])
                         })
                         names(subset) <- names.subset
                       }
                       
                     }else{ #  if(operator.withinR != "none" && operator.betweenR != "none")
                       
                       if(is.character(subset.save)){
                         
                         subset <- Reduce(operator.betweenR,
                                          lapply(subset.save, function(x){
                                            Reduce(operator.withinR, 
                                                   selectRegion(object, region = x)[[1]]
                                            )
                                          }))
                         
                       }else{
                         
                         for(iter_region in 1:n.subset){
                   
                           if(is.null(subset.save[[iter_region]])){
                             subset.save[[iter_region]] <- selectRegion(object, region = names.subset[iter_region], type = "names")
                             }
                           
                           if(iter_region == 1){
                             subset <- Reduce(operator.withinR, 
                                              selectRegion(object, region = names.subset[1], region.value = subset.save[[1]])
                             )
                           }else{
                             
                             subset <- do.call(operator.betweenR, 
                                               list(subset, 
                                               Reduce(operator.withinR, 
                                                      selectRegion(object, region = names.subset[iter_region], region.value = subset.save[[iter_region]])
                                               ))
                             )
                           }
                           
                         }
                         
                       }
                         
                       
                     } 
                     
                     #### export 
                     return(subset)
                   }
)
