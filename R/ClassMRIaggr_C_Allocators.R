#**********************************************************************
#**********************************************************************
#*************         0C Class MRIaggr             *******************
#**********************************************************************
#**********************************************************************

##### C) Allocators ############################################
# allocClinic         : ok 
# allocContrast       : ok
# allocDescStats      : to be checked
# allocHemisphere
# allocNormalization
# allocTable
# allocW
# supprContrast       : ok
# supprDescStats      : ok


####>>> allocClinic ####

methods::setReplaceMethod(f = "allocClinic", 
                          signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, add = FALSE, 
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"),
                                                value){
                            
                            #### test
                            if(checkArguments){
                              
                              validLogical(value = add, validLength = 1, method = "allocClinic")
                              validLogical(value = overwrite, validLength = 1, method = "allocClinic")
                              
                              if(!is.data.frame(value)){
                                stop("allocClinic[MRIaggr] : wrong specification of \'value\' \n", 
                                     "\'value\' must be a data.frame \n", 
                                     "is(value) : ", paste(is(value), collapse = " "), " \n")
                              }
                              
                            }
                            
                            #### main
                            if(ncol(object@clinic) == 0){
                              
                              ## update
                              object@clinic <- value
                              validObject(object)
                              
                              ## display
                              if(verbose){
                                cat("allocClinic[MRIaggr] : @clinic has been allocated \n", sep = "")                                       
                              }
                              
                            }else{
                              if(verbose){sauveClinic <- object@clinic}
                              
                              names_alloc <- names(value)[names(value) %in% names(object@clinic) == FALSE]
                              names_replace <- names(value)[names(value) %in% names(object@clinic)]  
                              
                              if(length(names_alloc) > 0){     
                                ## update
                                object@clinic <- data.frame(object@clinic, value)                       
                              }
                              if(length(names_replace) > 0){
                                
                                if(overwrite == FALSE){
                                  stop("allocClinic[MRIaggr] : \'value\' contains elements already present in \'@clinic\' \n", 
                                       "clinical parameter", if(length(names_replace) == 1){"s"}, " already present in @clinic : ", paste(names_replace, collapse = " "), "\n", 
                                       "set \'overwrite\' to TRUE to overwrite them \n")
                                }
                              }else{
                                ## update
                                object@clinic[,names_replace] <- value[,names_replace]
                              }
                              
                              validObject(object)
                              
                              ## display
                              if(length(names_alloc) > 0){
                                cat("allocClinic[MRIaggr] : parameter", if(length(names_alloc) > 1){"s"}, 
                                    " \"", paste(names_alloc, collapse = "\" \""), "\" \n", 
                                    "                       ", if(length(names_alloc) == 1){"has"}else{"have"}, " been added to @clinic \n", sep = "")
                              }
                              if(length(names_replace) > 0){
                                cat("allocClinic[MRIaggr] : parameter", if(length(names_replace) > 1){"s"}, 
                                    " \"", paste(names_replace, collapse = "\" \""), "\" \n", 
                                    "                       ", if(length(names_replace) == 1){"has"}else{"have"}, " been updated in @clinic \n", sep = "")
                              }
                              
                            }
                            
                            return(object)
                          }
)

####>>> allocContrast ####

methods::setReplaceMethod(f  = "allocContrast", 
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, param = NULL, ls.MergeParam = NULL, default_value = NULL, refValue = 0,
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            #### preparation
                            if(is.null(param)){
                              param.value <- names(value)
                            }else{
                              param.value <- param
                            }
                            
                            if(!is.data.table(value)){
                              value <- as.data.table(value)
                            }
                            coords.value <- key(value)
                            param.object <- selectParameter(object, type = "contrast")
                            
                            #### tests
                            
                            if(checkArguments){
                              
                              if(length(coords.value) > 0){ ## optionnal checking of the coordinates
                                coords.object <- selectParameter(object, type = "coords")
                                
                                if(verbose > 1){cat("* check matching between coordinates of \'value\' and \'object@contrast\'")}
                                
                                if( any(coords.value != coords.object) ){
                                  stop("allocContrast[MRIaggr] : keys do not match between \'value\' and \'object@contrast\' \n", 
                                       "names of the coordinates in object:",paste(coords.object, collapse = " "),"\n", 
                                       "key(value):",paste(coords.value, collapse = " "),"\n")
                                }
                                
                                if( nrow(value) != selectN(object) ){
                                  stop("allocContrast[MRIaggr] : number of observations do not match between \'value\' and \'object@contrast\' \n", 
                                       "number of observations in object:", selectN(object),"\n", 
                                       "nrow(value):",nrow(value),"\n")
                                }
                                
                                if(any(value[,coords.value, with = FALSE] - selectCoords(object)) ){
                                  diff <- value[,coords.value, with = FALSE] - selectCoords(object)
                                  stop("allocContrast[MRIaggr] : coordinates do not match between \'value\' and \'object@contrast\' \n", 
                                       "number of mismatch by coordinate (", paste(coords.value, collapse = " ; "),") : ", paste( apply(abs(diff)>0,2, sum), collapse = " ; "), "\n")
                                }
                                
                                # remove coords
                                value[,index_coords := NULL, drop = FALSE]
                              }else{
                                validDimension(value1 = value, validDimension = selectN(object), name1 = "value", name2 = "@contrast", type = "nrow", method = "allocContrast[MRIaggr]")
                              }
                              
                              if(length(param.value) != NCOL(value)){
                                stop("allocContrast[MRIaggr] : mismatch between \'param\' and the number of parameters in \'value\' \n", 
                                     "length(param) : ", length(param), " (", paste(param, collapse = " "), ") \n", 
                                     "ncol(value) : ", ncol(value), " (", paste(names(value), collapse = " "), ") \n")
                              }
                              
                            }
                            
                            
                            #### preparation 
                            if(!is.null(param)){
                              setnames(value, param.value)
                            }
                            
                            ## merge parameters
                            res <- initMergeParam(ls.MergeParam = ls.MergeParam, param = param.value,  refValue = refValue,
                                                  checkArguments = checkArguments, init = TRUE, method = "constMRIaggr")
                            name.merge <-  res$name.merge
                            name2merge <-  res$name2merge
                            n.merge <- res$n.merge
                            ls.indexMerge <- res$ls.indexMerge
                            
                            ## handle new and update
                            param.update <- setdiff(param.value[param.value %in% param.object],name2merge)
                            param.new <- setdiff(param.value[param.value %in% param.object == FALSE],name2merge)
                            
                            if(length(name.merge)>0){
                              paramList.update <- name.merge[name.merge %in% param.object]
                              n.paramList.update <- length(paramList.update)
                            }else{
                              paramList.update <- NULL
                              n.paramList.update <- 0
                            }
                            
                            if(length(name.merge)>0){
                              paramList.new <- name.merge[name.merge %in% param.object == FALSE]
                              n.paramList.new <- length(paramList.new)
                            }else{
                              paramList.new <- NULL
                              n.paramList.new <- 0
                            }
                            
                            #### tests 2
                            if(checkArguments){
                              
                              if("rowNames" %in% param.value){
                                stop("allocContrast[MRIaggr] : cannot allocate a parameter with name \"rowNames\" \n", 
                                     "\"rowNames\" is a reserved name \n", 
                                     "column of concern : ", which(param.value == "rowNames"), "\n")   
                              }
                              
                              if("hemisphere" %in% param.value){  
                                validCharacter(value = value[,"hemisphere"], name = "hemisphere", validLength = NULL, validValues = c("left", "right", "undefined"), method = "allocContrast")  
                              }
                              
                              if( (length(param.update) > 0 || length(paramList.update) > 0) && overwrite == FALSE){ 
                                stop("allocContrast[MRIaggr] : -  some names of \'param\' are already present in \'object@contrast\' \n", 
                                     "redundant names : ", paste(c(param.update,paramList.update), collapse = " "), "\n", 
                                     "set \'overwrite\' to TRUE to perform this operation \n")
                              }
                              
                              if(n.paramList.update > 0){
                                test.names <- lapply(paramList.update, function(x){ls.MergeParam[[x]] %in% selectRegion(object, region = x, type = "names")})
                                
                                if(any( unlist(test.names) )){
                                  stop("allocContrast[MRIaggr] : - wrong specification of \'ls.MergeParam\' \n", 
                                       "the allocator cannot update existing levels of a list parameter  \n", 
                                       "list parameter: \"",paste(name.merge[which(unlist(lapply(test.names,any)))], collapse = "\" \""),"\" \n",
                                       "existings levels: \"",paste(unlist(lapply(1:n.merge,function(x){ls.MergeParam[[x]][which(test.names[[x]])]})), collapse = "\" \""),"\" \n")
                                }
                                
                              }
                              if(is.null(default_value)){
                                default_value <- data.frame(matrix(NA, ncol = length(c(param.new,paramList.new)), nrow = 1))
                                names(default_value) <- c(param.new,paramList.new)
                              }
                              
                              if( (is.data.frame(default_value) == FALSE || sum(names(default_value) != c(param.new,paramList.new)) > 0 ) ){
                                stop("allocContrast[MRIaggr] : wrong specification of \'default_value\' \n", 
                                     "\'default_value\' must be a data.frame with names : \"", paste(c(param.new,paramList.new), collapse = "\" \""), "\" \n", 
                                     "is(default_value): ", paste(is(default_value), collpase = " "), "\n", 
                                     "names(default_value) : \"", paste(names(default_value), collapse = "\" \""), "\"\n")
                              }
                              
                            }
                            
                            #### main
                            if(n.paramList.new > 0){
                              object@contrast[,paramList.new := .(list()), with = FALSE]
                            }
                            
                            
                            
                            
                            for(iter_param in 1:length(param.value)){
                              
                              param_tempo <- param.value[iter_param]
                              
                              if(param_tempo %in% name2merge) { ## merge
                                index.Merge <- which(unlist(lapply(ls.MergeParam, function(x){param_tempo %in% x})))
                                ls.indexMerge[[index.Merge]] <- c(ls.indexMerge[[index.Merge]], 
                                                                  setNames(list(which(value[[param_tempo]]!=refValue)),  param[iter_param]))
                                
#                                 object@contrast[ls.indexMerge[[index.Merge]][[param_tempo]], name.merge[index.Merge] :=  list(lapply(.SD[[1]], FUN = function(x){c(x, param_tempo)})),
#                                                 .SDcol = name.merge[index.Merge], with = TRUE] 
                                
                                object@contrast[ls.indexMerge[[index.Merge]][[param_tempo]], name.merge[index.Merge] :=  list(lapply(1:.N, FUN = function(x){
                                  eval(parse(text = paste(
                                    "c(.SD[[1]][[x]],list(",param_tempo,"= value[ls.indexMerge[[index.Merge]][[param_tempo]][x]][[param_tempo]] ))"
                                  )))
                                })),
                                .SDcol = name.merge[index.Merge], with = TRUE] 
                                
                                
                                
                                if(name.merge[index.Merge] %in% paramList.new){ ## new default value
                                  
                                  if(name.merge[index.Merge] %in% names(object@default_value) == FALSE){ # if not already updated
                                    object@default_value <- data.frame(object@default_value, 
                                                                       default_value[,name.merge[index.Merge], drop = FALSE], stringsAsFactors = FALSE)
                                  }
                                  
                                }else{
                                  
                                  if(name.merge[index.Merge] %in% names(default_value)){ 
                                    object@default_value[,name.merge[index.Merge]] <- default_value[,name.merge[index.Merge], drop = FALSE]
                                  }
                                  
                                }
                                
                              }else{ ## standard parameter
                                object@contrast[,param_tempo := value[[param_tempo]], with = FALSE]
                                
                                if(param_tempo %in% param.new){ ## update default value
                                  object@default_value <- data.frame(object@default_value, 
                                                                     default_value[,param_tempo, drop = FALSE], stringsAsFactors = FALSE)
                                }else{
                                  if(param_tempo %in% names(default_value)){ 
                                    object@default_value[,param_tempo] <- default_value[,param_tempo, drop = FALSE]
                                  }
                                }
                              }
                              
                            }
                            
                            #### update index region
                            if(n.merge>0){
                              region <- selectRegion(object)
                              
                              if(n.paramList.update > 0){
                                
                                for(iter_update in 1:n.paramList.update){
                                  region[[paramList.update[iter_update]]] <- c(region[[paramList.update[iter_update]]],
                                                                               ls.indexMerge[[paramList.update[iter_update]]])
                                }
                                
                              }
                              if(n.paramList.new > 0){
                                region <- c(region,
                                            ls.indexMerge[paramList.new])  
                              }
                              
                              object@region$contrast <- region
                            }
                            
                            #### check object
                            validObject(object)
                            
                            #### display
                            if(verbose > 0){
                             
                              n.update <- length(c(param.update,paramList.update))
                              if(n.update > 0){                    
                                cat("allocContrast[MRIaggr] : Cartography \"", paste(c(param.update,paramList.update), collapse = "\" \""), "\" \n", 
                                    "                         ",if(n.update==1){"has"}else{"have"}," been updated \n", sep = "")                      
                                
                              }
                              
                              
                              n.new <- length(c(param.new,paramList.new))
                              if(n.new > 0){                    
                                cat("allocContrast[MRIaggr] : Cartography \"", paste(c(param.new,paramList.new), collapse = "\" \""), "\" \n", 
                                    "                         ",if(n.new==1){"has"}else{"have"}," been allocated \n", sep = "")                      
                                
                              }     
                              
                            }
                            
                            #### export
                            return(object)
                          }
)

####>>> allocDescStats ####

methods::setReplaceMethod(f  = "allocDescStats", # $contrast must be protected against the user
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, name, 
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value)
                          { #### preparation
                            test.overwrite <- name %in% selectParameter(object, type = "ls_descStats")
                            
                            #### tests 
                            if(checkArguments){
                              
                              if(length(name) > 1){
                                stop("allocDescStats[MRIaggr] : Only one element can be allocated at a time \n",            
                                     "length(name) : ", name, "\n")
                              }                   
                              
                              if(sum(test.overwrite) > 0 && overwrite == FALSE){
                                stop("allocDescStats[MRIaggr] : The requested field(s) already exist in object@ls_descStats \n", 
                                     "Set \'overwrite\' to TRUE to replace the field \n", 
                                     "already existing fields : ", name[test.overwrite], "\n")
                              }
                              
                            }
                            
                            #### initialisation
                            param.update <- names(object@ls_descStats)[(names(object@ls_descStats) %in% name)]
                            param.new <- name[(name %in% names(object@ls_descStats) == FALSE )]
                            
                            #### main
                            if(length(param.update) > 0){
                              object@ls_descStats[[param.update]] <- value 
                            }
                            
                            if(length(param.new) > 0){                    
                              object@ls_descStats <- c(object@ls_descStats, list(value))
                              names(object@ls_descStats)[length(object@ls_descStats)] <- param.new
                            }
                            
                            
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              
                              if(length(param.update) == 1){
                                cat("allocDescStats[MRIaggr] : Element \"", param.update, "\" \n", 
                                    "                          has been updated \n", sep = "")
                              }
                              if(length(param.update) > 1){
                                cat("allocDescStats[MRIaggr] : Elements \"", paste(param.update, collapse = "\" \""), "\" \n", 
                                    "                          have been \n", sep = "")
                              }
                              
                              if(length(param.new) == 1){
                                cat("allocDescStats[MRIaggr] : Element \"", param.new, "\" \n", 
                                    "                          has been allocated \n", sep = "")
                              }
                              if(length(param.new) > 1){
                                cat("allocDescStats[MRIaggr] : Elements \"", paste(param.new, collapse = "\" \""), "\" \n", 
                                    "                          have been allocated \n", sep = "")
                              }                     
                            }
                            
                            #### export
                            return(object)
                          }
)

####>>> allocHemisphere ####

methods::setReplaceMethod(f = "allocHemisphere",
                          signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, 
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            if("midplane" %in% names(value) || verbose){
                              midplane_tempo <- selectMidplane(object)
                            }
                            if("hemispheres" %in% names(value) || verbose){
                              hemispheres_tempo <- selectHemispheres(object)
                            }
                            
                            #### tests
                            if(checkArguments){
                              if(!is.list(value) || is.null(names(value)) || any(names(value) %in% c("midplane", "hemispheres", "data") == FALSE)){
                                stop("allocHemisphere[MRIaggr] : wrong specification of \'value\' \n", 
                                     "\'value\' must be a list containing some of the following elements \"midplane\" \"hemispheres\" \"data\" \n", 
                                     "incorrect names : ", paste(names(value)[names(value) %in% c("midplane", "hemispheres", "data") == FALSE], collapse = " "), "\n", 
                                     "is(value)  : ", paste(is(value), collapse = " "), "\n"
                                )
                              }
                            }
                            
                            if("midplane" %in% names(value)){
                              
                              if(!is.null(midplane_tempo) && any(!is.na(midplane_tempo)) && overwrite == FALSE){
                                stop("allocHemisphere[MRIaggr] : midplane already existing in \'object\' \n",                          
                                     "set \'overwrite\' to TRUE to replace it \n")
                              }
                              object@region$midplane <- value$midplane
                            }
                            
                            if("hemispheres" %in% names(value)){
                             
                              if( (hemispheres_tempo$right != "undefined" || hemispheres_tempo$left != "undefined") && overwrite == FALSE){
                                stop("allocHemisphere[MRIaggr] : hemispheres already existing in \'object\' \n",                          
                                     "set \'overwrite\' to TRUE to replace it \n")
                              }
                              object@region@hemispheres <- value$hemispheres
                            }
                            
                            if("data" %in% names(value)){
                              
                              if(!is.data.frame(value$data)){
                                stop("allocHemisphere[MRIaggr] : wrong specification of the data element of \'value\' \n",                          
                                     "data must be a data.frame \n", 
                                     "type of data : ", paste(is(value$data), collapse = " "), "\n")
                              }							  
                              validDimension(value1 = value$data, validDimension = selectN(object), name1 = "data", name2 = NULL, type = "nrow", method = "allocHemisphere[MRIaggr]")
                              validNames(value = value$data, name = "data", validValues = c("hemisphere", "i_hemisphere", "j_hemisphere"), method = "allocHemisphere[MRIaggr]")
                              
                              allocContrast(object, overwrite = overwrite) <- value$data
                            }                 
                            
                            validObject(object)
                            
                            if(verbose){
                              if("midplane" %in% names(value)){
                                cat("allocHemisphere[MRIaggr] : @midplane has been ")
                                if(sum(!is.na(midplane_tempo)) == 0){cat("alllocated \n")}else{cat("updated\n ")}
                              }
                              
                              if("hemispheres" %in% names(value)){
                                cat("allocHemisphere[MRIaggr] : @hemispheres has been ")
                                if(sum(hemispheres_tempo != "undefined") == 0){cat("allocated \n")}else{cat("updated\n ")}
                              }                       
                            }
                            
                            return(object)
                          }
)

####>>> allocNormalization ####

methods::setReplaceMethod(f = "allocNormalization", 
                          signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, 
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            if(verbose){sauveNormalization <- object@normalization}
                            
                            #### tests
                            if(checkArguments){
                            if(overwrite == FALSE && !is.null(object@normalization) && length(object@normalization) > 0 && any(!is.na(object@normalization))){
                              stop("allocNormalization[MRIaggr] : normalization already existing in \'object\' \n",                          
                                   "set \'overwrite\' to TRUE to replace it \n")
                              
                            }
                            
                            valid_names <- c("norm_global", 
                                             "normMu_slice_both", "normSigma_slice_both", 
                                             "normMu_slice_left", "normSigma_slice_left", 
                                             "normMu_slice_right", "normSigma_slice_right", 
                                             "normMu_3slices_both", "normSigma_3slices_both", 
                                             "normMu_3slices_left", "normSigma_3slices_left", 
                                             "normMu_3slices_right", "normSigma_3slices_right") 
                            
                            validNames(value = value, validValues = valid_names,  method = "allocNormalization[MRIaggr]")
                            }
                            
                            #### main
                            object@normalization <- value
                            
                            #### check
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              cat("allocNormalization[MRIaggr] : @normalization has been ")
                              if(length(sauveNormalization) == 0){cat("allocated \n")}else{cat("updated\n ")}                    
                            }
                            
                            #### export
                            return(object)
                          }
)

####>>> allocTable ####

methods::setReplaceMethod(f = "allocTable", 
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, type,
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            
                            #### tests
                            if(checkArguments){
                              
                              validCharacter(value = type, validLength = 1, validValues = c("lesion", "reperfusion", "hypoperfusion"), method = "allocTable[MRIaggr]")
                              
                              if(!is.null(selectTable(object, type = type)) && any(!is.na(selectTable(object, type = type))) && overwrite == FALSE){
                                stop("allocTable[MRIaggr] : Table already existing in \'object\' \n",                          
                                     "set \'overwrite\' to TRUE to replace it \n")
                              }
                              
                            }
                            
                            #### main
                            if(type == "lesion"){
                              if(verbose){saveTable<- selectTable(object, type = "lesion")}
                              object@table$lesion <- value
                            }
                            
                            if(type == "reperfusion"){
                              if(verbose){saveTable<- selectTable(object, type = "reperfusion")}
                              object@table$reperfusion <- value
                            }                 
                            
                            if(type == "hypoperfusion"){
                              if(verbose){saveTable<- selectTable(object, type = "hypoperfusion")}
                              object@table$hypoperfusion <- value
                            }
                            
                            #### check
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              cat(paste("allocTable[MRIaggr] : @table$", type, " has been ", sep = ""))
                              if(length(saveTable) == 0){cat(" allocated \n")}else{cat(" updated \n ")}
                            }
                            
                            #### export
                            return(object)
                          }              
)

####>>> allocW ####

methods::setReplaceMethod(f  = "allocW", 
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, type, 
                                                overwrite = FALSE, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            #### tests
                            if(checkArguments){
                            validCharacter(value = type, validLength = 1:3, validValues = c("Wmatrix", "Wblocks", "upper"), method = "allocW[MRIaggr]")
                            
                            if(!is.list(value)){
                              stop("allocW[MRIaggr] : wrong specification of \'value\' \n",            
                                   "\'value\' must be a list \n", 
                                   "is(value) : ", paste(is(value), collapse = " "), "\n")
                            }
                            
                            if(any(type %in% names(value) == FALSE)){
                              stop("allocW[MRIaggr] : wrong specification of \'value\' \n",            
                                   "\'value\' does not contains elements announced by \'type\' \n", 
                                   "missing elements : ", paste(type[type %in% names(value) == FALSE], collapse = " "), "\n")
                            }
                            }
                            
                            #### main
                            if("Wmatrix" %in% type){
                              object@W$Wmatrix <- Matrix::drop0(value$Wmatrix, is.Csparse = TRUE)
                              object@W$upper <- NULL # indicate the presence of a neighborhood matrix
                            }
                            
                            if("Wblocks" %in% type){
                              object@W$Wblocks <- value$Wblocks
                            }
                            
                            if("upper" %in% type){
                              object@W$upper <- value$upper
                            }
                            
                            #### check
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              new <- type[type %in% names(value)]
                              
                              cat("allocW[MRIaggr] : Element", if(length(new) > 1){"s"}, " \"", paste(new, collapse = "\" \""), "\" ", 
                                  if(length(new) > 1){"have"}else{"has"}, " been updated \n", sep = "")
                            }
                            
                            #### export
                            return(object)
                          }
)

####>>> supprContrast ####

methods::setReplaceMethod(f  = "supprContrast", 
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            #### preparation
                            nom <- initParameter(object = object, param = value, checkArguments = checkArguments, init = TRUE, accept.coords = FALSE,
                                                 arg_name = "value", method = "supprContrast")  
                            
                            #### main
                            object@contrast[,value := NULL, with = FALSE]
                            object@default_value <- object@default_value[,-which( names(object@default_value) == nom ), drop = FALSE]
                            
                            region <- selectRegion(object)
                            if(length(region) > 0 &&  value %in% names(region)){
                              object@region$contrast <- region[-which(names(region) == value)]
                            }
                            
                            #### check
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              cat("supprContrast[MRIaggr] : ")
                              if(length(nom) == 1){
                                cat("Cartography \"", nom, "\" \n", 
                                    "                         has been removed \n", sep = "")
                              }else{
                                cat("Cartographies \"", paste(nom, collapse = "\" \""), "\" \n", 
                                    "                         have been removed \n", sep = "")                      
                              }                  
                            }
                            
                            #### export
                            return(object)
                          }
)

####>>> supprDescStats ####

methods::setReplaceMethod(f  = "supprDescStats", 
                          signature  = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                          definition = function(object, 
                                                checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"), 
                                                value){
                            
                            #### tests
                            if(checkArguments){
                            validCharacter(value = value, validLength = 1, validValues = selectParameter(object, "ls_descStats"), 
                                           method = "supprDescStats[MRIaggr]")
                            }
                            
                            #### main
                            for(iter_l in value){
                              object@ls_descStats[[iter_l]] <- NULL
                            }
                            
                            #### check
                            validObject(object)
                            
                            #### display
                            if(verbose){
                              cat("supprDescStats[MRIaggr] : ")
                              if(length(value) == 1){
                                cat("Element \"", value, "\" \n", 
                                    "                          has been removed \n", sep = "")
                              }else{
                                cat("Elements \"", paste(value, collapse = "\" \""), "\" \n", 
                                    "                          have been \n", sep = "")                      
                              }                     
                            }
                            
                            return(object)
                          }
)

