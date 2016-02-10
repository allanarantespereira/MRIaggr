#**********************************************************************
#**********************************************************************
#*************         0D3 Class MRIaggr            *******************
#**********************************************************************
#**********************************************************************
#

#####  D3) Methods const ############################################################
# constCompressMRIaggr: ok (mask has not be tested)
# constReduceMRIaggr:   ok

####>>> constCompressMRIaggr ####

methods::setMethod(f  = "constCompressMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, compression.factor, slice_var = "k", param = TRUE, 
                                         mask = NULL, threshold = 0.49, 
                                         checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose")){
                     
                     #### preparation
                     var.coords <- selectParameter(object, type = "coords")
                     compression_var <- setdiff(var.coords, slice_var)
                     n.inital <- selectFieldDim(MRIaggr.Pa1, coords = compression_var, format = "vector")
                     n.slices <- selectFieldDim(MRIaggr.Pa1, coords = slice_var, format = "vector")
                     
                     param <- initParameter(object, param = param, init = TRUE, checkArguments = checkArguments, 
                                            accept.coords = FALSE, method = "constCompressMRIaggr")
                     n.param <- length(param)
                     
                     data.initial <- selectContrast(object, param = param, coords = TRUE)
                     
                     #### tests ####
                     if(checkArguments){
                       
                       validNumeric(value = compression.factor, name = "compression.factor", validLength = 1, 
                                    min = 1, method = "constCompressMRIaggr[MRIaggr]")
                       
                       validCharacter(value = slice_var, name = "slice_var", validLength = 1, method = "constCompressMRIaggr[MRIaggr]")
                       
                       if(n.inital[1] %% compression.factor > 0 || n.inital[2] %% compression.factor > 0){
                         stop("constCompressMRIaggr[MRIaggr] : wrong specification of \'compression.factor\' \n", 
                              "@fieldDim$",compression_var[1]," or @fieldDim$",compression_var[2]," is not a multiple of \'compression.factor\' \n", 
                              "@fieldDim$",compression_var[1]," / compression.factor : ", n.inital[1] / compression.factor, "\n", 
                              "@fieldDim$",compression_var[2]," / compression.factor : ", n.inital[2] / compression.factor, "\n")
                       }
                       
                       test.character <- data.initial[,lapply(.SD, is.numeric), .SDcols = param]
                       if(any(test.character == FALSE)){              
                         stop("constCompressMRIaggr[MRIaggr] : wrong specification of \'param\' \n", 
                              "param must only correspond to numeric parameters \n", 
                              "non numeric parameters: \"", paste(param[test.character == FALSE], collapse = "\" \""), "\" \n")
                       }
                       
                     }
                     
                     #### initialization ####
                     
                     n.final <- n.inital / compression.factor
                     n_px.px <- compression.factor^2
                     
                     data.final <- data.table(matrix(integer(0), nrow = prod(n.final) * n.slices, ncol = length(var.coords) ),
                                              matrix(numeric(0), nrow = prod(n.final) * n.slices, ncol = n.param)
                     )
                     names(data.final) <- names(data.initial)
                     
                     #### mise en place ####
                     seq_coord1 <- matrix(NA, ncol = n.inital[1], nrow = compression.factor)
                     for(iter_k in 1:compression.factor){
                       seq_coord1[iter_k,] <- seq(iter_k, by = compression.factor, length.out = n.final[1]) 
                     }
                     
                     seq_coord2 <- matrix(NA, ncol = n.inital[2], nrow = compression.factor)
                     for(iter_k in 1:compression.factor){
                       seq_coord2[iter_k,] <- seq(iter_k, by = compression.factor, length.out = n.final[2]) 
                     }
                    
                     index_carto <- 0
                     
                     #### loop ####
                     if(verbose > 1){cat("slice ",slice_var," = ")}
                     
                     for(iter_slice in 1:n.slices){
                       if(verbose > 1){cat(iter_slice, " ", sep = "")}
                       
                       index_slice <- which(data.initial[[slice_var]] == iter_slice)
                   
                       data.slice <- dt2array(contrast = data.initial[index_slice, param, with = FALSE],
                                              coords = data.initial[index_slice, compression_var, with = FALSE],
                                              range.coords = n.inital)$contrast
                       
                       for(iter_param in 1:n.param){              
                         nom_param <- param[iter_param]
                         matrix.init_tempo <- data.slice[[iter_param]]
                         matrixNA.init_tempo <- is.na(matrix.init_tempo)
                         
                         matrix.col_tempo <- matrix(NA, nrow = n.final[1], ncol = n.inital[2])
                         matrixNA.col_tempo <- matrix(NA, nrow = n.final[1], ncol = n.inital[2])
                         matrix.final_tempo <- matrix(NA, nrow = n.final[1], ncol = n.final[2])
                         matrixNA.final_tempo <- matrix(NA, nrow = n.final[1], ncol = n.final[2])
                         
                         for(iter_row in 1:n.final[1]){
                           matrix.col_tempo[iter_row,] <- colSums(matrix.init_tempo[seq_coord1[,iter_row],])
                           matrixNA.col_tempo[iter_row,] <- colSums(matrixNA.init_tempo[seq_coord1[,iter_row],])
                         }
                         
                         for(iter_col in 1:n.final[2]){
                           matrix.final_tempo[,iter_col] <- rowSums(matrix.col_tempo[,seq_coord2[,iter_col]])
                           matrixNA.final_tempo[,iter_col] <- rowSums(matrixNA.col_tempo[,seq_coord2[,iter_col]])
                         }
                         
                         index_NA <- matrixNA.final_tempo > n_px.px / 2
                         
                         matrix.final_tempo <- matrix.final_tempo / (n_px.px - matrixNA.final_tempo)
                         
                         if(sum(index_NA) > 0){matrix.final_tempo[index_NA] <- NA}
                         
                        
                         new_data <- array2dt(matrix.final_tempo, name_newparam = nom_param, names_coords = compression_var, na.rm = FALSE)
                         
                         if(iter_param == 1){
                           index_carto <- max(index_carto) + 1:nrow(new_data)
                         }     
                         
                         
                         ## update
                         data.final[index_carto, nom_param := new_data[,nom_param, with = FALSE], with =  FALSE ] # update parameter
                         
                         if(iter_param == 1){ # update coordinates
                           data.final[index_carto, compression_var := new_data[,compression_var, with = FALSE], with =  FALSE ]
                           data.final[index_carto, slice_var := iter_slice, with =  FALSE ]
                         }
                         
                       }
                     }
                     if(verbose > 1){cat("\n")}
                     
                     # Reglement des parametres binaires :
                     if(!is.null(mask)){
                       data.final[,mask := lapply(data.final[,mask, with = FALSE], function(x){as.numeric(mask > threshold)}), with = FALSE]
                     }
                     
                     #### new object
                     data.table::setkeyv(data.final, key(data.initial))
                     
                     fieldDim <- selectFieldDim(object)
                     fieldDim[compression_var] <- n.final
                     voxelDim <- selectVoxelDim(object)
                     voxelDim[compression_var] <- unlist(sapply(compression_var,
                                                         function(x){ if(is.numeric(voxelDim[x])){voxelDim[x]/compression.factor}else{voxelDim[x]} }
                                                           ))
                     
                     region <- object@region
                     if("midplane" %in% names(region) && all( is.na(region$midplane) ) == FALSE){
                      region$midplane[compression_var] <- region$midplane[compression_var] / compression.factor
                     }
                     if(length(region$contrast)>0){
                       if(verbose){
                       cat("constCompressMRIaggr[MRIaggr] : the index of the regions stored in @region$contrast has been removed from the new object \n")
                       }
                       region$contrast <- list()
                     }
                     
                     ##
                     y <- new(Class =  "MRIaggr", 
                              clinic = object@clinic,                     
                              contrast = data.final, 
                              default_value = object@default_value[,param, drop = FALSE], 
                              fieldDim = fieldDim, 
                              history = c(object@history, 
                                          list(constCompressMRIaggr = list(call = match.call(call = sys.call(sys.parent())), date = date()))
                              ), 
                              identifier = object@identifier, 
                              ls_descStats = object@ls_descStats,
                              normalization = object@normalization, 
                              region = region,
                              table = object@table,
                              voxelDim = voxelDim, 
                              W = object@W
                     )
                     
                     if(verbose){
                       cat("constCompressMRIaggr[MRIaggr] : MRIaggr has been compressed from (",paste(compression_var, collapse = ","),") = (",paste(n.inital, collapse = ","), ") \n", 
                           "                                                              to (",paste(compression_var, collapse = ","),") = (",paste(n.final, collapse = ","), ") \n", sep = "")
                     }
                     
                     #### export
                     return(y)           
                     
                   }
)

####>>> constReduceMRIaggr ####

methods::setMethod(f  = "constReduceMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask, numeric2logical = FALSE, keep.index = TRUE,
                                         checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose"))
                   { 
                     #### preparation 
                     if(length(mask) == 1 && is.character(mask)){                            
                       mask <- selectContrast(object, param = mask, format = "vector")              
                       supprContrast(object) <- "mask"      
                       test.loadData <- TRUE
                     }else{
                       test.loadData <- FALSE
                     }
                     
                     if(numeric2logical == TRUE){
                       mask <- as.logical(mask)
                     }
                     
                     #### tests
                     if(checkArguments){
                       if(test.loadData == FALSE){
                       validDimension(value1 = mask, validDimension = selectN(object), name1 = "mask", name2 = NULL, type = "length", method = "constReduceMRIaggr[MRIaggr]")
                       }
                       if(numeric2logical == FALSE){  
                        validLogical(value = mask, validLength = NULL, method = "constReduceMRIaggr[MRIaggr]")
                        }
                      }
                     
                     #### initialisation
                     if(keep.index == TRUE){
                       allocDescStats(object, name = "index_sauve") <- selectContrast(object, rowNumber = TRUE, param = FALSE, format = "vector")
                     }
                     index.mask <- which(mask)
                     
                     region <- object@region
                     if(length(region$contrast)>0){ ## update the index
                       region$contrast <- lapply(region$contrast, function(x){
                         lapply(x, function(y){which(index.mask %in% y)})}
                       )
                     }
                     
                     #### main
                     y <- new(Class =  "MRIaggr", 
                              clinic = object@clinic,                     
                              contrast = object@contrast[index.mask], 
                              default_value = object@default_value, 
                              fieldDim = object@fieldDim, 
                              history = c(object@history, 
                                          list(constReduceMRIaggr = list(call = match.call(call = sys.call(sys.parent())), date = date()))
                              ), 
                              identifier = object@identifier, 
                              ls_descStats = object@ls_descStats,
                              normalization = object@normalization, 
                              region =region,
                              table = object@table,
                              voxelDim = object@voxelDim, 
                              W = object@W
                     )
                     
                     cat("constReduceMRIaggr[MRIaggr] : MRIaggr_red has been created \n")
                     
                     return(y)
                   }
)

methods::setMethod(f  = "writeMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, num = NULL, norm_mu = FALSE, norm_sigma = FALSE, range.coords = NULL, default_value = NA, 
                                         filename, format, gzipped = TRUE, verbose = optionsMRIaggr("verbose"), size = "NA_integer_")
                   { 
                     data <- selectContrast(object, param = param, num = num, norm_mu = norm_mu, norm_sigma = norm_sigma)
                     coords <- selectCoords(object, num = num, format = "matrix")
                     
                     array <- dt2array(data, coords = coords, format = "any", 
                                       default_value = default_value, range.coords = range.coords)$contrast[[1]]
                     
                     writeMRI(data = array, filename = filename, format = format, gzipped = gzipped, verbose = verbose, size = size)              
                     
                   }
)

