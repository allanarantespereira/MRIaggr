#**********************************************************************
#**********************************************************************
#*************         0D1 Class MRIaggr            *******************
#**********************************************************************
#**********************************************************************
#


##### D1) Method calc ############################################
# calcBrainMask
# calcContralateral
# calcDistMask
# calcDistTissues
# calcFilter
# calcGroupsMask
# calcHemisphere
# calcNormalization
# calcRegionalContrast
# calcROCthreshold
# calcSmoothMask
# calcTableHypoReperf
# calcROCthreshold
# calcSmoothMask
# calcTableHypoReperf
# calcTableLesion
# calcThresholdMRIaggr
# calcTissueType
# calcW

####>>> calcBrainMask ####

methods::setMethod(f  = "calcBrainMask", 
                   signature = "MRIaggr", 
                   definition = function(object, param, type = "kmeans", 
                                         th.breaks = 100, th.smoothing = TRUE, th.select_optima = 1, th.upper = TRUE, plot = TRUE, 
                                         kmeans.n_groups = 2:4, kmeans.Neighborhood = 3, 
                                         skull.param = NULL, skull.n_groups = 3, 
                                         filename = paste("calcBrainMask", type, object@identifier, sep = "_"), 
                                         update.object = FALSE, overwrite = FALSE, ...){
                     
                     if(plot == TRUE){ #### graphical options
                       ## get graphical arguments
                       optionsMRIaggr.eg <- optionsMRIaggr()					 
                       dots.arguments <- list(...)
                       names_dots.arguments <- names(dots.arguments)
                       
                       validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                      validValues = c("height", "path", "res", "unit", "verbose", "width", "window"), 
                                      refuse.NULL = FALSE, method = "calcBrainMask[MRIaggr]")
                       
                       ## set specific display
                       if(length(names_dots.arguments) > 0){
                         optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                       }
                       
                       ## create all variables
                       height <- optionsMRIaggr.eg$height
                       path <- optionsMRIaggr.eg$path
                       res <- optionsMRIaggr.eg$res
                       unit <- optionsMRIaggr.eg$unit
                       verbose <- optionsMRIaggr.eg$verbose
                       width <- optionsMRIaggr.eg$width
                       window <- optionsMRIaggr.eg$window
                     }
                     
                     #### preparation ####
                     initParameter(object = object, param = param, checkArguments = TRUE, init = FALSE, 
                                   accept.coords = FALSE, accept.index = FALSE, accept.mask = FALSE, method = "calcBrainMask")
                     carto <- selectContrast(object, param = param, format = "data.frame")
                     coords <- selectCoords(object)
                     
                     p <- length(param)
                     validCharacter(value = type, validLength = 1, validValues = c("kmeans", "threshold"), method = "calcBrainMask[MRIaggr]")
                     
                     if(plot == TRUE){
                       scale <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                           n.plot = 1, mfrow = c(1,1), xlim = NULL, ylim = NULL, method = "calcBrainMask[MRIaggr]")$scale
                     }
                     
                     if(type == "threshold"){
                       
                       #### test ####
                       if(length(param) != 1){
                         stop("calcBrainMask[MRIaggr] : wrong specification of  \'param\' \n", 
                              "must have length 1 \n", 
                              "length(param) : ", length(param), "\n")
                       }
                       
                       #### initialization ####
                       res <- list()
                       
                       if(length(th.breaks) == 1){
                         breaks <- seq(min(carto[,param]), max(carto[,param]), length.out = th.breaks)  
                         
                         #                 test.breaks <- sapply(breaks, function(x){sum(carto[,param] > x)})
                         #                 seqY <- seq(test.breaks[1], 0, length.out = th.breaks)
                         #                 breaks <- unique(sapply(seqY, function(x){breaks[which.min(abs(test.breaks-x))]}))
                         
                         breaks <- c(seq(breaks[1], breaks[2], length.out = 12), utils::tail(breaks, length(breaks) - 2))
                         breaks <- c(utils::head(breaks, length(breaks) - 2), seq(breaks[length(breaks) - 1], breaks[length(breaks)], length.out = 12))
                         
                         unvalid.breaks <- c(1:8, (length(breaks) - 7):length(breaks))
                       }else{
                         breaks <- th.breaks              
                         unvalid.breaks <- NULL
                       }
                       n.breaks <-  length(breaks)
                       
                       res$analysis <- matrix(NA, nrow = n.breaks, ncol = 5)
                       colnames(res$analysis) <- c("threshold", "Nb", "dNb", "dNb.filtered", "optima")
                       
                       ####
                       
                       res$analysis[,"threshold"] <- breaks
                       res$analysis[,"Nb"] <- sapply(res$analysis[,"threshold"], 
                                                     function(x){sum(carto[,param] > x)}
                       )
                       
                       res$analysis[,"dNb"] <- c(NA, res$analysis[-n.breaks,"Nb"] - res$analysis[-1,"Nb"])
                       
                       if(th.smoothing > 0 && identical(th.smoothing, TRUE)){
                         res$analysis[,"dNb.filtered"] <- res$analysis[,"dNb"]                
                       }
                       if(th.smoothing > 0 && !identical(th.smoothing, TRUE)){
                         res$analysis[,"dNb.filtered"] <- as.numeric(stats::filter(res$analysis[,"dNb"], stats::dbinom(0:th.smoothing, th.smoothing, 0.5), sides = 2))                
                       }
                       
                       #### smoothing ####
                       cumul <- 0
                       bandwidth <- 2
                       
                       res$analysis[,"dNb"]
                       
                       
                       while(length(cumul) == 1 && bandwidth <= 15){
                         
                         test.post <- c(NA, res$analysis[-n.breaks,"dNb.filtered"] > res$analysis[-1,"dNb.filtered"])
                         unvalid.post <- c(unvalid.breaks, which(is.na(test.post)), which(is.na(test.post)) - 1, length(test.post))
                         res$analysis[,"optima"] <- 0
                         
                         # optima
                         for(iter_th in setdiff(1:length(test.post), unvalid.post)){
                           if(test.post[iter_th] == test.post[iter_th + 1]){
                             cumul[1] <- cumul[1] + 1
                           }else{
                             res$analysis[iter_th, "optima"] <- 1
                             cumul <- c(0, cumul)
                           }
                         }
                         
                         if(cumul[1] <= 1){cumul <- cumul[-1]} # first may be truncated
                         if(cumul[length(cumul)] <= 1){cumul <- cumul[-length(cumul)]} # last may be truncated
                         
                         # smoothing
                         if(identical(th.smoothing, TRUE) && any(cumul <= 1)){     
                           bandwidth <- bandwidth + 1
                           res$analysis[,"dNb.filtered"] <- as.numeric(stats::filter(res$analysis[,"dNb"], stats::dbinom(0:bandwidth, bandwidth, 0.5), sides = 2))                               
                           cumul <- 0
                         }               
                         
                       }
                       
                       #### optima ####
                       optima <- which(res$analysis[,"optima"] == 1)
                       if(length(optima) < th.select_optima){
                         stop("calcBrainMask[MRIaggr] : wrong specification of \'th.select_optima\' \n", 
                              "only ", length(optima), " optima are available \n", 
                              "(", paste(round(res$analysis[optima, "threshold"], optionsMRIaggr("digit.result")), collapse = " "), ") \n", 
                              "requested \'th.select_optima\' : ", th.select_optima, "\n")
                       } 
                       
                       res$th_opt <- data.frame(matrix(NA, ncol = length(optima), nrow = 2))
                       names(res$th_opt) <- 1:length(optima)
                       rownames(res$th_opt) <- c("Th", "derivative")
                       
                       res$th_opt["Th",] <- res$analysis[optima, "threshold"]
                       res$th_opt["derivative",] <- res$analysis[optima, "dNb.filtered"]
                       
                       if(verbose == TRUE){
                         
                         if(th.smoothing > 0){
                           cat("Derivative has been smoothed with a Gaussian kernel of width ", bandwidth, " breaks \n", sep = "")
                         }
                         
                         traceSeuil <- as.numeric(res$th_opt["Th",])
                         traceDerivative <- as.numeric(res$th_opt["derivative",])
                         traceUnit <- nchar(round(traceSeuil))
                         
                         if(verbose == TRUE){cat("Threshold : derivative (selected)\n")}
                         for(iter_optima in 1:length(optima)){   
                           
                           diff <- 6 - nchar(round(traceSeuil[iter_optima], digits = 5-traceUnit[iter_optima]))
                           
                           cat(rep(" ", max(traceUnit) - traceUnit[iter_optima]), 
                               round(traceSeuil[iter_optima], digits = 5 - traceUnit[iter_optima]), rep(" ", diff + traceUnit[iter_optima]), " : ", traceDerivative[iter_optima], sep = "")
                           if(iter_optima == th.select_optima){cat(rep(" ", max(0,10 - nchar(as.numeric(traceDerivative[iter_optima]))) ), " (*)", sep = "")}
                           cat("\n")
                         }
                         
                       }
                       
                       #### select observations ####
                       
                       if(th.upper){
                         res$best_group <- selectContrast(object, param = param, format = "vector") > as.numeric(res$th_opt["Th",th.select_optima])
                       }else{
                         res$best_group <-  selectContrast(object, param = param, format = "vector") < as.numeric(res$th_opt["Th",th.select_optima])
                       }
                       
                       res$mask_name <- paste(param, "threshold", th.select_optima, sep = ".")
                       
                       #### plot ####
                       
                       if(plot == TRUE){
                         
                         initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                           mfrow = c(p, 2), bg = NULL, pty = NULL, mar = rep(3, 4), mgp = c(2, 0.5, 0))
                         
                         graphics::plot(res$analysis[,"threshold"], res$analysis[,"Nb"], xlab = param, ylab = "Nb", main = "Number of voxels", type = "o")
                         graphics::points(breaks[optima[-th.select_optima]], res$analysis[optima[-th.select_optima], "Nb"], col = grDevices::rainbow(length(optima))[-th.select_optima], pch = 15)
                         graphics::points(breaks[optima[th.select_optima]], res$analysis[optima[th.select_optima], "Nb"], col = grDevices::rainbow(length(optima))[th.select_optima], pch = 8, cex = 2)
                         graphics::abline(v = breaks[optima[th.select_optima]], col = grDevices::rainbow(length(optima))[th.select_optima])
                         graphics::legend("topright", legend = 1:length(optima), col = grDevices::rainbow(length(optima)), bty = "n", pch = 20)
                         
                         graphics::plot(res$analysis[,"threshold"], res$analysis[,"dNb.filtered"], xlab = param, ylab = if(bandwidth > 0){paste("dNb.filtered - width = ", bandwidth, sep = "")}else{"dNb"}, main = "Derivative", type = "o")
                         graphics::points(breaks[optima[-th.select_optima]], res$analysis[optima[-th.select_optima],"dNb.filtered"], col = grDevices::rainbow(length(optima))[-th.select_optima], pch = 15)
                         graphics::points(breaks[optima[th.select_optima]], res$analysis[optima[th.select_optima],"dNb.filtered"], col = grDevices::rainbow(length(optima))[th.select_optima], pch = 8, cex = 2)
                         graphics::abline(v = breaks[optima[th.select_optima]], col = grDevices::rainbow(length(optima))[th.select_optima])
                         
                         if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                           grDevices::dev.off()
                         }
                         
                       }
                       
                     }
                     
                     
                     if(type == "kmeans"){
                       
                       #### combinaisons possibles
                       names_combin <- NULL
                       for(iter_group in kmeans.n_groups){
                         names_combin <- c(names_combin, 
                                           paste("G", iter_group, ".=.", 2:iter_group, sep = ""), 
                                           if(iter_group == 3){"G3.>=.2"}, 
                                           if(iter_group > 3){paste("G", iter_group, ".>=.", seq(3, iter_group - 1), sep = "")}, 
                                           if(iter_group > 3){paste("G", iter_group, ".<=.", seq(3, iter_group - 1), sep = "")}
                         )
                       }
                       
                       #### stockage des resultats 
                       res <- list()
                       res$kmeans <- list()
                       res$potential <- data.frame(matrix(NA, ncol = 2, 
                                                          nrow = length(names_combin)))
                       names(res$potential) <- c("nb_groups", "V") 
                       rownames(res$potential) <- names_combin
                       res$best_V <- -1
                       
                       #### potentiel
                       if(is.null(kmeans.Neighborhood)){
                         stop("calcBrainMask[MRIaggr] : \'kmeans.Neighborhood\' must be defined \n", 
                              "in order to choose the right kmeans partition \n")                
                       }
                       
                       for(iter_group in 1:length(kmeans.n_groups)){
                         
                         if(verbose){cat(iter_group, " brain groups : ", sep = "")}
                         res$kmeans[[iter_group]] <- stats::kmeans(carto[,param], centers = kmeans.n_groups[iter_group])
                         
                         conversion <- rank(res$kmeans[[iter_group]]$center)
                         res$kmeans[[iter_group]]$cluster <- conversion[res$kmeans[[iter_group]]$cluster]
                         res$kmeans[[iter_group]]$size <- table(res$kmeans[[iter_group]]$cluster)
                         res$kmeans[[iter_group]]$center <- tapply(carto[,param], res$kmeans[[iter_group]]$cluster, mean)
                         
                         test.combin <- lapply(strsplit(names_combin, split = ".", fixed = TRUE), function(x){x[1] == paste("G", kmeans.n_groups[iter_group], sep = "")})
                         combin_group <- strsplit(names_combin, split = ".", fixed = TRUE)[which(unlist(test.combin))]
                         for(iter_combin in 1:length(combin_group)){
                           combin_tempo <- combin_group[[iter_combin]]
                           
                           if(combin_tempo[2] == "="){
                             kmeans_logical <- res$kmeans[[iter_group]]$cluster == as.numeric(combin_tempo[3])
                           }
                           if(combin_tempo[2] == ">="){
                             kmeans_logical <- res$kmeans[[iter_group]]$cluster >= as.numeric(combin_tempo[3])
                           }
                           if(combin_tempo[2] == "<="){
                             kmeans_logical <- as.logical((res$kmeans[[iter_group]]$cluster != 1) * as.numeric(res$kmeans[[iter_group]]$cluster <= as.numeric(combin_tempo[3])))
                           }
                           
                           if(verbose){cat(paste(combin_tempo, collapse = "."), " ", sep = "")}
                           
                           res$potential[paste(combin_tempo, collapse = "."), "nb_groups"] <- kmeans.n_groups[iter_group]
                           
                           V <- calcFilter(dt2array(contrast = as.logical(kmeans_logical), coords = coords)$contrast[[1]], 
                                           filter = paste("3D_I", kmeans.Neighborhood, sep = ""), norm.filter = FALSE)$res
                           
                           res$potential[paste(combin_tempo, collapse = "."), "V"] <- mean(V[kmeans_logical == 1])
                           
                           
                           #### retain the best
                           if(res$potential[paste(combin_tempo, collapse = "."), "V"] > res$best_V){
                             res$best_V <- res$potential[paste(combin_tempo, collapse = "."), "V"]
                             res$best_group <- kmeans_logical
                             res$mask_name <- paste(param, "kmeans", paste(combin_tempo, collapse = "."), sep = ".")
                           }
                           
                         }
                         
                         if(verbose){cat("\n")}
                       }
                       
                     }
                     
                     if(!is.null(skull.param)){
                       
                       if(verbose){cat("skull groups : ")}
                       
                       initParameter(object = object, param = skull.param, checkArguments = TRUE, init = FALSE, 
                                     accept.coords = FALSE, accept.index = FALSE, accept.mask = FALSE, method = "calcBrainMask", 
                                     arg_name = "skull.param", long_name = "skull parameter")
                       
                       carto <- selectContrast(object, param = skull.param, format = "data.frame")
                       res$potential_skull <- data.frame(matrix(NA, ncol = 2, 
                                                                nrow = length(skull.n_groups) + 1))
                       names(res$potential_skull) <- c("nb_groups", "V") 
                       rownames(res$potential_skull) <- c("row", skull.n_groups)
                       res$potential_skull["row",] <- c(NA, res$best_V)
                       index_Skull <- NULL
                       
                       # identification de la boite cranienne
                       for(iter_kmeans in 1:length(skull.n_groups)){
                         if(verbose){cat(skull.n_groups[iter_kmeans], " ", sep = "")}  
                         # kmeans
                         res_kmeans <- stats::kmeans(carto[,skull.param], centers = skull.n_groups[iter_kmeans])
                         
                         group <- (1:skull.n_groups[iter_kmeans])[-which.max(res_kmeans$size)]
                         recouvrement <- sapply(group, 
                                                function(x){mean(res$best_group[res_kmeans$cluster == x])}
                         )              
                         index_Skull_test <- which(res_kmeans$cluster == group[which.min(recouvrement)])
                         kmeans_logical <- res$best_group
                         kmeans_logical[index_Skull_test] <- FALSE
                         
                         # evaluation
                         V <- calcFilter(dt2array(contrast = kmeans_logical, coords = coords)$contrast[[1]], 
                                         filter = paste("3D_I", kmeans.Neighborhood, sep = ""), norm.filter = FALSE)$res
                         
                         res$potential_skull[iter_kmeans + 1,] <- c(skull.n_groups[iter_kmeans], mean(V[kmeans_logical == 1]))
                         
                         if(which.max(res$potential_skull$V) == (iter_kmeans + 1)){
                           index_Skull <- index_Skull_test
                         }
                       }
                       if(verbose){cat("\n")}
                       
                       if(verbose){cat("intersect mask : ", sum(res$best_group[index_Skull]), "\n", sep = "")}
                       
                       # boite cranienne retiree du masque
                       if(length(index_Skull > 0)){
                         res$best_group[index_Skull] <- FALSE
                       }
                     }
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                   }
)

####>>> calcContralateral ####

methods::setMethod(f  = "calcContralateral", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, num = NULL, type = "mean", param.ref = NULL, distband = 1, lambda = 1, 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     if(is.null(num)){num <- 1:object@fieldDim$k}
                     
                     if(any(c("i_hemisphere", "j_hemisphere", "hemisphere") %in% selectParameter(object) == FALSE)){
                       stop("calcContralateral[MRIaggr] : missing elements in \'object@contrast\' \n", 
                            "missing parameters : \"i_hemisphere\" \"j_hemisphere\" \"hemisphere\" \n", 
                            "use calcHemisphere method with update.objet = TRUE to compute and allocate these elements \n")
                     }
                     
                     validCharacter(value = type, validLength = 1, validValues = c("mean", "median", "1NN_penalised"), method = "calcContralateral[MRIaggr]")
                     
                     if(type == "1NN_penalised" && is.null(param.ref)){
                       stop("calcContralateral[MRIaggr] : argument \'param.ref\' must be specified if \'type\' is \"1NN_penalised\" \n")
                     }
                     
                     if(type == "1NN_penalised"){
                       sd <- stats::sd(selectContrast(object, param = param.ref, num = num, format = "vector"), na.rm = TRUE) 
                     }else{
                       sd <- 0
                       param.ref <- param[1]
                     }
                     
                     if(is.null(param))
                     {param <- selectParameter(object)}
                     
                     data <- selectContrast(object, param = unique(c("index", param.ref, param)), num = num, coords = TRUE)
                     n.px <- nrow(data)
                     p <- length(param)
                     
                     data_miroir <- data.frame(matrix(NA, nrow = n.px, ncol = length(param) + 5))
                     names(data_miroir) <- c("index", "i_hemisphere", "j_hemisphere", "k", "hemisphere", paste(param, "_contro", sep = ""))
                     data_miroir$index <- data$index
                     data_miroir$k <- data$k
                     
                     # changement de repere
                     data_miroir[,c("i_hemisphere", "j_hemisphere", "hemisphere")] <- selectContrast(object, param = c("i_hemisphere", "j_hemisphere", "hemisphere"))
                     
                     ####### reperage des slices - hemi-left
                     index.plot_lesionL <- numeric(0)
                     index.plot_controL <- numeric(0)
                     
                     if(type == "mean"){type_moy <- TRUE}else{type_moy <- FALSE}
                     if(type == "median"){type_med <- TRUE}else{type_med <- FALSE}
                     if(type == "1NN_penalised"){type_NN <- TRUE}else{type_NN <- FALSE}
                     
                     if(verbose){cat(paste("Left slice : ", sep = ""))}
                     for(iter_slice in num)
                     { 
                       if(verbose){cat(paste(iter_slice, " ", sep = ""))}
                       # reperage des px de la slice
                       index_k <- which(data_miroir$k == iter_slice)
                       index_k_lesion <- stats::na.omit(which(data_miroir[index_k, "hemisphere"] == "left"))
                       index_k_contro <- stats::na.omit(which(data_miroir[index_k, "hemisphere"] == "right"))
                       
                       if(length(index_k_lesion) == 0 || length(index_k_contro) == 0){next}                
                       param_C <- unique(c(param, param.ref))
                       
                       res <- calcContro_cpp(contrast = as.matrix(data[index_k,param_C, drop = FALSE]), 
                                             coords_px = as.matrix(data_miroir[index_k, c("i_hemisphere", "j_hemisphere")]), 
                                             index_k = index_k_lesion - 1, index_k_contro = index_k_contro - 1, 
                                             d_lim = distband, lambda = lambda, param_ref = which(names(data[,param_C, drop = FALSE]) == param.ref) - 1, var_ref = sd, 
                                             type_moy = type_moy, type_med = type_med, type_NN = type_NN, verbose = verbose)
                       data_miroir[index_k[index_k_lesion], paste(param, "_contro", sep = "")] <- res$valeur_NormContro[,param_C %in% param]
                       index.plot_lesionL <- c(index.plot_lesionL, index_k[index_k_lesion[which(res$index_plot_k)]])
                       index.plot_controL <- c(index.plot_controL, index_k[index_k_contro[which(res$index_plot_k_contro)]])                
                     }
                     
                     if(verbose){cat("\n")}
                     
                     
                     ######## reperage des slices - lesion right
                     index.plot_lesionR <- numeric(0)
                     index.plot_controR <- numeric(0)
                     
                     if(type == "mean"){type_moy <- TRUE}else{type_moy <- FALSE}
                     if(type == "median"){type_med <- TRUE}else{type_med <- FALSE}
                     if(type == "1NN_penalised"){type_NN <- TRUE}else{type_NN <- FALSE}
                     
                     if(verbose){cat(paste("Right slice : ", sep = ""))}
                     for(iter_slice in num)
                     { 
                       if(verbose){cat(paste(iter_slice, " ", sep = ""))}
                       # reperage des px de la slice
                       index_k <- which(data_miroir$k == iter_slice)
                       index_k_lesion <- stats::na.omit(which(data_miroir[index_k, "hemisphere"] == "right"))
                       index_k_contro <- stats::na.omit(which(data_miroir[index_k, "hemisphere"] == "left"))
                       
                       if(length(index_k_lesion) == 0 || length(index_k_contro) == 0){next}
                       
                       ##### C
                       param_C <- unique(c(param, param.ref))
                       
                       res <- calcContro_cpp(contrast = as.matrix(data[index_k,param_C, drop = FALSE]), 
                                             coords_px = as.matrix(data_miroir[index_k, c("i_hemisphere", "j_hemisphere")]), 
                                             index_k = index_k_lesion - 1, index_k_contro = index_k_contro - 1, 
                                             d_lim = distband, lambda = lambda, param_ref = which(names(data[,param_C, drop = FALSE]) == param.ref) - 1, var_ref = sd, 
                                             type_moy = type_moy, type_med = type_med, type_NN = type_NN, verbose = verbose)
                       
                       data_miroir[index_k[index_k_lesion], paste(param, "_contro", sep = "")] <- res$valeur_NormContro[,param_C %in% param]
                       index.plot_lesionR <- c(index.plot_lesionR, index_k[index_k_lesion[which(res$index_plot_k)]])
                       index.plot_controR <- c(index.plot_controR, index_k[index_k_contro[which(res$index_plot_k_contro)]])                
                     }
                     if(verbose){cat("\n")}
                     
                     
                     return(list(data = cbind(data[,c("i", "j")], data_miroir), 
                                 index_plot = list(index.plot_lesionR = index.plot_lesionR, 
                                                   index.plot_lesionL = index.plot_lesionL, 
                                                   index.plot_controR = index.plot_controR, 
                                                   index.plot_controL = index.plot_controL), 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite
                     )
                     )
                   }
)

####>>> calcDistMask ####

methods::setMethod(f  = "calcDistMask", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask, name_newparam = paste("dist", mask, sep = "_"), 
                                         spatial_res = c(1, 1, 1), numeric2logical = FALSE, Neighborhood = "3D_N10", 
                                         verbose =  optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     #### initialisation ####
                     # initPackage("RANN", method = "calcDistMask[MRIaggr]")
                     
                     initParameter(object = object, param = mask, checkArguments = TRUE, init = FALSE, 
                                   accept.coords = FALSE, accept.index = FALSE, method = "calcDistMask", 
                                   arg_name = "mask", long_name = "mask")
                     
                     p <- length(mask)
                     
                     validCharacter(value = name_newparam, validLength = p, method = "calcDistMask[MRIaggr]")
                     validNumeric(value = spatial_res, validLength = 3, min = 0, method = "calcDistMask[MRIaggr]")
                     
                     data <- selectContrast(object, param = mask, coords = TRUE)
                     data$i_scaled <- data$i * spatial_res[1]
                     data$j_scaled <- data$j * spatial_res[2]
                     data$k_scaled<- data$k * spatial_res[3]
                     
                     res <- data.frame(matrix(NA, nrow = selectN(object), ncol = p))
                     names(res) <- name_newparam
                     
                     #### loop ####
                     
                     for(iter_mask in 1:p){
                       if(verbose){cat(iter_mask, " ", sep = "")}
                       
                       # test logical 
                       data[,mask[iter_mask]] <- initMask(object, mask[iter_mask], checkArguments = TRUE, init = numeric2logical, 
                                                          arg_name = "mask", long_name = "mask", method = "calcDistMask", format = "vector")
                       
                       index_mask <- which(data[,mask[iter_mask]] == TRUE)
                       
                       if(length(index_mask) > 0){
                         
                         mask_outline <- pointsOutline(data[index_mask, c("i", "j", "k")], filter = Neighborhood)
                         mask_outline$i <- mask_outline$i * spatial_res[1]
                         mask_outline$j <- mask_outline$j * spatial_res[2]
                         mask_outline$k <- mask_outline$k * spatial_res[3]
                         
                         res[index_mask, iter_mask] <- 0
                         res[-index_mask, iter_mask] <- RANN::nn2(data = mask_outline, 
                                                                  query = data[-index_mask, c("i_scaled", "j_scaled", "k_scaled")], 
                                                                  k = 1)$nn.dists              
                       }
                       
                     }
                     if(verbose){cat("\n")}
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }
)

####>>> calcDistTissues ####


methods::setMethod(f = "calcDistTissues", 
                   signature = "MRIaggr", 
                   definition = function(object, param, param.membership, num = NULL, hemisphere = "both"){
                     
                     ### import
                     data_membership <- selectContrast(object, num = num, param = param.membership, hemisphere = hemisphere, format = "data.frame")
                     data_membership <- apply(data_membership, 2, as.numeric)
                     p <- length(param)
                     n.membership <- length(param.membership)
                     
                     test_range <- apply(data_membership, 2, range)
                     if(any(test_range < 0) || any(test_range > 1)){
                       stop("calcDistTissues[MRIaggr] : wrong specification of \'param.membership\' \n", 
                            "param.membership values must be in [0;1] \n", 
                            "proposed parameters range in : ", paste(range(test_range), collapse = " "), "\n", 
                            "problematic parameters : ", paste(param.membership[colSums( (test_range < 0) + (test_range > 1) ) > 0], collapse = " "), "\n")
                     }
                     
                     data_param <-  selectContrast(object, num = num, param = param, hemisphere = hemisphere, format = "data.frame")
                     
                     ### calcul des moments
                     moments <- data.frame(matrix(NA, nrow = n.membership, ncol = 4*p + 1))
                     names(moments) <- c(paste(param, "mu", sep = "_"), paste(param, "sigma", sep = "_"), paste(param, "skewness", sep = "_"), paste(param, "kurtosis", sep = "_"), "nb.vx")
                     rownames(moments) <- param.membership
                     
                     for(iter_membership in 1:n.membership){
                       moments[param.membership[iter_membership], "nb.vx"] <- sum(data_membership[,param.membership[iter_membership]])
                       
                       for(iter_p in 1:p){
                         
                         moments[param.membership[iter_membership], paste(param[iter_p], "mu", sep = "_")] <- stats::weighted.mean(data_param[,param[iter_p]], w = data_membership[,param.membership[iter_membership]], na.rm = TRUE)  
                         diff_mu <- data_param[,param[iter_p]]-moments[param.membership[iter_membership], paste(param[iter_p], "mu", sep = "_")]
                         
                         moments[param.membership[iter_membership], paste(param[iter_p], "sigma", sep = "_")] <- sqrt(stats::weighted.mean(diff_mu^2, w = data_membership[,param.membership[iter_membership]], na.rm = TRUE))  
                         diff_sigma <- diff_mu / moments[param.membership[iter_membership], paste(param[iter_p], "sigma", sep = "_")]
                         
                         moments[param.membership[iter_membership], paste(param[iter_p], "skewness", sep = "_")] <- stats::weighted.mean(diff_sigma^3, w = data_membership[,param.membership[iter_membership]], na.rm = TRUE)
                         
                         moments[param.membership[iter_membership], paste(param[iter_p], "kurtosis", sep = "_")] <- stats::weighted.mean(diff_sigma^4, w = data_membership[,param.membership[iter_membership]], na.rm = TRUE)
                         
                       }
                       
                     }
                     
                     return(moments)   
                   }
)

####>>> calcFilter ####

methods::setMethod(f  = "calcFilter", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, filter, norm.filter = TRUE, bilateral = FALSE, na.rm = FALSE, name_newparam = NULL, 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
                   { 
                     # preparation
                     initParameter(object = object, param = param, checkArguments = TRUE, init = FALSE, 
                                   accept.coords = FALSE, accept.index = FALSE, method = "calcFilter")
                     p <- length(param)
                     
                     if(is.null(name_newparam)){
                       if(is.character(filter)){
                         name_newparam <- paste(param, filter, sep = "_")
                       }else{
                         name_newparam <- paste(param, "filtered", sep = "_")
                       }
                     }
                     
                     validDimension(value1 = param, value2 = name_newparam, name1 = "param", name2 = "name_newparam", type = "length", method = "calcFilter[MRIaggr]")
                     
                     carto <- selectContrast(object, param = param, coords = TRUE, format = "data.frame")
                     carto <- dt2array(contrast = carto[,param], 
                                       range.coords = object@fieldDim, 
                                       coords = carto[,c("i", "j", "k")])
                     
                     res <- data.frame(matrix(NA, ncol = p + 3, nrow = selectN(object)))
                     names(res) <- c("i", "j", "k", name_newparam)
                     res[,c("i", "j", "k")] <- selectCoords(object)
                     indexData <- res$i + object@fieldDim$i * (res$j - 1) + object@fieldDim$i * object@fieldDim$j * (res$k - 1)
                     
                     for(iter_param in 1:p){
                       if(verbose){cat(param[iter_param], " ", sep = "")}
                       tempo <- calcFilter(carto$contrast[[iter_param]], filter = filter, norm.filter = norm.filter, 
                                           bilateral = bilateral, na.rm = na.rm)
                       Mfilter <- tempo$filter
                       index_NNA <- which(!is.na(tempo$res))
                       tempo <- array2dt(array = tempo$res, 
                                         name_newparam = name_newparam[iter_param], names_coords = c("i", "j", "k"))
                       
                       res[which(indexData %in% index_NNA), name_newparam[iter_param]] <- tempo[,name_newparam[iter_param]] 
                     }
                     if(verbose){cat("\n")}
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 filter = Mfilter, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }
)

####>>> calcGroupsMask ####

methods::setMethod(f  = "calcGroupsMask", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask, numeric2logical = FALSE, 
                                         W = "ifany", W.range, W.spatial_res = c(1, 1, 1), 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = TRUE){
                     
                     #### preliminaries 
                     # test parameter
                     initParameter(object = object, param = mask, checkArguments = TRUE, init = FALSE, 
                                   accept.coords = FALSE, accept.index = FALSE, method = "calcGroupsMask", 
                                   arg_name = "mask", long_name = "parameters")
                     p <- length(mask)
                     
                     # coords
                     coords <- selectCoords(object)
                     validNumeric(value = W.spatial_res, validLength = 3, min = 0, method = "calcGroupsMask[MRIaggr]")
                     
                     # data lesion
                     carto <- selectContrast(object, param = mask, format = "data.frame")
                     
                     
                     carto[,mask] <- initMask(object, mask, checkArguments = TRUE, init = numeric2logical, 
                                              arg_name = "mask", long_name = "mask", method = "calcGroupsMask", format = "matrix")
                     
                     if(identical(W, "ifany") && any(!is.na(object@W$Wmatrix))){
                       W <- selectW(object)
                     }else if(identical(W, "ifany")){
                       W <- NULL
                     } 
                     
                     
                     if(!is.null(W)){
                       
                       if("dgCMatrix" %in% class(W) == FALSE){
                         stop("calcGroupsMask[MRIaggr] : wrong specification of \'W\' \n", 
                              "\'W\' must be of class dgCMatrix \n", 
                              "proposed class of \'W\' : ", paste(class(W), collapse = " "), "\n")
                       }
                       
                       validDimension(value1 = W, validDimension = rep(selectN(object),2), name1 = "W", name2 = NULL, type = c("nrow", "ncol"), method = "calcGroupsMask[MRIaggr]")
                       
                     }
                     
                     #### identification des groupes
                     if(verbose){
                       ncharMax <- max(nchar(mask), 15)
                       cat("mask parameter", rep(" ", ncharMax - 14), " : number of observations per spatial group (total) \n", sep = "")
                     }
                     
                     res <- list()
                     
                     for(iter_param in 1:p){
                       if(verbose){cat(mask[iter_param], rep(" ", ncharMax-nchar(mask[iter_param])), " : ", sep = "")}
                       
                       param_tempo <- mask[iter_param]
                       index_N <- which(carto[,param_tempo] == TRUE)                         
                       
                       if(length(index_N) > 1){
                         if(is.null(W)){
                           W_lesion <- calcW(coords[index_N,], method = "euclidean", range = W.range, 
                                             upper = NULL, format = "dgCMatrix", row.norm = FALSE, spatial_res = W.spatial_res)$W
                         }else{              
                           W_lesion <- W[index_N, index_N]
                         }
                         
                         res[[iter_param]] <-  calcGroupsW(W_lesion, max_groups = 10000, verbose = FALSE)               
                         carto[index_N, iter_param] <- res[[iter_param]]$group
                         res[[iter_param]]$group <- carto[,iter_param]
                         
                       }else{
                         res[[iter_param]] <- list()
                         res[[iter_param]]$group <- rep(0, selectN(object))
                         res[[iter_param]]$group_size <- if(length(iter_param) == 1){1}else{NA}
                         res[[iter_param]]$group_number <- length(iter_param)
                         res[[iter_param]]$group_max <- if(length(iter_param) == 1){1}else{NA}
                         
                         if(length(index_N) == 1){
                           res[[iter_param]]$group[index_N] <- 1  
                         }
                       }
                       
                       if(verbose){
                         text_group <- paste(res[[iter_param]]$group_size, collapse = " ")
                         cat(text_group, rep(" ", max(0, 40 - nchar(text_group))), " (", length(index_N), ") \n", sep = "")
                       } 
                       
                     }
                     names(res) <- mask
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite)
                     )
                   }
)

####>>> calcHemisphere ####

methods::setMethod(f  = "calcHemisphere", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, num = NULL, p = 1, subset = NULL, penalty = "symmetry", mask = NULL, numeric2logical = FALSE, n.points = 100, 
                                         gridSearch = TRUE, i_test = seq(-20, 20, by = 5), angle.test = seq(-30, 30, by = 5), unit_angle = "degree", 
                                         NelderMead = TRUE, maxit = 100, reltol = 0.001, 
                                         plot = TRUE, filename = paste(object@identifier, "_calcHemisphere", sep = ""), 
                                         update.object = FALSE, overwrite = FALSE, ...)
                   { 
                     if(plot == TRUE){ #### graphical options
                       ## get graphical arguments
                       optionsMRIaggr.eg <- optionsMRIaggr(c("height", "path", "res", "unit", "width", "window"))					 
                       dots.arguments <- list(...)
                       names_dots.arguments <- names(dots.arguments)
                       
                       validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                      validValues = c("height", "path", "res", "unit", "verbose", "width", "window"), 
                                      refuse.NULL = FALSE, method = "calcHemisphere[MRIaggr]")
                       
                       ## set specific display
                       if(length(names_dots.arguments) > 0){
                         optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                       }
                       
                       ## create all variables
                       height <- optionsMRIaggr.eg$height
                       path <- optionsMRIaggr.eg$path
                       unit <- optionsMRIaggr.eg$unit
                       res <- optionsMRIaggr.eg$res
                       width <- optionsMRIaggr.eg$width
                       window <- optionsMRIaggr.eg$window
                     }
                     
                     verbose <- optionsMRIaggr("verbose")
                     
                     #### preliminary tests ####
                     if(gridSearch == FALSE && NelderMead == FALSE){
                       stop("calcHemisphere[MRIaggr] : arguments gridSearch and NelderMead should not be simultaneously FALSE \n")
                     }
                     
                     validNumeric(value = p, validLength = 1, min = 0, method = "calcHemisphere[MRIaggr]")
                     
                     if(length(param) != 1){
                       stop("calcHemisphere[MRIaggr] : \'param\' must have lenght 1 \n", 
                            "proposed \'param\' : ", paste(param, collapse = " "), "\n")
                     }
                     
                     validCharacter(value = penalty, validLength = 1, validValues = c("symmetry", "asymmetry"), method = "calcHemisphere[MRIaggr]")
                     validCharacter(value = unit_angle, validLength = 1, validValues = c("radian", "degree"), method = "calcHemisphere[MRIaggr]")
                     
                     if(!is.null(mask)){
                       initParameter(object = object, param = mask, checkArguments = TRUE, init = FALSE, 
                                     accept.coords = FALSE, accept.index = FALSE, method = "calcHemisphere", 
                                     arg_name = "mask", long_name = "mask")
                       
                       mask <- selectContrast(object, param = mask, format = "matrix")
                       if(numeric2logical == TRUE){mask <- apply(mask, 2, as.logical)}
                       
                       if(!is.logical(mask)){
                         stop("calcHemisphere[MRIaggr] : type of \'mask\' is not logical \n", 
                              "proposed type : ", paste(is(mask), collapse = " "), "\n", 
                              "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
                       }
                       mask <- rowSums(mask) > 0
                     }
                     
                     if(plot == TRUE){
                       scale <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                           n.plot = 1, mfrow = c(1, 1), xlim = NULL, ylim = NULL, method = "calcHemisphere[MRIaggr]")$scale
                     }
                     
                     #### initialisation ####
                     
                     #### data
                     num <- initNum(object, num = num, method = "calcHemisphere")
                     
                     data <- selectContrast(object, param = param, num = num, subset = subset, coord = TRUE)
                     data <- data[is.na(data[,param]) == FALSE,]
                     n <- nrow(data)
                     n.num <- length(num)
                     ls.indexK <- lapply(num, function(x){which(data$k == x) - 1})
                     
                     sd_data <- stats::sd(data[,param], na.rm = TRUE)
                     sdp_data <- sd_data^p
                     
                     if(n < 2){
                       stop("calcHemisphere[MRIaggr] : there is not enougth data to distinguish 2 hemispheres \n", 
                            "nrow(data) : ", n, "\n")
                     }
                     
                     #### seed
                     deg2rad <- 2 * pi / 360
                     rad2deg <- 360 / 2 * pi
                     if(penalty == "symmetry"){penaltyNA <- 3}else{penaltyNA <- 1} # penalised data with no controlateral correspondant
                     
                     i_median <- stats::median(unique(data$i))
                     j_median <- stats::median(unique(data$j))
                     angle_median <- 0
                     
                     df.optimum <- data.frame(position_i = i_median, 
                                              position_j = j_median, 
                                              angle_rad = angle_median, 
                                              asymetrie = NA)  
                     
                     #### grid search ####
                     if(gridSearch == TRUE){
                       
                       # initialisation
                       if(unit_angle == "degree"){
                         angle.test <- deg2rad * angle.test
                       }
                       
                       grid <- expand.grid(i = i_median + i_test, j = df.optimum$position_j, angle = angle_median + angle.test)
                       grid$rank_i <- as.numeric(as.factor(rank(grid$i - i_median)))
                       grid$rank_angle <- as.numeric(as.factor(rank(grid$angle)))
                       grid$rank <- grid$rank_i + grid$rank_angle - 2
                       
                       n.i_test <- length(i_test)
                       n.angle.test <- length(angle.test)
                       n.grid <- nrow(grid)
                       
                       grid$penalty <- 0
                       grid$nb <- 0
                       grid$moy <- 0
                       
                       if(penalty == "symmetry"){optimum <- -Inf}else{optimum <- + Inf}  
                       order_optimum <- + Inf
                       
                       if(verbose == TRUE){
                         index_trace <- round(c(1, seq(1, n.grid, length.out = 10)[c(-1,-10)], n.grid))
                       }
                       
                       #### computation
                       if(verbose == TRUE){
                         cat("Grid Search : ", n.grid, " parametrisations \n", 
                             "i     : ", paste(i_test, collapse = " "), " (in voxels) \n", 
                             # "j     : ", paste(df.optimum$position_j, collapse = " "), "\n", 
                             "angle : ", paste(round(360 * angle.test / (2 * pi), 1), collapse = " "), " (in degrees) \n", sep = "")
                       }
                       
                       for(iter_grid in 1:n.grid){ 
                         
                         if(verbose == TRUE && iter_grid %in% index_trace){cat("*")}
                         
                         res_cpp <- calcHemi_cpp(coordsI = data$i, coordsJ = data$j, ls_indexK = ls.indexK, n_num = n.num, value = data[,param], n = n, 
                                                 i_pos = grid[iter_grid, "i"], j_pos = grid[iter_grid, "j"], angle_pos = grid[iter_grid, "angle"], 
                                                 penaltyNA = penaltyNA, sd_data = sd_data, p = p, symetrie = (penalty == "symmetry"))
                         #                 cat(res_cpp$numberAssociated, " ", res_cpp$pcNA, " ", res_cpp$criteria / res_cpp$numberAssociated, " ", res_cpp$compromise, "\n")
                         grid[iter_grid, "penalty"] <- res_cpp$criteria
                         grid[iter_grid, "nb"] <- res_cpp$numberAssociated
                         grid[iter_grid, "moy"] <- res_cpp$compromise
                         
                         if(penalty == "symmetry"){
                           if(res_cpp$compromise > optimum || (res_cpp$compromise == optimum) && (grid$rank[iter_grid] < order_optimum) ){
                             optimum <- res_cpp$compromise
                             order_optimum <- grid$rank[iter_grid]
                           }
                         }else{
                           if(res_cpp$compromise < optimum || (res_cpp$compromise == optimum) && (grid$rank[iter_grid] < order_optimum) ){
                             optimum <- res_cpp$compromise
                             order_optimum <- grid$rank[iter_grid]
                           }
                         }
                       }    
                       
                       index_optimum <- which.min(abs(grid$moy-optimum))
                       
                       df.optimum$position_i <- grid[index_optimum, "i"]
                       df.optimum$angle_rad <- grid[index_optimum, "angle"]
                       df.optimum$asymetrie <- optimum    
                       
                       test.bordI <- df.optimum$position_i %in% c(i_test[1], utils::tail(i_test, 1))
                       test.bordAngle <- df.optimum$angle_rad %in% c(angle.test[1], utils::tail(angle.test, 1))
                       
                       if(test.bordI || test.bordAngle){
                         gridSearch.cv <- FALSE
                         if(verbose){cat("\n")}
                       }else{
                         gridSearch.cv <- TRUE
                         if(verbose){cat(" cv \n")}
                       }
                       
                     }else{
                       gridSearch.cv <- NA
                     }
                     
                     #### optimization ####
                     
                     #optim.control <- list(fnscale = -1, maxit = maxit, reltol = reltol, verbose = if(verbose > 1){verbose}else{0})
                     
                     optim.control <- list(fnscale = if(penalty == "symmetry"){-1}else{1}, 
                                           maxit = maxit, reltol = reltol, 
                                           trace = if(verbose > 1){verbose}else{0})
                     
                     
                     if(NelderMead == TRUE){
                       
                       if(verbose == TRUE){cat("Nelder-Mead with optim \n")}
                       dim_carto <- selectFieldDim(object)
                       
                       res <- stats::optim(par = list(df.optimum$position_i, df.optimum$angle_rad), 
                                           fn = function(value){
                                             
                                             if(abs(value[2]) > pi){if(penalty == "symmetry"){return(-Inf)}else{return(Inf)}}
                                             
                                             calcHemi_cpp(coordsI = data$i, coordsJ = data$j, ls_indexK = ls.indexK, n_num = n.num, value = data[,param], n = n, 
                                                          i_pos = value[1], j_pos = df.optimum$position_j, angle_pos = value[2], 
                                                          penaltyNA = penaltyNA, sd_data = sd_data, p = p, symetrie = (penalty == "symmetry"))$compromise                             
                                           }, 
                                           method = "Nelder-Mead", control = optim.control
                       )
                       
                       df.optimum$position_i <- res$par[1]
                       df.optimum$angle_rad <- res$par[2]
                       df.optimum$asymetrie <- res$value         
                     }
                     
                     
                     #### midplane ####
                     j <- seq(min(data$j), max(data$j), length.out = n.points)
                     midplane <-  data.frame(i = df.optimum$position_i + sin(df.optimum$angle_rad) * (j-df.optimum$position_j), j = j )
                     
                     #### hemisphere caracterization ####      
                     data_miroir <- data.frame(matrix(NA, nrow = selectN(object), ncol = 8))
                     names(data_miroir) <- c("index", "i", "j", "k", "i_hemisphere", "j_hemisphere", "hemisphere")
                     data_miroir[,c("i", "j", "k")] <- selectCoords(object)
                     data_miroir$index <- selectContrast(object, "index", format = "vector")
                     if(!is.null(mask)){data_miroir$mask <- mask}
                     
                     # changement de repere
                     data_miroir$i_hemisphere <- cos(df.optimum$angle_rad) * (data_miroir$i-df.optimum$position_i) - sin(df.optimum$angle_rad) * (data_miroir$j-df.optimum$position_j)
                     data_miroir$j_hemisphere <- sin(df.optimum$angle_rad) * (data_miroir$i-df.optimum$position_i) + cos(df.optimum$angle_rad) * (data_miroir$j-df.optimum$position_j)
                     
                     data_miroir$hemisphere <- "undefined"
                     data_miroir$hemisphere[data_miroir$i_hemisphere > 0] <- "left"
                     data_miroir$hemisphere[data_miroir$i_hemisphere < 0] <- "right"
                     
                     hemispheres <- data.frame(left = "defined", right = "defined", stringsAsFactors = FALSE)
                     if(!is.null(mask)){              
                       
                       index_lesion <- which(data_miroir$mask == TRUE)
                       table_lesion <- table(data_miroir$hemisphere[index_lesion])
                       
                       if("left" %in% names(table_lesion)){
                         hemispheres$left <- "lesion"
                       }else{
                         hemispheres$left <- "contralateral"
                       }
                       
                       if("right" %in% names(table_lesion)){
                         hemispheres$right <- "lesion"
                       }else{
                         hemispheres$right <- "contralateral"
                       }
                     }
                     
                     #### display ####
                     if(plot == TRUE && gridSearch == TRUE){
                       
                       initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                         mfrow = c(1, 2), bg = NULL, pty = NULL, mar = rep(3, 4), mgp = c(2, 0.5, 0))
                       plot.seq_moy <- seq(min(grid$moy), max(c(grid$moy, df.optimum$asymetrie)), length.out = 5)
                       plot.seq_i <- unique(grid$i)
                       plot.seq_angle <- unique(grid$angle)
                       
                       # optimum selon i
                       palette <- grDevices::rainbow(n.angle.test)
                       graphics::plot(NA, 
                                      xlim = range(c(plot.seq_i, df.optimum$position_i)), xlab = "i", 
                                      ylim = c(plot.seq_moy[1], plot.seq_moy[5]), ylab = penalty, 
                                      axes = FALSE, main = "color = angles")
                       graphics::box()
                       
                       graphics::axis(1, at = plot.seq_i, labels = round(plot.seq_i, optionsMRIaggr("digit.result")))
                       graphics::axis(2, at = plot.seq_moy, labels = round(plot.seq_moy, optionsMRIaggr("digit.result")))
                       
                       for(iter_angle in unique(grid$rank_angle)){
                         index_angle <- which(grid$rank_angle == iter_angle)
                         col <- palette[iter_angle]
                         graphics::points(grid[index_angle, c("i", "moy")], lty = 1, col = col, type = "l")
                       }
                       graphics::points(df.optimum[,c("position_i", "asymetrie")], col = "black", cex = 2, pch = 15)
                       
                       # optimum selon angle
                       palette <- grDevices::rainbow(n.i_test)
                       graphics::plot(NA, 
                                      xlim = range(c(plot.seq_angle, df.optimum$angle_rad)), xlab = "angle", 
                                      ylim = c(plot.seq_moy[1], plot.seq_moy[5]), ylab = penalty, 
                                      axes = FALSE, main = "color = i")
                       graphics::box()
                       graphics::axis(1, at = plot.seq_angle, labels = round(plot.seq_angle, optionsMRIaggr("digit.result")))
                       graphics::axis(2, at = plot.seq_moy, labels = round(plot.seq_moy, optionsMRIaggr("digit.result")))
                       
                       for(iter_i in unique(grid$rank_i)){
                         index_i <- which(grid$rank_i == iter_i)
                         col <- palette[iter_i]                                
                         graphics::points(grid[index_i, c("angle", "moy")], lty = 1, col = col, type = "l")               
                       }
                       graphics::points(df.optimum[,c("angle_rad", "asymetrie")], col = "black", cex = 2, pch = 15)
                       
                       if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){grDevices::dev.off()}
                       
                     }
                     
                     #### export ####
                     res <- list()
                     
                     res$data <- data_miroir[,c("i_hemisphere", "j_hemisphere", "hemisphere")]
                     res$hemispheres <- hemispheres
                     res$midplane <- midplane
                     res$optimum <- df.optimum
                     res$grid <- grid
                     res$cv <- gridSearch.cv  
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }  
)

####>>> calcNormalization ####

methods::setMethod(f  = "calcNormalization", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, mu_type = "mean", sigma_type = "sd", rm.CSF = FALSE, rm.WM = FALSE, rm.GM = FALSE, 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
                   {
                     # controle preliminaire
                     if(rm.GM == TRUE && rm.WM == TRUE){
                       stop("calcNormalization[MRIaggr] : \'rm.GM\' and \'rm.WM\' cannot be simultaneously TRUE \n", 
                            "set at least one of them to FALSE \n")
                     }
                     
                     if(is.character(mu_type)){
                       
                       validCharacter(value = mu_type, validLength = 1, validValues = c("min", "1Q", "median", "3Q", "max", "mean"), method = "calcNormalization[MRIaggr]")
                       
                       mu_type <- switch(mu_type, 
                                         "min" = 0, 
                                         "1Q" = 0.25, 
                                         "median" = 0.5, 
                                         "3Q" = 0.75, 
                                         "max" = 1)
                     }
                     
                     if(is.numeric(mu_type)){
                       
                       if(mu_type > 1 || mu_type < 0){
                         stop("calcNormalization[MRIaggr] : wrong specification of \'mu_type\' \n", 
                              "\'mu_type\' must be between  0 and 1 \n", 
                              "proposed \'mu_type\' : ", mu_type, "\n")
                       }
                       
                       mu_norm <- function(x, w){return(stats::quantile(x[w > 0.5], probs = mu_type, na.rm = TRUE))}
                     }else{
                       mu_norm <- function(x, w){stats::weighted.mean(x, w = w, na.rm = TRUE)}
                     }
                     
                     validCharacter(value = sigma_type, validLength = 1, validValues = c("sd", "mad"),, method = "calcNormalization[MRIaggr]")
                     
                     if(sigma_type == "sd"){
                       sigma_norm <- function(x, w, center){stats::sd(x[w > 0.5], na.rm = TRUE)}
                     }else{
                       sigma_norm <- function(x, w, center){stats::mad(x[w > 0.5], center = center, na.rm = TRUE)}
                     }
                     
                     if("hemisphere" %in% selectParameter(object) == FALSE){
                       warning("calcNormalization[MRIaggr] : missing \"hemisphere\" parameter in \'object\' \n", 
                               "the computations will be incomplete \n", 
                               "use the calcHemisphere and calcContralateral function to compute and allocate it \n")
                       test.hemi <- FALSE
                     }else{
                       test.hemi <- TRUE
                     }
                     
                     # extraction des coordonnees
                     coords_both <- selectCoords(object, hemisphere = "both")
                     
                     # extraction des parametres
                     carto_both <- selectContrast(object, param = param, hemisphere = "both", format = "matrix")
                     p <- length(param)
                     
                     # extraction des hemispheres
                     if(test.hemi){
                       index_hemiL <- which(object@contrast$hemisphere == "left")
                       index_hemiR <- which(object@contrast$hemisphere == "right")
                     }
                     
                     # suppression des parties majoritairement CSF
                     if(rm.CSF == TRUE || rm.WM == TRUE || rm.GM == TRUE){
                       
                       if(any(c("WM", "GM", "CSF") %in% selectParameter(object) == FALSE)){
                         stop("calcNormalization[MRIaggr] : impossible to remove the CSF \n", 
                              c("WM", "GM", "CSF")[c("WM", "GM", "CSF") %in% selectParameter(object) == FALSE], " is not available in \'x\' \n", 
                              "use the calcTissueType function to obtain compute and allocate these parameters \n")}            
                       
                       param_tissue <- c("CSF", "WM", "GM")[c(rm.CSF == FALSE, rm.WM == FALSE, rm.GM == FALSE)]
                       w.tissue <- rowSums(selectContrast(object, param = param_tissue, hemisphere = "both", format = "matrix"))
                     }else{
                       w.tissue <- rep(1, selectN(object))
                     }
                     
                     #### global
                     if(verbose){cat("global \n")}
                     
                     norm_global <- data.frame(matrix(NA, ncol = length(param), nrow = 6))
                     names(norm_global) <- param
                     row.names(norm_global) <- c("mu_both", "mu_left", "mu_right", "sigma_both", "sigma_left", "sigma_right")
                     
                     norm_global["mu_both",] <- apply(carto_both, 2, function(x){mu_norm(x, w.tissue)})
                     if(test.hemi){
                       norm_global["mu_left",] <- apply(carto_both[index_hemiL,, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[index_hemiL])})
                       norm_global["mu_right",] <- apply(carto_both[index_hemiR,, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[index_hemiR])})
                     }
                     norm_global["sigma_both",] <- sapply(1:p, function(x){sigma_norm(x = carto_both[,x], w = w.tissue, center = norm_global["mu_both", x])})
                     if(test.hemi){
                       norm_global["sigma_left",] <- sapply(1:p, function(x){sigma_norm(x = carto_both[index_hemiL, x], w = w.tissue[index_hemiL], center = norm_global["mu_left", x])}) 
                       norm_global["sigma_right",] <- sapply(1:p, function(x){sigma_norm(x = carto_both[index_hemiR, x], w = w.tissue[index_hemiR], center = norm_global["mu_right", x])}) 
                     }
                     
                     #### par coupe
                     if(verbose){cat("slice : ")}
                     
                     normMu_slice_both <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_slice_both) <- param
                     normMu_slice_left <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_slice_left) <- param
                     normMu_slice_right <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_slice_right) <- param
                     normSigma_slice_both <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_slice_both) <- param
                     normSigma_slice_left <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_slice_left) <- param
                     normSigma_slice_right <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_slice_right) <- param
                     
                     normMu_3slices_both <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_3slices_both) <- param
                     normMu_3slices_left <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_3slices_left) <- param
                     normMu_3slices_right <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normMu_3slices_right) <- param
                     normSigma_3slices_both <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_3slices_both) <- param
                     normSigma_3slices_left <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_3slices_left) <- param
                     normSigma_3slices_right <- data.frame(matrix(NA, ncol = length(param), nrow = object@fieldDim$k))
                     names(normSigma_3slices_right) <- param
                     
                     
                     for(iter_k in 1:object@fieldDim$k){
                       if(verbose){cat(iter_k, " ", sep = "")}
                       
                       index_k <- which(coords_both$k == iter_k)
                       
                       normMu_slice_both[iter_k,] <- apply(carto_both[index_k,, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[index_k])})
                       normSigma_slice_both[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[index_k, x], w = w.tissue[index_k], center = normMu_slice_both[iter_k, x])})               
                       
                       if(test.hemi){
                         normMu_slice_left[iter_k,] <- apply(carto_both[intersect(index_k, index_hemiL),, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[intersect(index_k, index_hemiL)])})
                         normSigma_slice_left[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[intersect(index_k, index_hemiL), x], w = w.tissue[intersect(index_k, index_hemiL)], center = normMu_slice_left[iter_k, x])}) 
                         
                         normMu_slice_right[iter_k,] <- apply(carto_both[intersect(index_k, index_hemiR),, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[intersect(index_k, index_hemiR)])})
                         normSigma_slice_right[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[intersect(index_k, index_hemiR), x], w = w.tissue[intersect(index_k, index_hemiR)], center = normMu_slice_right[iter_k, x])}) 
                       }
                       
                       index_k <- which(coords_both$k %in% seq(iter_k - 1, iter_k + 1))
                       
                       normMu_3slices_both[iter_k,] <- apply(carto_both[index_k,, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[index_k])})
                       normSigma_3slices_both[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[index_k, x], w = w.tissue[index_k], center = normMu_3slices_both[iter_k, x])})               
                       
                       if(test.hemi){
                         normMu_3slices_left[iter_k,] <- apply(carto_both[intersect(index_k, index_hemiL),, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[intersect(index_k, index_hemiL)])})
                         normSigma_3slices_left[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[intersect(index_k, index_hemiL), x], w = w.tissue[intersect(index_k, index_hemiL)], center = normMu_3slices_left[iter_k, x])})
                         
                         normMu_3slices_right[iter_k,] <- apply(carto_both[intersect(index_k, index_hemiR),, drop = FALSE], 2, function(x){mu_norm(x, w.tissue[intersect(index_k, index_hemiR)])})
                         normSigma_3slices_right[iter_k,] <- sapply(1:p, function(x){sigma_norm(x = carto_both[intersect(index_k, index_hemiR), x], w = w.tissue[intersect(index_k, index_hemiR)], center = normMu_3slices_right[iter_k, x])}) 
                       }
                     }
                     cat("\n")
                     
                     res <- list(norm_global = norm_global, 
                                 normMu_slice_both = normMu_slice_both, 
                                 normSigma_slice_both = normSigma_slice_both, 
                                 normMu_slice_left = normMu_slice_left, 
                                 normSigma_slice_left = normSigma_slice_left, 
                                 normMu_slice_right = normMu_slice_right, 
                                 normSigma_slice_right = normSigma_slice_right, 
                                 normMu_3slices_both = normMu_3slices_both, 
                                 normSigma_3slices_both = normSigma_3slices_both, 
                                 normMu_3slices_left = normMu_3slices_left, 
                                 normSigma_3slices_left = normSigma_3slices_left, 
                                 normMu_3slices_right = normMu_3slices_right, 
                                 normSigma_3slices_right = normSigma_3slices_right)
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite)
                     )
                   }
)

####>>> calcSummaryRegion ####
methods::setMethod(f  = "calcSummaryRegion", 
                   signature  = "MRIaggr", 
                   definition = function(object, param = NULL, region, region.value = NULL, fct.summary,
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, na.rm = FALSE,
                                         long2wide = TRUE, sep = ".",
                                         format = "data.table",
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     #### tests
                     validCharacter(value = format, name = "format", validLength = 1, validValues = c("data.table","data.frame"), method = "calcSummaryRegion")
                     
                     #### initialisation
                     if(is.null(region.value)){
                       region.value <- selectRegion(object, region = region, type = "names")
                     }
                     n.regions <- length(region.value)
                     
                     if(is.null(param)){
                       param <- selectParameter(object, type = "contrastOnly")
                     }
                     n.param <- length(param)
                     
                     if(!is.list(fct.summary)){
                       fct.summary <- list(fct.summary)
                     }
                     n.fct_summary <- length(fct.summary)
                     
                     if(is.null(names(fct.summary))){
                       names.fct_summary <- paste0("fct",1:n.fct_summary)
                     }else{
                       names.fct_summary <- names(fct.summary)
                     }
                     
                     ##### main
                    subset <- list(NA)
                    names(subset) <- region
                    df.res <- NULL
                    
                     for(iter_region in  1:n.regions){
                       subset[[1]] <- names.regions[iter_region]
                       dt.Rcontrast <- selectContrast(object, param = param, subset = subset, coords = TRUE)
                       #dt.Rcontrast[, lapply(param,summary.fct) ,with = FALSE]
                       
                       res_tempo <- NULL
                       res_tempo <- lapply(1:n.fct_summary, function(x){dt.Rcontrast[, lapply(.SD, fct.summary[[x]]),.SDcols = param]})
                       res_tempo <- data.frame(matrix(unlist(res_tempo), ncol = n.param, byrow = TRUE))
                       colnames(res_tempo) <- param
                       
                       df.res <- rbind(df.res,
                                       cbind( subset, fct = names.fct_summary, res_tempo)
                       )
                     }
#                    
                    #### long to wide
                    if(long2wide){
                      df.res <- reshape(df.res, idvar = region, direction = "wide", timevar = "fct", sep = sep)
                      
                      if(format == "data.table"){
                        df.res <- as.data.table(df.res, key = c(region))  
                      }
                    }else{
                      if(format == "data.table"){
                        df.res <- as.data.table(df.res, key = c(region,"fct"))  
                      }
                    }
                    
                    #### export
                    return(list(res = df.res, 
                                verbose = verbose, 
                                update.object = update.object, 
                                overwrite = overwrite))
                    
                   }
)

####>>> calcRegionalContrast ####

methods::setMethod(f  = "calcRegionalContrast", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, bandwidth, power = 2, diagonal = FALSE,                                
                                         W = "ifany", W.range, W.spatial_res = c(1, 1, 1), 
                                         num = NULL, hemisphere = "both", name_newparam = paste(param, "regional", sep = "_"),
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     # definition de la matrice de voisinage
                     if(identical(W, "ifany") && any(!is.na(object@W$Wmatrix))){
                       if(verbose){cat(" loading W ... ")}
                       W <- selectW(object, num = num, hemisphere = hemisphere, upper = NULL)
                       
                     }else if(is.null(W) || identical(W, "ifany")){
                       if(verbose){cat(" computing W ... ")}
                       
                       W <- calcW(selectCoords(object, num = num, hemisphere = hemisphere), spatial_res = W.spatial_res, 
                                  range = W.range, upper = NULL, format = "dgCMatrix", row.norm = FALSE)$W               
                     }
                     
                     validClass(value = W, validClass = "dgCMatrix", superClasses = TRUE, method = "calcRegionalContrast[MRIaggr]")
                     
                     # convert distance to weights
                     if(diagonal == TRUE){
                       diag(W) <- -1
                       W@x[W@x > 0] <- EDK(W@x[W@x >= 0], bandwidth = bandwidth, power = power)  
                       W@x[W@x == -1] <- EDK(0, bandwidth = bandwidth, power = power)    
                     }else{
                       W@x <- EDK(W@x, bandwidth = bandwidth, power = power)                    
                     }
                     
                     # normalization
                     Rsum <- spam::rowSums(W)   
                     Rsum[is.na(Rsum)] <- -1                
                     W <- W / Rsum
                     
                     if(verbose){cat("W ready \n ")}
                     
                     carto <- selectContrast(object, param = param, num = num, hemisphere = hemisphere, format = "matrix")
                     p <- length(param)
                     
                     validDimension(value1 = W, validDimension = rep(nrow(carto),2), name1 = "W", name2 = NULL, type = c("nrow", "ncol"), method = "calcRegionalContrast[MRIaggr]")
                     validCharacter(value = name_newparam, validLength = p, method = "calcRegionalContrast[MRIaggr]")
                     
                     # calcul des valeurs regionales
                     carto_regional <- data.frame(matrix(NA, nrow = nrow(carto), ncol = p))
                     names(carto_regional) <- name_newparam
                     
                     if(verbose){cat("param : ")}
                     for(iter_param in 1:p){
                       if(verbose){cat(param[iter_param], " ", sep = "")}
                       carto_regional[,iter_param] <- as.matrix( get("%*%", envir = asNamespace("spam"))(W, carto[,iter_param, drop = FALSE]) )
                     }
                     if(verbose){cat("\n")}
                     
                     return(list(res = carto_regional, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite)
                     )
                     
                   }
)

####>>> calcROCthreshold ####

methods::setMethod(f  = "calcROCthreshold", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, mask, plot = "ROC_Youden", digit = 10,
                                         filename = paste(object@identifier, "calcROCthreshold", plot, sep = "_"), 
                                         update.object = FALSE, overwrite = FALSE, ...){
                     
                     # initPackage("ROCR", method = "calcROCthreshold[MRIaggr]")
                     
                     if(plot != FALSE){ #### graphical options
                       ## get graphical arguments
                       optionsMRIaggr.eg <- optionsMRIaggr()					 
                       dots.arguments <- list(...)
                       names_dots.arguments <- names(dots.arguments)
                       
                       validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                      validValues = c("height", "numeric2logical", "path", "res", "unit", "verbose", "width", "window"), 
                                      refuse.NULL = FALSE, method = "calcROCthreshold[MRIaggr]")
                       
                       ## set specific display
                       if(length(names_dots.arguments) > 0){
                         optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                       }
                       
                       ## create all variables
                       digit.legend <- optionsMRIaggr.eg$digit.legend
                       height <- optionsMRIaggr.eg$height
                       numeric2logical <- optionsMRIaggr.eg$numeric2logical
                       path <- optionsMRIaggr.eg$path
                       res <- optionsMRIaggr.eg$res
                       unit <- optionsMRIaggr.eg$unit
                       verbose <- optionsMRIaggr.eg$verbose
                       width <- optionsMRIaggr.eg$width
                       window <- optionsMRIaggr.eg$window
                     }
                     
                     initParameter(object = object, param = c(mask, param), checkArguments = TRUE, init = FALSE, 
                                   accept.coords = FALSE, accept.index = FALSE, 
                                   arg_name = "mask / param", method = "calcROCthreshold")     
                     
                     validDimension(value1 = mask, value2 = param, name1 = "mask", name2 = "param", type = "length", method = "calcROCthreshold[MRIaggr]")
                     
                     validCharacter(value = plot, validLength = 1, validValues = c(FALSE, "ROC_Youden", "ROC_prev", "boxplot_Youden", "boxplot_prev"), method = "calcROCthreshold[MRIaggr]")
                     
                     p <- length(mask)
                     
                     if(plot != FALSE){                       
                       res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                              n.plot = p, mfrow = NULL, xlim = NULL, ylim = NULL, method = "calcROCthreshold[MRIaggr]")
                       scale <- res.init$scale
                       mfrow <- res.init$mfrow
                     }
                     
                     #### mise en place du JDD ####
                     data <- selectContrast(object, param = c(mask, param))
                     
                     data[,mask] <- initMask(object, mask, checkArguments = TRUE, init = numeric2logical, 
                                             arg_name = "mask", long_name = "mask", method = "calcROCthreshold", format = "matrix")
                     
                     if(!is.null(digit)){
                       for(iter_param in param){data[,iter_param] <- round(data[,iter_param], digits = digit)}
                     }
                     
                     
                     #### etude des thresholds ####
                     res.ROC <- data.frame(matrix(NA, ncol = 10, nrow = p))
                     names(res.ROC) <- c("mask", "param", "AUC", "AUPRC", "OptTh_Youden", "Se_Youden", "Sp_Youden", "OptTh_prev", "Se_prev", "Sp_prev")
                     res.ROC[,"mask"] <- mask
                     res.ROC[,"param"] <- param
                     if(!is.null(plot) && plot %in% c("ROC_Youden", "ROC_prev")){
                       data_plot <- list()
                     }
                     
                     for(iter_param in 1:p){
                       
                       if(sum(data[,mask[iter_param]] == TRUE) == 0){
                         warning("calcROCthreshold[MRIaggr] : only FALSE were found for ", mask[iter_param], "\n")
                         next
                       }
                       
                       if(sum(data[,mask[iter_param]] == FALSE) == 0){
                         warning("calcROCthreshold[MRIaggr] : only TRUE were found for ", mask[iter_param], "\n")
                         next
                       }
                       
                       #### calcul du ROC ####
                       roc_tempo <- list()
                       prediction_tempo <- ROCR::prediction(data[,param[iter_param]], data[,mask[iter_param]])
                       performance_tempo <- ROCR::performance(prediction_tempo, x.measure = "spec", measure = "sens")
                       
                       roc_tempo$Specificity <- performance_tempo@x.values[[1]]
                       roc_tempo$Sensitivity <- performance_tempo@y.values[[1]]
                       roc_tempo$Threshold <- performance_tempo@alpha.values[[1]]
                       if(!is.null(plot) && plot %in% c("ROC_Youden", "ROC_prev")){
                         data_plot[[iter_param]] <- data.frame(Specificity = roc_tempo$Specificity, 
                                                               Sensitivity = roc_tempo$Sensitivity, 
                                                               Threshold = roc_tempo$Threshold)                
                       }
                       
                       performance_tempo <- ROCR::performance(prediction_tempo, x.measure = "rec", measure = "prec")
                       
                       roc_tempo$Recall <- performance_tempo@x.values[[1]]
                       roc_tempo$Precision <- performance_tempo@y.values[[1]]
                       
                       #### summary statistics ####
                       res.ROC[iter_param, "AUC"] <- ROCR::performance(prediction_tempo, measure = "auc")@y.values[[1]]
                       res.ROC[iter_param, "AUPRC"] <- calcAUPRC(x = NULL, y = NULL, performance = performance_tempo)["AUPRC"]
                       
                       #### seuil optimal : Youden ####
                       OptTh <- which.max(roc_tempo$Sensitivity + roc_tempo$Specificity)
                       res.ROC$OptTh_Youden[iter_param] <- roc_tempo$Threshold[OptTh]
                       res.ROC$Se_Youden[iter_param] <- roc_tempo$Sensitivity[OptTh]
                       res.ROC$Sp_Youden[iter_param] <- roc_tempo$Specificity[OptTh]
                       
                       #### seuil optimal : prevalence ####
                       prevalence <- mean(data[,mask[iter_param]])
                       OptTh <- which.max(prevalence * roc_tempo$Sensitivity + (1 - prevalence) * roc_tempo$Specificity)
                       res.ROC$OptTh_prev[iter_param] <- roc_tempo$Threshold[OptTh]
                       res.ROC$Se_prev[iter_param] <- roc_tempo$Sensitivity[OptTh]
                       res.ROC$Sp_prev[iter_param] <- roc_tempo$Specificity[OptTh]
                       
                     }
                     
                     #### display ####     
                     if(plot != FALSE){
                       
                       if(p > 4){
                         stop("calcROCthreshold[MRIaggr] : Too many parameters to be ploted \n", 
                              "maximum number of parameters : 4 \n ", 
                              "length(mask) : ", length(mask), "\n")
                       }
                       
                       initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                         mfrow = mfrow, bg = NULL, pty = NULL, mar = rep(3, 4), mgp = c(1.5, 0.5, 0))  
                       
                       for(iter_param in 1:p){             
                         
                         if(plot %in% c("ROC_Youden", "ROC_prev")){
                           prevalence <- mean(data[,mask[iter_param]])
                           
                           graphics::plot(1 - data_plot[[iter_param]]$Specificity, data_plot[[iter_param]]$Sensitivity, 
                                          main = paste(mask[iter_param], " ~ ", param[iter_param], " - patient ", object@identifier), 
                                          ylab = "sensitivity", xlab = "1 - specificity", type = "o")
                           graphics::points(c(0,1), c(0,1), type = "l", col = "grey")
                           
                           if(plot == "ROC_prev"){
                             graphics::points(1 - res.ROC$Sp_prev[iter_param], res.ROC$Se_prev[iter_param], col = "red", pch = 15, cex = 1.5)
                             graphics::legend(x = 0.2, y = 0.6, xjust = 0, yjust = 1, col = "red", pch = 15, 
                                              legend = paste("Th_Youden : ", round(res.ROC$OptTh_prev[iter_param], digit.legend), "\n", 
                                                             "Se : ", round(res.ROC$Se_prev[iter_param], digit.legend), "\n", 
                                                             "Sp : ", round(res.ROC$Sp_prev[iter_param], digit.legend), "\n"), 
                                              bty = "n")
                           }
                           
                           if(plot == "ROC_Youden"){
                             graphics::points(1 - res.ROC$Sp_Youden[iter_param], res.ROC$Se_Youden[iter_param], col = "blue", pch = 15, cex = 1.5)
                             graphics::legend(x = 0.2, y = 0.6, xjust = 0, yjust = 1, col = "blue", pch = 15, 
                                              legend = paste("Th[Y / prev] : ", round(res.ROC$OptTh_Youden[iter_param], digit.legend), "\n", 
                                                             "Se : ", round(res.ROC$Se_Youden[iter_param], digit.legend), "\n", 
                                                             "Sp : ", round(res.ROC$Sp_Youden[iter_param], digit.legend), "\n", 
                                                             "prevalence : ", round(prevalence, digit.legend)), 
                                              bty = "n")
                           }
                         }
                         
                         if(plot %in% c("boxplot_Youden", "boxplot_prev")){
                           graphics::boxplot(data[,param[iter_param]] ~ data[,mask[iter_param]], 
                                             ylab = param[iter_param], main = paste(mask[iter_param], " - ", object@identifier))  
                           
                           if(plot == "boxplot_Youden"){
                             graphics::abline(h = res.ROC$OptTh_Youden[iter_param], col = "red")
                             graphics::text(x = 1.5, y = 1.1 * res.ROC$OptTh_Youden[iter_param], col = "red", labels = "Youden", bty = "n")
                           }
                           if(plot == "boxplot_prev"){
                             graphics::abline(h = res.ROC$OptTh_prev[iter_param], col = "blue")
                             graphics::text(x = 1.5, y = 1.1 * res.ROC$OptTh_prev[iter_param], col = "blue", labels = "prev", bty = "n")
                           }
                         }
                         
                       }
                       
                       if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                         grDevices::dev.off()
                       }
                       
                     }
                     
                     return(list(res = res.ROC, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                   }
)

####>>> calcSmoothMask ####

methods::setMethod(f  = "calcSmoothMask", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask = "mask", numeric2logical = FALSE, 
                                         size_2Dgroup = 50, Neighborhood_2D = "3D_N8", rm.2Dhole = FALSE, 
                                         size_3Dgroup = "unique", Neighborhood_3D = "3D_N10", rm.3Dhole = TRUE, erosion.th = 0.75,  
                                         Vmask_min = 0.25, Vbackground_max = 0.75, Neighborhood_V = "3D_N10", 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
                   { #### preparation
                     if(is.null(size_2Dgroup) || (size_2Dgroup != FALSE && size_2Dgroup != "unique" && is.numeric(size_2Dgroup) == FALSE)){
                       stop("calcSmoothMask[MRIaggr] : wrong specification of \'size_2Dgroup\' \n", 
                            "\'size_2Dgroup\' must be FALSE, numeric or be \"unique\" \n", 
                            "proposed value  : ", size_2Dgroup, "\n")
                     }
                     
                     if(is.null(size_3Dgroup) || (size_3Dgroup != FALSE && size_3Dgroup != "unique" && is.numeric(size_3Dgroup) == FALSE)){
                       stop("calcSmoothMask[MRIaggr] : wrong specification of \'size_3Dgroup\' \n", 
                            "\'size_3Dgroup\' must be FALSE, numeric or be \"unique\" \n", 
                            "proposed value  : ", size_3Dgroup, "\n")
                     }
                     
                     if(is.null(erosion.th) || (erosion.th != FALSE &&  is.numeric(erosion.th) == FALSE)){
                       stop("calcSmoothMask[MRIaggr] : wrong specification of \'erosion.th\' \n", 
                            "\'erosion.th\' must be FALSE or numeric  \n", 
                            "proposed value  : ", erosion.th, "\n")
                     }
                     
                     coords <- selectCoords(object)            
                     n <- selectN(object)            
                     num <- unique(coords$k)
                     n.slices <- length(num)
                     
                     index_tous <- 1:n            
                     
                     #### Application du mask            
                     index_fond <- NULL
                     index_mask <- index_tous
                     
                     if(is.character(mask)){
                       mask <- selectContrast(object, param = mask, format = "matrix", coord = FALSE)
                     }
                     if(numeric2logical == TRUE){
                       mask <- as.logical(mask)
                     }
                     if(!is.logical(mask)){
                       stop("calcSmoothMask[MRIaggr] : \'mask\' is not of type logical \n", 
                            "proposed type : ", paste(is(mask), collapse = " "), "\n", 
                            "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
                     }    
                     
                     index_fond <- index_tous[which(mask == 0)]
                     index_mask.ref <- index_tous[which(mask == 1)] 
                     index_mask <- index_mask.ref
                     
                     #### Exclusion des petits groupes 2D
                     if( identical(size_2Dgroup, FALSE) == FALSE && length(index_mask) > 0)
                     { if(verbose){cat("rm small2D : ")} 
                       
                       group2D <- calcGroupsCoords(coords = coords[index_mask,], 
                                                   Neighborhood = Neighborhood_2D, verbose = FALSE)
                       
                       valid_group2D <- numeric()
                       for(iter_num in 1:n.slices){
                         index_num <- which(group2D$df.group[,"k"] == iter_num)
                         group2D_num <- unique(group2D$df.group[index_num, "group"])
                         
                         if(size_2Dgroup == "unique"){
                           valid_group2D <- c(valid_group2D, 
                                              group2D_num[which.max(group2D$group_size[group2D_num])])
                         }else{
                           valid_group2D <- c(valid_group2D, 
                                              group2D_num[group2D$group_size[group2D_num] > size_2Dgroup])
                         }
                       }
                       index_unvalid <- which(group2D$df.group$group %in%  valid_group2D == FALSE)
                       
                       if(length(index_unvalid > 0)){
                         index_fond <- sort(union(index_fond, index_mask[index_unvalid]))
                         index_mask <- sort(index_mask[-index_unvalid])
                       }
                       
                       if(verbose){cat(" ", sum(group2D$group_size[-valid_group2D]), 
                                       " vx ; ", length(group2D$group_size[-valid_group2D]), " groups \n", sep = "")}              
                     }
                     
                     # Rebouchage des trous des trous 2D
                     if( identical(rm.2Dhole, FALSE) == FALSE  && length(index_fond) > 0)
                     {  if(verbose){cat("add hole2D : ")}
                       
                       group2D <- calcGroupsCoords(coords = coords[index_fond,], 
                                                   Neighborhood = Neighborhood_2D, verbose = FALSE)
                       
                       valid_group2D <- numeric()
                       for(iter_num in 1:n.slices){
                         index_num <- which(group2D$df.group[,"k"] == iter_num)
                         group2D_num <- unique(group2D$df.group[index_num, "group"])
                         
                         valid_group2D <- c(valid_group2D, 
                                            group2D_num[which.max(group2D$group_size[group2D_num])])                 
                       }
                       index_unvalid <- which(group2D$df.group$group %in%  valid_group2D == FALSE)
                       
                       if(length(index_unvalid > 0)){
                         index_mask <- sort(union(index_mask, index_fond[index_unvalid]))
                         index_fond <- sort(index_fond[-index_unvalid])
                       }
                       
                       if(verbose){cat(" ", sum(group2D$group_size[-valid_group2D]), 
                                       " vx ; ", length(group2D$group_size[-valid_group2D]), " groups \n", sep = "")}
                     }  
                     
                     #### Exclusion des petits groupes 3D
                     if( identical(size_3Dgroup, FALSE) == FALSE && length(index_mask) > 0)
                     { if(verbose){cat("rm small3D : ")}
                       
                       group3D <- calcGroupsCoords(coords = coords[index_mask,], 
                                                   Neighborhood = Neighborhood_3D, verbose = FALSE)
                       
                       if(size_3Dgroup == "unique"){
                         valid_group3D <- which.max(group3D$group_size)
                       }else{
                         valid_group3D <- which(group3D$group_size > size_3Dgroup)                          
                       }         
                       index_unvalid <- which(group3D$df.group$group %in% valid_group3D == FALSE)
                       
                       if(length(index_unvalid > 0)){
                         index_fond <- sort(union(index_fond, index_mask[index_unvalid]))              
                         index_mask <- sort(index_mask[-index_unvalid])              
                       }
                       
                       if(verbose){cat(" ", sum(group3D$group_size[-valid_group3D]), 
                                       " vx ; ", length(group3D$group_size[-valid_group3D]), " groups \n", sep = "")}              
                     }
                     
                     #### erosion
                     if( identical(erosion.th, FALSE) == FALSE && size_3Dgroup == "unique" && !is.null(size_3Dgroup)  && length(index_mask) > 0){
                       
                       if(verbose){cat("erosion : ")}
                       
                       Amask <- dt2array(rep(0, n), 
                                         coords = coords)$contrast[[1]]              
                       
                       # identification des observations a eroder
                       Amask[index_fond] <- FALSE
                       Amask[index_mask] <- TRUE
                       Vlocal <- calcFilter(Amask, filter = Neighborhood_V, norm.filter = TRUE)$res
                       
                       index_erosion <- intersect(which(Vlocal < erosion.th), index_mask)
                       
                       index_fond <- sort(union(index_fond, index_erosion))
                       index_mask <- sort(setdiff(index_mask, index_erosion))
                       
                       # identification des nouveaux groupes spatiaux
                       Amask[index_fond] <- NA
                       Amask[index_mask] <- TRUE
                       
                       group3D <- calcGroupsCoords(array = Amask, 
                                                   Neighborhood = Neighborhood_3D, verbose = FALSE)
                       
                       valid_group3D <- which.max(group3D$group_size)
                       
                       index_unvalid <- which(group3D$df.group$group %in%  valid_group3D == FALSE)
                       
                       if(length(index_unvalid > 0)){
                         index_fond <- sort(union(index_fond, index_mask[index_unvalid]))
                         index_mask <- sort(index_mask[-index_unvalid])              
                       }
                       
                       # gestion des observations erodes
                       Amask[index_fond] <- FALSE
                       Amask[index_mask] <- TRUE
                       Vlocal <- calcFilter(Amask, filter = "2D_N8", norm.filter = TRUE)$res
                       
                       index_fond <- setdiff(index_fond, index_erosion[which(Vlocal[index_erosion] > 0)])
                       index_fond <- sort(index_fond)              
                       index_mask <- union(index_mask, index_erosion[which(Vlocal[index_erosion] > 0)])
                       index_mask <- sort(index_mask)
                       
                       if(verbose){cat(" ", sum(group3D$group_size[-valid_group3D]) - length(intersect(index_fond, index_erosion[which(Vlocal[index_erosion] > 0)])), 
                                       " vx ; ", length(group3D$group_size[-valid_group3D]), " groups \n", sep = "")} 
                     }
                     
                     #### Rebouchage des trous des trous 3D
                     if( identical(rm.3Dhole, FALSE) == FALSE  && length(index_fond) > 0)
                     {  if(verbose){cat("add hole3D : ")}
                       
                       group3D <- calcGroupsCoords(coords = coords[index_fond,], 
                                                   Neighborhood = Neighborhood_3D, verbose = FALSE)
                       
                       valid_group3D <- which.max(group3D$group_size)                         
                       index_unvalid <-which(group3D$df.group$group %in%  valid_group3D == FALSE)
                       
                       if(length(index_unvalid > 0)){
                         index_mask <- sort(union(index_mask, index_fond[index_unvalid]))
                         index_fond <- sort(index_fond[-index_unvalid])          
                       }
                       
                       if(verbose){cat(" ", sum(group3D$group_size[-valid_group3D]), 
                                       " vx ; ", length(group3D$group_size[-valid_group3D]), " groups \n", sep = "")}
                     }  
                     
                     # lissage du masque
                     if( identical(Vmask_min, FALSE) == FALSE  || identical(Vbackground_max, FALSE) == FALSE){
                       
                       if(verbose){cat("smoothing (add / rm) : ")}
                       
                       Amask <- dt2array(rep(0, n), 
                                         coords = coords)$contrast[[1]]
                       iter_max <- 20
                       iter <- 1
                       n.modif <- 1
                       
                       while(n.modif > 0 && iter < iter_max){
                         Amask[index_fond] <- FALSE
                         Amask[index_mask] <- TRUE
                         n.modif <- 0
                         
                         Vlocal <- calcFilter(Amask, filter = Neighborhood_V, norm.filter = TRUE)$res
                         
                         if(!is.null(Vmask_min)){
                           index_rm <-  intersect(index_mask, which(Vlocal < Vmask_min))
                           n.modif <- n.modif + length(index_rm)
                           if(length(index_rm) > 0){
                             index_fond <- union(index_fond, index_rm)
                             index_mask <- setdiff(index_mask, index_rm)                
                           }
                         }
                         
                         if(!is.null(Vbackground_max)){
                           index_add <-  intersect(index_fond, which(Vlocal > Vbackground_max))
                           n.modif <- n.modif + length(index_add)
                           if(length(index_add) > 0){
                             index_fond <- setdiff(index_fond, index_add)
                             index_mask <- union(index_mask, index_add)                
                           } 
                         }
                         if(verbose){cat("(", iter, ") ", length(index_add), "/", length(index_rm), "  ", sep = "")}
                         
                         iter <- iter + 1
                       }
                       if(verbose){cat("\n")}
                     }
                     
                     
                     
                     # mise en forme pour l export            
                     mask <- rep(FALSE, n)
                     mask[index_mask] <- TRUE
                     
                     return(list(res = data.frame(mask = mask, coords), 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                   }
)

####>>> calcTableHypoReperf ####

methods::setMethod(f  = "calcTableHypoReperf", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, timepoint, threshold = 1:10, sep = "_", norm_mu = FALSE, norm_sigma = FALSE, 
                                         mask = NULL, numeric2logical = FALSE, param.update = "reperf", 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
                   {
                     #### preparation ####
                     if(length(timepoint) %in% c(1, 2) == FALSE){
                       stop("calcTableHypoReperf[MRIaggr] : wrong specification of \'timepoint\'\n", 
                            "timepoint must have length 1 (hypoperfusion only) or 2 (hypoperfusion and reperfusion) \n", 
                            "length(timepoint) : ", length(timepoint), "\n")
                     }
                     
                     validCharacter(value = param.update, validLength = 1, validValues = c("shift", "reperf", "reperf_pc", "deperf", "deperf_pc"), method = "calcTableHypoReperf[MRIaggr]")
                     
                     #### hypoperfusion ####
                     
                     ## mise en place 
                     param_time1 <- paste(param, timepoint[1], sep = sep)
                     initParameter(object = object, param = param_time1, checkArguments = TRUE, init = FALSE, accept.coords = FALSE, 
                                   arg_name = "param_time1", method = "calcTableHypoReperf")
                     if(length(timepoint) == 2){
                       param_time2 <- paste(param, timepoint[2], sep = sep)
                       initParameter(object = object, param = param_time2, checkArguments = TRUE, init = FALSE, accept.coords = FALSE, 
                                     arg_name = "param_time2", method = "calcTableHypoReperf")
                     }
                     if(!is.null(mask)){
                       initParameter(object = object, param = mask, checkArguments = TRUE, init = FALSE, accept.coords = FALSE, 
                                     arg_name = "mask", method = "calcTableHypoReperf")
                     }
                     
                     n.param <- length(param)
                     n.threshold <- length(threshold)
                     
                     index_mask <- selectContrast(object, param = "index", format = "vector")
                     n.mask <- length(index_mask)
                     
                     res <- list()
                     if(length(timepoint) == 1){
                       res$volume_hypo <- data.frame(matrix(NA, nrow = n.threshold, ncol = 3 * n.param + 1))
                       names(res$volume_hypo ) <- c("threshold", paste("Vhypo.", param_time1, sep = ""), 
                                                    paste("Vmismatch.", param, sep = ""), paste("PCmismatch.", param, sep = ""))
                     }else{
                       res$volume_hypo <- data.frame(matrix(NA, nrow = n.threshold, ncol = 4 * n.param + 1))
                       names(res$volume_hypo ) <- c("threshold", paste("Vhypo.", param_time1, sep = ""), paste("Vhypo.", param_time2, sep = ""), 
                                                    paste("Vmismatch.", param, sep = ""), paste("PCmismatch.", param, sep = ""))
                     }
                     row.names(res$volume_hypo) <- c(as.character(threshold))
                     res$volume_hypo$threshold <- threshold
                     
                     ## data hypoperfusion
                     dataRaw.time1 <- selectContrast(object, param = param_time1, norm_mu = norm_mu, norm_sigma = norm_sigma)
                     
                     data.time1 <- apply(dataRaw.time1, 2, function(x){
                       tempo <- cut(x, breaks = c(-Inf, threshold - 10^{-12}, Inf))
                       levels(tempo) <- c(threshold[1] - 1, threshold, threshold[n.threshold] + 1)
                       return(as.numeric(as.character(tempo)))
                     })
                     dataRaw.time1[dataRaw.time1 < 0] <- 0
                     dataRaw.time1[dataRaw.time1 > threshold[n.threshold]] <- threshold[n.threshold]
                     
                     if(!is.null(mask)){
                       data.mask <- initMask(object, mask, test = TRUE, checkArguments = numeric2logical, 
                                             arg_name = "mask", long_name = "mask", method = "calcTableHypoReperf")
                       
                       test.mask <- (data.mask == FALSE)
                     }
                     
                     #### reperfusion ####
                     if(length(timepoint) == 2){
                       ## mise en place
                       res$pixel <- data.frame(matrix(NA, nrow = selectN(object), ncol = 5 * n.param + 3))
                       res$pixel[,1:3] <- selectCoords(object)
                       names(res$pixel) <- c("i", "j", "k", paste(param, "_shift", sep = ""), 
                                             paste(param, "_reperf", sep = ""), paste(param, "_reperf_pc", sep = ""), 
                                             paste(param, "_deperf", sep = ""), paste(param, "_deperf_pc", sep = ""))
                       
                       res$volume_reperf <- data.frame(matrix(NA, nrow = n.threshold, ncol = 12 * n.param + 1))
                       names(res$volume_reperf) <- c("threshold", 
                                                     paste("Vreperf.", param, sep = ""), paste("PCreperf.", param, sep = ""), paste("VreperfW.", param, sep = ""), paste("PCreperfW.", param, sep = ""), 
                                                     paste("Vdeperf.", param, sep = ""), paste("PCdeperf.", param, sep = ""), paste("VdeperfW.", param, sep = ""), paste("PCdeperfW.", param, sep = ""), 
                                                     paste("Vshift_reperf.", param, sep = ""), paste("PCshift_reperf.", param, sep = ""), 
                                                     paste("Vshift_deperf.", param, sep = ""), paste("PCshift_deperf.", param, sep = ""))
                       row.names(res$volume_reperf) <- c(as.character(threshold))
                       res$volume_reperf$threshold <- threshold
                       
                       ## data              
                       dataRaw.time2 <- selectContrast(object, param = param_time2, norm_mu = norm_mu, norm_sigma = norm_sigma)
                       
                       data.time2 <- apply(dataRaw.time2, 2, function(x){
                         tempo <- cut(x, breaks = c(-Inf, threshold - 10^{-12}, Inf))
                         levels(tempo) <- c(threshold[1] - 1, threshold, threshold[n.threshold] + 1)
                         return(as.numeric(as.character(tempo)))
                       })
                       dataRaw.time2[dataRaw.time2 < 0] <- 0
                       dataRaw.time2[dataRaw.time2 > threshold[n.threshold]] <- threshold[n.threshold]
                       
                       res$pixel[index_mask, paste(param, "_shift", sep = "")] <- dataRaw.time2 - dataRaw.time1
                     }
                     
                     #### calcul ####
                     
                     for(iter_param in param){
                       if(verbose){cat(iter_param, " ", sep = "")}
                       
                       iter_param1 <- paste(iter_param, sep, timepoint[1], sep = "")
                       index_H0pos <- which(dataRaw.time1[,iter_param1] > 0)
                       index_H0lim <- which(dataRaw.time1[,iter_param1] < threshold[n.threshold])
                       
                       if(length(timepoint) == 2){
                         iter_param2 <- paste(iter_param, sep, timepoint[2], sep = "")
                         iter_paramShift <- paste(iter_param, "_shift", sep = "")
                         iter_paramReperf <- paste(iter_param, "_reperf", sep = "")
                         iter_paramReperfPC <- paste(iter_param, "_reperf_pc", sep = "")
                         iter_paramDeperf <- paste(iter_param, "_deperf", sep = "")
                         iter_paramDeperfPC <- paste(iter_param, "_deperf_pc", sep = "")
                         
                         # reperf px              
                         index_reperf <- which(res$pixel[index_mask, iter_paramShift] < 0)
                         res$pixel[index_mask[index_reperf], iter_paramReperf] <- -data.time1[index_reperf, iter_param1]
                         res$pixel[index_mask[index_reperf], iter_paramReperfPC] <- res$pixel[index_mask[index_reperf], iter_paramShift] / (data.time1[index_reperf, iter_param1])
                         
                         # deperf px
                         index_deperf <- which(res$pixel[index_mask, iter_paramShift] > 0)
                         res$pixel[index_mask[index_deperf], iter_paramDeperf] <- data.time1[index_deperf, iter_param1]
                         res$pixel[index_mask[index_deperf], iter_paramDeperfPC] <- res$pixel[index_mask[index_deperf], iter_paramShift] / (threshold[n.threshold]-data.time1[index_deperf, iter_param1])
                       }
                       
                       for(iter_threshold in 1:n.threshold){
                         test.hypoH0 <- data.time1[,iter_param1] >= threshold[iter_threshold]
                         res$volume_hypo[iter_threshold, paste("Vhypo.", iter_param1, sep = "")] <- sum(test.hypoH0, na.rm = TRUE)
                         
                         if(!is.null(mask)){
                           res$volume_hypo[iter_threshold, paste("Vmismatch.", iter_param, sep = "")] <- sum( test.hypoH0*test.mask, na.rm = TRUE)
                         }
                         
                         if(length(timepoint) == 2){
                           test.hypoH2 <- data.time2[,iter_param2] >= threshold[iter_threshold]
                           res$volume_hypo[iter_threshold, paste("Vhypo.", iter_param2, sep = "")] <- sum(test.hypoH2, na.rm = TRUE)
                           
                           test.hyperH0 <- data.time1[,iter_param1] < threshold[iter_threshold]
                           test.hyperH2 <- data.time2[,iter_param2] < threshold[iter_threshold]
                           index_n0Reperf <- intersect(index_H0pos, index_reperf)
                           index_n0Deperf <- intersect(index_H0lim, index_deperf)
                           w.reperf <- rep(0, n.mask)
                           w.reperf[index_n0Reperf] <- -res$pixel[index_mask, iter_paramShift][index_n0Reperf]/dataRaw.time1[index_n0Reperf, iter_param1]
                           w.deperf <- rep(0, n.mask)
                           w.deperf[index_n0Deperf] <- res$pixel[index_mask, iter_paramShift][index_n0Deperf] / (threshold[n.threshold]-dataRaw.time1[index_n0Deperf, iter_param1])
                           
                           # Vreperf et Vdeperf
                           res$volume_reperf[iter_threshold, paste("Vreperf.", iter_param, sep = "")] <- sum( test.hypoH0 * test.hyperH2, na.rm = T)
                           res$volume_reperf[iter_threshold, paste("VreperfW.", iter_param, sep = "")] <- sum( test.hypoH0 * test.hyperH2 * w.reperf, na.rm = T)
                           
                           #                   if(!is.null(mask)){
                           res$volume_reperf[iter_threshold, paste("Vdeperf.", iter_param, sep = "")] <- sum( test.hyperH0 * test.hypoH2, na.rm = T)
                           res$volume_reperf[iter_threshold, paste("VdeperfW.", iter_param, sep = "")] <- sum( test.hyperH0 * test.hypoH2 * w.deperf, na.rm = T)
                           #                   }
                           
                           # Vshift
                           res$volume_reperf[iter_threshold, paste("Vshift_reperf.", iter_param, sep = "")] <- sum(res$pixel[index_mask, iter_paramShift] <= -threshold[iter_threshold], na.rm = TRUE)
                           res$volume_reperf[iter_threshold, paste("Vshift_deperf.", iter_param, sep = "")] <- sum(res$pixel[index_mask, iter_paramShift] >= threshold[iter_threshold], na.rm = TRUE)
                         }    
                       }
                       
                       if(!is.null(mask)){
                         res$volume_hypo[,paste("PCmismatch.", iter_param, sep = "")] <- res$volume_hypo[,paste("Vmismatch.", iter_param, sep = "")]/sum(data.mask == TRUE)
                       }
                       if(length(timepoint) == 2){
                         res$volume_reperf[,paste("PCreperf.", iter_param, sep = "")] <- res$volume_reperf[,paste("Vreperf.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]
                         res$volume_reperf[,paste("PCdeperf.", iter_param, sep = "")] <- res$volume_reperf[,paste("Vdeperf.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]
                         res$volume_reperf[,paste("PCreperfW.", iter_param, sep = "")] <- res$volume_reperf[,paste("VreperfW.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]
                         res$volume_reperf[,paste("PCdeperfW.", iter_param, sep = "")] <- res$volume_reperf[,paste("VdeperfW.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]
                         res$volume_reperf[,paste("PCshift_reperf.", iter_param, sep = "")] <- res$volume_reperf[,paste("Vshift_reperf.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]
                         res$volume_reperf[,paste("PCshift_deperf.", iter_param, sep = "")] <- res$volume_reperf[,paste("Vshift_deperf.", iter_param, sep = "")]/res$volume_hypo[,paste("Vhypo.", iter_param1, sep = "")]              
                       }
                     }
                     if(verbose){cat("\n")}
                     
                     
                     
                     #### export ####
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 param.update = param.update, 
                                 overwrite = overwrite)
                     )
                   }
)

####>>> calcTableLesion ####

methods::setMethod(f  = "calcTableLesion", 
                   signature  = "MRIaggr", 
                   definition = function(object, maskN, mask = NULL, numeric2logical = FALSE, 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     fieldDim <- object@fieldDim
                     params <- names(object@contrast)
                     
                     #### tests preliminaires : maskN
                     initParameter(object = object, param = maskN, checkArguments = TRUE, init = FALSE, accept.coords = TRUE, 
                                   arg_name = "maskN", long_name = "parameters", method = "calcTableLesion")          
                     n.maskN <- length(maskN)
                     
                     data <- selectContrast(object, param = maskN, coords = "k", format = "data.frame")
                     
                     data[,maskN] <- initMask(object, maskN, checkArguments = TRUE, init = numeric2logical, 
                                              arg_name = "mask", long_name = "mask", method = "calcTableLesion", format = "matrix")
                     
                     #### tests preliminaires : mask
                     if(length(mask) > 0){
                       if(is.character(mask)){
                         mask <- selectContrast(object, param = mask, format = "vector")
                       }else{
                         mask <- as.vector(mask)
                         if(!is.null(mask) && length(mask) != selectN(object)){
                           stop("calcTableLesion[MRIaggr] : length of \'mask\' incompatible with \'object\' \n", 
                                "number of observations in \'object\' : ", selectN(object), "\n", 
                                "length(mask) : ", length(mask), "\n")
                         }
                       }
                       if(numeric2logical == TRUE){mask <- as.logical(mask)}
                       if(!is.logical(mask)){
                         stop("calcTableLesion[MRIaggr] : type of \'mask\' is not logical \n", 
                              "proposed type : ", paste(is(mask), collapse = " "), "\n", 
                              "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
                       }
                       data <- data.frame(data, mask = mask)
                     }
                     
                     matlevels <- matrix(levels(interaction(maskN, maskN, sep = "_outside_")), 
                                         nrow = n.maskN)
                     
                     names.res <- c(maskN, 
                                    matlevels[lower.tri(matlevels)], 
                                    matlevels[upper.tri(matlevels)]
                     )
                     
                     if(!is.null(mask)){
                       names.res <- c("mask", names.res, paste(maskN, "_outside_mask", sep = ""))
                     }
                     
                     res <- data.frame(matrix(NA, ncol = length(names.res), nrow = fieldDim$k + 1))
                     names(res) <- names.res
                     rownames(res)[fieldDim$k + 1] <- c("total")
                     
                     
                     #### calcul des tables : sans interaction
                     for(iter_param in c(if(!is.null(mask)){"mask"}, maskN)){
                       table_tempo <- tapply(data[,iter_param], data[,"k"], 
                                             function(x){table_tempo <- table(x) ; if("TRUE" %in% names(table_tempo)){table_tempo["TRUE"]}else{0}})  
                       res[,iter_param] <- c(table_tempo, sum(table_tempo))              
                     }
                     
                     #### calcul des tables : avec interaction
                     for(iter_param in c(matlevels[lower.tri(matlevels)], matlevels[upper.tri(matlevels)])){
                       param_tempo <- unlist(strsplit(iter_param, split = "_outside_"))
                       
                       table_tempo <- tapply(data[,param_tempo[1]] > data[,param_tempo[2]], data[,"k"], function(x){table(x)["TRUE"]})  
                       table_tempo[is.na(table_tempo)] <- 0
                       res[,iter_param] <- c(table_tempo, sum(table_tempo))              
                     }
                     
                     #### calcul des tables : avec le mask
                     if(!is.null(mask)){
                       for(iter_param in maskN){
                         param_tempo <- c(iter_param, "mask")                
                         table_tempo <- tapply(data[,param_tempo[1]]-data[,param_tempo[2]], data[,"k"], function(x){table(x)["1"]})  
                         table_tempo[is.na(table_tempo)] <- 0
                         res[,paste(iter_param, "mask", sep = "_outside_")] <- c(table_tempo, sum(table_tempo))              
                       }
                     }
                     
                     
                     #### export
                     return(list(res = res, 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }
)

####>>> calcThresholdMRIaggr ####

methods::setMethod(f  = "calcThresholdMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, hemisphere = "both", rm.CSF = FALSE, threshold = 1:10, decreasing = FALSE, 
                                         GRalgo = FALSE, W = "ifany", seed = NULL, numeric2logical = FALSE, W.range, W.spatial_res = rep(1, 3), 
                                         name_newparam = paste(param, "Th", sep = "_"), 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     #### pre initialization 
                     data <- selectContrast(object, param = c("index", param), hemisphere = hemisphere)
                     n <- nrow(data)
                     
                     if(is.numeric(seed)){
                       seed_tempo <- seed
                       seed <- rep(FALSE,n)
                       seed[seed_tempo] <- TRUE
                     }
                     
                     #### test
                     if (optionsMRIaggr("checkArguments")) {
                       
                       validDimension(value1 = param, value2 = name_newparam, type = "length", method = "calcThresholdMRIaggr[MRIaggr]")
                       
                       if(GRalgo == TRUE){
                         if (is.null(seed)) {
                           stop("calcThresholdMRIaggr : wrong specification of argument \'seed\' \n",
                                "argument \'seed\' must not be NULL if \'GRalgo\' is set to TRUE \n")
                         } else if (is.character(seed)) {
                           validCharacter(value = seed, validLength = NULL, validValues = selectParameter(object), refuse.NULL = FALSE, method = "calcThresholdMRIaggr[MRIaggr]")
                         } else {
                           validInteger(value = seed, validLength = NULL, min = 1, max = n, refuse.duplicates = TRUE, method = "calcThresholdMRIaggr[MRIaggr]")
                         }
                       }
                       
                       validLogical(value = update.object, validLength = 1, method = "calcThresholdMRIaggr[MRIaggr]")
                       validLogical(value = overwrite, validLength = 1, method = "calcThresholdMRIaggr[MRIaggr]")
                       
                     }
                     
                     #### initialization 
                     if (GRalgo == TRUE) {
                       
                       if (is.character(seed)) {
                         if(length(seed)==1){
                           data$seed <- selectContrast(object, param = seed, hemisphere = hemisphere, format = "vector")
                         }else{
                           data$seed <- as.logical(apply(selectContrast(object, param = seed, hemisphere = hemisphere, format = "matrix"),1,prod))
                         }
                       } else{
                         data$seed <- FALSE
                         data$seed[seed] <- TRUE
                       }
                       
                       if (identical(W, "ifany") &&  identical(selectW(object, type = "upper"), TRUE)) {
                         
                         W <- selectW(object, subset_W = data$index, upper = FALSE)
                         
                       }else if(is.null(W) || identical(W, "ifany")){
                         
                         if(verbose == TRUE){cat("computing W ... \n")}
                         coords <- selectCoords(object)[data$index,]                
                         W <- calcW(coords, range = W.range, method = "euclidean", upper = NULL, format = "dgCMatrix", row.norm = TRUE, 
                                    spatial_res = W.spatial_res)$W                
                         
                       }
                       
                     }
                     
                     if(rm.CSF == TRUE){
                       index_CSF <- which(apply(selectContrast(object, param = c("CSF", "GM", "WM"), hemisphere = hemisphere), 1, which.max) == 1)             
                       if(length(index_CSF)>0){
                         data <- data[-index_CSF,]
                         if(GRalgo == TRUE){
                           W <- W[-index_CSF, -index_CSF]
                         }
                       }
                     }
                     
                     #### computation 
                     res.th <- calcThreshold(contrast = cbind(data, CSF = 0, hemisphere = hemisphere), param = param, 
                                             hemisphere = if(hemisphere == "both"){NULL}else{hemisphere}, rm.CSF = rm.CSF, 
                                             threshold = threshold, decreasing = decreasing, 
                                             W = W, GRalgo = GRalgo, numeric2logical = numeric2logical, seed = if(GRalgo == TRUE){"seed"}else{FALSE}, verbose = verbose)
                     
                     names(res.th) <- name_newparam
                     
                     #### export
                     res <- data.frame(matrix(0, nrow = selectN(object), ncol = ncol(res.th)))
                     names(res) <- name_newparam
                     res[data$index,] <- res.th
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 name_newparam = name_newparam, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }
)

####>>> calcTissueType ####

methods::setMethod(f  = "calcTissueType", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, niter = 100, nnei = 6, 
                                         beta = if(sub == TRUE){0.3}else{0.7}, sub = TRUE, digit = 0, verbose = TRUE, 
                                         name_newparam = c("CSF","GM","WM"), update.object = FALSE, overwrite = FALSE){
                     
                     initPackage("mritc", method = "calcTissueType[MRIaggr]")
                     
                     carto <- selectContrast(object, param = param, format = "vector")
                     coords <- selectCoords(object)
                     
                     init  <- mritc::initOtsu(round(carto, digits = digit), m = length(name_newparam) - 1)
                     
                     mask <- dt2array(rep(1, selectN(object)), 
                                      coords, default_value = 0)$contrast[[1]]
                     
                     W <-  mritc::makeMRIspatial(mask, nnei = nnei, sub = sub)
                     
                     res <- mritc::mritc.bayes(y = carto, 
                                               neighbors = W$neighbors, 
                                               blocks = W$blocks, 
                                               mu = init$mu, 
                                               sigma = init$sigma, 
                                               sub = sub, 
                                               niter = niter, 
                                               subvox = W$subvox, 
                                               verbose = verbose)
                     
                     order <- order(res$mu, decreasing = FALSE)
                     res$prob <- res$prob[,order]
                     res$sigma <- res$sigma[order]
                     res$mu <- res$mu[order]
                     
                     return(list(res = res, 
                                 verbose = verbose, 
                                 name_newparam = name_newparam, 
                                 update.object = update.object, 
                                 overwrite = overwrite)
                     )
                   }
)

####>>> calcW ####

methods::setMethod(f  = "calcW", 
                   signature  = "MRIaggr", 
                   definition = function(object, range, spatial_res = c(1, 1, 1), num = NULL, hemisphere = "both", subset = NULL, 
                                         upper = TRUE, format = "dgCMatrix", row.norm = FALSE, calcBlockW = FALSE, 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     #### test package W
                     # initPackage("spam", method = "calcW[MRIaggr]")
                     # initPackage("Matrix", method = "calcW[MRIaggr]")
                     
                     resW <- calcW(object = selectCoords(object, num = num, hemisphere = hemisphere, subset = subset), 
                                   spatial_res = spatial_res, range = range, upper = upper, format = format, row.norm = row.norm, 
                                   calcBlockW = calcBlockW)
                     
                     W <- Matrix::drop0(resW$W, is.Csparse = TRUE)       
                     
                     return(list(res = list(W = W, blocks = resW$blocks, upper = upper), 
                                 verbose = verbose, 
                                 update.object = update.object, 
                                 overwrite = overwrite))
                     
                   }
)
