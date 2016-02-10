#**********************************************************************
#**********************************************************************
#*************         A Functions MRIaggr          *******************
#**********************************************************************
#**********************************************************************
#

#####  A) Functions calc #############################################################
# calcBlockW
# calcGroupsCoords
# calcThreshold
# calcGroupsW
# calcW
# calcFilter


####>>> calcBlockW ####
calcBlockW <- function(W, site_order = NULL, dist.center = NULL, dist.max = Inf, verbose = optionsMRIaggr("verbose")){
  
  ## test
  # initPackage(package = "spam", method = "calcBlockW")
  # initPackage(package = "Matrix", method = "calcBlockW")
  
  
  if(!is.null(dist.center)){
    site_order <- order(dist.center) - 1
  }else{
    dist.center <- rep(0, length(W@p) - 1)
    if(dist.max < Inf){
      stop("calcBlockW : \'dist.max\' cannot be specified if argument \'dist.center\' is NULL \n")
    }
  }
  
  if(is.null(site_order)){
    site_order <- -1
  }else{
    if(length(site_order) != (length(W@p) - 1)){
      stop("calcBlockW : wrong specification of \'site_order\' \n", 
           "length(site_order) : ", length(site_order), "\n", 
           "length(W@p) - 1 : ", length(W@p) - 1, "\n") 
    }
    
    if(min(site_order) != 0 || max(site_order) != (length(W@p) - 2) ){
      stop("calcBlockW : wrong specification of \'site_order\' \n", 
           "must contains values between 0 and ", (length(W@p) - 2), " \n", 
           "range(site_order) = ", paste(range(site_order), collapse = " "), "\n") 
    }
  }
  
  ## conversion des indices depuic C a R
  res <- calcBlockW_cpp(W_i = W@i, W_p = W@p, site_order = site_order, dist_center = dist.center, dist_max = dist.max, verbose = verbose)
  res$ls_groups <- lapply(res$ls_groups, function(x){x + 1})
  
  ## export
  return(res)
  
}

####>>> calcGroupsCoords ####
calcGroupsCoords <- function(coords, array = NULL, Neighborhood, max_groups = 10000, verbose = optionsMRIaggr("verbose")){
  
  if(is.null(array)){
    if(any(coords < 0)){
      stop("calcGroupsCoords : wrong specification of  \'coords\' \n", 
           "\'coords must\' be positve \n", 
           "min(coords) : ", min(coords), "\n")
    }
    
    if(nrow(coords) == 0){
      return(list(indexArray.group = NULL,
                  indexCoords.group = NULL, 
                  df.group = NULL, 
                  group_size = 0)
      )
    }
    
    if(nrow(coords) == 1){
      return(list(indexArray.group = NULL,
                  indexCoords.group = list(1),  
                  df.group = data.frame(coords, index = 1, group = 1), 
                  group_size = 1)
      )
    }
    
    coords_NNA <- apply(coords, 2, function(x){x-min(x)})
    dim_coordsNNA <- apply(coords_NNA, 2, max) + 1
    p <- length(dim_coordsNNA)
    index_NNA <- rowSums(sweep(coords_NNA, 2, c(1, cumprod(dim_coordsNNA)[-p]), FUN = "*"))
  }else{
    coords_NNA <- which(!is.na(array), arr.ind = TRUE) - 1
    
    if(nrow(coords_NNA) == 0){
      return(list(indexArray.group = NULL,
                  indexCoords.group = NULL, 
                  df.group = NULL, 
                  group_size = 0)
      )
    }
    
    if(nrow(coords_NNA) == 1){
      return(list(ls.group = list(1), 
                  df.group = data.frame(coords_NNA, index = 1, group = 1), 
                  group_size = 1)
      )
    }
    
    dim_coordsNNA <- dim(array)
    p <- length(dim_coordsNNA)
    index_NNA <- which(!is.na(array)) - 1
  }
  
  #### definition de Neighborhood
  if(length(Neighborhood) == 1 && is.character(Neighborhood)){
    Neighborhood <- initNeighborhood(Neighborhood, method = "calcGroupsCoords")
  }
  
  validDimension(value1 = Neighborhood, validDimension = p, name1 = "Neighborhood", names2 = NULL, type = "ncol", method = "calcGroupsCoords")
  
  #### Rcpp
  res_cpp <-  calcGroupsCoords_cpp(coords_NNA = coords_NNA, 
                                   index_NNA = index_NNA, 
                                   min_index_NNA = 0,#index_NNA[1],
                                   max_index_NNA = index_NNA[length(index_NNA)],
                                   Neighborhood = Neighborhood, 
                                   coords_max = dim_coordsNNA, 
                                   max_groups = max_groups, 
                                   verbose = verbose)
  
  if(res_cpp$cv == FALSE){
    warning("calcGroupsCoords : maximum number of groups reached \n", 
            "the identification of the spatial groups may not be complet \n", 
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  #### export  
  indexArray.group <- list()
  for(iter_group in 1:length(res_cpp$group_size)){
    indexArray.group[[iter_group]] <- index_NNA[res_cpp$group == iter_group]
  }
  
  indexCoords.group <- list()
  for(iter_group in 1:length(res_cpp$group_size)){
    indexCoords.group[[iter_group]] <- which(res_cpp$group == iter_group)
  }
  
  if(is.null(array)){
    df.group <- data.frame(coords, group = res_cpp$group)
  }else{
    df.group <- NULL
  }
  
  return(list(indexArray.group = indexArray.group,
              indexCoords.group = indexCoords.group, 
              df.group = df.group, 
              group_size = res_cpp$group_size)
  )
}

####>>> calcThreshold ####
calcThreshold <- function(contrast, param, hemisphere = NULL, rm.CSF = FALSE, threshold = 1:10, decreasing = FALSE, 
                          GRalgo = FALSE, W = NULL, seed = NULL, numeric2logical = FALSE, verbose = optionsMRIaggr("verbose")){
  
  
  #### pre test
  if( optionsMRIaggr("checkArguments")) {
    
    validClass(value = contrast, validClass = c("matrix","data.frame"), superClasses = FALSE, method = "calcThreshold")
    validLogical(value = numeric2logical, validLength = 1, method = "calcThreshold")
    
    if (is.character(seed) == TRUE) {  
      
      validCharacter(seed, validLength = NULL, validValues = names(contrast), method = "calcThreshold")
      sapply(seed,function(x){
        if(numeric2logical == TRUE){
          validNumeric(contrast[,x], name = "seed", validLength = NULL, method = "calcThreshold")
        }else{
          validLogical(contrast[,x], name = "seed", validLength = NULL, method = "calcThreshold")
        }
      })
      
    }
    
  }
  
  #### pre initialization
  #if(GRalgo == TRUE){
  # initPackage(package = "spam", argument = "GRalgo = \"TRUE\"", method = "calcThreshold")
  # initPackage(package = "Matrix", argument = "GRalgo = \"TRUE\"", method = "calcThreshold")
  #}
  
  p <- length(param)
  n <- nrow(contrast)
  if(decreasing == TRUE){
    threshold <- - threshold
    contrast[,param] <- - contrast[,param]
  }
  threshold <- sort(threshold)
  n.threshold <- length(threshold)
  
  if(is.character(rm.CSF) == TRUE){
    param.CSF <- rm.CSF
    rm.CSF <- TRUE
  }else{ 
    param.CSF <- "CSF"
  }
  
  if (is.character(seed) == TRUE) {
    seed <- contrast[,seed, drop = FALSE]
    if(numeric2logical == TRUE){ seed <- apply(seed, 2, as.logical)}
    seed <-  which(rowSums(seed) > 0)
  }
  
  #### test
  if( optionsMRIaggr("checkArguments")) {
    
    validNames(value = param, validValues = names(contrast), method = "calcThreshold")
    
    if (rm.CSF == TRUE && (param.CSF %in% names(contrast) == FALSE)) {
      stop("calcThreshold : wrong specification of \'contrast\' \n", 
           "contrast must be contains a column \"", param.CSF, "\" if rm.CSF = TRUE \n", 
           "column names of \'contrast\' : ", paste(names(contrast), collapse = " "), " \n")
    }
    
    if (!is.null(hemisphere)) {
      
      if ("hemisphere" %in% names(contrast) == FALSE) {
        stop("calcThreshold : wrong specification of \'contrast\' \n", 
             "contrast must be contains a columns \"hemisphere\" if \'hemisphere\' is not NULL \n", 
             "column names of \'contrast\' : ", paste(names(contrast), collapse = " "), " \n")
      }
      
      if (hemisphere %in% unique(contrast$hemisphere) == FALSE) {
        stop("calcThreshold : wrong specification of \'hemisphere\' \n", 
             "\'hemisphere\' does not match element of the hemisphere column in contrast \n", 
             "unique(contrast$hemisphere) : ", unique(contrast$hemisphere), "\n", 
             "proposed \'hemisphere\' : ", hemisphere, " \n")
      }
      
    }
    
    validLogical(value = rm.CSF, validLength = 1, method = "calcThreshold")
    validNumeric(value = threshold, validLength = NULL, refuse.duplicates = TRUE, method = "calcThreshold")
    validLogical(value = decreasing, validLength = 1, method = "calcThreshold")
    validLogical(value = GRalgo, validLength = 1, method = "calcThreshold")
    
    if(GRalgo == TRUE){
      
      validDimension(value1 = W, validDimension = c(n,n), name1 = "W", name2 = "contrast", type = c("nrow", "ncol"), method = "calcThreshold")
      validInteger(value = seed, validLength = NULL, min = 1 , max = n, refuse.duplicates = TRUE, method = "calcThreshold")
      
    }
    
    validLogical(value = verbose, validLength = 1, method = "calcThreshold")
    
  }
  
  #### step 1 - CSF et hemi
  if(verbose == TRUE){cat("Step 1 : ")}
  index.perf <- 1:n
  
  if(!is.null(hemisphere)){
    if(verbose == TRUE){cat("keep only the \"", hemisphere, "\" hemisphere", sep = "")}
    index.perf <- intersect(index.perf, which(contrast$hemisphere %in% hemisphere))    
  } else{
    if(verbose == TRUE){cat("keep both hemipheres")}
  }
  
  if(rm.CSF == TRUE){
    if(verbose == TRUE){cat(", remove CSF ")}
    index.perf <- intersect(index.perf, which(contrast[,param.CSF] < 0.5))
  } else{
    if(verbose == TRUE){cat(", keep CSF ")}
  }
  
  contrast <- contrast[,param, drop = FALSE]
  min_perf <- threshold[1] - 1
  contrast[-index.perf,] <- min_perf
  if(verbose == TRUE){cat("\n")}
  
  #### etape 2 - seuillage
  if(verbose == TRUE){cat("Step 2 : thresholding ")}
  tempo <- matrix(min_perf, nrow = n, ncol = p)
  colnames(tempo) <- param
  
  
  for(iter_threshold in threshold){
    if(verbose == TRUE){cat("*")}
    for(iter_param in 1:p){
      tempo[contrast[,param[iter_param]] >= iter_threshold, iter_param] <- iter_threshold      
    }
  }
  
  contrast <- tempo
  rm(tempo)
  gc()
  if(verbose == TRUE){cat("\n")}
  
  #### etape 3 - GR
  if(GRalgo == TRUE){
    if(verbose == TRUE){cat("Step 3 : Growing Region \n", sep = "")}
    
    for(iter_param in param){
      if(verbose == TRUE){cat(iter_param, " ", sep = "")}
      
      for(iter_threshold in 1:n.threshold){
        if(verbose == TRUE){cat("*")}
        
        index_threshold <- which(contrast[,iter_param] >= threshold[iter_threshold])
        
        contrastBis <- rep(0, n)
        contrastBis[index_threshold] <- 100
        
        resGR <- GRalgo(contrastBis, W = W, seed = seed[seed >= threshold[iter_threshold]], 
                        sigma_max = 0.0001, range = c(threshold[iter_threshold], Inf), breaks = seq(-1, 109, by = 10), step = 10, 
                        operator = "sd", iter_max = 1000, keep.lower = FALSE, keep.upper = FALSE,
                        history.sigma = FALSE, history.step = FALSE, history.front = FALSE)
        
        if(length(setdiff(index_threshold, resGR$GR)) == 0){next}
        
        if(iter_threshold == 1){
          contrast[setdiff(index_threshold, resGR$GR), iter_param] <- min_perf
        }else{
          contrast[setdiff(index_threshold, resGR$GR), iter_param] <- threshold[iter_threshold - 1]
        }
        
      }
      if(verbose == TRUE){cat("\n")}
    }
    if(verbose == TRUE){cat("\n")}
  }
  
  if(decreasing == TRUE){
    contrast <- - contrast  
  }
  names(contrast) <- param
  
  return(as.data.frame(contrast))
}

####>>> calcGroupsW ####
calcGroupsW <- function(W, subset = NULL, max_groups = 10000, verbose = optionsMRIaggr("verbose")){ 
  
  #### tests
  # initPackage(package = "spam", method = "calcGroupsW")
  # initPackage(package = "Matrix", method = "calcGroupsW")
  
  validClass(value = W, validClass = "dgCMatrix", superClasses = TRUE, method = "calcGroupsW")
  
  
  #### initialization
  if(is.null(subset)){
    subset <- -1
  }else{
    subset <- subset - 1
  }
  
  #### call C++ function
  resCpp <- calcGroupsW_cpp(W_i = W@i, W_p = W@p, subset = subset, max_groups = max_groups, verbose = verbose)
  
  if(resCpp$cv == FALSE){
    warning("calcGroupsW : maximum number of groups reached \n", 
            "the identification of the spatial groups may not be complet \n", 
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  #### export
  res <- list(group = resCpp$group, 
              group_subset = resCpp$group_subset, 
              group_size = resCpp$group_size, 
              group_number = resCpp$nb_groups, 
              group_max = which.max(resCpp$group_size)              
  )
  
  return(res)
}

####>>> calcW ####
methods::setMethod(f  = "calcW", 
                   signature  = "data.frame", 
                   definition = function(object, range, method = "euclidean", upper = NULL, format = "dgCMatrix", row.norm = FALSE, 
                                         spatial_res = rep(1, ncol(object)), calcBlockW = FALSE)
                   { 
                     p <- ncol(object)
                     
                     #### tests
                     # initPackage(package = "spam", method = "calcW")
                     # initPackage(package = "Matrix", method = "calcW")
                     
                     validCharacter(value = format, validLength = 1, validValues = c("spam", "dgCMatrix"), refuse.NULL = TRUE, method = "calcW")
                     validNumeric(value = range, validLength = 1, min = 0, refuse.NA = TRUE, refuse.NULL = TRUE, method = "calcW")                     
                     validDimension(value1 = spatial_res, validDimension = p, name2 = NULL, type = "length", method = "calcW")
                     
                     #### initialisation
                     for(iter_c in 1:p){
                       object[,iter_c] <-  object[,iter_c]*spatial_res[iter_c]
                     }
                     
                     #### computing            
                     W <- spam::nearest.dist(object, method = "euclidean", delta = range, upper = upper)
                     
                     if(format == "dgCMatrix"){
                       W <- spam::as.dgCMatrix.spam(W)
                       if(row.norm){
                         pSum <- spam::rowSums(W)
                         pSum[pSum == 0] <- -1
                         W <- W / pSum
                       }
                       W <- Matrix::drop0(W, is.Csparse = TRUE)
                     }      
                     
                     
                     if(calcBlockW == TRUE){
                       dist.center <- sqrt(rowSums(sweep(object, MARGIN = 2, STATS = apply(object, 2, median), FUN = "-")^2))
                       blocks <- calcBlockW(W = W, site_order = NULL, dist.center = dist.center, dist.max = range, verbose = FALSE)
                     }else{
                       blocks <- NULL
                     }
                     
                     #### export
                     return(list(res = list(W = W, blocks = blocks, upper = upper), 
                                 update.object = FALSE, 
                                 overwrite = FALSE)
                     )
                   }
)


####>>> calcFilter ####
methods::setMethod(f  = "calcFilter", 
                   signature  = "array", 
                   definition = function(object, filter, norm.filter = TRUE, bilateral = FALSE, na.rm = FALSE){
                     
                     #### test
                     validLogical(value = norm.filter, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     validLogical(value = bilateral, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     validLogical(value = na.rm, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     
                     if(length(dim(object)) %in% c(2, 3) == FALSE){
                       stop("calcFilter[array] : wrong specification of \'object\' \n", 
                            "\'object\' must be an array of dimension 2 or 3 \n", 
                            "dim(object) : ", paste(dim(object), collapse = " "), "\n")
                     }
                     
                     #### Filtres disponibles 
                     if(length(filter) == 1 && is.character(filter)){
                       
                       filter_split <- strsplit(filter, split = "")[[1]]
                       
                       if(filter_split[[4]] == "N"){
                         
                         res <- initNeighborhood(filter, method = "calcFilter")
                         filter_split <- list(ncol(res), "D", "_", filter_split[[4]], nrow(res))
                         
                         if(filter_split[[1]] == 2){
                           filter <- matrix(NA, nrow = 3, ncol = 3)
                           for(iter in 1:filter_split[[5]]){
                             filter[res[iter, 1] + 2, res[iter, 2] + 2] <- 1  
                           }
                         }else{
                           filter <- array(NA, dim = c(3, 3, 3))
                           for(iter in 1:filter_split[[5]]){
                             filter[res[iter, 1] + 2, res[iter, 2] + 2, res[iter, 3] + 2] <- 1  
                           }
                         }
                         
                       }else{
                         
                         res <- initFilter(filter, method = "calcFilter")
                         filter <- res$filter
                         filter_split <- res$filter_split
                         
                       }
                     }else{ # Perso
                       filter_split <- c(length(dim(filter)), "D", "_", "P", dim(filter)[1])
                     }
                     
                     if(filter_split[[4]] == "S" && na.rm == FALSE){
                       warning("calcFilter[array] : gradient values may be incorrect at the edges \n", 
                               "set \'na.rm\' to TRUE to remove these values \n")
                     }
                     if(filter_split[[4]] == "M" && bilateral == TRUE){
                       warning("calcFilter[array] : there is no edge preserving median filtering \n", 
                               "\'bilateral\' will be ignored \n")
                     }
                     
                     ### preparation
                     p <- dim(filter)
                     p.dim <- length(p)
                     p_ref <- sapply(1:p.dim, function(x){stats::median(1:p[x])})
                     if(p.dim > 3){
                       stop("calcFilter[array] : wrong specification of \'filter\' \n", 
                            "can only handle 1D 2D and 3D filters \n", 
                            "dimension of the proposed filter : ", p.dim, "\n")
                     }
                     
                     M.dim <- dim(object)
                     Mres <- array(NA, dim = dim(object))
                     
                     Ind.operateur_compr <- which(as.vector(filter) != 0)
                     Vec.operateur_compr <- as.vector(filter)[Ind.operateur_compr]
                     
                     # filtrage 2D
                     if(filter_split[[4]] %in% c("N", "G", "I", "P", "S") && length(M.dim) == 2){
                       
                       resCpp <- filtrage2D_cpp(M_data = object, 
                                                M_operateur = filter, 
                                                index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                bilateral = bilateral, 
                                                na_rm = na.rm)
                       
                       if(norm.filter == TRUE){
                         Mres <- resCpp$Mres / resCpp$Wres
                       }else{
                         Mres <- resCpp$Mres
                       }
                       
                     }
                     
                     
                     # filtrage 3D
                     if(filter_split[[4]] %in% c("N", "G", "I", "P", "S") && length(M.dim) == 3){
                       if(filter_split[[1]] == 2){
                         for(iter_k in 1:(M.dim[3])){
                           
                           resCpp <- filtrage2D_cpp(M_data = as.matrix(object[,,iter_k, drop = TRUE]), 
                                                    M_operateur = filter, 
                                                    index_data = which(!is.na(object[,,iter_k, drop = FALSE]), arr.ind = TRUE) - 1, 
                                                    bilateral = bilateral, 
                                                    na_rm = na.rm)
                           
                           if(norm.filter == TRUE){
                             Mres[,, iter_k] <- resCpp$Mres / resCpp$Wres
                           }else{
                             Mres[,, iter_k] <- resCpp$Mres
                           }      
                         }
                       }else{  
                         
                         resCpp <- filtrage3D_cpp(Vec_data = as.vector(object), p_data = dim(object), 
                                                  Vec_operateur = as.vector(filter), p_operateur = dim(filter), 
                                                  index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                  bilateral = bilateral, 
                                                  na_rm = na.rm)                
                         
                         if(norm.filter == TRUE){
                           Mres <- resCpp$Mres / resCpp$Wres
                         }else{
                           Mres <- resCpp$Mres
                         }
                       }
                       
                     }           
                     
                     # filtrage median 2D
                     if(filter_split[[4]] == "M" && length(M.dim) == 2){
                       Mres <- filtrage2Dmed_cpp(M_data = object, 
                                                 M_operateur = filter, 
                                                 index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                 na_rm = na.rm)
                     }
                     
                     # filtrage median 3D
                     if(filter_split[[4]] == "M" && length(M.dim) == 3){
                       if(filter_split[[1]] == 2){
                         for(iter_k in 1:(M.dim[3])){
                           Mres[,, iter_k] <- filtrage2Dmed_cpp(M_data = as.matrix(object[,,iter_k, drop = TRUE]), 
                                                                M_operateur = filter, 
                                                                index_data = which(!is.na(object[,,iter_k, drop = FALSE]), arr.ind = TRUE) - 1, 
                                                                na_rm = na.rm)
                         }
                       }else{
                         
                         Mres <- filtrage3Dmed_cpp(Vec_data = as.vector(object), p_data = dim(object), 
                                                   Vec_operateur = as.vector(filter), p_operateur = dim(filter), 
                                                   index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                   na_rm = na.rm) 
                       }
                     }
                     
                     
                     ### export
                     return(list(res = Mres, 
                                 filter = filter, 
                                 update.object = FALSE, 
                                 overwrite = FALSE)
                     )
                   }
)