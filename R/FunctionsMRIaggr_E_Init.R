#**********************************************************************
#**********************************************************************
#*************         A Functions MRIaggr          *******************
#**********************************************************************
#**********************************************************************
#

#####  A) Functions calc #############################################################
# initCol
# initDisplayWindow
# initFilter
# initIndex
# initNeighborhood
# initPackage
# initPoints
# initWindow

####>>> initCol ####
initCol <- function(object, names.coords, names.param,
                    pch, col, palette, breaks, legend, type.breaks, 
                    checkArguments, method){
  
  # > image is used in all but three cases
  # col has been specifed by the user
  # there are more than one parameter (color have to be defined with a multiparametric palette)
  # user has specified a pch value
  
  n.coords <- length(names.coords)
  n.param <- length(names.param)
  n.px <- nrow(object)
  available_palette <- c("grey.colors", "gray.colors", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors", "tim.colors")
  Object.modif <- FALSE
  
  #### tests ####
  if(checkArguments){
    
    if((n.param %in% 1:3) == FALSE){
      stop(method, " : wrong specification of \'object\' \n", 
           "\'object\' must have between 1 and 3 parameters \n", 
           "number of parameters: ", n.param," \n",
           "names of the parameters: ", paste(names.param, collapse = " ")," \n")    
    }
    
    if(n.param %in% c(2:3)){
      validCharacter(value = palette, validLength = 1, validValues = c("rgb","hsv"), refuse.NULL = TRUE, method = method)
      
      if( object[,any(.SD>1 || .SD<0 ), .SDcols = "TTP_t0"] ){
        
        res_tempo <- object[,lapply(.SD, function(x){round(range(x),1)}), .SDcols = names.param]
        res_tempo <- lapply(1:n.param,
                            function(x){paste0("(",names.param[x],") ",paste(res_tempo[[x]], collapse = " ; "))})
        
        stop(method, " : wrong specification of \'object\' \n", 
             "when using several parameters, they must take values in [0;1] (e.g. be membership probabilities) \n",                   
             "current range (param): ", paste(unlist(res_tempo), collapse = " | "), "\n")
      }
      
    }else{
      
      if(length(palette) == 1 && palette %in% available_palette == FALSE){
        stop(method, " : wrong specification of \'palette\' \n", 
             "available palette for 1 parameter : \"", paste(available_palette, collapse = "\" \""), "\" \n", 
             "proposd palette : ", palette, "\n")
      }
      
      if(length(palette) > 500){
        stop(method, " : wrong specification of \'palette\' \n", 
             "too many colors in palette: it must be bellow 500 \n", 
             "length(palette): ", palette, "\n")
      }
      
      if(length(breaks) == 1){
        
        if(breaks > 500 || breaks < 2){
          stop(method, " : wrong specification of \'breaks\' \n", 
               "if the number of breaks is specified it must be between 2 and 500 \n", 
               "proposed breaks: ", breaks, "\n")
        }
        
        if( any(duplicated(breaks)) || any( abs(sort(breaks) - breaks) > 10^{-16} ) ){
          stop(method, " : the elements of \'breaks\' must be distinct \n", 
               "length(breaks) : ", length(breaks), "\n", 
               "length(unique(breaks)) : ", length(unique(breaks)), "\n")
        }
        
        if(length(palette) != 1 && length(breaks) != (length(palette) + 1)){
          stop(method, " : \'breaks\' and \'palette\' are incompatible \n", 
               "length(palette) must be equal to length(breaks) + 1 \n", 
               "length(palette) : ", length(palette), " \n", 
               "length(breaks) : ", length(breaks), " \n")
        }
        
      } 
      
      
      
    }
    
  }
  
  
  #### dealing with multiple parameters ####
  index_Nduplicated <- NULL
  index_order <- NULL
  
  if(n.param %in% c(2:3)){
    
    Object.modif <- TRUE
    
    if(n.param == 2){      
      
      col <- do.call(palette, # should add grDevices 
                     list(object[[names.param[1]]],object[[names.param[2]]], 1))
      object[,names.param[2] := NULL, with = FALSE]
      
    }else{
      
      col <- do.call(palette, # should add grDevices
                     list(object[[names.param[1]]],object[[names.param[2]]], object[[names.param[3]]]))
      
      # ordering dataset by Tier
      order1 <- order(object[[names.param[1]]], decreasing = TRUE)[seq(1, n.px, length.out = round(n.px / 3))]            
      order2 <- setdiff(order(object[[names.param[2]]], decreasing = TRUE), order1)[seq(1, n.px - round(n.px / 3), length.out = round(n.px / 3))]
      order3 <- setdiff(order(object[[names.param[3]]], decreasing = TRUE), c(order1, order2))[seq(1, n.px - 2 * round(n.px / 3), length.out = n.px - 2 * round(n.px / 3))]
      
      object[order2, names.param[1] := .SD + 1, with = FALSE, .SDcols = names.param[2]]
      object[order3, names.param[1] := .SD + 2, with = FALSE, .SDcols = names.param[3]]
      object[,names.param[2] := NULL, with = FALSE]            
      object[,names.param[3] := NULL, with = FALSE]            
      
      index_Nduplicated <- which(duplicated(object[[names.param[1]]]) == FALSE)
      index_order <- order(object[index_duplicated][[names.param[1]]])
    }
    names.param <- names.param[1]
  }
  
  #### come back to one parameter
  if(is.null(col)){
    
    range.object <- range(object[[names.param]], na.rm = TRUE)
    
        ## only one value  
        if(range.object[1] == range.object[2]){
          if(length(palette == 1) && palette %in% available_palette){
            palette <- rep(do.call(palette,list(1)), 2)  
          }else{
            palette <- rep(palette[1], 2)
          }
          breaks <- c(unique_tempo - 10^{-12}, range.object[1], unique_tempo + 10^{-12})
        }
    
    
        ## infinte values
        if(any(is.infinite(range.object))){
          Object.modif <- TRUE
          
          index.infinite <- which(is.infinite(object[[names.param]]))
          range.object <- range(object[-index.infinite][[names.param]], na.rm = TRUE)
          
          
          ## maxInf
          maxInf <-  max(c(range.object[2], 99999))
          index.infiniteP <- which(object[["TTP_t0"]] > maxInf)
          
          if(length(index.infiniteP)>0){
            object[index.infiniteP, names.param := maxInf,with = FALSE]
            warning(method, " : \'object\' values contains +Inf values \n", 
                    "they are set to ", maxInf, "\n")
          }
          
          ## minInf
          minInf <-  min(c(range.object[2], -99999))
          index.infiniteP <- which(object[["TTP_t0"]] < minInf)
          
          if(length(index.infiniteM)>0){
            object[index.infiniteM, names.param := minInf,with = FALSE]
            warning(method, " : \'object\' values contains -Inf values \n", 
                    "they are set to ", minInf, "\n")
          }
          
          #range.object <- range(object[[names.param]], na.rm = TRUE)
        }
    
    ## breaks 
    if(length(breaks) == 1){
      
      if(length(palette) > 1){breaks <- length(palette) + 1}
      
      breaks_sauve <- switch(type.breaks, 
                             "range" = seq(range.object[1], range.object[2], length.out = breaks), 
                             "range_center" = seq(-max(abs(range.object)), max(abs(range.object)), length.out = breaks), 
                             "quantile" = stats::quantile(object[,1,with = FALSE], probs = seq(0, 1, length.out = breaks))
      )
      breaks_sauve <- unique(breaks_sauve)
      
      breaks <- breaks_sauve
      if(length(breaks) == 1){breaks <- c(breaks - 10^{-12}, breaks, breaks + 10^{-12})}
      breaks[1] <- breaks[1] - 10^{-12}
      breaks[length(breaks)] <- breaks[length(breaks)] + 10^{-12}
      
    }else{
      breaks <- sort(breaks)
      breaks_sauve <- breaks
    }
    
    ## palette
    if(length(palette) == 1 && palette %in% available_palette){
      if(palette %in% c("grey.colors", "gray.colors", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors")){
        palette <- do.call(palette, # should add grDevices::
                                  list(length(breaks) - 1)
        )
      }else if(palette == "tim.colors"){
        palette <- do.call(palette, # should add fields::
                           list(length(breaks) - 1)
        )
      }
    }
    palette_sauve <- palette
    
    ## integrate extreme contrast in range
    if(max(breaks) <= max(object[[names.param]], na.rm = TRUE)){
      breaks[length(breaks)] <- max(object[[names.param]], na.rm = TRUE) + 10^{-12}
    }
    if(min(breaks) >= min(object[[names.param]], na.rm = TRUE)){
      breaks[1] <- min(object[[names.param]], na.rm = TRUE) - 10^{-12}           
    }
    
    if(!is.null(pch) || any(object[, lapply(.SD, is.integer), .SDcols = names.coords] == FALSE) ){
      col <- palette[findInterval(object[[names.param]], breaks, all.inside = TRUE)]
      col[is.na(object[[names.param]])] <- NA
      if(is.null(pch)){pch <- 15}
    }
    
  }else{
    palette_sauve <- col[order(unique(object[[names.param]]))]
    breaks_sauve <- seq(min(object[[names.param]], na.rm = TRUE), 
                        max(object[[names.param]], na.rm = TRUE), 
                        length.out = length(palette_sauve) - 1)
    
    if(any(is.na(object[[names.param]]))){
      col[object[[names.param]]] <- NA
    }
    
    if(is.null(pch)){pch <- 15}
  }
  
  ##
  if(length(unique(breaks_sauve)) == 1){
    
    if(is.null(col)){
      breaks_sauve <- c(breaks_sauve - 10^{-12}, breaks_sauve, breaks_sauve + 10^{-12})
    }else{
      breaks_sauve <- seq(breaks_sauve - 10^{-12}, breaks_sauve + 10^{-12}, length.out = length(unique(col)))
    }
    
  }
  
  #### export ####
  res <- list()
  res$object <- if(Object.modif == TRUE){object}else{NULL}
  res$palette <- palette
  res$breaks <- breaks
  res$col <- col
  res$pch <- pch
  res$palette_sauve <- palette_sauve
  res$breaks_sauve <- breaks_sauve
  res$index_Nduplicated <- index_Nduplicated
  res$index_order <- index_order
  return(res)
  
}

####>>> initDisplayWindow ####
initDisplayWindow <- function(window, filename, path, width, height, scale, res, 
                              mfrow, bg, pty, mar, mgp, n.contrast = 1){
  
  if(window %in% c("png", "eps", "svg", "pdf")){
    switch(window, 
           "eps" = grDevices::postscript(file = paste(path, filename, ".eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
           "svg" = grDevices::svg(filename = paste(path, filename, ".svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
           "png" = grDevices::png(filename = paste(path, filename, ".png", sep = ""), width = width * scale, height = height * scale, res = res), 
           "pdf" = grDevices::pdf(file = paste(path, filename, ".pdf", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE, paper = "special")
    )
  }
  
  if(window == TRUE){grDevices::dev.new()}
  
  if(!is.null(mfrow)){
    
    if(n.contrast == 1){ # cas uniparametrique
      graphics::par(mfrow = mfrow)
      
    }else{ # cas multiparametrique avec legende
      M.layout <- matrix(NA, nrow = mfrow[1], ncol = mfrow[2]+n.contrast - 1)
      M.layout[,1:mfrow[2]] <- 1:prod(mfrow)
      M.layout[1:(mfrow[1] - 1),  - (1:mfrow[2])] <-  M.layout[1:(mfrow[1] - 1), mfrow[2]]    
      M.layout[mfrow[1],  - (1:mfrow[2])] <-  seq(prod(mfrow) + 1, prod(mfrow) + n.contrast - 1)  
      widths.layout <- c(rep(1 / mfrow[2], mfrow[2] - 1), 
                         rep(1 / (mfrow[2] * n.contrast), n.contrast)
      )
      
      graphics::layout(M.layout, widths = widths.layout)
    }
  }
  
  if(!is.null(bg)){graphics::par(bg = bg)}
  if(!is.null(pty)){graphics::par(pty = pty)}
  if(!is.null(mar)){graphics::par(mar = mar)}   
  if(!is.null(mgp)){graphics::par(mgp = mgp)}     
  
}

####>>> initFilter ####
initFilter <- function(filter, method){
  
  filter_name <- as.character(filter)
  filter_split <- as.list(strsplit(filter_name, split = "")[[1]])
  
  #### tests ####
  # 1
  if(filter_split[[1]] %in% c("2", "3") == FALSE){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[1] : \"2\" \"3\" \n", 
         "proposed filter[1] :  ", filter_split[[1]], "\n")
  }
  filter_split[[1]] <- as.numeric(filter_split[[1]])
  
  # 2-3
  if(filter_split[[2]] != "D" || filter_split[[3]] != "_" ){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[2:3] : \"D_\" \n", 
         "proposed filter[2:3] :  ", filter_split[[2]], " ", filter_split[[3]], "\n")
  }
  
  # 4
  if(filter_split[[4]] %in% c("G", "M", "S", "I") == FALSE){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[4] : \"G\", \"M\", \"S\" or \"I\" \n", 
         "proposed filter[4] :  ", filter_split[[4]], "\n")
  }
  
  # 5-6
  if(filter_split[[4]] == "S"){  
    
    if(filter_split[[1]] == 2 && filter_split[[5]] %in% c("x", "y") == FALSE){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\" or \"y\" \n", 
           "proposed letter: ", filter_split[[5]], "\n")
    }
    if(filter_split[[1]] == 3 && filter_split[[5]] %in% c("x", "y", "z") == FALSE){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\", \"y\" or \"z\" \n", 
           "proposed letter: ", filter_split[[5]], "\n")
    }
    
    if(length(filter_split) != 5){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then it must contains only five letters \n", 
           "proposed nomber of letter: ", length(filter_split), "\n")
    }
    
  }else{
    
    if(length(filter_split) == 5){
      if(filter_split[[5]] %in% c("3", "5", "7", "9") == FALSE){
        stop(method, " : wrong specification of \'filter\' \n", 
             "the fifth letter of \'filter\' must correspond to \"3\", \"5\", \"7\", or \"9\" \n", 
             "proposed letter: ", filter_split[[5]], "\n")
      }
      filter_split[[5]] <- as.numeric(filter_split[[5]])
    }else if(length(filter_split) == 6){
      if(filter_split[[5]] %in% c("1", "3", "5", "7", "9") == FALSE){
        stop(method, " : wrong specification of \'filter\' \n", 
             "the fifth letter of \'filter\' be odd \n", 
             "proposed letter: ", filter_split[[5]], "\n")
      }
      
      if(filter_split[[6]] %in% as.character(1:9) == FALSE){
        stop(method, " : wrong specification of \'filter\' \n",            
             "the six letter of \'filter\' must correspond to \"1\", \"2\" ... \"9\" \n", 
             "proposed letter: ", filter_split[[6]], "\n")
      }
      
      filter_split[[5]] <- as.numeric(filter_split[[6]]) + 10 * as.numeric(filter_split[[5]])
      filter_split[[6]] <- NULL
    }else{
      stop(method, " : wrong specification of \'filter\' \n", 
           "\'filter\' must contains 5 or 6 letters \n", 
           "nb of letters proposed : ", length(filter_split), "\n")
    }
    
  }
  
  
  #### creation of the filters ####  
  if(filter_split[[1]] == 2){
    if(filter_split[[4]] == "I"){ # immediate neighborhood
      filter <- matrix(0, nrow = filter_split[[5]], ncol = filter_split[[5]])
      index_1 <- rowSums((which(filter == 0, arr.ind = TRUE) - stats::median(1:filter_split[[5]]))^2) <= filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]] == "G"){ # gaussian
      val_Filter <- stats::dbinom(0:(filter_split[[5]] - 1), filter_split[[5]] - 1, 0.5)
      filter <- matrix(val_Filter, nrow = filter_split[[5]]) %*% matrix(val_Filter, ncol = filter_split[[5]])
    }
    if(filter_split[[4]] == "S"){ # sobel
      if(filter_split[[5]] == "x"){filter <- c(1,2,1) %*% t(c(1,0,-1))}
      if(filter_split[[5]] == "y"){filter <- c(1,0,-1) %*% t(c(1,2,1))}      
    }
    if(filter_split[[4]] == "M"){ # median
      filter <- matrix(1, nrow = filter_split[[5]], ncol = filter_split[[5]])
    }     
  }
  if(filter_split[[1]] == 3){
    if(filter_split[[4]] == "I"){
      filter <- array(0, dim = rep(filter_split[[5]], 3))     
      index_1 <- rowSums((which(filter == 0, arr.ind = TRUE) - stats::median(1:filter_split[[5]]))^2) <= filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]] == "G"){
      val_Filter <- stats::dbinom(0:(filter_split[[5]] - 1), filter_split[[5]] - 1, 0.5)
      filter <- array(val_Filter, dim = filter_split[[5]]) %o% array(val_Filter, dim = filter_split[[5]]) %o% array(val_Filter, dim = filter_split[[5]])
    }
    if(filter_split[[4]] == "S"){      
      if(filter_split[[5]] == "x"){filter <- c(1,2,1) %o% c(1,0,-1) %o% c(0,1,0)}
      if(filter_split[[5]] == "y"){filter <- c(1,0,-1) %o% c(1,2,1) %o% c(0,1,0)}
      if(filter_split[[5]] == "z"){filter <- c(0,1,0) %o% c(1,2,1) %o% c(1,0,-1)}
    }  
    if(filter_split[[4]] == "M"){
      filter <- array(1, dim = rep(filter_split[[5]], 3))  
    }      
  }
  
  return(list(filter = filter, 
              filter_split = filter_split)
  )
}

####>>> initIndex2 ####
initIndex2 <- function(index, slices, names.coords,
                       indexNum = NULL, 
                       outline.default, cex.default, pch.default, col.default, lwd.default, lty.default, filter.default, method){
  
  #### initialization ####
  n.coords <- length(names.coords)
  
  if(is.matrix(index) || is.data.frame(index)){
    index <- list(coords = as.data.table(index))
  }
  if(is.data.table(index)){
    index <- list(coords = index)
  }
  
  #### tests ####
  if(!is.list(index) || "coords" %in% names(index) == FALSE){
    stop(method, " : wrong specification of \'index", num, "\' \n", 
         "\'index", num, "\' it must be a list containing an element named \"coords\" containing the coordinates of the points to display \n", 
         "is(index", num, "): ", paste(is(index), collapse = ""), "\n",
         "names(index", num, ")$coords : ", paste(names(index), collapse = ""), "\n")
  }

  validNames(index$coords, name = paste0("index",indexNum,"$coords"), 
             validLength = length(names.coords), validValues = names.coords, method = method)
  
  #### main
 
  ## subset to the relevant slices 
  if(!is.null(slices)){
    index$coords <- index$coords[, .SD[.SD[[names.coords[n.coords]]] %in% slices] ]
  }
  
  #### outline
  test.IndexOutlineT <- ("outline" %in% names(index) == TRUE) && (index$outline == TRUE)
  test.IndexOutlineF <- ("outline" %in% names(index) == TRUE) && (index$outline == FALSE)
  
  if(test.IndexOutlineT || (outline.default == TRUE && !test.IndexOutlineF)){
    if("filter" %in% names(index) == TRUE){
      filter <- index$filter
    }else{
      filter <- filter.default
    }
    index$coords <-  pointsOutline(index$coords, filter = filter)
  }
  
  ## prepare for plot
  setkeyv(index$coords, names.coords[n.coords])
  
  ## display options
  if( (indexNum != "Line")*("cex" %in% names(index) == FALSE) ){index$cex <- cex.default}
  if( (indexNum != "Line")*("pch" %in% names(index) == FALSE) ){index$pch <- pch.default}
  if( ("col" %in% names(index) == FALSE) ){index$col <- col.default}
  if( (indexNum == "Line")*("lwd" %in% names(index) == FALSE) ){index$lwd <- lwd.default}
  if( (indexNum == "Line")*("lty" %in% names(index) == FALSE) ){index$lty <- lty.default}

  #### export ####
  return(index)
  
}

####>>> initNeighborhood ####
initNeighborhood <- function(Neighborhood, method){
  
  filter_name <- as.character(Neighborhood)
  valid_names <- c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26")
  
  validCharacter(value = Neighborhood, validLength = 1, 
                 validValues = c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26"), 
                 refuse.NULL = TRUE, method = method)
  
  Neighborhood_split <- unlist(strsplit(Neighborhood, split = ""))
  
  p.Neighborhood <- as.numeric(Neighborhood_split[[1]])
  n.Neighborhood <- if(length(Neighborhood_split) == 5){
    as.numeric(Neighborhood_split[5])
  }else{
    sum(as.numeric(Neighborhood_split[5:6]) * c(10, 1))
  }
  
  Neighborhood <- matrix(0, nrow = n.Neighborhood, ncol = p.Neighborhood)
  
  Neighborhood[1:4,1:2] <- rbind(c(-1,0), 
                                 c(0,-1), 
                                 c(1,0), 
                                 c(0,1))
  
  if(n.Neighborhood %in% c(8,10,18,26)){
    Neighborhood[5:8,1:2] <- rbind(c(-1,-1), 
                                   c(1,1), 
                                   c(-1,1), 
                                   c(1,-1)
    )
  }
  
  if(n.Neighborhood %in% c(6,10,18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 1),] <- rbind(c(0,0,1), 
                                                          c(0,0,-1)
    )
  }
  
  if(n.Neighborhood %in% c(18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 7),] <- rbind(c(1,0,1), 
                                                          c(0,1,1), 
                                                          c(-1,0,1), 
                                                          c(0,-1,1), 
                                                          c(1,0,-1), 
                                                          c(0,1,-1), 
                                                          c(-1,0,-1), 
                                                          c(0,-1,-1)
    )
  }
  
  if(n.Neighborhood == 26){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 7),] <- rbind(c(1,1,1), 
                                                          c(-1,1,1), 
                                                          c(-1,-1,1), 
                                                          c(1,-1,1), 
                                                          c(1,1,-1), 
                                                          c(-1,1,-1), 
                                                          c(-1,-1,-1), 
                                                          c(1,-1,-1)
    )
  }
  
  return(Neighborhood)
}

####>>> initPackage ####
initPackage <- function(package, argument = NULL, tryAttach = FALSE, method){
  
  test.package <- requireNamespace(package, quietly = TRUE)
  if(test.package == FALSE){
    stop(method, " : this function ", if(!is.null(argument)){paste("with argument ", argument, " ", sep = "")}, "requires to have installed the ", package, " package to work \n")
  }
  if(tryAttach == TRUE && (paste("package:", package, sep = "") %in% search() == FALSE) ){
    try(attachNamespace(package))
  }
}

####>>> initPoints ####

initPoints <- function(subset = NULL, validValues = NULL, refusedValues = NULL, 
                       checkArguments = optionsMRIaggr("checkArguments"), 
                       method, ...){
  
  ## preparation   
  dots.arguments <- list(...)
  names_dots.arguments <- names(dots.arguments)
  n.dots.arguments <- length(dots.arguments)
  
  if(!is.null(subset)){
    optionsMRIaggr.eg <- optionsMRIaggr(subset, format = "list") 
    dots.arguments <- dots.arguments[names_dots.arguments[names_dots.arguments %in% subset]]
    names_dots.arguments <- names(dots.arguments)
  }else{
    optionsMRIaggr.eg <- optionsMRIaggr() 
  }
  
  ## tests
  if(checkArguments && n.dots.arguments > 0){
    
    if(is.null(validValues)){
      validValues <- names(optionsMRIaggr.eg)
    }
    
    validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                   validValues = setdiff(validValues, refusedValues),
                   refuse.NULL = FALSE, method = method)
    
    validOptionsMRIaggr(dots.arguments, method = method)
  }
  
  ## set specific display
  if(n.dots.arguments > 0){
    optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
  }
  
  ## export
  return(optionsMRIaggr.eg)
  
}

####>>> initWindow ####
initWindow <- function(window, filename, path, width, height, unit, res, 
                       n.plot, mfrow, xlim, ylim, checkArguments, method){
  
  #### tests ####
  if(checkArguments){
    validCharacter(window, validLength = 1, validValues = c(TRUE,FALSE,"png","eps","svg","pdf"), refuse.NULL = FALSE, method = method)
    
    validCharacter(filename, validLength = 1, method = method)
    
    validNumeric(width, validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
    
    validNumeric(height, validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
    
    validCharacter(path, validLength = 1, refuse.NULL = FALSE, method = method)
    
    validPath(path, method = method)
    
    validCharacter(unit, validLength = 1, validValues = c("px","in","cm","mm"), refuse.NULL = TRUE, method = method)
  }
  
  #### initialization ####
  scale <- switch(unit, 
                  "px" = 1, 
                  "in" = 90, 
                  "cm" = 35.43, 
                  "mm" = 3.543)
  
  if(is.null(mfrow)){
    if(n.plot <= 9){
      mfrow <- c(ceiling(n.plot / ceiling(sqrt(n.plot))), ceiling(sqrt(n.plot)))
    }else{
      mfrow <- c(3, 3)
    }              
  }else{
    mfrow <- mfrow
  }
  
  n.graph_par_window <- prod(mfrow)
  
  if(!is.null(xlim)){xlim.plot <- xlim}else{xlim.plot <- NULL}
  if(!is.null(ylim)){ylim.plot <- ylim}else{ylim.plot <- NULL}
  
  #### export ####
  res <- list()
  res$scale <- scale
  res$mfrow <- mfrow
  res$n.graph_par_window <- n.graph_par_window
  res$xlim.plot <- xlim.plot
  res$ylim.plot <- ylim.plot
  
  return(res)  
  
}
