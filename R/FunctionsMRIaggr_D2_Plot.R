#**********************************************************************
#**********************************************************************
#*************         B Plot MRIaggr               *******************
#**********************************************************************
#**********************************************************************
#

#####  B) Functions Plot #############################################################
# legendMRI
# legendMRI2
# multiplot   ### initIndex2 a mettre a jour
# outline
# plotMRI
# pointsOutline

####>>> legendMRI ####
legendMRI <- function(breaks, palette, mar, cex, cex.main, main, quantiles, digit){
  
  ### definition de scale 
  if(length(breaks) > 8){
    at_legend <- seq(min(breaks), max(breaks), length.out = 8)
  }else{
    at_legend <- breaks
  }
  at_legend <- signif(at_legend, digit)
  n.legend <- length(at_legend)
  n.breaks <- length(breaks)
  seq_scale <- 1:(n.breaks - 1)
  
  if(!is.null(mar)){graphics::par(mar = mar)}
  graphics::image(1, (breaks[-n.breaks] + breaks[-1]) / 2, rbind(seq_scale), col = palette, 
                  axes = FALSE, ylab = "", xlab = "")
  
  graphics::title(main, cex.main = cex.main)
  graphics::axis(2, cex.axis = cex, at = c(at_legend[-n.legend], 0.99 * at_legend[n.legend]), labels = at_legend, las = 2)            
  
  if(!identical(quantiles, FALSE)){
    graphics::abline(h = quantiles, col = c("red", "blue", "green", "blue", "red"), lwd = 3, lty = c(1, 2, 3, 2, 1))
  }
  
  # export
  return(invisible(TRUE))
}

####>>> legendMRI2 ####
legendMRI2 <- function(param, palette, mar, cex, cex.main){
  
  ### definition de scale 
  if(!is.null(mar)){graphics::par(mar = mar)}
  n.param <- length(param)
  seq_FALSE <- 0
  seq_TRUE <- seq(0, 1, 0.1)
  
  for(iter_param in 1:n.param){
    eval(parse(text = paste(
      "graphics::image(1, seq(0, 1, 0.1), rbind(1:10), col = ", palette, "(seq_", iter_param == 1, ", seq_", iter_param == 2, ", seq_", iter_param == 3, "), 
      axes = FALSE, ylab = \"\", xlab = \"\")", 
      sep = "")))
    graphics::title(param[iter_param], cex.main = cex.main)
    graphics::axis(2, cex.axis = cex, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), las = 2)            
  }
  
  # export
  return(invisible(TRUE))
}

####>>> multiplot ####
methods::setMethod(f  = "multiplot", 
                   signature  = "data.frame", 
                   definition = function(object, contrast = NULL, param = NULL, slices = NULL, 
                                         index1 = NULL, index2 = NULL, index3 = NULL, indexLine = NULL,
                                         col = NULL, pch = NULL, xlim = NULL, ylim = NULL, 
                                         main = NULL, main.legend = NULL, filename = "multiplotDF", ...){
                    
                     optionsMRIaggr.eg <- initPoints(subset = "checkArguments", ..., method = "multiplot[data.frame]")
                     
                     if(is.vector(contrast) || is.matrix(contrast)){
                       contrast <- as.data.frame(contrast)
                     }

                     #### test ####
                     if(optionsMRIaggr.eg$checkArguments){
                       
                       if(!is.null(contrast) && is.data.frame(contrast) == FALSE){
                         stop("multiplot[data.frame] : wrong specification of \'contrast\' \n", 
                              "must be either NULL or data.frame or matrix or vector \n", 
                              "is(contrast) : ", paste(is(contrast), collapse = " "), "\n")
                       }
                       
                       validDimension(value1 = contrast, value2 = object, name1 = "contrast", name2 = "object", type = "nrow", method ="multiplot[data.frame]")
                        
                     }
                    
                     
                     #### conversion
                     names.coords <- names(object)
                     object <- as.data.table(cbind(object,contrast))
                     setkeyv(object,names.coords)
                     
                     #### main
                     res <- multiplot(object, param = param, slices = slices, 
                                      index1 = index1, index2 = index2, index3 = index3,
                                      col = col, pch = pch, xlim = xlim, ylim = ylim, 
                                      main = main, main.legend = main.legend, filename = filename, 
                                      ...)
                     
                     #### export 
                     return(invisible(res))
                     
                   }
)

methods::setMethod(f  = "multiplot", 
                   signature  = "data.table", 
                   definition = function(object, param = NULL, slices = NULL, 
                                         index1 = NULL, index2 = NULL, index3 = NULL, indexLine = NULL,
                                         col = NULL, pch = NULL, xlim = NULL, ylim = NULL, 
                                         main = NULL, main.legend = NULL, noPlot = FALSE, filename = "multiplotDT", ...){
                     
                     names.coords <- key(object)
                     if(is.null(param)){
                       names.param <- setdiff(names(object), names.coords)
                     }else{
                       names.param <- param
                     }
                     n.param <- length(names.param)
                     
                     if(!is.null(index3)){
                       test.intersection1 <- length(index3) == 1 && index3 == "intersection"
                       test.intersection2 <- is.list(index3) && "coords" %in% names(index3) && length(index3$coords) == 1 && index3$coords == "intersection"
                     }
                     
                     #### graphical options
                     optionsMRIaggr.eg <- initPoints(..., method = "multiplot[data.table]")
                     checkArguments <- optionsMRIaggr.eg$checkArguments
                     legend <- optionsMRIaggr.eg$legend
                     window <- optionsMRIaggr.eg$window
                     
                     #### test ####
                     if(checkArguments){
                      
                       if(!is.null(slices) && any(slices %in% unique(object[[names.coords[3]]]) == FALSE)){
                         stop("multiplot[data.frame] : \'slices\' inconsistent with \'object\' \n", 
                              "requested slices : \"", paste(slices, collapse = "\" \""), "\" \n", 
                              "slices missing in \'object\' : \"", paste(slices[slices %in% unique(object[[names.coords[3]]]) == FALSE], collapse = " "), "\" \n")
                       }
                       
                       if(!is.null(param) && any(param %in% names(object) == FALSE)){
                         stop("multiplot[data.frame] : \'param\' inconsistent with \'object\' \n", 
                              "requested parameters no found in \'object\': \"", paste(param[param %in% setdiff(names(object), names.coords) == FALSE], collapse = "\" \""), "\" \n", 
                              "parameters in \'object\': \"", paste(setdiff(names(object), names.coords), collapse = "\" \""), "\" \n",
                              "coordinates in \'object\': \"", paste(names.coords, collapse = "\" \""), "\" \n")
                       }
                       
                       if(length(names.param) > 3){
                         stop("multiplot[data.frame] : \'object\' must contains at most 3 parameters \n", 
                              "number of parameters : ", length(names.param), "\n",
                              "names of the parameters : ", paste( names.param, collapse = " "), "\n")
                       }
                       
                       if(length(names.coords) %in% 2:3 == FALSE){
                         stop("multiplot[data.frame] : multiplot can only handle 2 or 3 dimensional datasets \n", 
                              "number of coordinates : ", length(names.coords), "\n",
                              "names of the coordinates : ", paste( names.coords, collapse = " "), "\n")
                       }
                       
                       if(any(c(names.coords,names.param) %in% c("MRIaggr.colXX","MRIaggr.coordsXX"))){
                         stop("multiplot[data.frame] : \'object\' must not contain a variable named \"MRIaggr.colXX\" or \"MRIaggr.coordsXX\" \n")
                       }
                       
                       if(!is.null(index3) && (test.intersection1 || test.intersection2)){
                           
                           if(is.null(index1)){
                             stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index1\' is missing  \n")
                           }
                           if(is.null(index2)){
                             stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index2\' is missing  \n")
                           }
                       }
                     }
                     
                     #### initialisation
                     
                     # deal with 2D dataset
                     if(length(names.coords) == 2){
                       object[,MRIaggr.coordsXX := 1]
                       names.coords <- c(names.coords,"MRIaggr.coordsXX")
                       # setkeyv(object, "MRIaggr.coordsXX")
                     }
                    
                     # order the coordinates
                     # setkeyv(object, names.coords[3])
                     
                     # slices
                     if(is.null(slices)){
                       slices <- unique(object[[names.coords[3]]])
                     }
                     n.slices <- length(slices)
                     
                     if(is.null(main.legend)){
                       main.legend <- paste(names.param, collapse = "-")
                     }
                     
                     # legend
                     if(!is.null(legend)){
                       if(legend == FALSE){n.plot <- n.slices}
                       if(legend == TRUE){n.plot <- n.slices + 1}
                       if(legend == "only"){n.plot <- 1}
                     }else{
                       n.plot <- n.slices
                     }

                     ## windows
                     par.init <- graphics::par()
                     
                     res.init <- initWindow(window = optionsMRIaggr.eg$window, filename = filename, 
                                            path = optionsMRIaggr.eg$path, width = optionsMRIaggr.eg$width, height = optionsMRIaggr.eg$height, unit = optionsMRIaggr.eg$unit, res = optionsMRIaggr.eg$res, 
                                            n.plot = n.plot, mfrow = optionsMRIaggr.eg$mfrow, xlim = xlim, ylim = ylim, 
                                            checkArguments = checkArguments, method = "multiplot[data.table]")
                     scale <- res.init$scale
                     mfrow <- res.init$mfrow
                     n.graph_par_window <- res.init$n.graph_par_window
                     xlim.plot <- res.init$xlim.plot
                     ylim.plot <- res.init$ylim.plot
                     
                     ## color and breaks 
                     res.init <- initCol(object = object, names.coords = names.coords,  names.param = names.param,
                                         pch = pch, col = col, 
                                         palette = optionsMRIaggr.eg$palette, breaks = optionsMRIaggr.eg$breaks, legend = legend, type.breaks = optionsMRIaggr.eg$type.breaks, 
                                         checkArguments = checkArguments, method = "multiplot[data.table]")
                     
                     if(!is.null(res.init$object)){
                     object <- res.init$object
                     }
                     palette <- res.init$palette 
                     breaks <- res.init$breaks
                     
                     if(!is.null(res.init$object)){
                       object[,MRIaggr.colXX := res.init$col]
                       col <- TRUE  
                     }else{
                       col <- FALSE  
                     }
                     
                     pch <- res.init$pch
                     palette_sauve <- res.init$palette_sauve
                     breaks_sauve <- res.init$breaks_sauve
                     index_duplicated <- res.init$index_duplicated
                     index_order <- res.init$index_order
                     
                     ## index            
                     if(!is.null(index1)){
                       index1 <- initIndex2(index = index1, slices = slices, names.coords = names.coords,
                                             indexNum = 1, 
                                             outline.default = optionsMRIaggr.eg$outline.index, cex.default = optionsMRIaggr.eg$cex.index[1], pch.default = optionsMRIaggr.eg$pch.index[1], 
                                             col.default = optionsMRIaggr.eg$col.index[1], filter.default = optionsMRIaggr.eg$filter.index,
                                             method = "multiplot[data.frame]")
                     }
                     
                     if(!is.null(index2)){
                       index2 <- initIndex2(index = index2, slices = slices, names.coords = names.coords,
                                             indexNum = 2, 
                                             outline.default = optionsMRIaggr.eg$outline.index, cex.default = optionsMRIaggr.eg$cex.index[2], pch.default = optionsMRIaggr.eg$pch.index[2], 
                                             col.default = optionsMRIaggr.eg$col.index[2], filter.default = optionsMRIaggr.eg$filter.index, 
                                             method = "multiplot[data.frame]")
                     }
                     
                     if(!is.null(index3)){
                       
                       if(test.intersection1 || test.intersection2){
                         index3$coords <- merge(index1$coords, index2$coords, all = TRUE)
                       }
                       
                       index3 <- initIndex2(index = index3, slices = slices, names.coords = names.coords,
                                             indexNum = 3, 
                                             outline.default = optionsMRIaggr.eg$outline.index, cex.default = optionsMRIaggr.eg$cex.index[3], pch.default = optionsMRIaggr.eg$pch.index[3], 
                                             col.default = optionsMRIaggr.eg$col.index[3], filter.default = optionsMRIaggr.eg$filter.index, 
                                             method = "multiplot[data.frame]")
                     }
                     
                     if(!is.null(indexLine)){
                       indexLine <- initIndex2(index = indexLine, slices = slices, names.coords = names.coords,
                                               indexNum = "Line", 
                                               outline.default = FALSE, lty.default = 2, lwd.default = optionsMRIaggr.eg$lwd.midplane, 
                                               method = "multiplot[data.frame]")
                     }
                     
                     #### NO PLOT ####
                     if(noPlot == TRUE){
                       return(list(names.param = names.param,
                                   names.coords = names.coords,
                                   mainPlot = list(object = object,
                                                   breaks.plot = breaks, 
                                                   palette.plot = palette,
                                                   optionsMRIaggr.eg = optionsMRIaggr.eg,
                                                   xlim = xlim.plot,
                                                   ylim = ylim.plot,
                                                   col = col,
                                                   pch = pch, 
                                                   main = main,
                                                   legend = legend),
                                   legendPlot = list(breaks_sauve = breaks_sauve, 
                                                     palette_sauve = palette_sauve, 
                                                     n.param = n.param),
                                   indexLine =  if(!is.null(indexLine)){
                                     indexLine
                                   }else{NULL},
                                   index1 =  if(!is.null(index1)){
                                     index1
                                   }else{NULL},
                                   index2 =  if(!is.null(index2)){
                                     index2
                                   }else{NULL},
                                   index3 =  if(!is.null(index3)){
                                     index3
                                   }else{NULL}
                       ))
                     }
                     
                        
                     
                     #### display plot ####
                     compteur <- 1
                     
                     if(!identical(legend, "only")){
                       for(iter_slice in 1:n.slices){
                         
                         # device
                         if(!is.null(window) && compteur == 1){
                           
                           if(window %in% c("png", "eps", "svg", "pdf") && iter_num > 1){grDevices::dev.off()}
                           slice_max_tempo <- min(max(slices), slices[iter_slice + n.graph_par_window - 1], na.rm = TRUE)
                           filename_all <- paste(filename, "_", paste(names.param, collapse = "-"),
                                                 "(",names.coords[3],"=", slices[iter_slice], "-", slice_max_tempo, ")", sep = "")                  
                           
                           initDisplayWindow(window = window, filename = filename_all, 
                                             path = optionsMRIaggr.eg$path, width = optionsMRIaggr.eg$width, height = optionsMRIaggr.eg$height, scale = optionsMRIaggr.eg$scale, res = optionsMRIaggr.eg$res, 
                                             mfrow = mfrow, bg = optionsMRIaggr.eg$bg, pty = optionsMRIaggr.eg$pty, mar = optionsMRIaggr.eg$mar, mgp = optionsMRIaggr.eg$mgp, 
                                             n.contrast = if(identical(legend, TRUE) && ((n.slices-iter_slice) < prod(mfrow) ) ){n.param}else{1})
                         }
                         
                         # data
                         indexK <- object[,.I[.SD == slices[iter_slice]] , .SDcols = names.coords[3]]
                         
                         if(optionsMRIaggr.eg$num.main == TRUE){
                           mainK <- paste0(main, 
                                          paste0(" ",names.coords[3],"=",slices[iter_slice]), 
                                          paste0(" (", paste(names.param, collapse = "-"), ")")
                                          )
                         }else{
                           mainK <- main
                         }
                         
                         
                         plot.test <- plotMRI(contrast = object[indexK, names.param[1], with = FALSE], 
                                              coords = object[indexK, names.coords[1:2], with = FALSE], 
                                              breaks = breaks, col = if(col){object[indexK][["MRIaggr.colXX"]]}else{NULL}, palette = palette, 
                                              asp = optionsMRIaggr.eg$asp, 
                                              xlim = xlim, ylim = ylim, pch = pch, cex = optionsMRIaggr.eg$cex, axes = optionsMRIaggr.eg$axes, col.NA = optionsMRIaggr.eg$col.NA, pch.NA = optionsMRIaggr.eg$pch.NA, 
                                              xlab = optionsMRIaggr.eg$xlab, ylab = optionsMRIaggr.eg$ylab, main = mainK, cex.main = optionsMRIaggr.eg$cex.main)
                         
                         if(!is.null(indexLine) && indexLine$coords[J(slices[iter_slice]),.N]>0 ){
                           graphics::points(indexLine$coords[J(slices[iter_slice])], 
                                            type = "l", col = col_indexLine, lwd = lwd_indexLine, lty = lty_indexLine)
                         }
                         
                         if(!is.null(index1) && index1$coords[J(slices[iter_slice]),.N]>0  ){
                           graphics::points(index1$coords[J(slices[iter_slice]), .SD, .SDcols = names.coords[1:2]], 
                                            pch = index1$pch, cex = index1$cex, col = index1$col)
                         }
                         
                         if(!is.null(index2) && index2$coords[J(slices[iter_slice]),.N]>0 ){
                           graphics::points(index2$coords[J(slices[iter_slice]), .SD, .SDcols = names.coords[1:2]], 
                                            pch = index2$pch, cex = index2$cex, col = index2$col)
                         }
                         
                         if(!is.null(index3) && index3$coords[J(slices[iter_slice]),.N]>0 ){
                           graphics::points(index3$coords[J(slices[iter_slice]), .SD, .SDcols = names.coords[1:2]], 
                                            pch = index3$pch, cex = index3$cex, col = index3$col)
                         }
                         
                         compteur <- compteur + 1
                         if(compteur > n.graph_par_window)
                         {compteur <- 1}
                         
                       }            
                     }
                     
                     #### legend ####
                     if(is.null(legend)){
                       compteur <- 1
                       mfrow <- c(1, 1)
                       legend <- TRUE
                     }
                     
                     if(legend %in% c(TRUE, "only")){
                       
                       if(!is.null(window) && compteur == 1){
                         
                         if(window %in% c("png", "eps", "svg", "pdf") && iter_slice > 1){grDevices::dev.off()}
                         filename_all <- paste(filename, "_", paste(names.param, collapse = "-"), 
                                               "(",names.coords[3],"=", min(slices), "-", max(slices), ") - legend", sep = "")
                         
                         initDisplayWindow(window = window, filename = filename_all, 
                                           path = optionsMRIaggr.eg$path, width = optionsMRIaggr.eg$width, height = optionsMRIaggr.eg$height, scale = optionsMRIaggr.eg$scale, res = optionsMRIaggr.eg$res,
                                           mfrow = mfrow, bg = optionsMRIaggr.eg$bg, pty = optionsMRIaggr.eg$pty, mar = NULL, mgp = optionsMRIaggr.eg$mgp, n.contrast = n.param)
                         
                       }
                       
                       if(n.param == 1){
                         if(optionsMRIaggr.eg$quantiles.legend == TRUE){
                           quantiles.legend <- stats::quantile(object[[names.param]], na.rm = TRUE)
                         }else{
                           quantiles.legend <- quantiles.legend
                         }   
                         
                         plot.test <- legendMRI(breaks = breaks_sauve, palette = palette_sauve, 
                                                mar = optionsMRIaggr.eg$mar.legend, cex = optionsMRIaggr.eg$cex.legend, 
                                                main = paste(main.legend, collapse = " "), 
                                                cex.main = optionsMRIaggr.eg$cex.main, quantiles = quantiles.legend, digit = optionsMRIaggr.eg$digit.legend)
                         
                       }else{
                         plot.test <- legendMRI2(param = names.param, palette = palette, 
                                                 mar = optionsMRIaggr.eg$mar.legend, cex = optionsMRIaggr.eg$cex.legend, cex.main = optionsMRIaggr.eg$cex.main)
                         
                       }              
                     }
                     
                     
                     
                     if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                       grDevices::dev.off()
                     }
                     graphics::par(bg = par.init$bg, pty = par.init$pty, mar = par.init$mar, mgp = par.init$mgp)
                     
                     #### export 
                     return(invisible(list(breaks.plot = breaks, 
                                           palette.plot = palette, 
                                           breaks.legend = breaks_sauve, 
                                           palette.legend = palette_sauve, 
                                           quantiles.legend = if(legend != FALSE){quantiles.legend}else{NULL})))
                   }
) 

####>>> outline ####
outline <- function(n = 50, sequential = TRUE, min_dist = 1, 
                    col = c("blue", "red", "grey"), pch = 20, cex = c(0.75, 1, 0.75)){
  
  # initPackage(package = "RANN", method = "outline")
  
  #### one shot outline 
  if(sequential == FALSE){
    res_locator <- graphics::locator(n = n)
    res_locator$x <- c(res_locator$x, res_locator$x[1])
    res_locator$y <- c(res_locator$y, res_locator$y[1])
    n <- length(res_locator$x)
    
    # round
    res_locator$x <- round(res_locator$x)
    res_locator$y <- round(res_locator$y)  
    
    
    df_points <- matrix(NA, nrow = 0, ncol = 4)
    
    for(iter_point in 1:(n - 1)){
      n_tempo <- max(1 + abs(res_locator$x[iter_point]-res_locator$x[iter_point + 1]), 
                     1 + abs(res_locator$y[iter_point]-res_locator$y[iter_point + 1]))
      n_tempo <- ceiling(n_tempo)
      i_tempo <- seq(res_locator$x[iter_point], res_locator$x[iter_point + 1], length.out = n_tempo)
      j_tempo <- seq(res_locator$y[iter_point], res_locator$y[iter_point + 1], length.out = n_tempo)
      
      i_tempo <- round(i_tempo)
      j_tempo <- round(j_tempo)    
      
      if(n_tempo > 1){
        df_points <- rbind(df_points, 
                           cbind(i = i_tempo[-n_tempo], j = j_tempo[-n_tempo], edge = iter_point, points = c(1, rep(0, n_tempo-2))))
      }
    }
    
  }
  
  
  #### sequential outline
  if(sequential == TRUE){
    
    iter <- 1
    dist <- min_dist + 1
    df_points <- matrix(NA, nrow = 0, ncol = 4)
    
    while(iter <= n && dist > min_dist){
      res_locator <- graphics::locator(n = 1)
      
      # if user enter echap directly
      if(is.null(res_locator)){
        return(list(edge = NULL, 
                    surface = NULL)
        )
      }
      
      if(iter > 1){
        dist <- sqrt( (df_points[1, "i"] - res_locator$x)^2 + (df_points[1, "j"] - res_locator$y)^2)
        
        res_locator$x <- c(df_points[nrow(df_points), "i"], round(res_locator$x))
        res_locator$y <- c(df_points[nrow(df_points), "j"], round(res_locator$y))
        
        n_tempo <- max(1 + abs(res_locator$x[1] - res_locator$x[2]), 
                       1 + abs(res_locator$y[1] - res_locator$y[2]))
        n_tempo <- ceiling(n_tempo)
        if(n_tempo == 1){i_tempo <- c() ; next} # cas ou l on selectionne le meme points
        
        i_tempo <- round(seq(res_locator$x[1], res_locator$x[2], length.out = n_tempo)[-1])
        j_tempo <- round(seq(res_locator$y[1], res_locator$y[2], length.out = n_tempo)[-1])
        
        points <- c(rep(0, n_tempo - 2), 1)
      }else{
        
        # round
        res_locator$x <- round(res_locator$x)
        res_locator$y <- round(res_locator$y)  
        
        i_tempo <- res_locator$x
        j_tempo <- res_locator$y
        points <- 1  
      }
      
      if(dist < 1){
        if(length(i_tempo) > 1){
          df_points <- rbind(df_points, 
                             cbind(i = i_tempo[ - (n_tempo - 1)], j = j_tempo[ - (n_tempo - 1)], edge = iter, points = 0))
        }
      }else{
        df_points <- rbind(df_points, 
                           cbind(i = i_tempo, j = j_tempo, edge = iter, points = points))
      }
      
      graphics::points(i_tempo, j_tempo, 
                       col = col[points + 1], pch = pch, cex = cex[points + 1])
      
      
      
      iter <- iter + 1
    }
    
    # relier au premier point   
    n_tempo <- max(1 + abs(df_points[nrow(df_points),"i"] - df_points[1,"i"]), 
                   1 + abs(df_points[nrow(df_points),"j"] - df_points[1,"j"]))
    n_tempo <- ceiling(n_tempo)
    
    if(n_tempo > 2){
      i_tempo <- round(seq(df_points[nrow(df_points),"i"], df_points[1,"i"], length.out = n_tempo))
      j_tempo <- round(seq(df_points[nrow(df_points),"j"], df_points[1,"j"], length.out = n_tempo))
      
      df_points <- rbind(df_points, 
                         cbind(i = i_tempo[c(-1, -n_tempo)], j = j_tempo[c(-1, -n_tempo)], edge = iter, points = 0))
      
      graphics::points(i_tempo[c(-1, -n_tempo)], j_tempo[c(-1, -n_tempo)], 
                       col = col[points + 1], pch = pch, cex = cex[points + 1])
    }
  }
  df_points <- as.data.frame(df_points)
  
  # enlever d eventuels dupliques
  df_points <- df_points[duplicated(df_points[,c("i", "j")]) == FALSE,]
  
  
  #### filling the outline
  df_points <- cbind(df_points, valid = TRUE)
  j_level <- unique(df_points[,"j"])
  i_tempo <- c()
  j_tempo <- c()  
  
  for(iter_j in 1:length(j_level)){
    
    index_j <- which(df_points[,"j"] == j_level[iter_j])
    
    if(length(index_j) == 1){
      i_tempo <-  c(i_tempo, df_points[index_j, "i"])
      j_tempo <-  c(j_tempo, df_points[index_j, "j"])
    }else{
      index_j <- index_j[order(df_points[index_j, "i"])]
      
      # enlever les doublons d un meme trait
      diff_tempo <- c(FALSE, abs(diff(index_j)) %in% c(1, nrow(df_points) - 1))
      
      while(TRUE %in% diff_tempo && length(diff_tempo) > 2){   
        pos_FALSE <- which(diff_tempo == TRUE)[1]
        if(pos_FALSE %% 2 == 0){
          df_points$valid[index_j[pos_FALSE]] <- FALSE
          index_j <- index_j[-pos_FALSE]          
        }else{          
          df_points$valid[index_j[pos_FALSE - 1]] <- FALSE
          index_j <- index_j[ - (pos_FALSE - 1)]
        }
        diff_tempo <- diff_tempo[-pos_FALSE]
      }
      
      # enlever les V inverses
      df_Vpoints <- cbind(index = 1:nrow(df_points), df_points)[df_points$valid == TRUE,]
      index_Vj <- sapply(index_j, function(x){which(df_Vpoints$index == x)})
      index_m1 <- c(nrow(df_Vpoints), 1:nrow(df_Vpoints) - 1)
      index_p1 <- c(2:nrow(df_Vpoints), 1)
      
      diff_tempo <- rbind(df_Vpoints[,"i"]-df_Vpoints[index_m1,"i"], 
                          df_Vpoints[,"j"]-df_Vpoints[index_m1,"j"], 
                          df_Vpoints[,"i"]-df_Vpoints[index_p1,"i"], 
                          df_Vpoints[,"j"]-df_Vpoints[index_p1,"j"])[,index_Vj, drop = FALSE]
      
      test_V <- apply(diff_tempo, 2, function(x){
        test.H <- all(x == c(1,-1,-1,-1)) + all(x == c(0,-1,-1,-1)) + all(x == c(0,-1,0,-1)) + all(x == c(1,-1,0,-1))
        test.antiH <- all(x == c(-1,-1,1,-1)) + all(x == c(-1,-1,0,-1)) + all(x == c(0,-1,1,-1))
        return(test.H + test.antiH)       
      }) 
      
      test_Vinv <- apply(diff_tempo, 2, function(x){
        test.H <- all(x == c(1,1,-1,1)) + all(x == c(0,1,-1,1)) + all(x == c(0,1,0,1)) + all(x == c(1,1,0,1))
        test.antiH <- all(x == c(-1,1,1,1)) + all(x == c(-1,1,0,1)) + all(x == c(0,1,1,1))
        return(test.H + test.antiH)        
      }) 
      
      i_tempo <- c(i_tempo, df_points[index_j[test_V == 1], "i"])
      j_tempo <- c(j_tempo, df_points[index_j[test_V == 1], "j"])
      index_j <- index_j[c(test_Vinv + test_V) == 0]
      
      n.i <- length(index_j)
      
      # remplissage
      if(n.i > 0){
        for(iter_i in 1:floor(n.i / 2)){
          seq_i <- seq(df_points[index_j[2 * (iter_i - 1) + 1], "i"], df_points[index_j[2 * iter_i], "i"], by = 1)
          i_tempo <-  c(i_tempo, seq_i)
          j_tempo <-  c(j_tempo, rep(j_level[iter_j], length(seq_i)))
        }
      }
      
    }
    
  }
  
  df_fill <- data.frame(i = i_tempo, j = j_tempo)
  
  if(sequential == TRUE){
    test <- RANN::nn2(data = df_points[,c("i", "j")], query = df_fill[,c("i", "j")], k = 1)
    graphics::points(i_tempo[test$nn.dists > 0], j_tempo[test$nn.dists > 0], 
                     col = col[3], pch = pch, cex = cex[3])
  }
  
  
  return(list(edge = df_points, 
              surface = df_fill)
  )
}

####>>> plotMRI ####
plotMRI <- function(contrast, coords, breaks, palette, col, asp, 
                    xlim, ylim, pch, cex, axes, col.NA, pch.NA, xlab, ylab, main, cex.main){
  
  ### specific case
  if(NROW(contrast) == 0 || all(is.na(contrast)) || is.null(contrast)){
    
    graphics::plot(1, 1, col = col, xlab = "", ylab = "", axes = axes)
    graphics::legend("center", c("Missing", "data"), cex = 2, pch = 62, bty = "n")
    
    graphics::title(main = main)
    graphics::axis(1)
    graphics::axis(2)
    
    return(invisible(FALSE))
  }
  
  ####
  
  #### preparation
  ## xlim-ylim
  if(is.null(xlim)){                
    xlim <- c(min(coords[[1]]) - 0.5, max(coords[[1]]) + 0.5)
  }
  
  if(is.null(ylim)){
    ylim <- c(min(coords[[2]]) - 0.5, max(coords[[2]]) + 0.5)
  }
  
  if(is.null(col)){ #### reshape data
  
    array <- dt2array(coords = coords, contrast = contrast, format = "matrix", default_value = NA)
    
  }else { #### voxel size
    
    unique_x <- sort(unique(coords[,1]))
    n_x <- (unique_x[length(unique_x)] - unique_x[1]) # /  min(diff(unique_x))
    
    unique_y <- sort(unique(coords[,2]))
    n_y <- (unique_y[length(unique_y)] - unique_y[1]) # /  min(diff(unique_y))
    
    cex <- max(3 * graphics::par("pin") / (graphics::par("cin") * max(n_y,n_x))) * cex
    #     cex = 100 * min(graphics::par()$plt[2] - graphics::par()$plt[1], graphics::par()$plt[4] - graphics::par()$plt[3]) / max(xlim[2] - xlim[1], ylim[2] - ylim[1]) 
  }
  
  #### dislay #### 
  if(is.null(col)){
    
    graphics::image(x = array$unique_coords[[1]], y = array$unique_coords[[2]], z = array$contrast, 
                    breaks = breaks, col = palette, 
                    xlim = xlim, ylim = ylim, 
                    xlab = xlab, ylab = ylab, axes = axes, asp = asp)
     
  }else{
    
    graphics::plot(coords, 
                   pch = pch, cex = cex, col = col, xlim = xlim, ylim = ylim, 
                   axes = axes, asp = asp, xlab = xlab, ylab = ylab)
    
  }
  
  graphics::title(main = main, cex.main = cex.main)
  
  ## display missing values
  if(!is.null(col.NA)){
    
    test_NA <- which(is.na(contrast))
    
    if(sum(test_NA) > 0){
      graphics::points(coords[test_NA,], col = col.NA, pch = pch.NA, cex = cex[1] / 5)
    }
    
  }
  
  # export
  return(invisible(list(xlim = xlim,
                        ylim = ylim)))
  
}

####>>> pointsOutline ####
pointsOutline <- function(coords, filter = "2D_N4"){
  
  #     names_coords <- names(coords)
  #     if(is.null(names_coords)){names_coords <- letters[9:(8 + p)]}
  #     eval(parse(text = paste("coords <- coords[with(coords, order(", paste(names_coords[length(names_coords):1], collapse = ", "), ")),]", sep = "")))
  
  #### pathological case
  n.neighbors <- as.numeric(paste(strsplit(filter, split = "", fixed = TRUE)[[1]][ - (1:4)], collapse = ""))
  if(nrow(coords) <= n.neighbors){
    return(coords)
  }      
  
  #### conversion to array
  array <- dt2array(data.table(rep(1, nrow(coords))), 
                    coords = coords)$contrast[[1]]    
  
  #### order voxels by coordinates in the same order than array
  setkeyv(coords, rev(names(coords)))
  
  #### filter
  res.Filter <- calcFilter(array, filter = filter, norm.filter = FALSE)
  N.max <- max(calcFilter(array(1, dim = 2*dim(res.Filter$filter)), filter = filter, norm.filter = FALSE)$res)
  
  #### backconvert
#   res2 <- array2dt(array = res.Filter$res, coords = coords)[res < N.max][, res := NULL]
#  return(res2)
  
  index_array <- which(!is.na(array))
  index_outline <- which(res.Filter$res < N.max)

  return(coords[index_array %in% index_outline,])
}

