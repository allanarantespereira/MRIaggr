#**********************************************************************
#**********************************************************************
#*************         0D2 Class MRIaggr            *******************
#**********************************************************************
#**********************************************************************
#

#####  D2) Methods plot #############################################################
# boxplotMask
# heatmapMRIaggr
# multiplot
# plotDistClass
# pointsHemisphere TO BE REMOVED (ALSO FROM GENERIC)
# plotLesion3D
# plotTableLesion
# outlineMRIaggr
# show
# summary

####>>> boxplotMask ####

methods::setMethod(f = "boxplotMask", 
                   signature = "MRIaggr", 
                   definition = function(object, param, mask, num = NULL, rescale = TRUE,  
                                         ylim = NULL, col = c("white", "purple"), main = NULL, mgp = c(2, 0.5, 0), x.legend = "topright", y.legend = NULL, cex.legend = 0.8, 
                                         filename = paste(object@identifier, "boxplotMask", sep = "_"), ...){
                     
                     #### graphical options
                     ## get graphical arguments
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c("height", "hemisphere", "norm_mu", "norm_sigma", "numeric2logical", "path", "res", "unit", "width", "window"), 
                                    refuse.NULL = FALSE, method = "boxplotMask[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ## create all variables
                     height <- optionsMRIaggr.eg$height
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     norm_mu <- optionsMRIaggr.eg$norm_mu
                     norm_sigma <- optionsMRIaggr.eg$norm_sigma
                     numeric2logical <- optionsMRIaggr.eg$numeric2logical
                     path <- optionsMRIaggr.eg$path
                     res <- optionsMRIaggr.eg$res
                     unit <- optionsMRIaggr.eg$unit
                     width <- optionsMRIaggr.eg$width
                     window <- optionsMRIaggr.eg$window
                     
                     ### preparation
                     p <- length(param)
                     
                     if(length(mask) != 1){
                       stop("boxplotMask[MRIaggr] : argument \'mask\' must not contain several values \n", 
                            "number of mask requested : ", length(mask), "\n")
                     }
                     
                     carto <- selectContrast(object, param = param, num = num, hemisphere = hemisphere, norm_mu = norm_mu, norm_sigma = norm_sigma, format = "data.frame")
                     if(rescale == TRUE){carto <- scale(carto)}
                     
                     carto_mask <- selectContrast(object, param = mask, num = num, hemisphere = hemisphere, format = "data.frame")
                     
                     if(numeric2logical == TRUE){carto_mask[,1] <- as.logical(carto_mask[,1])}
                     if(!is.logical(carto_mask[,1])){
                       stop("boxplotMask[MRIaggr] : type of \'mask\' is not logical \n", 
                            "proposed type : ", paste(is(carto_mask), collapse = " "), "\n", 
                            "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
                     }
                     
                     carto_ls <- list()
                     for(iter_param in 1:p){
                       carto_ls[[paste(param[iter_param])]] <- carto[carto_mask[,1] == FALSE, iter_param]
                       carto_ls[[paste(param[iter_param], 1)]] <- carto[carto_mask[,1] == TRUE, iter_param]
                     }
                     
                     ### test export
                     res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                            n.plot = 1, mfrow = 1, xlim = NULL, ylim = NULL, 
                                            method = "boxplotMask[MRIaggr]")
                     scale <- res.init$scale
                     
                     #### display                       
                     initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                       mfrow = NULL, bg = NULL, pty = NULL, mar = NULL, mgp = mgp)              
                     
                     if(is.null(ylim)){
                       ylim <- range(carto)
                     }
                     
                     graphics::boxplot(carto_ls, ylim = ylim, col = col, main = main, names = rep("", p * 2))  
                     graphics::mtext(text = param, side = 1, line = 1, at = seq(from = 1.5, by = 2, length.out = p))
                     graphics::abline(v = seq(from = 2.5, by = 2, length.out = p - 1), lty = 2)
                     
                     graphics::legend(x = x.legend, y = y.legend, bty = "n", legend = c("other", mask), pch = c(21, 21), col = rep("black", 2), 
                                      cex = cex.legend)
                     graphics::legend(x = x.legend, y = y.legend, bty = "n", legend = c("other", mask), pch = c(20, 20), col = col, 
                                      cex = cex.legend)
                     
                     if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                       grDevices::dev.off()
                     }
                     
                   }
)

####>>> heatmapMRIaggr ####

methods::setMethod(f  = "heatmapMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, num = NULL, 
                                         rescale = TRUE, method = "pearson", points.values = TRUE, type = "image", breaks = NULL, 
                                         col = cm.colors(256), main = NULL, mgp = c(2, 0.5, 0), mar = c(4, 4, 1, 6), las = 1, cex.axis = 1, 
                                         filename = paste(object@identifier, "heatmapMRIaggr", sep = "_"), ...){
                     
                     #### graphical options
                     ## get graphical arguments
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c( "height", "hemisphere", "path", "res", "unit", "width", "window"), 
                                    refuse.NULL = FALSE, method = "heatmapMRIaggr[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ## create all variables
                     height <- optionsMRIaggr.eg$height
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     path <- optionsMRIaggr.eg$path
                     res <- optionsMRIaggr.eg$res
                     unit <- optionsMRIaggr.eg$unit
                     width <- optionsMRIaggr.eg$width
                     window <- optionsMRIaggr.eg$window
                     
                     ### preparation
                     p <- length(param)
                     carto <- selectContrast(object, param = param, num = num, hemisphere = hemisphere, norm_mu = FALSE, norm_sigma = FALSE, format = "matrix")
                     if(rescale == TRUE){
                       carto <- scale(carto)                     
                       attr(carto, "scaled:center") <- NULL
                       attr(carto, "scaled:scale") <- NULL 
                     }
                     
                     carto.corr <- stats::cor(carto, method = method)
                     rownames(carto.corr) <- param
                     colnames(carto.corr) <- param
                     grid <- expand.grid(1:p, 1:p)
                     
                     if(is.null(breaks)){breaks <- seq(min(carto.corr), max(carto.corr), length.out = length(col) + 1)}
                     
                     
                     ### test 
                     validCharacter(value = type, validLength = 1, validValues = c(FALSE, "image", "image.plot"), refuse.NULL = FALSE, method = "heatmapMRIaggr[MRIaggr]")
                     
                     res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                            n.plot = 1, mfrow = 1, xlim = NULL, ylim = NULL, 
                                            method = "boxplotMask[MRIaggr]")
                     scale <- res.init$scale
                     
                     
                     #### display
                     initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                       mfrow = NULL, bg = NULL, pty = NULL, mar = NULL, mgp = mgp)
                     
                     if(type == "image"){
                       graphics::image(1:p, 1:p, carto.corr, col = col, breaks = breaks, 
                                       main = main, axes = FALSE, xlab = "", ylab = "")
                       graphics::axis(1, at = 0:(p + 1), labels = c("", param, ""), las = las, cex.axis = cex.axis)
                       graphics::axis(2, at = 0:(p + 1), labels = c("", param, ""), las = las, cex.axis = cex.axis)
                       
                       if(points.values == TRUE){
                         graphics::text(x = grid[,1], y = grid[,2], 
                                        round(as.vector(carto.corr), digits = optionsMRIaggr("digit.legend"))
                         )
                       }   
                     }
                     
                     if(type == "image.plot"){
                       
                       graphics::par(mar = mar)
                       graphics::image(1:p, 1:p, carto.corr, col = col, breaks = breaks, 
                                       main = main, axes = FALSE, xlab = "", ylab = "")
                       graphics::axis(1, at = 0:(p + 1), labels = c("", param, ""), las = las, cex.axis = cex.axis)
                       graphics::axis(2, at = 0:(p + 1), labels = c("", param, ""), las = las, cex.axis = cex.axis)
                       
                       if(points.values == TRUE){
                         graphics::text(x = grid[,1], y = grid[,2], 
                                        round(as.vector(carto.corr), digits = optionsMRIaggr("digit.legend"))
                         )
                       }
                       
                       initPackage("fields", argument = "type = \"image.plot\"", method = "heatmapMRIaggr[MRIaggr]")
                       
                       fields::image.plot(1:p, 1:p, carto.corr, col = col, legend.only = TRUE)
                       
                     }
                     
                     
                     if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                       grDevices::dev.off()
                     }
                     
                     return(invisible(carto.corr))
                     
                     
                   }
)

####>>> multiplot ####

methods::setMethod(f  = "multiplot", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, 
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, subset = NULL, na.rm = FALSE,
                                         index1 = NULL, index2 = NULL, index3 = NULL, midplane = FALSE, 
                                         col = NULL, pch = NULL, xlim = NULL, ylim = NULL, 
                                         main = NULL, main.legend = NULL, noPlot = FALSE, filename = paste0(selectIdentifier(object),"_multiplotMRIaggr"), ...){
                     
                     optionsMRIaggr.eg <- initPoints(subset = c("norm_mu", "norm_sigma", "hemisphere", "slice_var",
                                                                "col.midplane", "lwd.midplane",
                                                                "checkArguments"), 
                                                     method = "multiplot[MRIaggr]",
                                                     ...)
                     
                     slice_var <- optionsMRIaggr.eg$slice_var
                     checkArguments <- optionsMRIaggr.eg$checkArguments
                     
                     #### data ####            
                     contrast <- selectContrast(object, param = param, norm_mu = optionsMRIaggr.eg$norm_mu, norm_sigma = optionsMRIaggr.eg$norm_sigma,
                                            slice_i = slice_i, slice_j = slice_j, slice_k = slice_k, subset = subset, na.rm = na.rm, hemisphere = optionsMRIaggr.eg$hemisphere, 
                                            coords = TRUE, format = "data.table")
              
                     n.px <- nrow(contrast)
                     
                     if( !identical(slice_var, key(contrast) ) ){
                       setkeyv(contrast, slice_var[slice_var %in% names(contrast)])
                     }
                     
                     #### tests
                     if(checkArguments){
                       
                       if(n.px == 0){
                         stop("multiplot[MRIaggr] : \'data\' contains no observation when restricted to \'slices\' and \'hemisphere\' \n", 
                              "requested slices : ", paste(slices, collapse = " "), " \n", 
                              "requested hemisphere : ", hemisphere, " \n")
                       }
                       
                       validLogical(midplane, name = "midplane", validLength = TRUE, method = multiplot[MRIaggr])
                       
                     }
                     
                     #### index  definition
                     if(!is.null(index1)){
                       index1 <- initIndex(object = object, index = index1, 
                                           slice_i = slice_i, slice_j = slice_j, slice_k = slice_k, hemisphere = optionsMRIaggr.eg$hemisphere, 
                                           arg_name = "index1", method = "multiplot")
                       if(!identical(slice_var,key(contrast))){
                         setkey(index1$coords, slice_var)
                       }
                       
                     }
                     
                     if(!is.null(index2)){
                       index2 <- initIndex(object = object, index = index2, 
                                           slice_i = slice_i, slice_j = slice_j, slice_k = slice_k, hemisphere = optionsMRIaggr.eg$hemisphere, 
                                           arg_name = "index2", method = "multiplot")
                       if(!identical(slice_var,key(contrast))){
                         setkey(index2$coords, slice_var)
                       }
                       
                     }
                     
                     if(!is.null(index3)){
                       index3 <- initIndex(object = object, index = index3, 
                                           slice_i = slice_i, slice_j = slice_j, slice_k = slice_k, hemisphere = optionsMRIaggr.eg$hemisphere, 
                                           arg_name = "index3", method = "multiplot")
                       if(!identical(slice_var,key(contrast))){
                         setkey(index3$coords, slice_var)
                       }
                       
                     }
                     
                     if(midplane){
                       indexLine <- list(coords = selectMidplane(object),
                                         col = optionsMRIaggr.eg$col.midplane,
                                         lwd = optionsMRIaggr.eg$lwd.midplane)
                       if(!identical(slice_var,key(contrast))){
                         setkey(indexLine$coords, slice_var)
                       }
                       
                     }else{
                       indexLine <- NULL
                     }
                     
                     #### preparation
                     if(is.null(main)){
                       main <- selectIdentifier(object)
                     }
                     if(is.null(main.legend)){
                       main.legend <- param
                     }
                     
                     #### main
                     res <- multiplot(object = contrast, index1 = index1, index2 = index2, index3 = index3, indexLine = indexLine,
                                      col = col, pch = pch, xlim = xlim, ylim = ylim, 
                                      main = main, main.legend = main.legend, noPlot = noPlot, filename = filename, 
                                      ...)
                     
                     #### export 
                     return(invisible(res))
                   }
)

####>>> plotDistClass ####

methods::setMethod(f  = "plotDistClass", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, param.membership, num = NULL,
                                         bw.adjust = 1, kernel = "gaussian", from = NULL, to = NULL, ylim = NULL, 
                                         col = 1:6, main = NULL, mgp = c(2, 0.5, 0), type = "l", pch = 20, lwd = 1, x.legend = "topright", y.legend = NULL, cex.legend = 0.8, 
                                         filename = paste(object@identifier, "plotDistClass", sep = "_"), ...){
                     
                     #### graphical options
                     ## get graphical arguments
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c( "height", "hemisphere", "norm_mu", "norm_sigma", "path", "res", "unit", "width", "window"), 
                                    refuse.NULL = FALSE, method = "plotDistClass[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ## create all variables
                     height <- optionsMRIaggr.eg$height
                     hemisphere <- optionsMRIaggr.eg$hemisphere
                     norm_mu <- optionsMRIaggr.eg$norm_mu
                     norm_sigma <- optionsMRIaggr.eg$norm_sigma
                     path <- optionsMRIaggr.eg$path
                     res <- optionsMRIaggr.eg$res
                     unit <- optionsMRIaggr.eg$unit
                     width <- optionsMRIaggr.eg$width
                     window <- optionsMRIaggr.eg$window
                     
                     ### preparation
                     n.membership <- length(param.membership)
                     
                     if(length(param) != 1){
                       stop("plotDistClass[MRIaggr] : wrong specification of \'param\'\n", 
                            "only one parameter can be specified \n", 
                            "length(param) : ", length(param), "\n")                  
                     }
                     
                     carto <- cbind(selectContrast(object, param = param, num = num, hemisphere = hemisphere, norm_mu = norm_mu, norm_sigma = norm_sigma, format = "data.frame"), 
                                    selectContrast(object, param = param.membership, num = num, hemisphere = hemisphere, format = "data.frame"))
                     if(is.null(from)){
                       from <- min(carto[,param])
                     }
                     if(is.null(to)){
                       to <- max(carto[,param])
                     }
                     
                     density <- list()
                     for(iter_membership in 1:n.membership){
                       
                       density[[iter_membership]] <- stats::density(carto[,param], adjust = bw.adjust, 
                                                                    kernel = kernel, 
                                                                    weights = carto[,param.membership[iter_membership]]/sum(carto[,param.membership[iter_membership]]), 
                                                                    from = from, to = to)
                     }
                     max_density <- max(unlist(lapply(density, function(x){max(x$y)})))
                     
                     ### test export
                     res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                            n.plot = 1, mfrow = 1, xlim = NULL, ylim = NULL, 
                                            method = "plotDistClass[MRIaggr]")
                     scale <- res.init$scale
                     
                     #### display
                     if(!is.null(window)){
                       
                       initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                         mfrow = NULL, bg = NULL, pty = NULL, mar = NULL, mgp = mgp)      
                       
                       if(is.null(ylim)){
                         ylim <- c(0, max_density)
                       }
                       
                       graphics::plot(NA, NA, xlim = c(from, to), xlab = param, ylab = "density", main = main, 
                                      ylim = ylim)  
                       
                       for(iter_membership in 1:n.membership){
                         graphics::points(density[[iter_membership]]$x, density[[iter_membership]]$y, 
                                          type = type, pch = pch, col = col[iter_membership], lwd = lwd)  
                       }
                       
                       graphics::legend(x = x.legend, y = y.legend, bty = "n", legend = param.membership, pch = pch, col = col, 
                                        cex = cex.legend)
                     }
                   }
)

####>>> orthplot ####

methods::setMethod(f  = "orthoplot", #### to be updated to enable data frame and data table object
                   signature  = "MRIaggr", 
                   definition = function(object, param, 
                                         slice_i = NULL, slice_j = NULL, slice_k = NULL, subset = NULL, na.rm = FALSE,
                                         index1 = NULL, index2 = NULL, index3 = NULL, midplane = FALSE, 
                                         col = NULL, pch = NULL, Ilim = NULL, Jlim = NULL, Klim = NULL, 
                                         main = NULL, main.legend = NULL,
                                         mfrow = c(2,2), col_line = "blue", lwd_line = 2, lty_line = 3, ...){

                     
                     unique_tempo <- lapply(selectCoords(object), unique)
                     levels_i <- unique_tempo[[1]]
                     levels_j <- unique_tempo[[2]]
                     levels_k <- unique_tempo[[3]]
                     range_i <- range(levels_i)
                     range_j <- range(levels_j)
                     range_k <- range(levels_k)
                    
                     if(!is.null(Ilim)){
                       slice_Ilim <- levels_i[which( (levels_i>=Ilim[1])*(levels_i<=Ilim[2]) == 1 )]
                       levels_i <- slice_Ilim
                     }else{
                       slice_Ilim <- NULL
                     }
                     if(!is.null(Ilim)){
                       slice_Jlim <- levels_j[which( (levels_j>=Jlim[1])*(levels_j<=Jlim[2]) == 1 )]
                       levels_j <- slice_Jlim
                     }else{
                       slice_Jlim <- NULL
                     }
                     if(!is.null(Ilim)){
                       slice_Klim <- levels_k[which( (levels_k>=Klim[1])*(levels_k<=Klim[2]) == 1 )]
                       levels_k <- slice_Klim
                     }else{
                       slice_Klim <- NULL
                     }
                     
                     res <- multiplot(object = object, param = param, 
                                      slice_i = slice_Ilim, slice_j = slice_Jlim, slice_k = slice_Klim, subset = subset, na.rm = na.rm,
                                      index1 = index1, index2 = index2, index3 = index3, midplane = midplane, 
                                      col = col, pch = pch, 
                                      main = main, main.legend = main.legend, noPlot = TRUE, ...)
                     
                     names.coords <- res$names.coords
                     names.param <- res$names.param
                     object <-  res$mainPlot$object
                     breaks <-  res$mainPlot$breaks
                     palette <- res$mainPlot$palette
                     optionsMRIaggr.eg <- res$mainPlot$optionsMRIaggr.eg
                     col <- res$mainPlot$col
                     pch <- res$mainPlot$pch
                     main <- res$mainPlot$main
                     legend <- res$mainPlot$legend
                     
                     breaks_sauve <- res$legendPlot$breaks_sauve
                     palette_sauve <- res$legendPlot$palette_sauve
                     n.param <- res$legendPlot$n.param
                    
                     #### tests
                     # say that height and filename and so on are useless
                     validLogical(optionsMRIaggr.eg$window, name = "window", validLength = 1, method = "orthoplot[MRIaggr]")
                     validLogical(legend, name = "legend", validLength = 1, method = "orthoplot[MRIaggr]")
                     validNumeric(slice_i, name = "slice_i", validLength = 1, validValues = levels_i, refuse.NULL = FALSE, method = "orthoplot[MRIaggr]")
                     validNumeric(slice_j, name = "slice_i", validLength = 1, validValues = levels_j, refuse.NULL = FALSE, method = "orthoplot[MRIaggr]")
                     validNumeric(slice_k, name = "slice_i", validLength = 1, validValues = levels_k, refuse.NULL = FALSE, method = "orthoplot[MRIaggr]")
                     # test mfrow = 4 if legend and 3 else (also be carefull for the number of paramter)
                     
                     #### prepare
                    
                     if(is.null(slice_i)){
                      slice_i <- levels_i[round(length(levels_i)/2)]
                     }
                     
                     if(is.null(slice_j)){
                       slice_j <- levels_j[round(length(levels_j)/2)]
                     }
                     
                     if(is.null(slice_k)){
                       slice_k <- levels_k[round(length(levels_k)/2)]
                     }
                     
                     xlimI <- NULL
                     ylimI <- NULL
                     xlimJ <- NULL
                     ylimJ <- NULL
                     xlimK <- NULL
                     ylimK <- NULL
                     
                     #### windows
                     if(!is.null(optionsMRIaggr.eg$bg)){graphics::par(bg = optionsMRIaggr.eg$bg)}
                     if(!is.null(optionsMRIaggr.eg$pty)){graphics::par(pty = optionsMRIaggr.eg$pty)}
                     if(!is.null(optionsMRIaggr.eg$mar)){graphics::par(mar = optionsMRIaggr.eg$mar)}   
                     if(!is.null(optionsMRIaggr.eg$mgp)){graphics::par(mgp = optionsMRIaggr.eg$mgp)}     
                     
#                      M.layout <- matrix(NA, nrow = mfrow[1], ncol = mfrow[2]+n.contrast - 1)
#                      M.layout[,1:mfrow[2]] <- 1:prod(mfrow)
#                      M.layout[1:(mfrow[1] - 1),  - (1:mfrow[2])] <-  M.layout[1:(mfrow[1] - 1), mfrow[2]]    
#                      M.layout[mfrow[1],  - (1:mfrow[2])] <-  seq(prod(mfrow) + 1, prod(mfrow) + n.contrast - 1)  
#                      widths.layout <- c(rep(1 / mfrow[2], mfrow[2] - 1), 
#                                         rep(1 / (mfrow[2] * n.contrast), n.contrast)
#                      )
                     
                     if(n.param == 1){ # cas uniparametrique
                       split.screen(figs = mfrow)    
                       
                     }else{ # cas multiparametrique avec legende
                       
                       stop("orthoplot[MRIaggr]: multiparametric not yet supported \n")
                       graphics::layout(M.layout, widths = widths.layout)
                     }
                     
                        # split display into two screens
                     
                     #### function
                     warperOrthoplot <- function(indexI, indexJ, indexK,
                                                 plot.i, plot.j, plot.k, plot.legend,
                                                 xlimI, ylimI, xlimJ, ylimJ, xlimK, ylimK){
                       
                       mainI <- paste0(main," (i=",slice_i,")")
                       mainJ <- paste0(main," (j=",slice_j,")")
                       mainK <- paste0(main," (k=",slice_k,")")
                       
                       ## projection on i (j,k)
                       screen(1)
                       
                       plotI.test <- plotMRI(contrast = object[indexI, names.param[1], with = FALSE], 
                                            coords = object[indexI, c("j","k"), with = FALSE], 
                                            breaks = breaks, col = if(col){object[indexI][["MRIaggr.colXX"]]}else{NULL}, palette = palette, 
                                            asp = optionsMRIaggr.eg$asp, 
                                            xlim = xlimI, ylim = ylimI, pch = pch, cex = optionsMRIaggr.eg$cex, axes = optionsMRIaggr.eg$axes, col.NA = optionsMRIaggr.eg$col.NA, pch.NA = optionsMRIaggr.eg$pch.NA, 
                                            xlab = optionsMRIaggr.eg$xlab, ylab = optionsMRIaggr.eg$ylab, main = mainI, cex.main = optionsMRIaggr.eg$cex.main)
                       
                       abline(v=slice_j, col = col_line, lwd = lwd_line, lty = lty_line)
                       abline(h=slice_k, col = col_line, lwd = lwd_line, lty = lty_line)
                       
                       ## projection on j (i,k)
                       screen(2)                     
                       plotJ.test <- plotMRI(contrast = object[indexJ, names.param[1], with = FALSE], 
                                            coords = object[indexJ, c("i","k") , with = FALSE], 
                                            breaks = breaks, col = if(col){object[indexJ][["MRIaggr.colXX"]]}else{NULL}, palette = palette, 
                                            asp = optionsMRIaggr.eg$asp, 
                                            xlim = xlimJ, ylim = ylimJ, pch = pch, cex = optionsMRIaggr.eg$cex, axes = optionsMRIaggr.eg$axes, col.NA = optionsMRIaggr.eg$col.NA, pch.NA = optionsMRIaggr.eg$pch.NA, 
                                            xlab = optionsMRIaggr.eg$xlab, ylab = optionsMRIaggr.eg$ylab, main = mainJ, cex.main = optionsMRIaggr.eg$cex.main)
                       
                       abline(v=slice_i, col = col_line, lwd = lwd_line, lty = lty_line)
                       abline(h=slice_k, col = col_line, lwd = lwd_line, lty = lty_line)
                       
                       ## projection on k (i,j)
                       screen(3)
                       plotK.test <- plotMRI(contrast = object[indexK, names.param[1], with = FALSE], 
                                            coords = object[indexK, c("i","j"), with = FALSE], 
                                            breaks = breaks, col = if(col){object[indexK][["MRIaggr.colXX"]]}else{NULL}, palette = palette, 
                                            asp = optionsMRIaggr.eg$asp, 
                                            xlim = xlimK, ylim = ylimK, pch = pch, cex = optionsMRIaggr.eg$cex, axes = optionsMRIaggr.eg$axes, col.NA = optionsMRIaggr.eg$col.NA, pch.NA = optionsMRIaggr.eg$pch.NA, 
                                            xlab = optionsMRIaggr.eg$xlab, ylab = optionsMRIaggr.eg$ylab, main = mainK, cex.main = optionsMRIaggr.eg$cex.main)
                       
                       abline(v=slice_i, col = col_line, lwd = lwd_line, lty = lty_line)
                       abline(h=slice_j, col = col_line, lwd = lwd_line, lty = lty_line)
                       
                       if(plot.legend){
                       screen(4)
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
                       
                       return(invisible(list(I = if(plot.i){plotI.test}else{NULL},
                                             J = if(plot.j){plotJ.test}else{NULL},
                                             K = if(plot.k){plotK.test}else{NULL})))
                     }
                     
                     
                     ### main 
                     continue <- TRUE
                     plot.i <- TRUE 
                     plot.j <- TRUE
                     plot.k <- TRUE
                     plot.legend <- TRUE
                     
                     while(continue){
                       
                       if(plot.i){
                         indexI <- object[ ,.I[.SD == slice_i], .SDcols = names.coords[1]]
                       }
                       if(plot.j){
                         indexJ <- object[ ,.I[.SD == slice_j], .SDcols = names.coords[2]]
                       }
                       if(plot.k){
                         indexK <- object[ ,.I[.SD == slice_k], .SDcols = names.coords[3]]
                       }
                      
                       ls.lim <- warperOrthoplot(indexI = indexI, indexJ = indexJ, indexK = indexK,
                                                 plot.i = plot.i, plot.j  = plot.j, plot.k = plot.k, plot.legend = plot.legend,
                                                 xlimI = xlimI, ylimI = ylimI, xlimJ = xlimJ, ylimJ = ylimJ, xlimK = xlimK, ylimK = ylimK)
                       
#                        res <- locator(n = 4)
#                        print(res)
                       xlimI <- ls.lim$xlimI
                       ylimI <- ls.lim$ylimI
                       xlimJ <- ls.lim$xlimJ
                       ylimJ <- ls.lim$ylimJ
                       xlimK <- ls.lim$xlimK
                       ylimK <- ls.lim$ylimK
                       plot.legend <- FALSE
                       
                       text <- readline("new coodinate value: ")
                     
                       if(is.null(text) || text %in%  c("Q","q","quit","Quit","QUIT","exit","Exit","EXIT")){
                         continue <- FALSE
                       } else {
                         
                         if(length(grep("=",x = text,fixed = TRUE)>0)){
                           text <- gsub(pattern = " ", replacement = "",
                                        strsplit(text, split = "=", fixed = TRUE)[[1]], 
                                        fixed = TRUE)
                         }else if(text == "+"){
                           text <- c(text.save[1],
                                     as.numeric(text.save[2])+1
                           )
                         }else if(text == "-"){
                           text <- c(text.save[1],
                                     as.numeric(text.save[2])+1
                           )
                         }else{
                           text <- c(text.save[1],
                                     text)
                         }
                         
                       
                      
                       if(text[1] == "i"){
                         slice_i <- as.numeric(text[2])
                         plot.i <- TRUE 
                         plot.j <- FALSE
                         plot.k <- FALSE
                       }else if(text[1] == "j"){
                         slice_j <- as.numeric(text[2])
                         plot.i <- FALSE 
                         plot.j <- TRUE
                         plot.k <- FALSE
                       }else if(text[1] == "k"){
                         slice_k <- as.numeric(text[2])
                         plot.i <- FALSE 
                         plot.j <- FALSE
                         plot.k <- TRUE
                       } else{
                         continue <- FALSE
                       }
                         text.save <- text
                       }
                       #                        res <- graphics::locator(n = 1)
                       #                        res$xI <- res$x 
                       #                        res$yI <- res$y 
                       #                        res$xJ <- res$x
                       #                        res$yJ <- res$y
                       #                        res$xK <- res$x 
                       #                        res$yK <- res$y
                       #                        
                       #                        test.deviceI <- (res$xI >= range_i[1])*(res$xI <= range_i[2])*(res$yI >= range_j[1])*(res$yI <= range_j[2])
                       #                        test.deviceJ <- (res$xJ >= range_i[1])*(res$xJ <= range_i[2])*(res$yJ >= range_j[1])*(res$yJ <= range_j[2])
                       #                        test.deviceK <- (res$xK >= range_i[1])*(res$xK <= range_i[2])*(res$yK >= range_j[1])*(res$yK <= range_j[2])
#                        
#                        if(is.null(res) ||  (test.deviceI==FALSE && test.deviceJ==FALSE && test.deviceK==FALSE )){
#                          continue <- FALSE
#                        }else if(test.deviceI){
#                          slice_j <- res$xK
#                          slice_k <- res$yK
#                        }else if(test.deviceI){
#                          slice_i <- res$yK
#                          slice_k <- res$xK
#                        }else {#if(test.deviceI){
#                          slice_i <- res$xK
#                          slice_j <- res$yK
#                        }
                     }
                     
                     
                     
                   }
)

####>>> plotLesion3D ####

methods::setMethod(f  = "plotLesion3D", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask, edge = FALSE, Neighborhood = "3D_N6", numeric2logical = FALSE, spatial_res = c(1, 1, 1), 
                                         xlim = NULL, ylim = NULL, zlim = NULL, type.plot = "shapelist3d", px_max = 10000, 
                                         radius = 1, type = "s", col = "red", col.edge = "black"){
                     
                     initPackage("rgl", argument = "type = \"image.plot\"", method = "plotLesion3D[MRIaggr]")
                     
                     # import
                     if(length(mask) != 1){
                       stop("plotLesion3D[MRIaggr] : argument \'mask\' must not contain several values \n", 
                            "number of mask requested : ", length(mask), "\n")
                     }					 
                     carto <- selectContrast(object, param = mask, coords = TRUE)
                     
                     validNumeric(value = spatial_res, validLength = 3, min = 0, method = "plotLesion3D[MRIaggr]")
                     
                     carto[,"i"] <- carto[,"i"]*spatial_res[1]
                     carto[,"j"] <- carto[,"j"]*spatial_res[2]
                     carto[,"k"] <- carto[,"k"]*spatial_res[3]
                     
                     if(numeric2logical == TRUE){carto[,mask] <- as.logical(carto[,mask])}
                     if(!is.logical(carto[,mask])){
                       stop("plotLesion3D[MRIaggr] : type of \'mask\' is not logical \n", 
                            "proposed type : ", paste(is(carto[,mask]), collapse = " "), "\n", 
                            "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
                     }
                     
                     index_mask <- which(carto[,mask] == 1)
                     n.mask <- length(index_mask)
                     
                     if(edge){
                       carto_Wmask <- calcFilter(object, param = mask, filter = Neighborhood, norm.filter = FALSE, verbose = FALSE)
                       nb_neighbor <- max(carto_Wmask[,paste(mask, Neighborhood, sep = "_")])
                       index_core <- which(carto_Wmask[,paste(mask, Neighborhood, sep = "_")] == nb_neighbor)
                       index_edge <- setdiff(index_mask, index_core)
                     }
                     
                     # test
                     validCharacter(value = type.plot, validLength = 1, validValues = c("plot3d", "shapelist3d"), method = "plotLesion3D[MRIaggr]")
                     
                     if(n.mask > px_max){
                       stop("plotLesion3D[MRIaggr] : the number of points to be displayed exceed the limit \'px_max\' \n", 
                            "number of points to be displayed : ", n.mask, "\n", 
                            "limit  : ", px_max, "\n")
                     }
                     
                     # preparation                      
                     if(is.null(xlim)){ xlim <- range(carto[index_mask, "i"])}
                     
                     if(is.null(ylim)){ ylim <- range(carto[index_mask, "j"])}
                     
                     if(is.null(zlim)){ zlim <- range(carto[index_mask, "k"])}
                     
                     # verbose
                     if(type.plot == "plot3d"){
                       
                       if(edge == FALSE){
                         rgl::plot3d(carto[index_mask, c("i", "j", "k")], xlim = xlim, ylim = ylim, zlim = zlim, 
                                     col = col, type = type, radius = radius)
                       }else{
                         rgl::plot3d(carto[index_core, c("i", "j", "k")], xlim = xlim, ylim = ylim, zlim = zlim, 
                                     col = col, type = type, radius = radius)
                         rgl::points3d(carto[index_edge, c("i", "j", "k")], 
                                       col = col.edge)
                       }
                     }
                     
                     if(type.plot == "shapelist3d"){
                       M <- diag(spatial_res)
                       if(edge == FALSE){
                         rgl::shapelist3d(rgl::cube3d(trans = M), x = carto[index_mask, "i"], y = carto[index_mask, "j"], z = carto[index_mask, "k"], 
                                          size = radius, xlim = xlim, ylim = ylim, zlim = zlim, col = col)                
                         rgl::axes3d(c('x', 'y', 'z'))
                         rgl::title3d(main = mask, xlab = "i", ylab = "j", zlab = "k")
                         
                       }else{
                         rgl::shapelist3d(rgl::cube3d(trans = M), x = carto[index_edge, "i"], y = carto[index_edge, "j"], z = carto[index_edge, "k"], 
                                          size = radius, col = col.edge)
                         rgl::axes3d(c('x', 'y', 'z'))
                         rgl::title3d(main = mask, xlab = "i", ylab = "j", zlab = "k")
                         
                       }
                     }
                     
                   }
)

####>>> plotTableLesion ####

methods::setMethod(f  = "plotTableLesion", 
                   signature  = "MRIaggr", 
                   definition = function(object, mask, num = NULL, type = "matplot", 
                                         col = 1:5, lty = 1:5, lwd = 1, mgp = c(2, 0.5, 0), mar = rep(3, 4), 
                                         main = paste("lesion - ", object@identifier, sep = ""), cex.legend = 1, cex.main = 1, cex.axis = 1, cex.lab = 1, 
                                         filename = paste(object@identifier, "plotTableLesion", sep = "_"), ...){                    
                     
                     #### graphical options
                     ## get graphical arguments
                     optionsMRIaggr.eg <- optionsMRIaggr()					 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = c("height", "path", "res", "unit", "width", "window"), 
                                    refuse.NULL = FALSE, method = "plotTableLesion[MRIaggr]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ## create all variables
                     height <- optionsMRIaggr.eg$height
                     path <- optionsMRIaggr.eg$path
                     res <- optionsMRIaggr.eg$res
                     unit <- optionsMRIaggr.eg$unit
                     window <- optionsMRIaggr.eg$window
                     width <- optionsMRIaggr.eg$width
                     
                     #### preparation ####
                     
                     num <- initNum(object, num = num, method = "plotTableLesion")
                     
                     table_lesion <- selectTable(object, type = "lesion")
                     n.num <- length(num)
                     p <- length(mask)
                     
                     validCharacter(value = mask, validLength = NULL, validValues = names(table_lesion), method = "plotTableLesion[MRIaggr]")
                     validCharacter(value = type, validLength = 1, validValues = c("matplot", "evolution"), method = "plotTableLesion[MRIaggr]")
                     
                     #### test export
                     res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                            n.plot = 1, mfrow = 1, xlim = NULL, ylim = NULL, 
                                            method = "boxplotMask[MRIaggr]")
                     scale <- res.init$scale
                     
                     #### display
                     initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                                       mfrow = NULL, bg = NULL, pty = NULL, mar = NULL, mgp = mgp)      
                     
                     if(type == "matplot"){
                       graphics::matplot(table_lesion[num, mask], 
                                         type = "l", xlab = "k", ylab = "observation number", lty = lty, lwd = lwd, col = col, main = main, axes = FALSE, 
                                         cex.main = cex.main, cex.lab = cex.lab)
                       graphics::axis(1, at = 1:n.num, labels = num, cex.axis = cex.axis)
                       graphics::axis(2, cex.axis = cex.axis)
                       
                       graphics::legend("topleft", legend = mask, lty = lty, col = col, bty = "n", cex = cex.legend)
                     }else{
                       if(!is.null(window)){
                         graphics::par(mfrow = c(1, p))
                         if(!is.null(mar)){graphics::par(mar = mar)}
                       }
                       
                       lwd.init <- lwd
                       col.init <- col
                       for(iter_p in 1:p){
                         mask_max <- max(table_lesion[num, mask[iter_p]])
                         col <- rep(col.init[1], n.num)
                         
                         if(iter_p > 1){
                           croissance <-  table_lesion[num, mask[iter_p]] - table_sauve
                           lwd <- lwd.init + 10 * abs(croissance) / max(abs(croissance))
                           col[croissance > 0] <- col.init[2]
                           col[croissance < 0] <- col.init[3]
                         }
                         
                         graphics::plot(NA, NA, xlim = c(0, mask_max), ylim = range(num), 
                                        xlab = mask[iter_p], ylab = "k", main = main, 
                                        cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab)
                         graphics::segments(y0=num, 
                                            y1=num, 
                                            x0 = 0, 
                                            x1=table_lesion[num, mask[iter_p]], 
                                            lwd = lwd, 
                                            col = col)         
                         
                         table_sauve <- table_lesion[num, mask[iter_p]]                
                       }
                     }
                     
                     if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                       grDevices::dev.off()
                     }
                     
                     
                   }
)

####>>> outlineMRIaggr ####

methods::setMethod(f  = "outlineMRIaggr", 
                   signature  = "MRIaggr", 
                   definition = function(object, param, index1 = NULL, num = NULL, hemisphere = "both", numeric2logical = FALSE, 
                                         xlim = NULL, ylim = NULL, legend = FALSE, palette = "terrain.colors", col = NULL, breaks = 25, 
                                         fill = TRUE, n = 50, sequential = TRUE, min_dist = 1, operator_index1 = "none", 
                                         col.outline = c("blue", "red", "grey"), pch = 20, cex = c(0.75, 1, 0.75), 
                                         name_newparam = "userMask", 
                                         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE){
                     
                     
                     num <- initNum(object = object, num = num, method = "outlineMRIaggr")
                     param <- initParameter(object = object, param = param, test = TRUE, init = TRUE, method = "outlineMRIaggr")      
                     
                     newparam <- selectCoords(object, coords = c("i", "j", "k", "index"))
                     edge <- matrix(NA, nrow = 0, ncol = 3)
                     surface <- matrix(NA, nrow = 0, ncol = 3)
                     
                     validLogical(value = legend, validLength = 1, refuse.NULL = FALSE, refuse.NA = TRUE, method = "outlineMRIaggr[MRIaggr]")
                     validCharacter(value = name_newparam, validLength = 1, method = "outlineMRIaggr[MRIaggr]")
                     
                     if(is.null(index1)){operator_index1 <- "none"}
                     validCharacter(value = operator_index1, validLength = 1, validValues = c("difference", "union", "intersection", "none"), method = "outlineMRIaggr[MRIaggr]")
                     
                     if(!is.null(index1)){
                       index1 <- initIndex(object = object, index = index1, num = num, hemisphere = hemisphere, numeric2logical = numeric2logical, indexNum = 1, 
                                           outline.default = optionsMRIaggr("outline.index"), cex.default = 1, pch.default = 20, col.default = "black", filter.default = "2D_N4", method = "outlineMRIaggr")
                     }
                     
                     grDevices::graphics.off()
                     for(iter_slice in 1:length(num)){
                       
                       # display 
                       repeat.slice <- 1
                       while(!is.na(repeat.slice) && repeat.slice == 1){
                         
                         multiplot(object, param = param, num = num[iter_slice], hemisphere = hemisphere, xlim = xlim, ylim = ylim, index1=index1, 
                                   col = col, breaks = breaks, palette = palette, legend = legend)
                         
                         res_outline <- outline(n = n, sequential = sequential, min_dist = min_dist, 
                                                col = col.outline, pch = pch, cex = cex)
                         
                         if(is.null(res_outline$edge)){
                           repeat.slice <- readline("Do you want to start again ? (answer 1)\n Otherwise the slice will be skipped  ")                
                         }else{
                           repeat.slice <- FALSE
                         }
                         
                       }
                       
                       if(!is.null(res_outline$edge)){
                         edge <- rbind(edge, 
                                       cbind(res_outline$edge[,c("i", "j")], k = num[iter_slice]))
                       }
                       
                       if(!is.null(res_outline$surface)){
                         surface <- rbind(surface, 
                                          cbind(res_outline$surface[,c("i", "j")], k = num[iter_slice]))
                       }
                     }
                     
                     newparam <- merge(x = newparam, y = data.frame(edge, edge = TRUE), by = c("i", "j", "k"), all.x = TRUE, all.y = FALSE)
                     newparam <- merge(x = newparam, y = data.frame(surface, surface = TRUE), by = c("i", "j", "k"), all.x = TRUE, all.y = FALSE)
                     
                     newparam$surface[is.na(newparam$surface)] <- FALSE
                     newparam$edge[is.na(newparam$edge)] <- FALSE
                     
                     if(operator_index1 != "none"){
                       newparam <- merge(x = newparam, y = data.frame(index1$coords, index1 = TRUE), by = c("i", "j", "k"), all.x = TRUE, all.y = FALSE)
                       newparam$index1[is.na(newparam$index1)] <- FALSE
                       
                       if(operator_index1 == "difference"){
                         newparam$edge <- as.logical(newparam$edge + newparam$index1 - 2 * newparam$edge * newparam$index1)
                         newparam$surface <- as.logical(newparam$surface + newparam$index1 - 2 * newparam$surface * newparam$index1)
                       }
                       if(operator_index1 == "intersection"){
                         newparam$edge <- newparam$edge * newparam$index1 
                         newparam$surface <- newparam$surface * newparam$index1 
                       }
                       if(operator_index1 == "union"){
                         newparam$edge <- as.logical(newparam$edge + newparam$index1)
                         newparam$surface <- as.logical(newparam$surface + newparam$index1)
                       }
                     }
                     
                     newparam <- newparam[order(newparam$index),]
                     
                     if(fill == TRUE){  
                       newparam <- eval(parse(text = paste("data.frame(newparam, ", name_newparam, " = newparam[,\"surface\"])", sep = "")))
                     }else{
                       newparam <- eval(parse(text = paste("data.frame(newparam, ", name_newparam, " = newparam[,\"edge\"])", sep = "")))
                     }
                     
                     return(list(res = newparam, 
                                 verbose = verbose, 
                                 name_newparam = name_newparam, 
                                 update.object = update.object, 
                                 overwrite = overwrite)
                     )
                   }           
)

####>>> show ####

methods::setMethod(f  = "show", 
                   signature  = "MRIaggr", 
                   definition = function(object){
                     
                     summary(object, verbose = FALSE)
                     
                   }
)

####>>> summary ####

methods::setMethod(f  = "summary", 
                   signature  = "MRIaggr", 
                   definition = function(object, param = FALSE, clinic = FALSE, descStats = FALSE, history = FALSE,
                                         checkArguments = optionsMRIaggr("checkArguments"), verbose = optionsMRIaggr("verbose")){
                     
                     #### test
                     if (checkArguments == TRUE) {
                       validLogical(value = param, validLength = 1, method = "summary[MRIaggr]")
                       validLogical(value = clinic, validLength = 1, method = "summary[MRIaggr]")
                       validLogical(value = descStats, validLength = 1, method = "summary[MRIaggr]")
                       validLogical(value = history, validLength = 1, method = "summary[MRIaggr]")
                       validLogical(value = verbose, validLength = 1, method = "summary[MRIaggr]")
                     }
                     
                     #### extract 
                     
                     # names and size
                     names.param <- selectParameter(object, type = "contrast")
                     n.param <- length(names.param)
                     names.coords <- selectParameter(object, type = "coords")
                     n.coords <- length(names.coords)
                     names.descStats <- selectParameter(object, type = "ls_descStats")
                     n.descStats <- length(names.descStats)
                     names.clinic <- selectParameter(object, type = "clinic")
                     n.clinic <- length(names.clinic)
                     n.data <- selectN(object)
                     
                     ncol.TableLesion <- ncol(selectTable(object, type = "lesion"))
                     ncol.TableReperfusion <- ncol(selectTable(object, type = "reperfusion"))
                     ncol.TableHypoperfusion <- ncol(selectTable(object, type = "hypoperfusion"))
                     
                     # test
                     test.normalization <- length(selectNormalization(object)) > 0
                     test.midplane <- all(is.na(selectMidplane(object))) == FALSE
                     test.voxelDim <- all(is.na(selectVoxelDim(object, unit = FALSE))) == FALSE
                     test.history <- length(selectHistory(object))
                     test.hemisphere <- any(is.na(selectHemispheres(object)) == FALSE)
                     
                     # slot
                     ls_descStats <- selectDescStats(object)
                     df.clinic <- selectClinic(object)
                     History <- selectHistory(object)
                     FieldDim <- selectFieldDim(object)
                     Identifier <- selectIdentifier(object)
                     Hemispheres <- selectHemispheres(object)
                     
                     #### print
                     cat("Object of class \'MRIaggr\' with identifier: ", Identifier , "\n", sep = "")
                     
                     cat("  # image dimensions (i, j, k): ", paste(FieldDim, collapse = " x "), " voxels \n", sep = "")
                     
                     if(test.voxelDim){
                     cat("  # voxel dimensions (i, j, k unit): ", paste(selectVoxelDim(object, unit = FALSE, format = "vector"), collapse = " x "), 
                         " ", selectVoxelDim(object, coords = FALSE, unit = TRUE, format = "vector"), " \n", sep = "")
                     }else{
                     cat("  # voxel dimensions (i, j, k unit): Unknown \n", sep = "")
                     }
                     
                     if(verbose == TRUE || test.history == TRUE){
                       cat("  # history: ", length(History), " calculations performed on the object \n", sep = "")
                       if(test.history == TRUE && history == TRUE){
                         cat(paste(names(History), collapse = " "), " \n", sep = " ")
                       }
                     }
                     
                     cat("  # number of observations: ", n.data, " \n", sep = "")
                     
                     cat("  # number of contrast parameters: ", n.param, " \n", sep = "")
                     if(param == TRUE){
                       utils::str(selectContrast(object), max.level = 1, give.attr = FALSE)
                     }
                     
                     if(verbose == TRUE || ncol.TableLesion > 0 || ncol.TableHypoperfusion > 0 || ncol.TableReperfusion > 0){
                       cat("  # tables: lesion ", if(ncol.TableLesion > 0){paste(ncol.TableLesion, " columns", sep = "")}else{"empty"}, 
                           ", hypoperfusion ", if(ncol.TableHypoperfusion > 0){paste(ncol.TableHypoperfusion, " columns", sep = "")}else{"empty"}, 
                           ", reperfusion ", if(ncol.TableReperfusion > 0){paste(ncol.TableReperfusion, " columns", sep = "")}else{"empty"}, "\n", sep = "")
                     }
                     
                     if(verbose == TRUE || test.normalization > 0){
                       cat("  # normalization values: ", if(test.normalization > 0){"present"}else{"empty"}, "\n", sep = "")
                     }
                     
                     if(verbose == TRUE || test.hemisphere == TRUE){
                       cat("  # hemisphere left: ", paste( unlist(Hemispheres$left), collapse = ""), 
                           " | hemisphere right: ", paste( unlist(Hemispheres$right), collapse = ""), 
                           "\n", sep = "")
                     }
                     
                     if(verbose == TRUE || test.midplane ){
                       cat("  # midplane coordinates: ", if(test.midplane){"present"}else{"empty"}, "\n", sep = "")
                     }
                     
                     if(verbose == TRUE || n.descStats > 0 ){
                       cat("  # ls_descStats: ", if(n.descStats == 0){"empty"}else{if(descStats == FALSE){paste(n.descStats, " elements", sep = "")}}, "\n", sep = "")
                       if(n.descStats > 0 && descStats == TRUE){
                         utils::str(ls_descStats, max.level = 1, give.attr = FALSE)              
                       }
                     }
                     
                     if(verbose == TRUE || n.clinic > 0 ){
                       cat("  # clinic : ", if(n.clinic == 0){"empty"}else{if(clinic == FALSE){paste(n.clinic, " parameters", sep = "")}}, "\n", sep = "")
                       if(n.clinic > 0 && clinic == TRUE){
                         utils::str(df.clinic, max.level = 1, give.attr = FALSE)              
                       }
                     }
                     
                   }
)
