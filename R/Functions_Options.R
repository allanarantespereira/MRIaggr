#**********************************************************************
#**********************************************************************
#*************         Default Options MRIaggr      *******************
#**********************************************************************
#**********************************************************************
#

#####  Default options #############################################################
# allocOptionsMRIaggr
# optionsMRIaggr
# selectOptionsMRIaggr
# validOptionsMRIaggr

MRIaggr.env <- new.env()

assign(".ls_optionsMRIaggr", 
       list(asp = 1, 
            axes = TRUE, 
            bg = "lightblue", 
            breaks = 50, 
            cex = 1, 
            cex.index = c(1, 1, 1), 
            cex.legend = 1.5, 
            cex.main = 1.5, 
            checkArguments = TRUE,
            col.index = c("red","purple","green"), 
            col.midplane = "red",
            col.NA = "lightyellow", 
            digit.legend = 2,
            digit.result = 2,
            digit.epsilon = 5,
            digit.percentage = 3,
            filter.index = "2D_N4",
            height = 500, 
            hemisphere = "both", 
            legend = TRUE, 
            lwd.midplane = 2,
            mar = rep(1.5, 4), 
            mar.legend = c(2, 7, 2, 2), 
            mfrow = NULL, 
            mgp = c(2, 0.5, 0), 
            norm_mu = FALSE, 
            norm_sigma = FALSE, 
            num.main = TRUE, 
            numeric2logical = FALSE, 
            operator.withinR = "union", 
            operator.betweenR = "union",
            outline.index = FALSE, 
            palette = "terrain.colors", 
            path = NULL, 
            pch.index = 20:22, 
            pch.NA = 8, 
            pty = NULL, 
            quantiles.legend = TRUE, 
            res = NA,
            slice_var = c("i","j","k"), 
            type.breaks = "range", 
            unit = "px", 
            verbose = TRUE,
            xlab = "", 
            ylab = "", 
            width = 1000, 
            window = FALSE
       ), 
       envir = MRIaggr.env)

####>>> allocOptionsMRIaggr ####
allocOptionsMRIaggr <- function(field, n.args){
  
  ## test
  validOptionsMRIaggr(field, method = "allocOptionsMRIaggr")
  
  ## main
  names.field <- names(field)
  .ls_optionsMRIaggr <- get(".ls_optionsMRIaggr", envir = MRIaggr.env)
  
  for(iter_field in 1:n.args){
    .ls_optionsMRIaggr[[names.field[iter_field]]] <- field[[iter_field]]
  }
  
  ## internal
  assign(".ls_optionsMRIaggr", .ls_optionsMRIaggr, envir = MRIaggr.env)
  
  ## for the user
  return(field)
}

####>>> optionsMRIaggr ####
optionsMRIaggr <- function(..., format = "reduce", reinit.options = FALSE){
  
  if(reinit.options == TRUE){
    
    assign(".ls_optionsMRIaggr", 
           list(asp = 1, 
                axes = TRUE, 
                bg = "lightblue", 
                breaks = 50, 
                cex = 1, 
                cex.index = c(1, 1, 1), 
                cex.legend = 1.5, 
                cex.main = 1.5, 
                checkArguments = TRUE,
                col.index = c("red","purple","green"), 
                col.midplane = "red",
                col.NA = "lightyellow", 
                digit.legend = 2,
                digit.result = 2,
                digit.epsilon = 5,
                digit.percentage = 3,
                filter.index = "2D_N4",
                height = 500, 
                hemisphere = "both", 
                legend = TRUE, 
                lwd.midplane = 2,
                mar = rep(1.5, 4), 
                mar.legend = c(2, 7, 2, 2), 
                mfrow = NULL, 
                mgp = c(2, 0.5, 0), 
                norm_mu = FALSE, 
                norm_sigma = FALSE, 
                num.main = TRUE, 
                numeric2logical = FALSE, 
                operator.withinR = "union", 
                operator.betweenR = "union",
                outline.index = FALSE, 
                palette = "terrain.colors", 
                path = NULL, 
                pch.index = 20:22, 
                pch.NA = 8, 
                pty = NULL, 
                quantiles.legend = TRUE, 
                res = NA,
                slice_var = c("i","j","k"), 
                type.breaks = "range", 
                unit = "px", 
                verbose = TRUE,
                xlab = "", 
                ylab = "", 
                width = 1000, 
                window = FALSE
           ), 
           envir = MRIaggr.env)
    
    return(invisible(get(".ls_optionsMRIaggr", envir = MRIaggr.env)))
  }
  
  args <- list(...)
  n.args <- length(args)
  validCharacter(format, name = "format", validLength = 1, validValues = c("list","reduce"), method = "optionsMRIaggr")
  
  #### si lecture 
  if(n.args == 0){ # retourne tout
    
    value <- selectOptionsMRIaggr()
    single <- FALSE
    
  }else if(all(unlist(lapply(args, is.character)))){ # retourne uniquement les arguments demandes
    
    args <- unlist(args)
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- selectOptionsMRIaggr(args)
    
  } else { #### si ecriture 
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- allocOptionsMRIaggr(args, n.args)
  }
  
  #### export
  if (single && format == "reduce"){ 
    value <- value[[1L]]
  }
  
  if(!is.null(names(args))){
    invisible(value)
  }else{
    value
  } 
}

####>>> selectOptionsMRIaggr ####

selectOptionsMRIaggr <- function(field = NULL){
  
  .ls_optionsMRIaggr <- get(".ls_optionsMRIaggr", envir = MRIaggr.env)
  
  if(is.null(field)){
    
    return(.ls_optionsMRIaggr)
    
  }else{
    
    validCharacter(value = field, validLength = NULL, validValues = names(.ls_optionsMRIaggr), method = "selectOptionsMRIaggr")
    
    return(.ls_optionsMRIaggr[field])
  }
  
}

####>>> validOptionsMRIaggr ####

validOptionsMRIaggr <- function(field, method){
  
  names.field <- names(field)
  
  ## global
  validCharacter(value = names.field, name = "field", validLength = NULL, validValues = names(selectOptionsMRIaggr()), refuse.NULL = TRUE, method = method)
  
  ## argument by argument
  if("asp" %in% names.field){
  validNumeric(value = field$asp, name = "asp", validLength = 1, min = 0, refuse.NULL = FALSE, method = method)
  }
  if("axes" %in% names.field){
  validLogical(value = field$axes, name = "axes", validLength = 1, method = method)
  }
  if("bg" %in% names.field){
    validCharacter(value = field$bg, name = "bg", validLength = 1, method = method)
  }
  if("breaks" %in% names.field){
    validNumeric(value = field$breaks, name = "breaks", refuse.duplicates = TRUE, validLength = NULL, method = method)
  }
  if("cex" %in% names.field){
    validNumeric(value = field$cex, name = "cex", validLength = 1, min = 0, method = method)
  }
  if("cex.index" %in% names.field){
    validNumeric(value = field$cex.index, name = "cex.index", validLength = 3, min = 0, method = method)
  }
  if("cex.legend" %in% names.field){
    validNumeric(value = field$cex.legend, name = "cex.legend", validLength = 1, min = 0, method = method)
  }
  if("cex.main" %in% names.field){
    validNumeric(value = field$cex.main, name = "cex.main", validLength = 1, min = 0, method = method)
  }
  if("checkArguments" %in% names.field){
    validLogical(value = field$checkArguments, name = "checkArguments", validLength = 1, method = method)
  }
  if("col.index" %in% names.field){
    validCharacter(value = field$col.index, name = "col.index", validLength = 3, method = method)
  }
  if("col.NA" %in% names.field){
    validCharacter(value = field$col.NA, name = "col.NA", validLength = 1, method = method)
  }
  if("col.midplane" %in% names.field){
    validCharacter(value = field$col.midplane, name = "col.midplane", validLength = 1, method = method)
  }
  if("digit.legend" %in% names.field){
    validInteger(value = field$digit.legend, name = "digit.legend", validLength = 1, min = 0, method = method)
  }
  if("digit.result" %in% names.field){
    validInteger(value = field$digit.result, name = "digit.result", validLength = 1, min = 0, method = method)
  }
  if("digit.epsilon" %in% names.field){
    validInteger(value = field$digit.epsilon, name = "digit.epsilon", validLength = 1, min = 0, method = method)
  }
  if("digit.percentage" %in% names.field){
    validInteger(value = field$digit.percentage, name = "digit.percentage", validLength = 1, min = 0, method = method)
  }
  if("filter.index" %in% names.field){
    validCharacter(value = field$filter.index, name = "filter.index", validLength = 1, validValues = c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26"), refuse.NULL = FALSE, method = method)
  }
  if("height" %in% names.field){
    validNumeric(value = field$height, name = "height", validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
  }
  if("hemisphere" %in% names.field){
    validCharacter(value = field$hemisphere, name = "hemisphere", validLength = 1, validValues = c("both", "left", "right", "lesion", "contralateral"), method = method)
  }
  if("legend" %in% names.field){
    validCharacter(value = field$legend, name = "legend", validLength = 1, validValues = c(TRUE, FALSE, "only"), refuse.NULL = FALSE, method = method)
  }
  if("lwd.midplane" %in% names.field){
    validNumeric(value = field$lwd.midplane, name = "lwd.midplane", validLength = 1, min = 0, method = method)
  }
  if("mar" %in% names.field){
    validNumeric(value = field$mar, name = "mar", validLength = 4, min = 0, refuse.NULL = TRUE, method = method)
  }
  if("mar.legend" %in% names.field){
    validNumeric(value = field$mar.legend, name = "mar.legend", validLength = 4, min = 0, method = method)
  }
  if("mfrow" %in% names.field){
    validInteger(value = field$mfrow, name = "mfrow", validLength = 2, min = 0, refuse.NULL = FALSE, method = method)
  }
  if("mgp" %in% names.field){
    validNumeric(value = field$mgp, name = "mgp", validLength = 3, min = 0, refuse.NULL = TRUE, method = method)
  }
  if("norm_mu" %in% names.field){
    validCharacter(value = field$norm_mu, name = "norm_mu", validLength = 1, 
                 validValues = c(FALSE, "global", "global_1slice", "global_3slices", "contralateral", "contralateral_1slice", "contralateral_3slices", "default_value"),
                 method = method)
  }
  if("norm_sigma" %in% names.field){
    validCharacter(value = field$norm_sigma, name = "norm_sigma", validLength = 1, 
                 validValues = c(FALSE, "global", "global_1slice", "global_3slices", "contralateral", "contralateral_1slice", "contralateral_3slices", "default_value"),
                 method = method)
  }
  if("num.main" %in% names.field){
    validLogical(value = field$num.main, name = "num.main", validLength = 1, method = method)
  }
  if("numeric2logical" %in% names.field){
    validLogical(value = field$numeric2logical, name = "numeric2logical",  validLength = 1, method = method)
  }
  if("outline.index" %in% names.field){
    validLogical(value = field$outline.index, name = "outline.index", validLength = 1, method = method)
  }
  if("operator.withinR" %in% names.field){
    validCharacter(field$operator.withinR, name = "operator.withinR", validLength = 1, validValues = c("union","intersect"), method = method)
  }
  if("operator.betweenR" %in% names.field){
    validCharacter(field$operator.betweenR, name = "operator.betweenR", validLength = 1, validValues = c(,"union","intersect"), method = method)
  }
  if("path" %in% names.field){
    validPath(value = field$path, name = "path", method = method)
  }
  if("pch.index" %in% names.field){
    validInteger(value = field$pch.index, name = "pch.index", validLength = 3, min = 0, method = method)
  }
  if("pch.NA" %in% names.field){
    validInteger(value = field$pch.NA, name = "pch.NA", validLength = 1, min = 0, method = method)
  }
  if("pty" %in% names.field){
    validCharacter(value = field$pty, name = "pty", validLength = 1, validValues = c("m", "s"), refuse.NULL = FALSE, method = method)
  }
  if("quantiles.legend" %in% names.field){
    validLogical(value = field$quantiles.legend, name = "quantiles.legend", validLength = 1, method = method)
  }
  if("res" %in% names.field){
    validNumeric(value = field$res, name = "res", validLength = 1, refuse.NA = FALSE, method = method)
  }
  if("slice_var" %in% names.field){
    validCharacter(value = field$slice_var, name = "slice_var", refuse.duplicates = TRUE, validLength = 3, validValues = c("i", "j", "k"), method = method)
  }
  if("type.breaks" %in% names.field){
    validCharacter(value = field$type.breaks, name = "type.breaks", validLength = 1, validValues = c("range", "range_center", "quantile"), method = method)
  }
  if("unit" %in% names.field){
    validCharacter(value = field$unit, name = "unit", validLength = 1, validValues = c("px", "in", "cm", "mm"), refuse.NULL = FALSE, method = method)  
  }
  if("verbose" %in% names.field){
    validLogical(value = field$verbose, name = "verbose",  validLength = 1, method = method)
  }
  if("xlab" %in% names.field){
    validCharacter(value = field$xlab, name = "xlab", validLength = 1, method = method)
  }
  if("ylab" %in% names.field){
    validCharacter(value = field$ylab, name = "ylab", validLength = 1, method = method)
  }
  if("width" %in% names.field){
    validNumeric(value = field$width, name = "width", validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
  }
  if("window" %in% names.field){
    validCharacter(value = field$window, name = "window", validLength = 1, validValues = c(TRUE, FALSE, "eps", "pdf", "png", "svg"), refuse.NULL = FALSE, method = method)
  }
  
  
}


