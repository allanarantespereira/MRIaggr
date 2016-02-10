#**********************************************************************
#**********************************************************************
#*************         validation functions         *******************
#**********************************************************************
#**********************************************************************
#

# validCharacter:   
# validClass:      
# validDimension:    
# validInteger:       
# validLogical:     
# validNames:     
# validNumeric:  
# validPath:    

#### validCharacter ####
validCharacter <- function(value, name = as.character(substitute(value)), validLength, 
                           validValues = "character", refuse.NULL = TRUE, refuse.duplicates = FALSE, method){
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", n.value, "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' contains duplicated values : \n",        
           "\"",paste(unique(value[duplicated(value)]), collapse = "\" \""), "\" \n")
    }
    
    #### check values
    if(identical(validValues,"character")){
      
      if(any(is.character(value) == FALSE)){
        stop(method, " : wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character"}else{"vector of characters"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(identical(validValues,"character_or_logical")){
      
      if(any( (is.character(value) == FALSE) * (is.logical(value) == FALSE) > 0 )){
        stop(method, " : wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character or logical"}else{"vector of characters or logicals"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = "\" \""), "\"\n")
      
    }
    
  }
  
}

#### validClass ####
validClass <- function(value, name = as.character(substitute(value)), validClass, 
                       superClasses = TRUE, method){
  
  if(superClasses == TRUE){
    
    if( all(is(value) %in% validClass == FALSE) ){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "superclass of \'", name, "\' must be one of the following \"", paste(validClass,collapse="\" \""), "\"  \n", 
           "proposed superclass : \"", paste(is(value),collapse="\" \""), "\" \n")
    }  
    
  }else{
 
    if( class(value) %in% validClass == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "class of \'", name, "\' must be \"", paste(validClass,collapse="\" \""),"\"  \n", 
           "proposed class : ", class(value)[[1]], "\n")
    }  
    
  }
  
  
  
}

#### validDimension ####
validDimension <- function(value1, value2 = NULL, name1 = as.character(substitute(value1)), name2 = as.character(substitute(value2)),
                           validDimension = NULL,
                           type = c("NROW","NCOL"), method){
  
  n.type <- length(type)
  
  #### dimension 1
  testDimension <- sapply(1:n.type, function(x){
    do.call(type[x], list(value1))
  })
  
  #### dimension 2
  
  
  if(is.null(validDimension)){
    
    validDimension <- sapply(1:n.type, function(x){
      do.call(type[x], list(value2))
    })
    test.validDimension <- TRUE
    
  }else if(is.null(name2)){
    
    test.validDimension <- FALSE
    
  }else{
    
    test.validDimension <- TRUE
    
  }
  
  #### main
  for(iter_type in 1:n.type){
    
    if(testDimension[iter_type] != validDimension[iter_type]){
      
      if(test.validDimension){
        stop(method, " : dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")  
      }else{
        stop(method, " : dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")
        
      }
      
    }
    
  }
    
  }
  
#### validInteger ####
validInteger <- function(value, name = as.character(substitute(value)), validLength, 
                         validValues = NULL, min = NULL, max = NULL, 
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method){
  
  validNumeric(value = value, name = name, validLength = validLength, min = min, max = max, 
               refuse.NA = refuse.NA, refuse.NULL = refuse.NULL, refuse.duplicates = refuse.duplicates, method = method)
  
  #### check integer
  if(any(value %% 1 > 0)){
    stop(method, " : wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must contain integers not doubles \n",        
         "invalid value(s) in ", name, " : ", paste(value[value %% 1 > 0], collapse = " "), "\n")
  }
  
}

#### validLogical ####
validLogical <- function(value, name = as.character(substitute(value)), validLength, 
                         refuse.NULL = TRUE, refuse.NA = TRUE, method){
  
  
  if(is.null(value)){
    
    #### NULL
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NA == FALSE){"or NA"}," and not NULL \n")
    }
    
  }else{ 
    
    #### Size
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    } 
    
    #### Type
    if(any(is.logical(value) == FALSE)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be ", if(refuse.NULL == FALSE){"NULL or "}, if(refuse.NA == FALSE){"NA or "},"TRUE or FALSE \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    if(refuse.NA == TRUE && any(is.na(value)) ){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NULL == FALSE){"or NULL"}," and not NA \n")
    }
    
  }
  
}

#### validNames ####
validNames <- function(value, name = as.character(substitute(value)), refuse.NULL = TRUE,
                       validLength = NULL, validValues = NULL, requiredValues = NULL, forbiddenValues = NULL,
                       method){
  
  ## type
  if(is.matrix(value)){
    value <- colnames(value)
  }
  
  if(is.data.table(value) || is.data.frame(value) || is.list(value)){
    value <- names(value)
  }
  
  ## tests
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
    stop(method, " : wrong specification of \'", name, "\' \n", 
         "names of \'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have ", paste(validLength, collapse = " or ")," names  \n", 
           "length(names(", name, ")) : ", n.value, "\n")
    }
    
    #### check content
    
    if(!is.null(requiredValues) && any(requiredValues %in% value == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must contains specific names \n",
           "missing names : \"",paste(requiredValues[requiredValues %in% value == FALSE], collapse = "\" \""),"\" \n", 
           "proposed names : \"", paste(value, collapse = "\" \""), "\"\n")  
      
    }
    
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid names for \'", name, "\' : \"",paste(validValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")  
      
    }
    
    if(!is.null(forbiddenValues) && any(value %in% forbiddenValues)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "forbidden names for \'", name, "\' : \"",paste(forbiddenValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value[value %in% forbiddenValues], collapse = " "), "\"\n")  
      
    }
    
    if(any(duplicated(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           name, " must not contain duplicated names \n", 
           "duplicated names : \"", paste(value[duplicated(value)], collapse = " "), "\"\n")  
    }
    
  }
  
}
#### validNumeric ####
validNumeric <- function(value, name = as.character(substitute(value)), validLength,
                         validValues = NULL , min = NULL, max = NULL,
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method){
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check length
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    }
    
    #### check NA
    if(refuse.NA == TRUE && any(is.na(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not contain NA \n", 
           "index of NA values : ", which(paste(is.na(value), collapse = " ")), "\n")
    }
    
    #### check numeric
    if(any( (is.numeric(value) == FALSE) * (is.na(value) == FALSE) )){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be a numeric \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' contains duplicated values : \n",        
           paste(unique(value[duplicated(value)]), collapse = " "), "\n")
    }
    
    #### check min value
    if(!is.null(min) && any(stats::na.omit(value) < min)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be bigger than ", min, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) < min], collapse = " "), "\n")
    }
    
    #### check max value
    if(!is.null(max) && any(stats::na.omit(value) > max)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be smaller than ", max, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) > max], collapse = " "), "\n")
    }
    
    #### check valid values
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
  }
}

#### validPath ####
validPath <- function(value, name = as.character(substitute(value)),
                      method){

  if(!is.null(value)){
    try_path <- dir.exists(value)
    
    if(any(is(try_path) == "try-error")){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "proposed ", name, " : \"", value, "\" \n", 
           "error : ", paste(try_path, collapse = " "), 
           "current ", name, " : ", getwd(), "\n")
    }
    
    
    if(substr(value, start = nchar(value), stop = nchar(value)) != "/"){
      warning(method, " : possible bad specification of \'", name, "\' \n", 
              "\'", name, "\' should end with a fsep (e.g. \"/\") \n", 
              "proposed ", name, " : ", value, "\n")
    }
  }
}
