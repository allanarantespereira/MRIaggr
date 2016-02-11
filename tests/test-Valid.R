#### validDimension ####

validDimension(value1 = 1:5, value2 = 1:5, method = "test")
validDimension(value1 = 1:5, value2 = 1:6, method = "test")
validDimension(value1 = 1:5, value2 = 1:5, type = "length" , method = "test")
validDimension(value1 = 1:5, value2 = 1:6, type = "length" , method = "test")

validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), type = "length" , method = "test")
validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), type = "nrow" , method = "test")
validDimension(value1 = matrix(1:5), value2 = matrix(1,2,2), method = "test")
validDimension(value1 = matrix(1:4,2,2), value2 = matrix(1,2,2), method = "test")