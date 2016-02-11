#### array2dt ####

#### matrix2array ####

set.seed(10)
n <- 4
Y <- rnorm(n^2)

## conversion
res1 <- df2array(expand.grid(1:n + 0.5, 1:n + 0.5), contrast = Y)
res2 <- df2array(expand.grid(1:n + 0.5, 1:n + 0.5), contrast = Y, format = "matrix")
res3 <- df2array(expand.grid(2 * (1:n), 2 * (1:n)), contrast = Y)
res4 <- df2array(expand.grid(2 * (1:n), 2 * (1:n)), contrast = cbind(Y ,Y, Y), range.coords = c(10,10))

## display
par(mfrow = c(2,2), mar = rep(2,4), mgp = c(1.5,0.5,0))
fields::image.plot(unique(res1$coords[[1]]), unique(res1$coords[[2]]), res1$contrast[[1]],
                   xlab = "", ylab = "")
fields::image.plot(unique(res2$coords[[1]]), unique(res2$coords[[2]]), res2$contrast,
                   xlab = "", ylab = "")
fields::image.plot(res3$contrast[[1]])
fields::image.plot(res4$contrast[[2]])



#### df2array ####

#### dt2array ####
