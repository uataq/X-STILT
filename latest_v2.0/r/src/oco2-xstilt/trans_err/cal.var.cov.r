# subroutine to deal with variance of a sum
# written by DW, 03/20/2018

# inputs are correlation length in meters or time scales in second,
# weighting functions w (default is NULL, can be no weighting func),
# and identity matrix for SD of different error sources,
# output is a matrix and a sum

# generalize and clear up codes, DW, 07/23/2018

cal.var.cov <- function(sd, L, w = NULL, x, type = c('hor', 'ver', 'time')){

  # create a diagonal grid, either spatial or temporal with ndim
  library(geosphere)
  num <- seq(1, nrow(x))
  grid <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))

  # fill grid with separation distance or time
  for (n in 1:length(num)){
    if (type == 'hor') {
      # calculate the horizontal distance between points first
      point1 <- data.frame(lon = x[n,'lon'], lat = x[n,'lat'])
      point2 <- data.frame(lon = x[ ,'lon'], lat = x[ ,'lat'])
      grid[n, ] <- distCosine(p1 = point1, p2 = point2) # in m
    }

    # calculate the vertical distance between levels first
    if (type == 'ver') grid[n, ] <- abs(x$hgt - x$hgt[n])

    # calculate the separation time between each two soundings
    if (type == 'time') grid[n, ] <- abs(x$dsec - x$dsec[n])
  }

  # calculate the correlation matrix, from 0 to 1, L in meters
  corr <- exp(-grid / L)

  # create identity matrix for error SD
  sdd <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))
  diag(sdd) <- sd  # fill diagonal elements with SD in ppm

  # matrix multiplications,
  # to calculate off-diagonal elements based on correlation and variances
  prod <- sdd %*% corr %*% t(sdd)

  # now create the weighting matrix
  if (length(w) > 0) {
    wgt <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))
    for(n in 1:length(num)){wgt[n,] <- w[n] * w}

    prod.wgt <- wgt * prod
    prod <- prod.wgt
  }

  tot.sd <- sqrt(sum(prod, na.rm=T))   # total weighted aggregated SD

  return(prod.wgt)
}
