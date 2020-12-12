# subroutine to deal with variance of a sum
# written by DW, 03/20/2018

# inputs are correlation length in meters or time scales in second,
# weighting functions w (default is NULL, can be no weighting func),
# and identity matrix for SD of different error sources,
# output is a matrix
# L should be in meters or second, can be NULL, meaning no correlation considered
# x should have column names of lon, lat (if type == 'hor');
#                             hgt (if for 'ver' type); or dsec (for 'time' type)
# generalize and clear up codes, DW, 07/23/2018

cal.var.cov <- function(sd, L = NULL, w = NULL, x,
                        type = c('hor', 'ver', 'time', NULL)[4]){

  # create a diagonal grid, either spatial or temporal with ndim
  library(geosphere)
  num <- seq(1, nrow(x))

  # create identity matrix for error SD
  sdd <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))
  diag(sdd) <- sd  # fill diagonal elements with SD in ppm

  # calculate the correlation matrix, from 0 to 1, L in meters
  if (!is.null(L) & !is.null(type)) {
    grid <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))

    # fill grid with separation distance or time
    for (n in 1:length(num)){

      if (type == 'hor') {
        # calculate the horizontal distance between points first
        point1 <- data.frame(lon = x[n, 'lon'], lat = x[n, 'lat'])
        point2 <- data.frame(lon = x[, 'lon'],  lat = x[, 'lat'])
        grid[n, ] <- distCosine(p1 = point1, p2 = point2) # in meters
      } # end if type == hor

      # calculate the vertical distance between levels first
      if (type == 'ver') grid[n, ] <- abs(x$hgt - x$hgt[n])

      # calculate the separation time between each two soundings
      if (type == 'time') grid[n, ] <- abs(x$dsec - x$dsec[n])
    } # end for n

    corr <- exp(- grid / L)  # calculate error correlation

    # matrix multiplications,
    # to calculate off-diagonal elements based on correlation and variances
    prod <- sdd %*% corr %*% t(sdd)

  } else {
    prod <- sdd %*% t(sdd)
  }

  # now create the weighting matrix
  if (!is.null(w)) {
    w2 <- array(0, dim = c(max(num), max(num)), dimnames = list(num, num))

    for(n in 1:length(num)) w2[n,] <- w[n] * w
    prod.wgt <- w2 * prod
    prod <- prod.wgt
  }

  #tot.sd <- sqrt(sum(prod, na.rm=T))   # total weighted aggregated SD
  return(prod)
}
