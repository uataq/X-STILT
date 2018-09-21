# calculate the areas of grid box,lon & lat, in unit of m2
# input the lat, lon (lower left) and spatial res from the main script
# written by Dien Wu

# better create a matrix with lon, lat that matches the matrix we used
# use third.dim for the additional dimension, mostly like time, or date.
# use raster layer instead of array()

area <- function(res, start.lon, start.lat, end.lon, end.lat,
	third.dim = NULL){

	library(geosphere)

	# input lat lon are lower left lat lon, thus, need one more lat,
	# for calculating the last area; but in this way, there will be one more area
	# (NA) for an additional lower left latitude
	end.lon <- end.lon
	end.lat <- end.lat + res
	lon <- seq(start.lon, end.lon, res)
	lat <- seq(start.lat, end.lat, res)

	dlat <- 0
	dlon <- 0	# initialize
	area <- array(0, dim = c(length(lon), length(lat), third.dim))
	colnames(area) <- lat
	rownames(area) <- lon

	# calculate the distance in latitudinal direction between two coordinates
	dlat <- distCosine(c(lon[1], lat[1]), c(lon[1], lat[2]))	#dlat all the same

	# calculate the distance in longitudinal direction between two coordinates,
	# which depends on latitude
	# dlon is the same along the same latitude
	# dlon is smaller as moving towards the pole
	for (i in 1:length(lat)) {
		dlon[i] <- distCosine(c(lon[1], lat[i]), c(lon[2], lat[i]))
	}

	# calculate the area of polygon
	for (i in 1:length(lat)){
		if (length(dim(area)) == 2) {
			area[,i] <- dlat * (dlon[i+1] + dlon[i]) / 2
		} else {
			area[,i,] <- dlat * (dlon[i+1] + dlon[i]) / 2
		}
	}

	# since we calculate one more area,
	# select the actual area according to the length of latitude
	if (length(dim(area)) == 2) {
		area <- area[,1:(length(lat) - 1)]
	} else {
		area <- area[,1:(length(lat) - 1),]
	}

  return(area)
} # end subroutine
