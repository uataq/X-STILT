#' script to define spatial domain
#' @author Dien Wu, 07/05/2018

#' @update:
#' use geocode and SpatialPoints to find lat/lon coordinates,
#' country and reg name, DW, DR, 08/15/2018

#' site can be a vector
get.lon.lat <- function(site, dlon, dlat, city.loc = NULL) {

  # spatial domains placing receptors and city center, help select OCO-2 data
  library(ggmap); library(rworldmap); library(sp)

  # location name to lon, lat coordinates
  if (is.null(city.loc)) 
    city.loc <- geocode(location = site, output = 'latlon', source = 'google',
      override_limit = T)

  # from https://stackoverflow.com/questions/21708488/
  # get-country-and-continent-from-longitude-and-latitude-point-in-r
  # use high res map from rworldxtra if you were concerned about detail
  countriesSP <- getMap(resolution = 'low')

  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP <- SpatialPoints(city.loc, proj4string = CRS(proj4string(countriesSP)))

  # use 'over' to get indices of the Polygons object containing each point
  indices <- over(pointsSP, countriesSP)
  lon.lat <- data.frame(
    minlon = city.loc$lon - dlon, maxlon = city.loc$lon + dlon,
    minlat = city.loc$lat - dlat, maxlat = city.loc$lat + dlat,
    citylon = city.loc$lon, citylat = city.loc$lat, cityid = site,
    countryid = indices$ADMIN, regid = indices$continent, iso3 = indices$ISO3)

  return(lon.lat)
}
