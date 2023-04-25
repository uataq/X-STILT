#' script to define spatial domain
#' @author Dien Wu, 07/05/2018

#' @update:
#' use geocode and SpatialPoints to find lat/lon coordinates,
#' country and reg name, DW, DR, 08/15/2018
#' site can be a vector

get.lon.lat = function(site, dlon, dlat, site.loc = NULL, api.key = NULL) {

  library(ggmap); library(rworldmap); library(sp); library(lutz)

  # location name to lon, lat coordinates
  if (is.null(site.loc)) {
  
    # Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
    # google API can be obtained from https://console.developers.google.com/
    # *** if the site lat/lon is known, no need for API
    if (is.null(api.key)) api.key = readLines('insert_ggAPI.csv')
    if (api.key == '') stop('Missing googleAPI for getting lat/lon of your site; put your API in insert_ggAPI.csv or mannually insert lat/lon, see @param site.loc...\n')
    register_google(key = api.key)
    site.loc = geocode(location = site, output = 'latlon', source = 'google',
                       override_limit = T)
  }

  # from https://stackoverflow.com/questions/21708488/
  # get-country-and-continent-from-longitude-and-latitude-point-in-r
  # use high res map from rworldxtra if you were concerned about detail
  countriesSP = getMap(resolution = 'low')

  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(site.loc, proj4string =CRS(proj4string(countriesSP)))
  indices = over(pointsSP, countriesSP)

  state_df = data.frame(state = NA, state_abb = NA)
  # if (as.character(indices$ISO3) == 'USA') {

  #   # get country, state, zipcode, and timezone
  #   state_df = site.loc %>% 
  #              tidygeocoder::reverse_geocode(lat = lat, long = lon, 
  #                                            address = site, 
  #                                            full_results = T) %>%
  #              dplyr::select(state) %>% 
  #              mutate(state_abb = state.abb[state.name %in% state])
  # }

  
  # get time.zone 
  tz = tz_lookup_coords(site.loc$lat, site.loc$lon)

  # convert indices from factors to characters
  lon.lat = data.frame(site_id = site, site_lon = site.loc$lon, 
                       site_lat = site.loc$lat, tz = tz, 
                       countryid = as.character(indices$ADMIN), 
                       iso3 = as.character(indices$ISO3),
                       regid = as.character(indices$continent), 
                       state = state_df$state, state_abb = state_df$state_abb, 
                       minlon = site.loc$lon - dlon, 
                       maxlon = site.loc$lon + dlon,
                       minlat = site.loc$lat - dlat, 
                       maxlat = site.loc$lat + dlat, 
                       stringsAsFactors = F)

  return(lon.lat)
}
