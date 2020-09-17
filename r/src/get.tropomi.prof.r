# function to find out the AK profiles from TROPOMI, 
# will find the most adjacent TROPOMI grid that receptor falls into 
# DW, 08/25/2020 

# add other species like NO2, DW, 09/02/2020
get.tropomi.prof <- function(tropomi.path, tropomi.speci, receptor) {

  # grabbing OCO-2 info
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  lon.lat <- data.frame(minlon = receptor$long - 1, maxlon = receptor$long + 1, 
                        minlat = receptor$lati - 1, maxlat = receptor$lati + 1)
  tropomi.info <- find.tropomi(tropomi.path, substr(timestr, 1, 8), lon.lat)
  if (length(tropomi.info$fn) == 0) stop('No TROPOMI file found for this timestr\n')



  # --------------------- load all TROPOMI grid center lat/lon
  tropomi.dat <- nc_open(file.path(tropomi.path, tropomi.info$fn))
  indx_along_track  <- ncvar_get(tropomi.dat, 'PRODUCT/scanline')
  indx_across_track <- ncvar_get(tropomi.dat, 'PRODUCT/ground_pixel')
  tropomi.lat <- ncvar_get(tropomi.dat, 'PRODUCT/latitude')
  tropomi.lon <- ncvar_get(tropomi.dat, 'PRODUCT/longitude')
  dimnames(tropomi.lat) <- dimnames(tropomi.lon) <- list(indx_across_track, indx_along_track)

  var_array <- array(  
    data = do.call(cbind, list(tropomi.lon, tropomi.lat)), dim = c(dim(tropomi.lon), 2), 
    dimnames = c(dimnames(tropomi.lon), list(c('lon', 'lat')))
  )   # assuming all matrices in the list have equal dimensions
  
  var_df <- dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
  colnames(var_df)[1:2] <- c('indx_across_track', 'indx_along_track')



  # ----------- locate the TROPOMI grid that is most close to X-STILT receptor 
  # average observations according to coordinates of limited receptors
  # calculate the spherical distance between a given obs and all receptors
  # and locate the one receptor that has the shortest distance
  tmp.dist <- geosphere::distCosine(var_df[, c('lon', 'lat')], 
                                    data.frame(lon = receptor$long, lat = receptor$lati)) # in meters
  cat(paste('get.tropomi.prof(): find a TROPOMI grid;', signif(min(tmp.dist) / 1E3, 3),
            'km away from receptor\n'))
  tmp.grid <- var_df[tmp.dist == min(tmp.dist), ]
  along.indx <- which(indx_along_track == tmp.grid$indx_along_track)
  across.indx <- which(indx_across_track == tmp.grid$indx_across_track)
  loc.indx <- c(across.indx, along.indx) 



  # ------------ grabbing 50 TROPOMI CO vertical levels or 34 NO2 levels
  time <- substr(ncvar_get(tropomi.dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
  timestr <- format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')[along.indx]
  
  # these matrics all have dims of [across, along]
  qa   <- ncvar_get(tropomi.dat, 'PRODUCT/qa_value')[across.indx, along.indx]  # quality assurances
  hsfc <- ncvar_get(tropomi.dat, 'INPUT_DATA/surface_altitude')[across.indx, along.indx]
  psfc <- ncvar_get(tropomi.dat, 'INPUT_DATA/surface_pressure')[across.indx, along.indx] / 100 # Pa to hPa
  sfc_class <- ncvar_get(tropomi.dat, 'INPUT_DATA/surface_classification')[across.indx, along.indx]

  ### combine all OCO-2 vertical profiles and other 1D variables
  find.lat <- tropomi.lat[across.indx, along.indx]
  find.lon <- tropomi.lon[across.indx, along.indx]



  # ------------- grabbing column concentrations of mol m-2 and AK profiles 
  if (tropomi.speci == 'CO') {
    
    xco <- ncvar_get(tropomi.dat, 'PRODUCT/carbonmonoxide_total_column')[across.indx, along.indx] # mol m-2
    xco.uncert <- ncvar_get(tropomi.dat, 'PRODUCT/carbonmonoxide_total_column_precision')[across.indx, along.indx]
    xh2o <- ncvar_get(tropomi.dat, 'DETAILED_RESULTS/water_total_column')[across.indx, along.indx]

    # Pressure at bottom of layer in hPa, from TOA to the surface
    lower.pres <- ncvar_get(tropomi.dat, 'DETAILED_RESULTS/pressure_levels')[, across.indx, along.indx] / 100 
    
    # compute pressure for upper level in hPa, press from TOA to sfc
    upper.pres <- c(lower.pres[1] * lower.pres[1] / lower.pres[2], 
                    lower.pres[1 : length(lower.pres) - 1] )	

    # zagl in meters
    zagl <- as.numeric(ncvar_get(tropomi.dat, 'PRODUCT/layer'))  

    # column AK in meter, 
    ak.co <- ncvar_get(tropomi.dat, 'DETAILED_RESULTS/column_averaging_kernel')[, across.indx, along.indx]

    # store as a list
    all.info <- list(tropomi.lat = find.lat, tropomi.lon = find.lon, ak = ak.co, 
                     lower.pres = lower.pres, upper.pres = upper.pres, agl = zagl, 
                     tropomi.zsfc = hsfc, tropomi.psfc = psfc, xco = xco, 
                     xco.uncert = xco.uncert, xh2o = xh2o)
  

  } else if (tropomi.speci == 'NO2') {  # grab NO2

    xno2.tropo <- ncvar_get(tropomi.dat, 'PRODUCT/nitrogendioxide_tropospheric_column')[across.indx, along.indx]
    xno2.tropo.uncert <- ncvar_get(tropomi.dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision')[across.indx, along.indx]
    
    # Water vapor and water liquid slant column density, in mol m-2
    xh2ov.slant <- ncvar_get(tropomi.dat, 'DETAILED_RESULTS/water_slant_column_density')[across.indx, along.indx]
    xh2ol.slant <- ncvar_get(tropomi.dat, 'DETAILED_RESULTS/water_liquid_slant_column_density')[across.indx, along.indx]

    ### compute 34 TROPOMI pressure levels using hybrid_sigma_pressure_coordinate
    # p(t, k, j, i, l) = ap(k, l) + b(k, l) *  ps(t, j, i); i, j for locations
    # ap: tm5_constant_a; b: tm5_constant_b; ps: surface air pressure
    # l=0 for base of layer, l=1 for top of layer.

    # TM5 hybrid A coefficient at upper and lower interface levels in Pa, dim of [34, 2]
    a <- ncvar_get(tropomi.dat, 'PRODUCT/tm5_constant_a') / 100   # convert to hPa now
    
    # TM5 hybrid B coefficient at upper and lower interface levels, unitless, dim of [34, 2]
    b <- ncvar_get(tropomi.dat, 'PRODUCT/tm5_constant_b')   
    
    # k from surface to top of atmosphere, 0 to 33 
    k <- ncvar_get(tropomi.dat, 'PRODUCT/layer')      
    
    # vertices: lower or upper pressure level boundaries, v = 0 for lower level
    v <- ncvar_get(tropomi.dat, 'PRODUCT/vertices')  

    # finally compute pressures for lower and upper level boundaries, in hPa now
    pres <- (a + b * psfc); dimnames(pres) <- list(c('lower', 'upper'), k)  
    pres.df <- as.data.frame(t(pres))

    # NO2 averaging kernel from TOA to surface, unitless, no need for any conversion
    ak.norm.no2 <- ncvar_get(tropomi.dat, 'PRODUCT/averaging_kernel')[, across.indx, along.indx]

    # store as a list
    # reverse TROPOMI NO2 pressure levels, now from TOA to sfc
    all.info <- list(tropomi.lat = find.lat, tropomi.lon = find.lon, 
                     ak = rev(ak.norm.no2), lower.pres = rev(pres.df$lower), 
                     upper.pres = rev(pres.df$upper), tropomi.zsfc = hsfc, 
                     tropomi.psfc = psfc, xno2.tropo = xno2.tropo, 
                     xno2.tropo.uncert = xno2.tropo.uncert, 
                     xh2ov.slant = xh2ov.slant, xh2ol.slant = xh2ol.slant)

  } else {
    stop('get.tropomi.prof(): please add codes for grabbing other TROPOMI species other than CO and NO2\n')
  }

  nc_close(tropomi.dat)
  all.info      # return both profiles and other retrivals
  
} # end of subroutine



