
# ---------------------------------------------------------------------------- #
# function to find out the AK*pwf, apriori profiles are,
# given the lat lon time from model receptors, Dien Wu, 01/11/2017
# ---------------------------------------------------------------------------- #
### Input variables--
#' @param oco.path oco path and file for searching satellite soundings;
#' @param timestr for finding correct oco.file
#' @param recp.lat @param recp.lon numeric numbers for receptor lat/lon, one at a time
#' @param diff.td allowable thredshold for difference in lat/lon between
#'                given receptors and all satellite soundings
# ---------------------------------------------------------------------------- #
### Updates--
# replace 'timestr', 'recp.lat', 'recp.lon' with 'receptor' from 'output', DW
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# ---------------------------------------------------------------------------- #

# must have one of oco.path or oco.fn
get.oco.info = function(oco.path, receptor, diff.td = 1E-4, oco.fn = NA){

  # grabbing OCO-2 info
  if (is.na(oco.fn)) {
    timestr  = strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
    oco.file = list.files(oco.path, pattern = paste0('_', substr(timestr, 3, 8), '_'))
    if (length(oco.file) == 0) stop('No OCO file found for this timestr\n')
    fn = file.path(oco.path, oco.file)
  } else fn = oco.fn 

  oco.dat = nc_open(fn)

  # grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  # may need to reverse later
  oco.lev = ncvar_get(oco.dat, 'levels')
  oco.lat = ncvar_get(oco.dat, 'latitude')
  oco.lon = ncvar_get(oco.dat, 'longitude')

  # add OCO-3 orbit number since mission start
  oco.orbit = ncvar_get(oco.dat, 'Sounding/orbit')

  # grabbing sounding ID for STILT receptors
  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  oco.id  = as.character(ncvar_get(oco.dat, 'sounding_id'))

  # locate the OCO2 data using lat, lon, when diff are both the smallest
  diff.lat = abs(oco.lat - receptor$lati)
  diff.lon = abs(oco.lon - receptor$long)

  # try to find the closest sounding lat/lon,
  # given receptor lat/lon and allowable difference thredshold, 'diff.td'
  # choose only if lat/lon indices are the same
  loc.index = intersect(which(diff.lat < diff.td), which(diff.lon < diff.td))

  # cannot find the sounding according to receptor lat/lon
  # if so, loose 'diff.td', or check OCO-2 version, or input lat/lon
  if(length(loc.index) != 1){

    cat('get.oco.info(): cannot find the receptor lat/lon from OCO file...')
    return()

  } else {

    # if one OCO-2 sounding found for a given receptor lat/lon --
    # return the oco lat, lon, ak, pwf, apriori, profiles
    find.lat = oco.lat[loc.index]
    find.lon = oco.lon[loc.index]
    find.id  = oco.id [loc.index]

    ## grab column co2, averaging kernel, pressure weight and prior CO2 profiles
    ## dimensions--[levels, soundingID]
    ap = ncvar_get(oco.dat, 'co2_profile_apriori')[, loc.index] # in ppm
    pwf = ncvar_get(oco.dat, 'pressure_weight')[, loc.index]  # pwf
    pres = ncvar_get(oco.dat, 'pressure_levels')[, loc.index] # press in hPa

    # normalized averaging kernel (unitless)
    ak.norm = ncvar_get(oco.dat, 'xco2_averaging_kernel')[, loc.index]

    ## dimensions--[soundingID]
    xco2 = ncvar_get(oco.dat, 'xco2')[loc.index]
    grdhgt = ncvar_get(oco.dat, 'Sounding/altitude')[loc.index] # mASL
    xco2.uncert = ncvar_get(oco.dat, 'xco2_uncertainty')[loc.index] # ret err

    # total column water vapor in kg m-2, convert to mol m-2
    xh2o = ncvar_get(oco.dat, 'Retrieval/tcwv')[loc.index] * 1E3 / 18  

    # satellite footprint, 1-8
    footprint = ncvar_get(oco.dat, 'Sounding/footprint')[loc.index]
    psfc  = ncvar_get(oco.dat, 'Retrieval/psurf')[loc.index] # sfc pressure

    # check whether is missing data
    ap[ap == -999999] = NA
    pwf[pwf == -999999] = NA
    xco2[xco2 == -999999] = NA
    pres[pres == -999999] = NA
    psfc[psfc == -999999] = NA
    grdhgt[grdhgt == -999999] = NA
    ak.norm[ak.norm == -999999] = NA
    footprint[footprint == -999999] = NA
    xco2.uncert[xco2.uncert == -999999] = NA

    # assign vertical dimnames
    attributes(ap)$names      = oco.lev
    attributes(pwf)$names     = oco.lev
    attributes(pres)$names    = oco.lev
    attributes(ak.norm)$names = oco.lev

    ### combine all OCO-2 vertical profiles and other 1D variables
    all.info = list(oco.id = find.id, oco.orbit = oco.orbit, oco.lat = find.lat, 
                     oco.lon = find.lon, ak.norm = ak.norm, pwf = pwf, 
                     pres = pres, ap = ap, oco.grdhgt = grdhgt, oco.psfc = psfc, 
                     oco.foot = footprint, oco.xco2 = xco2, 
                     oco.xco2.uncert = xco2.uncert, oco.xh2o = xh2o)
    
    nc_close(oco.dat)
    all.info      # return both profiles and other retrivals
  }
  
} # end of subroutine
