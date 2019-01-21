### subroutine to find out the AK*pwf, apriori profiles are,
# given the lat lon time from model receptors
# written by Dien Wu, 01/11/2017

### Input variables--
# oco2.path: oco2 path and file for searching satellite soundings;
# timestr: for finding correct oco2.file
# recp.lat, recp.lon: numeric numbers for receptor lat/lon, one at a time
# diff.td: allowable thredshold for difference in lat/lon between
#          given receptors and all satellite soundings

### Updates--
# replace 'timestr', 'recp.lat', 'recp.lon' with 'receptor' from 'output', DW

get.oco2.info <- function(oco2.path, receptor, diff.td = 1E-4){

  # grabbing OCO-2 info
  timestr   <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  oco2.file <- list.files(path = oco2.path, pattern = substr(timestr, 3, 8))
  oco2.dat  <- nc_open(file.path(oco2.path, oco2.file))

  # grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  # may need to reverse later
  oco2.lev <- ncvar_get(oco2.dat, 'levels')
  oco2.lat <- ncvar_get(oco2.dat, 'latitude')
  oco2.lon <- ncvar_get(oco2.dat, 'longitude')

  # grabbing sounding ID for STILT receptors
  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  oco2.id  <- as.character(ncvar_get(oco2.dat, 'sounding_id'))

  # locate the OCO2 data using lat, lon, when diff are both the smallest
  diff.lat <- abs(oco2.lat - receptor$lati)
  diff.lon <- abs(oco2.lon - receptor$long)

  # try to find the closest sounding lat/lon,
  # given receptor lat/lon and allowable difference thredshold, 'diff.td'
  # choose only if lat/lon indices are the same
  loc.index <- intersect(which(diff.lat < diff.td), which(diff.lon < diff.td))

  # cannot find the sounding according to receptor lat/lon
  # if so, loose 'diff.td', or check OCO-2 version, or input lat/lon
  if(length(loc.index) != 1){

    cat('get.oco2info(): cannot find the receptor lat/lon from OCO-2 file...')
    return()

  } else {

    # if one OCO-2 sounding found for a given receptor lat/lon --
    # return the oco2 lat, lon, ak, pwf, apriori, profiles
    find.lat <- oco2.lat[loc.index]
    find.lon <- oco2.lon[loc.index]
    find.id  <- oco2.id [loc.index]

    ## grab column co2, averaging kernel, pressure weight and prior CO2 profiles
    ## dimensions--[levels, soundingID]
    ap <- ncvar_get(oco2.dat, 'co2_profile_apriori')[, loc.index] # in ppm
    pwf <- ncvar_get(oco2.dat, 'pressure_weight')[, loc.index]  # pwf
    pres <- ncvar_get(oco2.dat, 'pressure_levels')[, loc.index] # press in hPa

    # normalized averaging kernel (unitless)
    ak.norm <- ncvar_get(oco2.dat, 'xco2_averaging_kernel')[, loc.index]

    ## dimensions--[soundingID]
    xco2 <- ncvar_get(oco2.dat, 'xco2')[loc.index]
    grdhgt <- ncvar_get(oco2.dat, 'Sounding/altitude')[loc.index] # mASL
    xco2.uncert <- ncvar_get(oco2.dat, 'xco2_uncertainty')[loc.index] # ret err

    # satellite footprint, 1-8
    footprint <- ncvar_get(oco2.dat, 'Sounding/footprint')[loc.index]
    psfc  <- ncvar_get(oco2.dat, 'Retrieval/psurf')[loc.index] # sfc pressure

    # check whether is missing data
    ap[ap == -999999] <- NA
    pwf [pwf == -999999]  <- NA
    xco2[xco2 == -999999] <- NA
    pres[pres == -999999] <- NA
    psfc[psfc == -999999] <- NA
    grdhgt[grdhgt == -999999] <- NA
    ak.norm[ak.norm == -999999] <- NA
    footprint[footprint == -999999] <- NA
    xco2.uncert[xco2.uncert == -999999] <- NA

    # assign vertical dimnames
    attributes(ap)$names      <- oco2.lev
    attributes(pwf)$names     <- oco2.lev
    attributes(pres)$names    <- oco2.lev
    attributes(ak.norm)$names <- oco2.lev

    ### combine all OCO-2 vertical profiles and other 1D variables
    all.info <- list(oco2.id = find.id, oco2.lat = find.lat,
                     oco2.lon = find.lon, ak.norm = ak.norm, pwf = pwf,
                     pres = pres, ap = ap, oco2.grdhgt = grdhgt,
                     oco2.psfc = psfc, oco2.foot = footprint,
                     oco2.xco2 = xco2, oco2.xco2.uncert = xco2.uncert)
    nc_close()
    all.info      # return both profiles and other retrivals
  }
  
} # end of subroutine
