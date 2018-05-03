### subroutine to find out the AK*PW, apriori profiles are,
# given the lat lon time from model receptors
# written by Dien Wu, 01/11/2017

### Input variables--
# ocopath, ocofile: oco2 path and file for searching satellite soundings;
# recp.lat, recp.lon: numeric numbers for receptor lat/lon
# diff.td: allowable thredshold for difference in lat/lon between
#          given receptors and all satellite soundings

### updates--
# allow for vector form of recp.lat and recp.lon, DW, 05/02/2018

get.oco2info <- function(ocopath, ocofile, recp.lat, recp.lon, diff.td=1E-4){

# load libraries
library(ncdf4)

# grabbing OCO-2 info
ocodat <- nc_open(file.path(ocopath,ocofile))

# grabbing OCO-2 levels, lat, lon
# level 1 to 20, for space-to-surface, level 20 is the bottom level
# may need to reverse later
oco.level <- ncvar_get(ocodat, "levels")
oco.lat <- ncvar_get(ocodat, "latitude")
oco.lon <- ncvar_get(ocodat, "longitude")

# grabbing warn levels
warnlevel <- ncvar_get(ocodat, "warn_level")

# grabbing time for STILT receptors
# YYYY MM DD HH mm ss m (millisecond) f (footprint)
id <- as.character(ncvar_get(ocodat, "sounding_id"))

all.prof <- list(); all.info <- NULL

# loop over all receptors
for(r in 1:length(recp.lat)){

  # locate the OCO2 data using lat, lon, when diff are both the smallest
  diff.lat <- abs(oco.lat - recp.lat[r])
  diff.lon <- abs(oco.lon - recp.lon[r])

  # try to find the closest sounding lat/lon,
  # given receptor lat/lon and allowable difference thredshold, "diff.td"
  lat.index <- which(diff.lat < diff.td)
  lon.index <- which(diff.lon < diff.td)

  # only if lat/lon indices are the same
  loc.index <- intersect(lat.index, lon.index)

  # cannot find the sounding according to receptor lat/lon
  # if so, loose "diff.td", or check OCO-2 version, or input lat/lon
  if(length(loc.index)!=1){
    cat("get.oco2info(): cannot find the receptor lat/lon from OCO-2 file...")
    next
  }

  # return the oco2 lat, lon, ak, pw, apriori, profiles
  find.lat <- oco.lat[loc.index]
  find.lon <- oco.lon[loc.index]
  find.id  <- id[loc.index]

  ## grab column co2, averaging kernel, pressure weight and prior CO2 profiles,

  # dimensions--[levels, soundingID]
  # normalized averaging kernel (unitless)
  ak.norm <- ncvar_get(ocodat, "xco2_averaging_kernel")[, loc.index]

  # pressure weighting (unitless)
  pw <- ncvar_get(ocodat, "pressure_weight")[, loc.index]

  # pressure in hPa
  pres <- ncvar_get(ocodat, "pressure_levels")[, loc.index]

  # CO2.apriori in ppm
  apriori <- ncvar_get(ocodat, "co2_profile_apriori")[, loc.index]

  # dimensions--[soundingID]
  # ground height measured in OCO-2 in meter ASL
  grdhgt <- ncvar_get(ocodat, "Sounding/altitude")[loc.index]

  # retrieved XCO2 and its uncertainty
  xco2 <- ncvar_get(ocodat, "xco2")[loc.index]
  xco2.uncert <- ncvar_get(ocodat, "xco2_uncertainty")[loc.index]

  # satellite footprint
  footprint<-ncvar_get(ocodat,"Sounding/footprint")[loc.index]

  #t_700 <- ncvar_get(ocodat, "Retrieval/T700")[loc.index] # temp at 700mb
  psfc<-ncvar_get(ocodat,"Retrieval/psurf")[loc.index] # retrieved sfc pressure

  # check whether is missing data
  pw[pw == -999999] <- NA
  xco2[xco2 == -999999] <- NA
  pres[pres == -999999] <- NA
  psfc[psfc == -999999] <- NA
  #t_700[t_700 == -999999] <- NA
  grdhgt[grdhgt == -999999] <- NA
  apriori[apriori == -999999] <- NA
  ak.norm[ak.norm == -999999] <- NA
  footprint[footprint==-999999]<-NA
  xco2.uncert[xco2.uncert == -999999]<-NA

  # assign vertical dimnames
  attributes(ak.norm)$names <- oco.level
  attributes(pw)$names <- oco.level
  attributes(pres)$names <- oco.level
  attributes(apriori)$names<-oco.level

  ### combine all OCO-2 vertical profiles
  tmp.prof <- data.frame(ak.norm, pw, pres, apriori)  # vertical profiles

  # combine all variables, just one number
  tmp.info <- data.frame(grdhgt, psfc, footprint, xco2, xco2.uncert, find.id,
                         find.lat, find.lon)

  all.prof[[r]] <- tmp.prof
  all.info <- rbind(all.info, tmp.info)
}

# return both profiles and variables
return(list(all.prof, all.info))

} # end of subroutine
