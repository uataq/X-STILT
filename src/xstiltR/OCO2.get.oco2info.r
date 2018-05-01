# subroutine to find out the AK*PW, apriori profiles are,
# given the lat lon time from model receptors
# written by Dien Wu, 01/11/2017

# Input variables--
# ocopath, ocofile: oco2 path and file for searching satellite soundings;
# recp.lat, recp.lon: numeric numbers for receptor lat/lon
# diff.td: allowable thredshold for difference in lat/lon between given receptors and all satellite soundings

get.oco2info<-function(ocopath,ocofile,recp.lat,recp.lon,diff.td=1E-4){

# load libraries
library(ncdf4)

# grabbing OCO-2 info
ocodat<-nc_open(paste(ocopath,ocofile,sep=""))

# grabbing OCO-2 levels, lat, lon
# level 1 to 20, for space-to-surface, level 20 is the bottom level
# may need to reverse later
oco.level<-ncvar_get(ocodat,"levels")
oco.lat<-ncvar_get(ocodat,"latitude")
oco.lon<-ncvar_get(ocodat,"longitude")

# grabbing warn levels
warnlevel<-ncvar_get(ocodat,"warn_level")

# grabbing time for STILT receptors
id<-as.character(ncvar_get(ocodat,"sounding_id"))	# YYYY MM DD HH mm ss m (millisecond) f (footprint)

# locate the OCO2 data using lat, lon, when diff are both the smallest
diff.lat<-abs(oco.lat-recp.lat)
diff.lon<-abs(oco.lon-recp.lon)

# try to find the closest sounding lat/lon, given receptor lat/lon and allowable difference thredshold, "diff.td"
lat.index<-which(diff.lat < diff.td)
lon.index<-which(diff.lon < diff.td)
loc.index<-intersect(lat.index, lon.index)  # only if lat/lon indices are the same
if(length(loc.index)!=1){cat("get.oco2info():cannot find the receptor lat lon from this OCO2 file...");next}

# return the oco2 lat, lon, ak, pw, apriori, profiles
find.lat<-oco.lat[loc.index]
find.lon<-oco.lon[loc.index]
find.id<-id[loc.index]

## grab column co2, averaging kernel, pressure weight and a priori CO2 profiles, only for TARGET REGION
# dimensions--[levels, soundingID]
# averaging kernel (non-dim)
ak.norm<-ncvar_get(ocodat,"xco2_averaging_kernel")[,loc.index];ak.norm[ak.norm==-999999]<-NA
attributes(ak.norm)$names<-oco.level

# pressure weighting (non-dim)
pw<-ncvar_get(ocodat,"pressure_weight")[,loc.index];pw[pw==-999999]<-NA
attributes(pw)$names<-oco.level

# pressure in hPa
pres<-ncvar_get(ocodat,"pressure_levels")[,loc.index];pres[pres==-999999]<-NA
attributes(pres)$names<-oco.level

# CO2.apriori in ppm
apriori<-ncvar_get(ocodat,"co2_profile_apriori")[,loc.index];apriori[apriori==-999999]<-NA
attributes(apriori)$names<-oco.level

### combine all OCO-2 vertical profiles
all.profiles<-data.frame(ak.norm, pw, pres, apriori)

# ground height measured in OCO-2 in meter ASL
grdhgt<-ncvar_get(ocodat,"Sounding/altitude")[loc.index];grdhgt[grdhgt==-999999]<-NA  # [sounding ID]

# retrieved XCO2 and its uncertainty
xco2.oco<-ncvar_get(ocodat,"xco2")[loc.index];xco2.oco[xco2.oco==-999999]<-NA
xco2.uncert.oco<-ncvar_get(ocodat,"xco2_uncertainty")[loc.index];xco2.uncert.oco[xco2.uncert.oco==-999999]<-NA

# temp at 700mb
#t_700<-ncvar_get(ocodat,"Retrieval/T700")[loc.index];t_700[t_700==-999999]<-NA

# satellite footprint
footprint<-ncvar_get(ocodat,"Sounding/footprint")[loc.index];footprint[footprint==-999999]<-NA

# retrieved surface pressure
psfc<-ncvar_get(ocodat,"Retrieval/psurf")[loc.index];psfc[psfc==-999999]<-NA

# combine all variables
other.info<-data.frame(grdhgt, psfc, footprint, xco2.oco, xco2.uncert.oco, find.id, find.lat, find.lon)

# return both profiles and variables
return(list(all.profiles, other.info))
} # end of subroutine
