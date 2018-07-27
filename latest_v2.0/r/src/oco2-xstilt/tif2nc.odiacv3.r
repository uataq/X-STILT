# subroutine to readin tiff format of 1kmx1km ODIAC emissions
# and return rds format with selected emissions given a lat/lon domain
# DW, update with ODIACv2017, 11/01/2017

# add missing power plant in CARVA, PP10 and PP14 in ODIAC, DW, 11/28/2017
# addTF  whether to add the missing PP in CARVA
# tiff.path, tiff format of ODIAC on lin group
# modify based on Ben's code, use rds format and flip emiss lat/lon, DW, 06/04/2018
# also, use centered lat/lon of emissions instead of lower left, match STILTv2
# use raster package that saves a lot of time, DW, 06/18/2018
# use raster to read tif file, instead of rgdal
# no need to store as RDS, directly read from tif file, DW, 06/19/2018
# gzTF, whether to unzip tiff file in the end

tif2nc.odiacv3 <- function(site, timestr, vname, workdir, foot.extent,
  store.path = file.path(workdir, 'in', 'ODIAC'), tiff.path, gzTF = T){

  library(Hmisc); library(raster)

  YYYYMM <- substr(timestr, 1, 6)
  mod <- Hmisc::monthDays(as.Date(paste0(substr(timestr, 1, 4), '-',
    substr(timestr, 5, 6), '-', substr(timestr, 7, 8))))

  # path and filename of the 1km GeoTiff ODIAC
  gzfile <- list.files(path = tiff.path, pattern = substr(timestr, 3, 6))
  cat(paste('tif2rds.odiac(): working on file', gzfile, '...\n'))

  if (grepl('.gz', gzfile)) {   # unzip gz file
    # always unzip to odiac directory, do not unzip to group_data
    cat('Unzipping up...\n')
    system(paste0('gunzip ', file.path(tiff.path, gzfile)))
    tiff.file <- substr(gzfile, 1, nchar(gzfile) - 3)

  } else {
    tiff.file <- gzfile
  } # end if gz

  gc()

  # after reading using readGDAL, the default dimension is (y,x)
  # 21600 rows and 43200 columns
  tiff.file <- file.path(tiff.path, tiff.file)
  emiss <- raster(tiff.file) # convert to raster
  print(emiss)
  cat('Done reading tiff file as raster.\n')

  # subset spatial domain
  print(foot.extent)
  sel.emiss <- crop(emiss, foot.extent)
  cat('Done subsetting emissions according to footprint domain.\n')

  # return the area map as [LON, LAT]
  cat('tif2nc.odiac(): calculating the area for each grid...\n')
  res <- res(sel.emiss)[1]
  ext <- sel.emiss@extent
  area.co2 <- area(res = res, start.lon = ext@xmin,
    start.lat = ext@ymin, end.lon = ext@xmax - res, end.lat = ext@ymax - res,
  	third.dim = NULL) # need lower left lat/lon, return dim [lon, lat]

  # fix bug to reformat area --> prepare for raster, DW, 06/22/2018
  area.co2 <- t(area.co2) # convert to [lat, lon]
  area.lat <- rownames(area.co2)
  area.co2 <- area.co2[length(area.lat):1, ] # flip lat, decreasing trend

  # convert area.co2 to raster form
  area.raster <- raster(area.co2, xmn = ext@xmin, xmx = ext@xmax,
    ymn = ext@ymin, ymx = ext@ymax, crs = crs(sel.emiss))
  plot(area.raster)

  # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
  sel.emiss <- sel.emiss * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
  sel.emiss <- sel.emiss / mod / 24 / 60 / 60	# convert per month to per second
  sel.emiss <- sel.emiss / area.raster		# convert per cell to per m2
  # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
  print(sel.emiss)
  
  # store as nc
  cat('Storing contribution map as nc format...\n')
  outname <- paste0('odiac', vname, '_1kmx1km_', YYYYMM, '_', site, '.nc')
  outfile <- file.path(store.path, outname)
  crs(sel.emiss) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(sel.emiss, outfile, overwrite = TRUE, format = "CDF",
    varname = "emiss", varunit = "micromole-CO2/m2/s",
    longname = "ODIAC emissions", xname = "lon", yname = "lat")

  if (gzTF) {
    # finally zip up the tif file
    cat('Zipping up.\n')
    system(paste0('gzip ', file.path(tiff.path, tiff.file)))
  }

  # finally, return nc filename
  return(outfile)
}
