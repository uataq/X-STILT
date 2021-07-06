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
# use area() in raster package to calculate grid area, DW, 01/20/2019

tif2nc.odiacv3 <- function(site, timestr, vname, workdir, foot.extent, 
                           tiff.path, gzTF = T) {

  library(Hmisc); library(raster)

  cat('Start reading and subsetting emissions that match foot...\n')
  emiss.path <- file.path(workdir, 'in')
  dir.create(emiss.path, showWarnings = F)  # create the dir to store emiss

  YYYYMM <- substr(timestr, 1, 6)
  mod <- Hmisc::monthDays(as.Date(paste0(substr(timestr, 1, 4), '-',
                                         substr(timestr, 5, 6), '-', 
                                         substr(timestr, 7, 8))))

  # path and filename of the 1km GeoTiff ODIAC
  gzfile <- list.files(path = tiff.path, pattern = substr(timestr, 3, 6), 
                       recursive = T, full.names = T)

  if (length(gzfile) == 0) {
    cat('NO tiff format of ODIAC emission found...\n'); return()}
    
  cat(paste('tif2nc.odiacv3(): working on file', gzfile, '...\n'))
  if (grepl('.gz', gzfile)) {   # unzip gz file
    # always unzip to odiac directory, do not unzip to group_data
    cat('Unzipping up...\n'); system(paste0('gunzip ', gzfile))
    tiff.file <- substr(gzfile, 1, nchar(gzfile) - 3)
  } else { tiff.file <- gzfile } # end if gz
  gc()

  # after reading using readGDAL, the default dimension is (y,x)
  # 21600 rows and 43200 columns
  emiss <- raster(tiff.file) # convert to raster
  cat('Done reading tiff file as raster..\n')

  # subset spatial domain
  sel.emiss <- crop(emiss, foot.extent)

  # Method 2 -- compute area using area() function in raster package
  area.raster <- raster::area(sel.emiss) * 1E6    # convert km2 to m2

  # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
  sel.emiss <- sel.emiss * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
  sel.emiss <- sel.emiss / mod / 24 / 60 / 60	# convert per month to per second
  sel.emiss <- sel.emiss / area.raster		# convert per cell to per m2
  # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint

  # store as nc
  cat('Storing cropped gridded emissions in nc files...\n')
  outname <- paste0('odiac', vname, '_1kmx1km_', YYYYMM, '_', site, '.nc')
  outfile <- file.path(emiss.path, outname)
  crs(sel.emiss) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  writeRaster(sel.emiss, outfile, overwrite = TRUE, format = 'CDF',
              varname = 'emiss', varunit = 'micromole-CO2/m2/s',
              longname = 'ODIAC emissions', xname = 'lon', yname = 'lat')

  # finally zip up the tif file
  if (gzTF) {cat('Zipping up.\n'); system(paste0('gzip ', tiff.file))}

  # finally, return nc filename
  return(outfile)
}
