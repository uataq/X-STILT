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

tif2nc.odiacv2 <- function(timestr, foot.path, foot.file, workdir,
  store.path = file.path(workdir, 'in', 'ODIAC'), tiff.path, vname, site){

  library(Hmisc); library(raster)

  YYYYMM <- substr(timestr, 1, 6)
  mod <- Hmisc::monthDays(as.Date(paste0(substr(timestr, 1, 4), '-',
    substr(timestr, 5, 6), '-', substr(timestr, 7, 8))))

  # path and filename of the 1km GeoTiff ODIAC
  gzfile <- list.files(path = tiff.path, pattern = substr(timestr, 3, 6))

  # check if RData file exist, no need to unzip and convert from tiff
  rds.file <- paste0('ODIACv', vname, '_', YYYYMM, '.rds')
  gzTF <- F # whether to unzip tiff file, initialize with F

  if (!file.exists(file.path(store.path, rds.file))) {
    cat(paste('tif2rds.odiac(): working on file', gzfile, '...\n'))

    if (grepl('.gz', gzfile)) {   # unzip gz file

      # always unzip to odiac directory, do not unzip to group_data
      cat('Unzipping up...\n')
      system(paste0('gunzip ', file.path(tiff.path, gzfile)))
      tiff.file <- substr(gzfile, 1, nchar(gzfile) - 3)
      gzTF <- T  # remember to unzip tiff file in the end

    } else {
      tiff.file <- gzfile
    } # end if gz

    # after reading using readGDAL, the default dimension is (y,x)
    # 21600 rows and 43200 columns
    cat('Reading tiff file.\n')
    emiss <- raster(file.path(tiff.path, tiff.file)) # convert to raster

    # assign as RData file
    cat('Saving emissions as RDS file, no need to unzip gz file next time.\n')
    saveRDS(object = emiss, file = file.path(store.path, rds.file))

  } else {
    # if RDS file found, read it
    cat('RDS format found...reading in emissions...\n')
    emiss <- readRDS(file.path(store.path, rds.file)) # in raster form

  } # end if !file.exists()

  # subset spatial domain
  foot.extent <- extent(raster(file.path(foot.path, foot.file[1])))
  sel.emiss <- crop(emiss, foot.extent)

  # return the area map as [LON, LAT]
  cat('tif2nc.odiac(): calculating the area...\n')
  res <- res(sel.emiss)[1]
  ext <- sel.emiss@extent
  area.co2 <- area(res = res, start.lon = ext@xmin,
    start.lat = ext@ymin, end.lon = ext@xmax - res, end.lat = ext@ymax - res,
  	third.dim = NULL) # need lower left lat/lon

  # convert area.co2 to raster form
  area.raster <- raster(area.co2, xmn = ext@xmin, xmx = ext@xmax,
    ymn = ext@ymin, ymx = ext@ymax, crs = crs(sel.emiss))

  # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
  sel.emiss <- sel.emiss * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
  sel.emiss <- sel.emiss / mod / 24 / 60 / 60	# convert per month to per second
  sel.emiss <- sel.emiss / area.raster		# convert per cell to per m2
  # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint

  # store as nc
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
  gc()

  # finally, return nc filename
  return(outfile)
}






# for converting tif plot
if (F) {

  timestr <- seq(as.Date('2014/01/01'), as.Date('2018/05/01'), by = 'month')
  timestr <- format(timestr, '%Y%m%d')
  YYYYMM <- substr(timestr, 1, 6)
  homedir <- '/uufs/chpc.utah.edu/common/home'
  workdir <- file.path(homedir, 'lin-group1/wde/github/stilt')
  source('r/dependencies.r') # source all functions
  tiff.path <- file.path(homedir, 'lin-group1/group_data/ODIAC/ODIAC2017',
    substr(timestr, 1,4))  # tif file from ODIAC website
  store.path <- file.path(workdir, 'in')

  for (t in 1:length(tiff.path)) {

    # path and filename of the 1km GeoTiff ODIAC
    gzfile <- list.files(path = tiff.path[t], pattern = substr(timestr[t], 3, 6))

    # check if RData file exist, no need to unzip and convert from tiff
    rds.file <- paste0('ODIACv2017_', YYYYMM[t], '.rds')
    cat(paste('tif2rds.odiac(): working on file', gzfile, '...\n'))

    if (grepl('.gz', gzfile)) {   # unzip gz file
      # always unzip to odiac directory, do not unzip to group_data
      cat('Unzipping up...\n')
      system(paste0('gunzip ', file.path(tiff.path[t], gzfile)))
      tiff.file <- substr(gzfile, 1, nchar(gzfile) - 3)

    } else {
      tiff.file <- gzfile
    } # end if gz

    # after reading using readGDAL, the default dimension is (y,x)
    # 21600 rows and 43200 columns
    cat('Reading tiff as raster.\n')
    emiss <- raster(file.path(tiff.path[t], tiff.file)) # convert to raster

    # assign as RData file
    cat('Saving emissions as RDS file, no need to unzip gz file next time.\n')
    save(emiss, file = file.path(store.path, 'test.RData'))
    save(emiss, file = file.path(store.path, rds.file))

    # finally zip up the tif file
    cat('Zipping up.\n')
    system(paste0('gzip ', file.path(tiff.path[t], tiff.file)))
    gc()
  }  # end for t

} # end if F
