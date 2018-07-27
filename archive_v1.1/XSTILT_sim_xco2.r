# create a namelist for running X-STILT by DW, 04/18/2018
# require trajectories and store somewhere, see 'XSTILT_run_trajec.r'

# clean codes, DW, 06/12/2018

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT')
setwd(workdir) # move to working directory

source(file.path(workdir, 'src/sourceall.r'))  # source all functions
library(ncdf4); library(ggplot2)

#------------------------------ STEP 1 --------------------------------------- #
#### CHOOSE CITIES, TRACKS AND OCO-2 LITE FILE VERSION ***
site  <- 'Riyadh'
met   <- c('1km', 'gdas1', 'gdas0p5')[3]    # customized WRF, 1 or 0.5deg GDAS
tt    <- 2                                # which track to model
oco2.version <- c('b7rb', 'b8r')[1]        # OCO-2 version

# input all track numbers to be modeled, can be YYYYMMDD OR YYYYMMDDHH ***
riyadh.timestr<- c('20141227', '20141229', '20150128', '20150817', '20151112',
  '20151216', '20160115', '20160216', '20160725', '20161031')

# final timestr, YYYYMMDD and file strings for trajec
track.timestr <- riyadh.timestr[tt]

filestr <- paste0(substr(track.timestr,1,4), 'x', substr(track.timestr,5,6), 'x',
  substr(track.timestr,7,8))

#### SET spatial domains for generating footprints or grabbing ODIAC emissions,
# these variables will determine the numpix.x, numpix.y, lon.ll, lat.ll;
# e.g., for Riyadh, 	# 0-50N, 0-60E,  ***

# !!!! these lat/lon ranges differ from those in 'XSTILT_run_trajec.r' !!!!
#if (site == 'Riyadh') {minlon <- 0; maxlon <- 60; minlat <- 0; maxlat <- 50}
if (site == 'Riyadh') {minlon <- 30; maxlon <- 50; minlat <- 15; maxlat <- 35}

#### CHOOSE THE CO2 SOURCES/SINKS TO BE MODELED ***
# total 5 components, no wildfire component for now, but can be easily added
ss  <- c(1, 2, 3, 4, 5)
sel <- c('anthro', 'bio', 'ocean', 'edp', 'apriori')[ss]
cat('# ----------- STEP 1: city, tracks selected ...\n')


#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
ppTF <- F       # whether use ODIAC with addtional PP10+PP14
selTF <- F      # true for only use 1-day trajs.
dmassTF <- T    # dmassTF when generating trajecfoot().
storeTF <- T    # whether store model results;
                # e.g., footprint, co2 contribution maps, weighted trajec.
nummodel <- 0	  # copy for running trajwind()
hourlyTF <- F   # true for using hourly ODIAC emissions;
                # false for using monthly mean emissions.
columnTF  <- T  # whether a column receptor or fixed receptor
forwardTF <- F  # forward or backward traj
uncertTF  <- F  # sim using trajec with transport error components
biascorrTF<- F  # whether use bias-corrected trajec

ak.weight <- T  # whether weighting footprint by averaging kernels
pw.weight <- T  # whether weighting footprint by pressure weighting functions
if (columnTF) pw.weight <- T  # if release from a column, always XCO2

# get string, based on desired simulations
# e.g., non-weighting, only PW, or both AK PW weighting
if (ak.weight * pw.weight) string <- 'akpw'
if (ak.weight == F & pw.weight == T) string <- 'pw'
if (ak.weight * pw.weight * columnTF == F) string <- 'orig'
cat('# ----------- STEP 2: flags turned on/off ...\n')


#------------------------------ STEP 3 --------------------------------------- #
#### INPUT/CHANGE PATHS and GET TRAJEC INFO ***
# intpath: path that stores input prior data, including ODIAC, CT, and OCO-2
inpath  <- file.path(homedir, 'lin-group5/wde/input_data')
oco2.path <- file.path(inpath, paste0('OCO-2/L2/OCO2_lite_', oco2.version))

# outpath: main output directory and model results (txtfile)
# trajpath: path that stores original non-weighted trajectories
outpath  <- file.path(workdir, 'output')
trajpath <- file.path(outpath, 'Traj')

orig.trajpath <- file.path(trajpath, 'orig_traj')
#orig.trajpath <- file.path(trajpath, 'orig_traj', site)
#orig.trajdir  <- dir(orig.trajpath, track.timestr)
#orig.trajpath <- file.path(orig.trajpath, orig.trajdir)
#orig.trajpath <- '/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/Riyadh/gdas0p5/multi_agl/orig_traj/2014122910'

### ncdfpath: path that stores model intermediate (e.g., weighted traj,
# footprint, contribution map), call find.create.dir() to create one if not found
ncdfpath <- find.create.dir(path = outpath, dir = 'NetCDF', workdir = workdir)

### paths that store 1) weighted traj (wgt.trajpath), 2) footprint (foot.path),
# 3) XCO2 contribution maps for each source/sink (xco2.path)
# *** if no directories, create them, NO NEED TO CHANGE THESE -->
# 1) path for storing weighted traj (wgt.trajpath)
# if no weighting, return original non-weighting trajpath
wgt.dir <- paste0(string, '_traj')
wgt.trajpath <- find.create.dir(path = trajpath, dir = wgt.dir, workdir = workdir)

# 2) path for storing time-integrated (weighted column) footprint
foot.dir <- paste0(string, '_intfoot')
foot.path  <- find.create.dir(path = ncdfpath, dir = foot.dir, workdir = workdir)

# 3) path for storing foot * emission results--> 'xco2.path'
xco2.dir <- paste0(string, '_xco2')
xco2.path  <- find.create.dir(path = ncdfpath, dir = xco2.dir, workdir = workdir)
cat('# ----------- STEP 3: input required paths ...\n')


#------------------------------ STEP 4 --------------------------------------- #
#### GET ORIGINAL TRAJEC (NO AK WEIGHTING YET) PATH and FILES ***
# get trajec files, names and info, using ident.to.info()
orig.outname <- list.files(pattern = filestr, path = orig.trajpath)

# get trajec info
info <- ident.to.info(ident = orig.outname)
recp.info <- info[[1]]	# for receptor info
agl.info  <- info[[2]]	  # for release level info
timestr   <- unique(recp.info$timestr)
recp.info$recp.ident <- orig.outname

#### call get.grdhgt() or read from existing txt file
# to interpolate modeled ground heights from metfile,
# where to run trajwind() for model grdhgt, where Copy lies
rundir <- file.path(homedir, 'lin-group5/wde/STILT_modeling/STILT_Exe')

#### get ground height info:

# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG ***
metpath <- file.path(homedir, 'u0947337', met)
met.format <- '%Y%m%d_gdas0p5'
metfile <- find.metfile(timestr = timestr, nhrs = -1, metpath = metpath,
                        met.format = met.format)

overwrite <- F   # true for regenerating ground heights
filenm <- paste0('trajwind_recp_info_', site, '_', timestr, '.txt')

# if found
if (length(list.files(path = outpath, pattern = filenm)) != 0 & overwrite == F) {

  filenm <- file.path(outpath, filenm)
  trajwind.info <- read.table(file = filenm, header = T, sep = ',', stringsAsFactors = F)

  # then merge with recp.info
  sel.colnm <- colnames(recp.info)[grep('recp.', colnames(recp.info))]
  recp.info <- merge(recp.info, trajwind.info, by = sel.colnm)

}else{  # if there's none, call get.grdhgt() to get info

  cat('Interpolating Ground Heights & winds @ receptors from met fields...\n')
  cat('It takes few minutes...\n\n')
  recp.info <- get.grdhgt(
    recp.info = recp.info, nummodel = nummodel, agl = 10, nhrs = -1,
    rundir = rundir, outpath = outpath, site = site, metpath = metpath,
    metfile = metfile)

  # store ground height and wind info etc. @ receptor for future use
  trajwind.info <- recp.info[, grep('recp.', colnames(recp.info))]
  filenm <- file.path(outpath, filenm)
  write.table(trajwind.info, file = filenm, sep = ',', row.names = F, quote = F)
}  # end if get grdhgt
cat('# ----------- STEP 4: get all receptor info ...\n')


#------------------------------ STEP 5 --------------------------------------- #
#### choose ODIAC version
odiac.ver   <- 3		# odiac version, 2 for v2016; 3 for v2017
odiac.vname <- c('2015a', '2016', '2017')[odiac.ver]
odiac.path  <- file.path(inpath, 'anthro_inventories/ODIAC')
odiac.domain<- paste0(minlat, '-', maxlat, 'N_', minlon, '-', maxlon, 'E')

## read ODIAC emissions to save time
if (hourlyTF) {
  odiac.path <- file.path(odiac.path, odiac.vname, 'hourly', site)
  odiac.file <- list.files(path = odiac.path, pattern = odiac.domain)
  odiac.file <- file.path(odiac.path, odiac.file)
} else {
  odiac.path <- file.path(odiac.path, odiac.vname)
  pattern <- paste0(substr(timestr,1,6), '_', odiac.domain, '.nc')
  if (ppTF) pattern <- paste0(substr(timestr,1,6), '_', odiac.domain, '_PP.nc')
  odiac.file <- list.files(path = odiac.path, pattern = pattern)
  odiac.file <- file.path(odiac.path, odiac.file)
} # end if hourly

# if cannot find the correct format of nc file for emissions given selected area
# and return ODIAC file name with path in front
if (length(odiac.file) == 0) {
  cat('NO ODIAC file found, creating one from tiff ODIAC...\n')
  odiac.file <- convert.tid2nc.odiac(
    track.timestr, odiac.vname, odiac.path, ppTF, minlat, maxlat, minlon, maxlon)
} # end if

## read in emissions, grab emission grid, [LAT, LON, (HR.back)]
# lat, lon should have already been converted to LOWER LEFT !!!
cat('take a while to readin ODIAC CO2 emissions...\n')
odiac.dat <- nc_open(odiac.file)
odiac.co2 <- ncvar_get(odiac.dat, 'odiac_co2_emiss') # [lat,lon] or [lat,lon,hr]
odiac.lat <- ncvar_get(odiac.dat, 'lat')
odiac.lon <- ncvar_get(odiac.dat, 'lon')

if(hourlyTF){
  odiac.hr <- ncvar_get(odiac.dat, 'hr')  # hours back are negative
  dimnames(odiac.co2) <- list(odiac.lat, odiac.lon, odiac.hr)
}else{
  dimnames(odiac.co2) <- list(odiac.lat, odiac.lon)
}

#### CHOOSE CT-NRT version
# no CT-NRTv2017 files before 2015/12/12
ct.version <- c('v2016-1', 'v2017')[2]
if (substr(timestr, 1, 8) < '20151212') ct.version <- 'v2016-1'
ctflux.path <- file.path(inpath, 'CT-NRT', ct.version, 'fluxes', 'optimized')
ctmole.path <- file.path(inpath, 'CT-NRT', ct.version, 'molefractions', 'co2_total')

## create txt files for results and
vv <- 'v5'  # diff version, avoid overwriting
txtpath <- outpath
txtname <- paste0('STILT_OCO2_XCO2_', site, '_', timestr, '_odiac', odiac.vname)
if (hourlyTF) {
  txtname <- paste0(txtname, '_hrly_',met, vv, '.txt')
} else {
  txtname <- paste0(txtname, '_', met, vv, '.txt')
}
cat('# ----------- STEP 5: get prior emissions/fluxes ...\n')


#------------------------------ STEP 6 --------------------------------------- #
#### INPUT X-STILT PARAMETERS FOR GENERATING FOOTPRINT
max.hour <- -72       # footprint hour

# horizontal resolution of ODIAC
lon.res.anthro <- 1/120
lat.res.anthro <- lon.res.anthro

# horizontal resolution of CT-NRT
lon.res.ct <- 1
lat.res.ct <- lon.res.ct
cat('# ----------- STEP 6: X-STILT footprint setup ...\n')


#------------------------------ STEP 7 --------------------------------------- #
#### create a namelist (data.frame) including all variables
# and  to subroutines
namelist <- list(
  met = met, site = site, timestr = timestr, filestr = filestr, sel = sel,
  minlat = minlat, maxlat = maxlat, minlon = minlon, maxlon = maxlon,
  oco2.path = oco2.path, ak.weight = ak.weight, pw.weight = pw.weight,
  orig.trajpath = orig.trajpath, wgt.trajpath = wgt.trajpath,
	foot.path = foot.path, xco2.path = xco2.path, ct.version = ct.version,
  ctflux.path = ctflux.path, ctmole.path = ctmole.path, selTF = selTF,
  dmassTF = dmassTF, columnTF = columnTF, uncertTF = uncertTF,
  forwardTF = forwardTF, biascorrTF = biascorrTF, hourlyTF = hourlyTF,
  storeTF = storeTF, txtpath = txtpath, txtname = txtname,
  max.hour = max.hour,lat.res.ct = lat.res.ct,lon.res.ct = lon.res.ct,
  lat.res.anthro = lat.res.anthro, lon.res.anthro = lon.res.anthro)

# store namelist and trajlist as RData file
filename <- paste0('simlist_', site, '_', track.timestr, '.txt')
filename <- file.path(outpath, filename)
write.table(x = t(namelist), file = filename, sep = '\n',
            col.names = F, row.names = F, quote = T)
cat('# ----------- STEP 7: Namelist created/stored ...\n')


#------------------------------ STEP 8 --------------------------------------- #
# start X-STILT runs, call oco2.get.xco2()
cat('\n\n# ----------- STEP 8: Simulation started ...\n')
sim.xco2(namelist = namelist, recp.info = recp.info, odiac.co2 = odiac.co2)


# end of script
