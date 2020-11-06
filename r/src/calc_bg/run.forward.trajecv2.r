### subroutine to define background region with forward-time run from a box
# see 'ggplot.forward.trajec.r' for more info on input variables
# by DW, 11/02/2017

# try box receptors/sources,'dxyp' in trajecmulti(), DW, 11/08/2017
# add time windows for releasing particles, DW, 11/09/2017
# clear things up and add more sites, DW, 04/16/2018

# *** default output directory will be automatically created and under 'xstilt_wd',
# in different Copies; make sure different nummodels are used for different
# overpass events, DW, 07/31/2018
# add customized data filtering, DW, 08/20/2018
# remove data filtering, always use QF = 0 for background estimates, DW, 10/31/2018
# rename output folder names, DW, 01/25/2019 


# ---------------------- function inside run.forward.trajec ------------------ #
convert.timestr2ident <- function(timestr, dtime, site.lon, site.lat, 
                                  agl, numpar, dxyp) {
  
  # get lat/lon for city center as well as time string
  clon <- signif(site.lon, 6)
  clat <- signif(site.lat, 6)
  lonstr <- ifelse(clon > 0, paste0(clon, 'E'), paste0(abs(clon), 'W'))
  latstr <- ifelse(clat > 0, paste0(clat, 'N'), paste0(abs(clat), 'S'))
  recp.date <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  rel.date <- recp.date + dtime * 60 * 60  # in second

  # calculate each release times when continuously releasing particles
  format.date <- strsplit.to.df(format(rel.date, '%Y-%m-%d-%H-%M-%S'), sep = '-') 
  format.date <- format.date %>% mutate_all(funs(as.numeric), colnames(format.date))
  colnames(format.date) <- c('yr', 'mon', 'day', 'hr', 'min', 'sec')
  format.date$lat <- clat 
  format.date$lon <- clon 

  # create traj names for each release time, in form of
  # YYYYxMMxDDxHHxmmxNxWxmAGLxPxdeg that works with Trajecmulti()
  ident <- paste0(paste(format(rel.date, '%Yx%mx%dx%Hx%M'), latstr, lonstr,
                        formatC(agl, width = 5, flag = 0), paste0(numpar, 'P'),
                        dxyp, sep = 'x'), '.rds')
  return(list(ident = ident, date.df = format.date))
} # end of convert.timestr2ident



# ---------------------- function inside run.forward.trajec ------------------ #
find.all.metfiles <- function(timestr, dtime, met.format, met.path, nhrs) {

  # get met files, run_time in form of 'yyyy-mm-dd HH:MM:SS (UTC)'
  recp.date <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  rel.date <- recp.date + dtime * 60 * 60  # in second

  met.files <- NULL
  for (r in 1:length(rel.date))
    met.files <- c(met.files, find_met_files(t_start = rel.date[r],
                                            met_file_format = met.format,
                                            n_hours = nhrs, met_path = met.path))
  met.files <- basename(unique(met.files))
  if (length(met.files) == 0) {
    cat('find.all.metfiles(): No meteo fields found...please check'); return()}

  return(met.files)
} # end of find.all.met.files



# ---------------------- run.forward.trajec ()  ----------------------------- #
run.forward.trajecv2 <- function(site, site.lon, site.lat, timestr, run_trajec, 
                                 run_hor_err = F, run_ver_err = F, xstilt_wd, 
                                 traj.path, box.len, dtime.from, dtime.to, 
                                 dtime.sep, nhrs, delt, agl, numpar, 
                               
                                 # met fields
                                 met, met.res, met.format, met.path, met.num = 1, 

                                 # horizontal and vertial trans error input
                                 raob.path, raob.format = 'fsl', siguverr = NULL, 
                                 err.path, overwrite = F, 
                                 nhrs.zisf = 24, zisf         # zi scaling factor
                                 ){
    
  setwd(xstilt_wd); source('r/dependencies.r', local = T) # source all functions
  cat(paste('\n\n## ----- Working on', site, 'on', timestr, '----- ##\n'))

  #------------------------------ STEP 3 --------------------------------- #
  ## get horizontal transport error component if run_hor_err = T
  # if overwrite = T, prepare RAOB data and compute model-data wind comparisons
  # *** if you do not want to run wind error analysis, set overwrite to FALSE 
  # and prescribe a horizontal wind error, e.g., siguverr = 3 (with unit of m/s)
  # for more info, please see get.uverr.r
  #lon.lat  <- get.lon.lat(site = site, dlon = 5, dlat = 5)  # near fields wind error
  lon.lat <- data.frame(citylon = site.lon, citylat = site.lat, 
                        minlon = site.lon - 5, maxlon = site.lon + 5, 
                        minlat = site.lat - 5, maxlat = site.lat + 5)

  hor.err <- get.uverr(run_hor_err, site, timestr, xstilt_wd, overwrite, 
                       raob.path, raob.format, nhrs, met, met.path, met.format, 
                       lon.lat, agl = c(0, 100), err.path, nfTF = T, siguverr)

  ## get vertical transport error component if run_ver_err = T
  pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)
  cat('run.forward.trajec(): Done with choosing met & inputting wind errors...\n')

  # determine dxpy using the desired box size and grid length of met fields, 
  # since final box length = 2 * dxyp * met.res, DW, 06/30/2020
  dxyp <- box.len / met.res / 2   

  # reformat trajec info that fits Trajecmulti.r
  dtime <- seq(dtime.from, dtime.to, dtime.sep)     # vector in hours
  results <- convert.timestr2ident(timestr, dtime, site.lon, site.lat, agl, numpar, dxyp)
  ident   <- results$ident
  date.df <- results$date.df

  # create out_forward dir for generating forward trajec
  traj.dir <- file.path(traj.path, paste0('out_forward_', timestr, '_', site, '_debug'))

  # if running trajec
  if (run_trajec) {

    # delete previous directories and then create new one
    system(paste0('rm -rf ', traj.dir), ignore.stderr = T)
    dir.create(traj.dir, showWarnings = FALSE, recursive = T)

    # linking AER_NOAA_branch's hymodelc and other executables to outpath
    exes <- list.files(file.path(xstilt_wd, 'exe'))
    file.symlink(file.path(xstilt_wd, 'exe', exes), traj.dir)
    
    # if using multiple receptors, or box of receptors or sources, turn it on,
    # then call updated Trajecmulti() instead of Trajec()
    varstrajec <- c('time', 'index', 'lat', 'lon', 'agl', 'grdht', 'foot',
                    'sampt', 'dmass', 'zi', 'pres')  # trajec output variables

    #### try box receptors/sources, DW, 11/08/2017
    # the updated Trajecmulti() and fortran codes will randomly place receptors
    # according to dxyp
    cat('run.forward.trajec(): Generating forward trajec...\n')
    met.files <- find.all.metfiles(timestr, dtime, met.format, met.path, nhrs)
    
    # Aggregate STILT/HYSPLIT namelist
    namelist <- list(capemin = -1, cmass = 0, conage = 48, cpack = 1, delt = 1,
                     dxf = 1, dyf = 1, dzf = 0.01, efile = '', emisshrs = 0.01,
                     frhmax = 3, frhs = 1, frme = 0.1, frmr = 0, frts = 0.1,
                     frvs = 0.1, hnf_plume = T, hscale = 10800, ichem = 8,
                     idsp = 2, initd = 0, k10m = 1, kagl = 1, kbls = 1, kblt = 5,
                     kdef = 0, khinp = 0, khmax = 9999, kmix0 = 250, kmixd = 3,
                     kmsl = 0, kpuff = 0, krand = 4, krnd = 6, kspl = 1, kwet = 1,
                     kzmix = 0, maxdim = 1, maxpar = numpar, mgmin = 10,
                     mhrs = 9999, ncycl = 0, ndump = 0, ninit = 1, nstr = 0,
                     nturb = 0, nver = 0, numpar = 1000, outdt = 0, p10f = 1,
                     pinbc = '', pinpf = '', poutf = '', qcycle = 0, rhb = 80,
                     rht = 60, splitf = 1, tkerd = 0.18, tkern = 0.18,
                     tlfrac = 0.1, tout = 0, tratio = 0.75, tvmix = 1,
                     varsiwant = c('time', 'indx', 'long', 'lati', 'zagl', 
                                   'foot', 'mlht', 'dens','samt', 'sigw', 'tlgr'),
                     veght = 0.5, vscale = 200, vscaleu = 200, vscales = -1,
                     wbbh = 0, wbwf = 0, wbwr = 0, wvert = FALSE, zicontroltf = F,
                     ziscale = 0) 

    calc_trajectory(namelist, rundir, emisshrs, hnf_plume = T, met_files, n_hours, 
                    output, rm_dat, timeout, w_option, z_top)
                    
    Trajecmulti(yr = date.df$yr - 2000, mon = date.df$mon, day = date.df$day, 
                hr = date.df$hr, mn  = date.df$min, outname = ident, 
                numpar = numpar, lat = date.df$lat, lon = date.df$lon,
                dxyp = rep(dxyp, nrow(date.df)), 
                dzp  = rep(0, nrow(date.df)),
                agl  = rep(agl, nrow(date.df)), 
                nhrs = rep(nhrs, nrow(date.df)),
                nummodel = timestr, 
                metd = c('fnl', 'awrf'), 
                outpath = traj.dir, 
                overwrite = run_trajec,
                metfile = met.files, 
                metlib = paste0(met.path, '/'),
                doublefiles = T, 
                rundir = dirname(traj.dir), 
                rundirname = basename(traj.dir), 
                varsout = varstrajec, 
                siguverr = hor.err$siguverr, 
                TLuverr = hor.err$TLuverr, 
                zcoruverr = hor.err$zcoruverr, 
                horcoruverr = hor.err$horcoruverr,
                hymodelc.exe = './hymodelc.aer', # use the AER version of hymodelc
                setup.list = list(DELT = delt, VEGHT = 0.5)) %>% invisible()
  
  } else {
    cat('run.forward.trajec(): do nothing, check run_trajec...\n')
  } # end if run_trajec
} # end of script


