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

run.forward.trajec <- function(site, site.lon, site.lat, timestr, 
                               run_trajec = T, run_hor_err = F, run_ver_err = F, 
                               
                               agl, delt, dtime, dxyp, dzp = 0, numpar,
                               nhrs, xstilt_wd, output.path, 
                               
                               # met fields
                               met, met.format, met.path, met.num = 1, 

                               # horizontal and vertial trans error input
                               raob.path, raob.format = 'fsl', siguverr = NULL, 
                               err.path, overwrite = F, 
                               nhrs.zisf = 24, zisf         # zi scaling factor
                               ){
  
  setwd(xstilt_wd); source('r/dependencies.r', local = T) # source all functions
  cat(paste('\n\n## ----- Working on', site, 'on', timestr, '----- ##\n'))
  if (is.list(dtime)) dtime <- dtime[[1]]

  #------------------------------ STEP 3 --------------------------------- #
  ## get horizontal transport error component if run_hor_err = T
  # if overwrite = T, prepare RAOB data and compute model-data wind comparisons
  # *** if you do not want to run wind error analysis, set overwrite to FALSE 
  # and prescribe a horizontal wind error, e.g., siguverr = 3 (with unit of m/s)
  # for more info, please see get.uverr.r
  #lon.lat  <- get.lon.lat(site = site, dlon = 5, dlat = 5)  # near fields wind error
  lon.lat <- data.frame(minlon = site.lon - 5, maxlon = site.lon + 5, 
                        minlat = site.lat - 5, maxlat = site.lat + 5)
  hor.err <- get.uverr(run_hor_err, site, timestr, xstilt_wd, overwrite,
                       raob.path, raob.format, nhrs, met, met.path, 
                       met.format, met.files, lon.lat, agl = c(0, 100), 
                       nfTF = T, siguverr = siguverr, err.path = err.path)

  ## get vertical transport error component if run_ver_err = T
  pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)
  cat('run.forward.trajec(): Done with choosing met & inputting wind errors...\n')

  # reformat trajec info that fits Trajecmulti.r
  results <- convert.timestr2ident(timestr, dtime, site.lon, site.lat, agl, 
                                   numpar, dxyp)
  ident <- results$ident
  date.df <- results$date.df

  # create out_forward dir for generating forward trajec
  rundirname <- paste0('out_forward_', timestr, '_', site)
  output.path <- file.path(output.path, rundirname)
  
  # if running trajec
  if (run_trajec) {

    # delete previous directories and then create new one
    system(paste0('rm -rf ', output.path), ignore.stderr = T)
    dir.create(output.path, showWarnings = FALSE, recursive = T)

    # linking AER_NOAA_branch's hymodelc and other executables to outpath
    exes <- list.files(file.path(xstilt_wd, 'exe'))
    file.symlink(file.path(xstilt_wd, 'exe', exes), output.path)
    
    # if using multiple receptors, or box of receptors or sources, turn it on,
    # then call updated Trajecmulti() instead of Trajec()
    varstrajec <- c('time', 'index', 'lat', 'lon', 'agl', 'grdht', 'foot',
                    'sampt', 'dmass', 'zi', 'pres')  # trajec output variables

    metfiles <- find.all.metfiles(timestr, dtime, met.format, met.path, nhrs)

    #### try box receptors/sources, DW, 11/08/2017
    # the updated Trajecmulti() and fortran codes will randomly place receptors
    # according to dxyp
    cat('run.forward.trajec(): Generating forward trajec...\n')
    Trajecmulti(yr = date.df$yr - 2000, mon = date.df$mon, day = date.df$day, 
                hr = date.df$hr, mn  = date.df$min, outname = ident, 
                numpar = numpar, lat = date.df$lat, lon = date.df$lon,
                dxyp = rep(dxyp, nrow(date.df)), 
                dzp  = rep(dzp, nrow(date.df)),
                agl  = rep(agl, nrow(date.df)), 
                nhrs = rep(nhrs, nrow(date.df)),
                nummodel = timestr, 
                metd = c('fnl', 'awrf'), 
                outpath = output.path, 
                overwrite = run_trajec,
                metfile = metfiles, 
                metlib = paste0(met.path, '/'),
                doublefiles = T, 
                rundir = dirname(output.path), 
                rundirname = basename(output.path), 
                varsout = varstrajec, 
                siguverr = hor.err$siguverr, 
                TLuverr = hor.err$TLuverr, 
                zcoruverr = hor.err$zcoruverr, 
                horcoruverr = hor.err$horcoruverr,
                hymodelc.exe = './hymodelc.aer', # use the AER version of hymodelc
                setup.list = list(DELT = delt, VEGHT = 0.5)) %>%
    invisible()
  } else {
    cat('run.forward.trajec(): do nothing, check run_trajec...\n')
  } # end if run_trajec
} # end of script





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

  metfiles <- NULL
  for (r in 1:length(rel.date))
    metfiles <- c(metfiles, find_met_files(t_start = rel.date[r],
                                            met_file_format = met.format,
                                            n_hours = nhrs, met_loc = met.path))
  metfiles <- basename(unique(metfiles))
  if (length(metfiles) == 0) {cat('No meteo fields found...please check'); return()}

  return(metfiles)
} # end of find.all.metfiles