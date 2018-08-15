### subroutine to define background region with forward-time run from a box
# by DW, 11/02/2017

# try box receptors/sources,'dxyp' in trajecmulti(), DW, 11/08/2017
# add time windows for releasing particles, DW, 11/09/2017
# clear things up and add more sites, DW, 04/16/2018

# read variables from output of 'create_namelist_trajec.r'
# remove namelist, use all variables instead, DW, 07/27/2018
# clear code for X-STILTv2, DW, 07/29/2018

# *** default output directory will be automatically created and under 'workdir',
# in different Copies; make sure different nummodels are used for different
# overpass events, DW, 07/31/2018

run.forward.trajec <- function(site, timestr, overwrite = F, nummodel = 0,
  lon.lat, delt, dxyp, dzp, dtime, agl, numpar, nhrs, workdir, outpath = NULL,
  siguverr = NULL, TLuverr = NULL, zcoruverr = NULL, horcoruverr = NULL,
  met.format, met.path, met.num = 1, plotTF = F, oco2.path, oco2.ver,
  zoom = 8, td = 0.05, clean.side = c('north','south', 'both')[3]){

  # get lat/lon for city center
  clon <- signif(lon.lat$citylon, 6)
  clat <- signif(lon.lat$citylon, 6)
  lonstr <- ifelse(clon > 0, paste0(clon, 'E'), paste0(abs(clon), 'W'))
  latstr <- ifelse(clat > 0, paste0(clat, 'N'), paste0(abs(clat), 'S'))
  cat(paste('run.forward.trajec(): Working on trajec at', latstr, '\n'))

  # last release time
  recp.date <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  recp.yr   <- as.numeric(substr(timestr, 1, 4))
  recp.mon  <- as.numeric(substr(timestr, 5, 6))
  recp.day  <- as.numeric(substr(timestr, 7, 8))
  recp.hr   <- as.numeric(substr(timestr, 9, 10))

  # calculate each release times when continuously releasing particles
  rel.date <- recp.date + dtime * 60 * 60  # in second
  format.date <- as.data.frame(matrix(as.numeric(unlist(
    strsplit(format(rel.date, '%Y-%m-%d-%H-%M-%S'), '-'))),
    byrow = T, ncol = 6), stringsAsFactors = F)
  colnames(format.date) <- c('yr', 'mon', 'day', 'hr', 'min', 'sec')

  # create traj names for each release time, in form of
  # YYYYxMMxDDxHHxmmxNxWxmAGLxPxdeg that works with Trajecmulti()
  ident <- paste(format(rel.date, '%Yx%mx%dx%Hx%M'), latstr, lonstr,
                 formatC(agl, width = 5, flag = 0), paste0(numpar, 'P'),
                 dxyp, sep = 'x')

  # if not outpath is inputted, here is the default directory
  if (is.null(outpath)) outpath <- paste0(workdir, '/Copy', nummodel, '/')

  # if running trajec
  if (overwrite) {

    # delete previous directories and then create new one
    system(paste0('rm -rf ', outpath), ignore.stderr = T)
    dir.create(outpath, showWarnings = FALSE)

    # linking AER_NOAA_branch's hymodelc and other executables to outpath
    exes <- list.files(file.path(workdir, 'exe/AER_NOAA_branch'))
    file.symlink(file.path(workdir, 'exe/AER_NOAA_branch', exes), outpath)

    # if using multiple receptors, or box of receptors or sources, turn it on,
    # then call updated Trajecmulti() instead of Trajec()
    varstrajec <- c('time', 'index', 'lat', 'lon', 'agl', 'grdht', 'foot',
                    'sampt', 'dmass', 'zi', 'pres')  # trajec output variables

    # get met files, run_time in form of 'yyyy-mm-dd HH:MM:SS (UTC)'
    metfiles <- NULL
    for (r in 1:length(rel.date))
      metfiles <- c(metfiles, find_met_files(t_start = rel.date[r],
                                             met_file_format = met.format,
                                             n_hours = nhrs,
                                             met_loc = met.path))
    metfiles <- basename(unique(metfiles))

    #### try box receptors/sources, DW, 11/08/2017
    # the updated Trajecmulti() and fortran codes will randomly place receptors
    # according to dxyp
    cat('Generating forward trajec...\n')
    Trajecmulti(yr  = format.date$yr - 2000, mon = format.date$mon,
                day = format.date$day, hr = format.date$hr,
                mn  = format.date$min, outname = ident, numpar = numpar,
                lat = rep(clat, nrow(format.date)),
                lon = rep(clon, nrow(format.date)),
                dxyp = rep(dxyp, nrow(format.date)),
                dzp  = rep(dzp, nrow(format.date)),
                agl  = rep(agl, nrow(format.date)),
                nhrs = rep(nhrs, nrow(format.date)),
                nummodel = nummodel, metd = c('fnl','awrf'),
                outpath = outpath, overwrite = overwrite,
                metfile = metfiles, metlib = paste0(met.path, '/'),
                doublefiles = T, rundir = paste0(workdir, '/'),
                varsout = varstrajec, siguverr = siguverr, TLuverr = TLuverr,
                zcoruverr = zcoruverr, horcoruverr = horcoruverr,
                setup.list = list(DELT = delt, VEGHT = 0.5)) %>%
    invisible()
  } # end if overwrite

  # if plotting...
  if (plotTF) {
    bg.info <- ggplot.forward.trajec(ident = ident, trajpath = outpath,
                                     site = site, timestr = timestr,
                                     oco2.path = oco2.path,
                                     oco2.ver = oco2.ver, zoom = zoom,
                                     lon.lat = lon.lat, td = td,
                                     font.size = rel(1.2),
                                     clean.side = clean.side)
    return(bg.info)
  } else {
    return()
  }
}

# end of script
