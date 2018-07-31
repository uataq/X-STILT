### subroutine to define background region with forward-time run from a box
# by DW, 11/02/2017

# try box receptors/sources,'dxyp' in trajecmulti(), DW, 11/08/2017
# add time windows for releasing particles, DW, 11/09/2017
# clear things up and add more sites, DW, 04/16/2018
# clean codes, DW, 06/12/2018

# read variables from output of 'create_namelist_trajec.r'
run.forward.trajec <- function(namelist, plotTF = F){

  library(ggplot2); library(reshape)

  # if using multiple receptors, or box of receptors or sources, turn it on,
  # then call updated Trajecmulti() instead of Trajec()
  nummodel <- namelist$nummodel
  varstrajec <- c('time', 'index', 'lat', 'lon', 'agl', 'grdht', 'foot',
                  'sampt', 'dmass', 'zi', 'pres')  # trajec output variables

  # get lat/lon for city center and agl, time info
  lat.lon  <- namelist$lat.lon
  city.lat <- lat.lon[5]
  city.lon <- lat.lon[6]
  dxyp     <- namelist$dxyp

  agl      <- namelist$agl  # fixed level near ground
  nlevel   <- length(agl)

  # get first releassing time
  rel.yr  <- as.numeric(substr(namelist$timestr, 1, 4))
  rel.mon <- as.numeric(substr(namelist$timestr, 5, 6))
  rel.day <- as.numeric(substr(namelist$timestr, 7, 8))
  rel.hr  <- namelist$oco2.hr

  #### try box receptors/sources, DW, 11/08/2017
  # fortran code from Thomas Nehrkornm that automatically places receptors
  npar <- 10000   # increase particle number, since we have box of receptors
  mn   <- 0         # minutes
  lat.str <- formatC(signif(city.lat, 6), width = 5, format = 'f', digits = 4, flag = 0)
  lon.str <- formatC(signif(city.lon, 6), width = 5, format = 'f', digits = 4, flag = 0)

  #create traj name
  ident <- paste0(rel.yr, 'x', formatC(rel.mon, width = 2, flag = 0), 'x',
                  formatC(rel.day, width = 2, flag = 0), 'x',
                  formatC(rel.hr, width = 2, flag = 0), 'x', lat.str, 'Nx',
                  lon.str, 'Ex', formatC(agl, width = 5, flag = 0), 'x',
                  formatC(npar, width = 5, flag = 0), 'Px', dxyp, 'deg')

  #### if continously release particles, need update receptor info
  if (namelist$windowTF) {

    npar  <- namelist$npar  # number of particles per run

    # offset hours, e.g., release every half an hour
    dtime <- namelist$dtime
    dhr   <- floor(dtime)       # in hours
    dmin  <- dtime * 60 - dhr * 60       # in minutes

    # update hours or days
    update.hr <- rel.hr + dhr  # same time zone as rel.hr, UTC
    update.day <- rep(rel.day, length(update.hr))

    forw.index <- update.hr >= 24
    update.hr [forw.index] <- update.hr [forw.index] - 24
    update.day[forw.index] <- update.day[forw.index] + 1

    back.index <- update.hr < 0
    update.hr [back.index] <- update.hr [back.index] + 24
    update.day[back.index] <- update.day[back.index] - 1

    #### update all receptor info, vectors for all input
    rel.yr  <- rep(rel.yr, length(dtime))
    rel.mon <- rep(rel.mon, length(dtime))
    rel.day <- update.day
    rel.hr  <- update.hr

    # mn in trajectmulti() has unit in minutes, simply is the receptor minutes...
    rel.min <- dmin
    rel.lat <- rep(city.lat, length(dtime))
    rel.lon <- rep(city.lon, length(dtime))
    rel.agl <- rep(agl, length(dtime))

    #### trajec names
    ident <- paste0(rel.yr, 'x', formatC(rel.mon, width = 2, flag = 0), 'x',
                   formatC(rel.day, width = 2, flag = 0), 'x',
                   formatC(rel.hr, width = 2, flag = 0), 'x',
                   formatC(rel.min, width = 2, flag = 0), 'x',
                   lat.str, 'Nx', lon.str, 'Ex',
                   formatC(rel.agl, width = 5, flag = 0), 'x',
                   formatC(npar, width = 5, flag = 0), 'Px', dxyp, 'deg')

  }  # end if windowTF

  #### the updated Trajecmulti() will randomly place receptors according to dxyp
  # call trajecmulti() to generate trajec, it takes time...
  nhrs <- namelist$nhrs
  nhrs <- rep(nhrs,length(rel.yr))

  cat(paste('run.forward.trajex(): Working on forward traj @', city.lat, 'N\n'))
  if (namelist$overwrite) {
    info <- Trajecmulti(yr = rel.yr - 2000, mon = rel.mon, day = rel.day,
                        hr = rel.hr, mn = rel.min, dxyp = dxyp, dzp = 0,
                        lat = rel.lat, lon = rel.lon, agl = rel.agl,
                        outname = ident, nhrs = nhrs, numpar = npar,
                        nummodel = nummodel, metd = c('fnl','awrf'), rundir = rundir,
                        metfile = namelist$metfile, metlib = namelist$metpath,
                        doublefiles = T, overwrite = namelist$overwrite,
                        outpath = namelist$outpath, varsout = varstrajec,
                        siguverr = siguverr, TLuverr = TLuverr,
                        zcoruverr = zcoruverr, horcoruverr = horcoruverr,
                        setup.list = list(DELT = namelist$delt, VEGHT = 0.5))
  } # end if overwrite

  # if plotting...
  if(plotTF){
     zoom <- 8
     pp <- ggplot.forward.trajec(ident = ident, trajpath = namelist$outpath,
                                 site = namelist$site, timestr = namelist$timestr,
                                 ocopath = namelist$oco2.path, lat.lon = lat.lon,
                                 zoom = zoom)
    return(pp)
  } # end if plotTF

}
# end of script
