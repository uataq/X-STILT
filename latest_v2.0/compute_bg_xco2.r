#' Compute background XCO2, given different methods:
#' M1. Trajec-endpoint (using CarbonTracker)
#' M2H. Regional daily median (based on Hakkareinen et al., 2016)
#' M2S. Localized normal statistics (based on Silva and Arellano, 2017)
#' M3. X-STILT overpass-specific (based on Wu et al., GMDD)
#' originated from 'create_namelist_oco2-xsilt.r'
#' @author Dien Wu, 04/18/2018

#' @updates:
#' add customized data filtering, DW, 08/20/2018


#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT/latest_v2.0')
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

#------------------------------ STEP 1 --------------------------------------- #
# insert target city
site <- 'LV'

# OCO-2 version, path
oco2.ver <- c('b7rb', 'b8r')[2]  # OCO-2 version
input.path <- file.path(homedir, 'lin-group5/wde/input_data')
output.path <-file.path(homedir, 'lin-group5/wde/github/result')

oco2.path <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_SIF_', oco2.ver))
txtpath <- file.path(output.path, 'oco2_overpass')

# date range for searching OCO-2 tracks, min, max YYYYMMDD
date.range <- c('20140101', '20181231')

# input dlat, dlon for spatial domain around city center
lon.lat <- get.lon.lat(site, dlon = 1, dlat = 3)

# 'thred.count' for at least how many soundings needed per 1deg lat range
# -> calculate a total thred on total # of soundings given 'lon.lat'
thred.count.per.deg <- 100  # number of soundings per degree
thred.count.per.deg.urban <- 50

# whether to re-search OCO-2 overpasses and output in txtfile
# if FALSE, read timestr from existing txt file;
# always TRUE, if doing first simulation for a new site
searchTF <- F

# whether search for overpasses over urban region,
# defined as city.lat +/- dlat, city.lon +/- dlon
urbanTF <- T; dlon.urban <- 0.5; dlat.urban <- 0.5

# call get.site.info() to get lon.lat and OCO2 overpasses info
# PLEASE add lat lon info in 'get.site.track'
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF,
  date.range, thred.count.per.deg, lon.lat, urbanTF, dlon.urban, dlat.urban,
  thred.count.per.deg.urban, txtpath)

# one can further subset 'oco2.track' based on sounding # over near city center
# one can further subset 'oco2.track' based on data quality
# see columns 'qf.count' or 'wl.count' in 'oco2.track'
# e.g., choose overpasses that have 100 soundings with QF == 0, & get reordered
if (urbanTF) oco2.track <- oco2.track %>% filter(tot.urban.count > 200)
oco2.track <- oco2.track %>% filter(qf.urban.count > 50)

# remove summtertime tracks
oco2.track <- oco2.track %>%
  filter(substr(timestr, 5, 6) < '05' | substr(timestr, 5, 6) > '08')

# finally narrow down and get timestr
#all.timestr <- oco2.track$timestr[c(2, 3, 5, 8, 9, 10, 13)]
#all.timestr <- oco2.track$timestr[c(2, 3, 5, 6, 7)]  # v7, Riyadh
all.timestr <- oco2.track$timestr[c(7, 8, 9, 10, 14)]  # v8, LV
print(all.timestr)

# once you have all timestr, you can choose whether to plot them on maps
# this helps you choose which overpass to simulate first, see 'tt' below
plotTF <- F
if (plotTF) {
  for(t in 1:length(all.timestr)){
  ggmap.obs.xco2(site, timestr = all.timestr[t], oco2.path, lon.lat, workdir,
    plotdir = workdir)
  ggmap.obs.sif(site, timestr = all.timestr[t], sif.path, lon.lat, workdir,
    plotdir = workdir)
  }
}

## choose which background method:
# M1. Trajec-endpoint (using CarbonTracker)
# M2H. Regional daily median (based on Hakkareinen et al., 2016)
# M2S. Localized normal statistics (based on Silva and Arellano, 2017)
# M3. X-STILT overpass-specific (based on Wu et al., GMDD)
method <- c('M1', 'M2H', 'M2S', 'M3')[4]

# ---------------------------- Trajec endpoint  --------------------------- #
# need to convolve footprint with diff fluxes and get endpoint CO2 as well
if (method == 'M1') {

  # grab footprint info
  foot.path <- file.path(workdir, 'out/by-id')
  foot.file <- file.path(foot.path, list.files(foot.path, pattern = 'foot.nc',
    recursive = T))

  #------------------------------ STEP 2.1 ----------------------------------- #
  # call function to convolve CT bio with footprints

}

# -------------------------- Regional daily median  -------------------------- #
if (method == 'M2H') {

  # minlon, maxlon, minlat, maxlat, used in Hakkareinen et al., 2016
  #reg.lon.lat <- c(-15, 60, 0, 60)  # For EU
  reg.lon.lat <- c(-130, -60, 0, 60)  # For US
  #reg.lon.lat <- c(60, 150, 0, 60)  # For Asia

  print(reg.lon.lat)
  mm <- ggplot.map(map = 'ggmap', center.lon = 25, center.lat = 20, zoom = 4)

  hakka.bg <- NULL
  for (t in 1:length(all.timestr)) {
    obs <- grab.oco2(oco2.path, all.timestr[t], reg.lon.lat) %>%
      #filter(wl <= 1) #%>%
      filter(qf == 0)
    m1 <- mm[[1]] + geom_point(data = obs, aes(lon, lat, colour = xco2),
      size = 0.4) + scale_colour_gradientn(colours = def.col(), name = 'XCO2')

    hakka.bg <- c(hakka.bg, median(obs$xco2))
  } # end for t
  hakka.bg <- data.frame(timestr = all.timestr, hakka.bg = hakka.bg)
  write.table(hakka.bg, file = paste0('M2H_bg_', site, '_', oco2.ver, '.txt'),
    quote = F, row.names = F, sep = ',')
} # end if M2H


# -------------------------- Normal statistics  -------------------------- #
if (method == 'M2S') {
  library(MASS)
  silva.bg <- NULL

  # 4x4 deg box around hotsplot, from Silva and Arellano, 2017
  lon.lat[1] <- lon.lat[5] - 2; lon.lat[2] <- lon.lat[5] + 2
  lon.lat[3] <- lon.lat[6] - 2; lon.lat[4] <- lon.lat[6] + 2
  print(lon.lat)

  for (t in 1:length(all.timestr)) {
    obs <- grab.oco2(oco2.path, all.timestr[t], lon.lat) %>% filter(qf == 0) #%>%
      #filter(wl <= 1)
    tmp.bg <- as.numeric(fitdistr(obs$xco2, 'normal')$estimate[1]) -
              as.numeric(fitdistr(obs$xco2, 'normal')$estimate[2])
    silva.bg <- c(silva.bg, tmp.bg)
  }
  silva.bg <- data.frame(timestr = all.timestr, silva.bg = silva.bg)
  write.table(silva.bg, file = paste0('M2S_bg_', site, '_', oco2.ver, '.txt'),
    quote = F, row.names = F, sep =',')
} # end if M2S


# ---------------------------- Overpass-specific ---------------------------- #
# need to run forward runs from a box around the city
if (method == 'M3') {

  bg.info <- NULL
  for (t in 1:length(all.timestr)) {
  #t = 1
    timestr <- all.timestr[t]
    print(timestr)

    #------------------------------ STEP 5.1 --------------------------------- #
    #### Whether forward/backward, release from a column or a box
    forwardTF  <- T  # forward or backward traj, if forward, release from a box
    delt       <- 2  # fixed timestep [min]; set = 0 for dynamic timestep
    run_trajec <- F  # whether to run forward traj, if T, will overwrite existing
    plotTF     <- T  # whether to calculate background and plot 2D density

    ### MUST-HAVE parameters about errors
    run_hor_err <- T  # run trajec with hor wind errors/calc trans error
    run_ver_err <- F  # run trajec with mixed layer height scaling

    # release particles from a box, +/-dxyp, +/-dxyp around the city center
    dxyp  <- 0.4               # degree around city center
    nhrs  <- 12                # since forward, nhrs should always be positive
    # release FROM ? hrs before obs time (-10), TO obs time (0), every 0.5 hour
    dtime <- seq(-10, 0, 0.5)  # time windows to release particles

    ### release particles from fixed levels
    agl    <- 10         # in mAGL
    numpar <- 1000       # particle # per each time window
    cat(paste('\n\nDone with receptor setup...\n'))

    #------------------------------ STEP 5.2 --------------------------------- #
    # path for the ARL format of WRF and GDAS
    # simulation_step() will find corresponding met files
    met.indx   <- 1
    met        <- c('hrrr', '1km', 'gdas0p5')[met.indx] # choose met fields
    met.path   <- file.path(homedir, 'u0947337', met)
    met.num    <- 1     # min number of met files needed

    # met file name convention
    met.format <- c('%Y%m%d.%Hz.hrrra', 'wrfout_', '%Y%m%d_gdas0p5')[met.indx]

    #### Whether obtaining wind errors, transport error component
    # require wind error comparisons stored in txtfile *****
    if (run_hor_err) {
      cat('+++ horizontal wind error component +++\n')

      # intput correlation lengthscale (in m) and timescales (in mins)
      # correlation timescale, horizontal and vertical lengthscales
      if (met == 'gdas0p5') {TLuverr <- 1*60; zcoruverr <- 600; horcoruverr <- 40}
      if (met == 'gdas1') {TLuverr <- 2.4*60; zcoruverr <- 700; horcoruverr <- 97}
      if (met == 'hrrr') {TLuverr <- 0.5*60; zcoruverr <- 400; horcoruverr <- 10}

      ## add errors, mainly siguverr, create a subroutine to compute siguverr
      # from model-data wind comparisons
      err.path <- file.path(homedir, 'lin-group5/wde/input_data/wind_err')
      err.path <- file.path(err.path, tolower(site), tolower(met))

      # call get.SIGUVERR() to interpolate most near-field wind errors
      #err.info <- get.siguverr(site, timestr, err.path, nfTF = F, forwardTF = T,
      #  lon.lat, nhrs)
      err.info <- NULL

      if (is.null(err.info)) {
        cat('no wind error found; make consevative assumption of siguverr...\n')
        # make a conservative assumption about the wind error, for the Middle East
        siguverr <- 2.5  # < 2 m/s for GDAS 1deg, based on Wu et al., GMDD

      } else {
        met.rad  <- err.info[[1]]
        siguverr <- as.numeric(err.info[[2]][1])    # SD in wind errors
        u.bias   <- as.numeric(err.info[[2]][2])
        v.bias   <- as.numeric(err.info[[2]][3])
        cat(paste('u.bias:', signif(u.bias,3), 'm/s; v.bias:',
          signif(v.bias,3), 'm/s\n'))
      }  # end if is.null(err.info)
      cat(paste('SIGUVERR:', signif(siguverr, 3), 'm/s..\n'))

    } else {  # if no wine error component used
      cat('NO horizontal wind error component for generating trajec...\n')
      siguverr    <- NA; TLuverr     <- NA
      horcoruverr <- NA; zcoruverr   <- NA
    }  # end if run_hor_err

    # no error covariance on PBL heights used for now
    # but one can assign values as below
    sigzierr <- NA; TLzierr <- NA; horcorzierr <- NA

    ### Besides horizontal wind error, do we want to account for PBL height?
    # add vertical trans error via ziscale *****
    if (run_ver_err) {
      zicontroltf <- 1              # 0 for FALSE; 1 for scaling, TRUE
      ziscale     <- rep(list(rep(0.8, 24)), nrecp)  # create as list
      # 1st # for scaling factor; 2nd # for # of hours (always use abs())
      cat(paste('+++ Mixed layer height scaling of', ziscale[[1]][1],
        'when generating trajec +++\n'))

    } else {
      cat('NO Mixed layer height scaling ...\n')
      zicontroltf <- 0
      ziscale     <- 1.0
    } # end if run_ver_err
    cat('Done with choosing met & inputting wind errors...\n')

    #------------------------------ STEP 5.3 -------------------------------- #
    # !!! need to add makefile for AER_NOAA_branch in fortran ;
    # link two hymodelcs in exe directory
    #clean.side <- c('south', 'both', 'north', 'south', 'north')[t]  # Riyadh
    clean.side <- 'south'  # all northern part for clean background for LV
    print(clean.side)

    if (run_trajec == T) outpath <- NULL  # will be in copies
    if (run_trajec == F)
      outpath <- file.path(homedir, 'lin-group5/wde/github/stilt', site,
        'out_forward/')

    # data filtering on observations
    data.filter <- c('QF', 0)
    tmp.info <- run.forward.trajec(site = site, timestr = timestr,
                                   overwrite = run_trajec, nummodel = t,
                                   lon.lat = lon.lat, delt = delt, dxyp = dxyp,
                                   dzp = 0, dtime = dtime, agl = agl,
                                   numpar = numpar, nhrs = nhrs,
                                   workdir = workdir, outpath = outpath,
                                   siguverr = siguverr, TLuverr = TLuverr,
                                   zcoruverr = zcoruverr,
                                   horcoruverr = horcoruverr,
                                   met.format = met.format, met.path = met.path,
                                   met.num = 1, plotTF = plotTF,
                                   oco2.path = oco2.path, oco2.ver = oco2.ver,
                                   zoom = 7, td = 0.05, perc = 0.05,
                                   clean.side = clean.side,
                                   data.filter = data.filter)
    print(tmp.info)
    bg.info <- rbind(bg.info, tmp.info)
  } # end for t

  write.table(bg.info, quote = F, row.names = F, sep = ',',
    file = file.path(outpath, paste0('M3_bg_', site, '_', oco2.ver, '.txt')))
} # end if method == 'M3'
