# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli
# edited by Dien Wu for X-STILT


# preserve original 'run_stilt.r' and modify 'run_stilt.r' as a subroutine
# add ziscale, DW, 07/25/2018
# add two before_* functions for getting ground height and AK PW weighting, DW, 01/24/2019
# allows for generate footprints with different resolutions, DW, 02/11/2019 
# add additional variables due to changes in simulation_step()'s argument, DW, 07/02/2020
# add/remove variables as migrating to STILT-HYSPLIT, DW, 07/03/2020 


run.xstilt <- function(namelist){

  # User inputs ----------------------------------------------------------------
  stilt_wd  <- namelist$xstilt_wd
  output_wd <- namelist$outdir
  lib.loc   <- .libPaths()[1]
  receptors <- namelist$recp.info

  # met fields 
  met_directory   <- namelist$met.path
  met_file_format <- namelist$met.format
  n_met_min <- namelist$met.num

  # modification for OCO-2 column simulations, DW, 05/24/2018
  # if NA, meaning no weighting nor OCO2 files for regular simualtions
  ak.wgt   <- namelist$ak.wgt     # whether weighted foot by averaging kernel
  pwf.wgt  <- namelist$pwf.wgt    # whether weighted foot by pres weighting
  oco.path <- namelist$oco.path
  
  # Model control
  rm_dat     <- T
  timeout    <- namelist$timeout  # in sec
  run_trajec <- namelist$run_trajec
  run_foot   <- namelist$run_foot
  n_hours    <- namelist$nhrs
  numpar     <- namelist$numpar
  varsiwant  <- namelist$varstrajec
  if (length(varsiwant) == 0)
    varsiwant  <- c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 'foot', 
                    'mlht', 'dens', 'samt', 'sigw', 'tlgr', 'temp', 'pres')

  # Transport error and PBL error input variables
  horcoruverr <- namelist$hor.err$horcoruverr
  siguverr    <- namelist$hor.err$siguverr
  tluverr     <- namelist$hor.err$TLuverr
  zcoruverr   <- namelist$hor.err$zcoruverr

  horcorzierr <- namelist$pbl.err$horcorzierr
  sigzierr    <- namelist$pbl.err$sigzierr
  tlzierr     <- namelist$pbl.err$TLzierr

  zicontroltf <- namelist$pbl.err$zicontroltf
  ziscale     <- namelist$pbl.err$ziscale   # prescribe PBL scaling, vector form

  # Footprint grid settings
  xmn   <- as.numeric(namelist$foot.info$xmn)
  xmx   <- as.numeric(namelist$foot.info$xmx)
  ymn   <- as.numeric(namelist$foot.info$ymn)
  ymx   <- as.numeric(namelist$foot.info$ymx)
  xres  <- as.numeric(namelist$foot.info$xres)
  yres  <- as.numeric(namelist$foot.info$yres)

  # footprint with diff resolution, can be NA or numbers
  xres2 <- namelist$foot.info$xres2
  yres2 <- namelist$foot.info$yres2
  foot.nhrs <- namelist$foot.info$foot.nhrs 
  hnf_plume      <- namelist$hnf_plume
  smooth_factor  <- namelist$smooth_factor
  time_integrate <- namelist$time_integrate
  projection     <- namelist$projection 

  # Transport and dispersion settings, use default setting in STILT-R v2
  capemin     <- -1
  cmass       <- 0
  conage      <- 48
  cpack       <- 1
  delt        <- namelist$delt
  dxf         <- 1
  dyf         <- 1
  dzf         <- 0.01
  efile       <- ''
  emisshrs    <- 0.01
  frhmax      <- 3
  frhs        <- 1
  frme        <- 0.1
  frmr        <- 0
  frts        <- 0.1
  frvs        <- 0.01
  hscale      <- 10800
  ichem       <- 8
  idsp        <- 2
  initd       <- 0
  k10m        <- 1
  kagl        <- 1
  kbls        <- 1
  kblt        <- 5
  kdef        <- 0
  khinp       <- 0
  khmax       <- 9999
  kmix0       <- 250
  kmixd       <- 3
  kmsl        <- 0
  kpuff       <- 0
  krand       <- 4
  krnd        <- 6
  kspl        <- 1
  kwet        <- 1
  kzmix       <- 0
  maxdim      <- 1
  maxpar      <- numpar
  mgmin       <- 10
  mhrs        <- 9999
  nbptyp      <- 1
  ncycl       <- 0
  ndump       <- 0
  ninit       <- 1
  nstr        <- 0
  nturb       <- 0
  nver        <- 0
  outdt       <- 0
  p10f        <- 1
  pinbc       <- ''
  pinpf       <- ''
  poutf       <- ''
  qcycle      <- 0
  rhb         <- 80
  rht         <- 60
  splitf      <- 1
  tkerd       <- 0.18
  tkern       <- 0.18
  tlfrac      <- 0.1
  tout        <- 0
  tratio      <- 0.75
  tvmix       <- 1
  veght       <- 0.5
  vscale      <- 200
  vscaleu     <- 200
  vscales     <- -1
  wbbh        <- 0
  wbwf        <- 0
  wbwr        <- 0
  wvert       <- FALSE
  w_option    <- 0
  zicontroltf <- 0
  ziscale     <- rep(list(rep(1, 24)), nrow(receptors))
  z_top       <- 25000

  # customized functions
  # before_*_xstilt() are two customized functions for OCO-2/XSTILT
  # they are loaded after sourcing the dependency, 
  # if columnTF is turned off, initial it with identity functions, DW, 01/24/2019
  if (columnTF == F) {
    before_footprint <- function() {output}
    before_trajec <- function() {output}
  } else {
    before_footprint <- before_footprint_xstilt
    before_trajec <- before_trajec_xstilt
  } # end if columnTF

  # Parallel simulation settings
  slurm_options <- namelist$slurm_options
  slurm         <- namelist$slurm
  n_cores       <- namelist$n_cores
  n_nodes       <- namelist$n_nodes
  jobname       <- namelist$jobname
  
  # Startup messages -----------------------------------------------------------
  message('Initializing STILT')
  message('Number of receptors: ', nrow(receptors))
  message('Number of parallel threads: ', n_nodes * n_cores)


  # Source dependencies --------------------------------------------------------
  setwd(stilt_wd)
  source('r/dependencies.r')

  # Structure out directory ----------------------------------------------------
  # Outputs are organized in three formats. by-id contains simulation files by
  # unique simulation identifier. particles and footprints contain symbolic
  # links to the particle trajectory and footprint files in by-id
  system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = T)
  if (run_trajec) {
    system(paste0('rm -r ', output_wd, '/by-id'), ignore.stderr = T)
    system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = T)
  }
  for (d in c('by-id', 'particles', 'footprints')) {
    d <- file.path(output_wd, d)
    if (!file.exists(d))
      dir.create(d, recursive = T)
  }

  # Met path symlink -----------------------------------------------------------
  # Auto symlink the meteorological data path to the working directory to
  # eliminate issues with long (>80 char) paths in fortran. Note that this
  # assumes that all meteorological data is found in the same directory.
  if ((nchar(paste0(met_directory, met_file_format)) + 2) > 80) {
    met_loc <- file.path(path.expand('~'), paste0('m', namelist$project))
    if (!file.exists(met_loc)) invisible(file.symlink(met_directory, met_loc))
  } else met_loc <- met_directory

  # removed 'isot', 'iconvert', 'outfrac', 'random' -- 
  # not used by the latest STILT
  output <- xstilt_apply(FUN = simulation_step,
                         slurm = slurm, 
                         slurm_options = slurm_options,
                         n_cores = n_cores,
                         n_nodes = n_nodes,
                         jobname = jobname, 
                         before_footprint = list(before_footprint),
                         before_trajec = list(before_trajec),
                         capemin = capemin,
                         cmass = cmass,
                         conage = conage,
                         cpack = cpack,
                         delt = delt, 
                         dxf = dxf,
                         dyf = dyf,
                         dzf = dzf,
                         efile = efile,
                         emisshrs = emisshrs,
                         frhmax = frhmax, 
                         frhs = frhs, 
                         frme = frme, 
                         frmr = frmr,
                         frts = frts, 
                         frvs = frvs, 
                         hnf_plume = hnf_plume,
                         horcoruverr = horcoruverr, 
                         horcorzierr = horcorzierr,
                         hscale = hscale,
                         ichem = ichem, 
                         idsp = idsp,
                         initd = initd,
                         k10m = k10m,
                         kagl = kagl, 
                         kbls = kbls, 
                         kblt = kblt, 
                         kdef = kdef,
                         khinp = khinp,
                         khmax = khmax, 
                         kmix0 = kmix0, 
                         kmixd = kmixd,
                         kmsl = kmsl, 
                         kpuff = kpuff, 
                         krand = krand,
                         krnd = krnd, 
                         kspl = kspl,
                         kwet = kwet,
                         kzmix = kzmix, 
                         lib.loc = lib.loc,
                         maxdim = maxdim, 
                         maxpar = maxpar,
                         met_file_format = met_file_format,
                         met_loc = met_loc,
                         mgmin = mgmin,
                         n_hours = n_hours,
                         n_met_min = n_met_min,
                         ncycl = ncycl,
                         ndump = ndump,
                         ninit = ninit,
                         nstr = nstr,
                         nturb = nturb,
                         numpar = numpar, 
                         nver = nver,
                         outdt = outdt,
                         output_wd = output_wd,
                         p10f = p10f,
                         pinbc = pinbc,
                         pinpf = pinpf,
                         poutf = poutf,
                         projection = projection,
                         qcycle = qcycle,
                         r_run_time = receptors$run_time,
                         r_lati = receptors$lati,
                         r_long = receptors$long,
                         r_zagl = receptors$zagl,
                         rhb = rhb,
                         rht = rht,
                         rm_dat = rm_dat,
                         run_foot = run_foot,
                         run_trajec = run_trajec, 
                         siguverr = siguverr,
                         sigzierr = sigzierr, 
                         smooth_factor = smooth_factor,
                         splitf = splitf,
                         stilt_wd = stilt_wd,
                         time_integrate = time_integrate,
                         timeout = timeout,
                         tkerd = tkerd,
                         tkern = tkern,
                         tlfrac = tlfrac,
                         tluverr = tluverr,
                         tlzierr = tlzierr,
                         tout = tout,
                         tratio = tratio,
                         tvmix = tvmix,
                         varsiwant = list(varsiwant),
                         veght = veght,
                         vscale = vscale,
                         vscaleu = vscaleu,
                         vscales = vscales,
                         w_option = w_option,
                         wbbh = wbbh,
                         wbwf = wbwf,
                         wbwr = wbwr,
                         wvert = wvert,
                         xmn = xmn,
                         xmx = xmx,
                         xres = xres,
                         ymn = ymn,
                         ymx = ymx,
                         yres = yres,
                         zicontroltf = zicontroltf,
                         ziscale = ziscale,
                         z_top = z_top,
                         zcoruverr = zcoruverr, 

  # pass additional variables to stilt_apply and then to simulation_step() 
  # needed for before_*_xstilt() for X-STILT, DW, 02/11/2019
                        oco.path = oco.path, 
                        ak.wgt = ak.wgt, 
                        pwf.wgt = pwf.wgt, 
                        overwrite_wgttraj = namelist$overwrite_wgttraj, 
                        run_hor_err = namelist$run_hor_err,
                        emiss.file = namelist$emiss.file, 
                        met = namelist$met, 
                        ct.ver = namelist$ct.ver, 
                        ctflux.path = namelist$ctflux.path, 
                        ctmole.path = namelist$ctmole.path, 
                        xres2 = list(xres2), 
                        yres2 = list(yres2), 
                        foot.nhrs = foot.nhrs)
}


# for debugging
if (F) {
  
  X = 1
  r_run_time = receptors$run_time[X]
  r_lati = receptors$lati[X]
  r_long = receptors$long[X]
  r_zagl = receptors$zagl[X]


  args <- list(oco.path = oco.path, ak.wgt = ak.wgt, pwf.wgt = pwf.wgt, 
               overwrite_wgttraj = namelist$overwrite_wgttraj, 
               run_hor_err = namelist$run_hor_err,
               emiss.file = namelist$emiss.file, met = namelist$met, 
               ct.ver = namelist$ct.ver, ctflux.path = namelist$ctflux.path, 
               ctmole.path = namelist$ctmole.path, xres2 = list(xres2), 
               yres2 = list(yres2), foot.nhrs = foot.nhrs)

}

