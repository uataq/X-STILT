# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli

# preserve original 'run_stilt.r' and modify 'run_stilt.r' as a subroutine
# add ziscale, DW, 07/25/2018
# add two before_* functions for getting ground height and AK PW weighting, DW, 01/24/2019

run.xstilt <- function(namelist){

  # User inputs ----------------------------------------------------------------
  stilt_wd  <- namelist$workdir
  output_wd <- namelist$outdir
  lib.loc   <- .libPaths()[1]
  receptors <- namelist$recp.info

  # met fields 
  met_directory   <- namelist$met.path
  met_file_format <- namelist$met.format
  n_met_min <- namelist$met.num

  # modification for OCO-2 column simulations, DW, 05/24/2018
  # if NA, meaning no weighting nor OCO2 files for regular simualtions
  ak.wgt    <- namelist$ak.wgt     # whether weighted foot by averaging kernel
  pwf.wgt   <- namelist$pwf.wgt    # whether weighted foot by pres weighting
  oco2.path <- namelist$oco2.path

  # Model control
  rm_dat     <- T
  timeout    <- namelist$timeout  # in sec
  run_trajec <- namelist$run_trajec
  run_foot   <- namelist$run_foot
  n_hours    <- namelist$nhrs
  numpar     <- namelist$numpar
  varsiwant  <- namelist$varstrajec
  if (length(varsiwant) == 0) {
    varsiwant  <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr',
                    'zsfc', 'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld',
                    'dmas', 'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc',
                    'dswf', 'wout', 'mlht', 'rain', 'crai', 'pres')
  }

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
  xmn            <- as.numeric(namelist$foot.info$xmn)
  xmx            <- as.numeric(namelist$foot.info$xmx)
  ymn            <- as.numeric(namelist$foot.info$ymn)
  ymx            <- as.numeric(namelist$foot.info$ymx)
  xres           <- as.numeric(namelist$foot.info$xres)
  yres           <- as.numeric(namelist$foot.info$yres)
  hnf_plume      <- namelist$hnf_plume
  smooth_factor  <- namelist$smooth_factor
  time_integrate <- namelist$time_integrate
  projection     <- namelist$projection 

  # Transport and dispersion settings
  conage      <- 48
  cpack       <- 1
  delt        <- namelist$delt
  dxf         <- 1
  dyf         <- 1
  dzf         <- 0.1
  emisshrs    <- 0.01
  frhmax      <- 3
  frhs        <- 1
  frme        <- 0.1
  frmr        <- 0
  frts        <- 0.1
  frvs        <- 0.1
  hscale      <- 10800
  ichem       <- 0
  iconvect    <- 0
  initd       <- 0
  isot        <- 0
  kbls        <- 1
  kblt        <- 1
  kdef        <- 1
  khmax       <- 9999
  kmix0       <- 250
  kmixd       <- 3
  kmsl        <- 0
  kpuff       <- 0
  krnd        <- 6
  kspl        <- 1
  kzmix       <- 1
  maxdim      <- 1
  maxpar      <- min(100000, numpar)
  mgmin       <- 2000
  ncycl       <- 0
  ndump       <- 0
  ninit       <- 1
  nturb       <- 0
  outdt       <- 0
  outfrac     <- 0.9
  p10f        <- 1
  qcycle      <- 0
  random      <- 1
  splitf      <- 1
  tkerd       <- 0.18
  tkern       <- 0.18
  tlfrac      <- 0.1
  tratio      <- 0.9
  tvmix       <- 1
  veght       <- 0.5
  vscale      <- 200
  w_option    <- 0
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

  # require list form for stilt_apply.r
  before_footprint <- list(before_footprint)
  before_trajec <- list(before_trajec)

  # Parallel simulation settings
  slurm_options <- namelist$slurm_options
  slurm         <- namelist$slurm
  n_cores       <- namelist$n_cores
  n_nodes       <- namelist$n_nodes

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
    met_loc <- file.path(path.expand('~'), paste0('m', project))
    if (!file.exists(met_loc)) invisible(file.symlink(met_directory, met_loc))
  } else met_loc <- met_directory

  output <- stilt_apply(FUN = simulation_step,
                        slurm = slurm, 
                        slurm_options = slurm_options,
                        n_cores = n_cores,
                        n_nodes = n_nodes,
                        before_footprint = list(before_footprint),
                        before_trajec = list(before_trajec),
                        conage = conage,
                        cpack = cpack,
                        delt = delt, 
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
                        ichem = ichem, 
                        iconvect = iconvect, 
                        initd = initd,
                        isot = isot, 
                        kbls = kbls, 
                        kblt = kblt, 
                        kdef = kdef,
                        khmax = khmax, 
                        kmix0 = kmix0, 
                        kmixd = kmixd,
                        kmsl = kmsl, 
                        kpuff = kpuff, 
                        krnd = krnd, 
                        kspl = kspl,
                        kzmix = kzmix, 
                        maxdim = maxdim, 
                        maxpar = maxpar,
                        lib.loc = lib.loc,
                        met_file_format = met_file_format,
                        met_loc = met_loc,
                        mgmin = mgmin,
                        n_hours = n_hours,
                        n_met_min = n_met_min,
                        ncycl = ncycl,
                        ndump = ndump,
                        ninit = ninit,
                        nturb = nturb,
                        numpar = numpar, 
                        outdt = outdt,
                        outfrac = outfrac,
                        output_wd = output_wd,
                        p10f = p10f,
                        projection = projection,
                        qcycle = qcycle,
                        r_run_time = receptors$run_time,
                        r_lati = receptors$lati,
                        r_long = receptors$long,
                        r_zagl = receptors$zagl,
                        random = random,
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
                        tratio = tratio,
                        tvmix = tvmix,
                        varsiwant = list(varsiwant),
                        veght = veght,
                        vscale = vscale,
                        w_option = w_option,
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
  # needed for X-STILT, DW, 01/25/2019
                        oco2.path = oco2.path, 
                        ak.wgt = ak.wgt, 
                        pwf.wgt = pwf.wgt)

}
