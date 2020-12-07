# STILT Dependency Loader
# For documentation, see https://github.com/uataq/stilt
# Ben Fasoli
# add Dien's OCO-2/XSTILT source codes, 09/18/2020, DW

####
if (!'xstilt_wd' %in% ls()) xstilt_wd <- getwd()
stilt_wd <- file.path(xstilt_wd, 'stilt')

# if fortran dependencies not exist, DW, 11/06/2020
permute_exe <- file.path(stilt_wd, 'r/src/permute.so')
hycs_std <- file.path(xstilt_wd, 'exe/hycs_std')

if ( !file.exists(permute_exe) | !file.exists(hycs_std) ) {
  setwd(file.path(xstilt_wd, 'stilt'))
  cat('need to setup STILT\n')
  system('chmod +x setup')
  system('./setup')

  if ( !file.exists(hycs_std) ) {
    cat('linking STILT executables & dependencies to correct X-STILT directory...\n')
    system(paste('ln -s', file.path(stilt_wd, 'exe/*'), 
                          file.path(xstilt_wd, 'exe/')))
  }

  setwd(xstilt_wd)
}



### Source stilt R functions
rsc <- dir(file.path(stilt_wd, 'r/src'), pattern = '.*\\.r$', full.names = T)

### Source R scripts for Dien's X-STILT, 05/23/2018, DW
rsc <- c(rsc, dir(file.path(xstilt_wd, 'r/src'), pattern = '.*\\.r$',
                  full.names = T, recursive = T))

invisible(lapply(rsc, source))

### Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL

# add few other packages for OCO-2/XSTILT, DW, 05/23/2018
load_libs('dplyr', 'ncdf4', 'parallel', 'raster', 'readr', 'ggplot2', 'ggmap', 
          'geosphere', 'reshape2', 'stringr', 'MASS', 'rworldmap', 'lutz', 
          'rslurm', lib.loc = lib.loc)

### Load permute fortran dll for footprint matrix permutation
if (!file.exists(permute_exe))
  stop('calc_footprint(): failed to find permute.so in stilt/r/src/')
dyn.load(permute_exe)

# Validate arguments and load dependencies if necessary
if ((!class(projection) == 'function') && ('projection' %in% ls()))
  validate_projection(projection)
if ('varsiwant' %in% ls())
  validate_varsiwant(varsiwant)
if (all(c('xmn', 'xmx', 'ymn', 'ymx') %in% ls()))
  validate_footprint_extent(xmn, xmx, ymn, ymx)