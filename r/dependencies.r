# STILT Dependency Loader
# For documentation, see https://github.com/uataq/stilt
# Ben Fasoli
# add Dien's OCO-2/XSTILT source codes, 05/23/2018, DW

####
if (!'xstilt_wd' %in% ls()) xstilt_wd <- getwd()

### Source stilt R functions
rsc <- dir(file.path(xstilt_wd, 'stilt_hysplit/r/src'), 
           pattern = '.*\\.r$', full.names = T)

### Source R scripts for Dien's X-STILT, 05/23/2018, DW
rsc <- c(rsc, dir(file.path(xstilt_wd, 'r/src'), pattern = '.*\\.r$',
                  full.names = T, recursive = T))

invisible(lapply(rsc, source))

### Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL

# add few other packages for OCO-2/XSTILT, DW, 05/23/2018
#devtools::install_github("SESYNC-ci/rslurm")

load_libs('dplyr', 'ncdf4', 'parallel', 'raster', 'readr', 'rslurm', 'ggplot2', 
          'ggmap', 'geosphere', 'reshape2', 'stringr', 'MASS', #'PBSmapping', 
          'ggpubr', 'rworldmap', 'lutz', lib.loc = lib.loc)

### Load permute fortran dll for footprint matrix permutation
permute_exe <- file.path(xstilt_wd, 'stilt_hysplit/r/src/permute.so')
if (!file.exists(permute_exe)) {
  cat('calc_footprint(): failed to find permute.so in stilt_hysplit/r/src/, compile now...\n')
  system(paste0('R CMD SHLIB ', xstilt_wd, '/stilt_hysplit/r/src/permute.f90'))
}
dyn.load(permute_exe)

# Validate arguments and load dependencies if necessary
if ((!class(projection) == 'function') && ('projection' %in% ls()))
  validate_projection(projection)
if ('varsiwant' %in% ls())
  validate_varsiwant(varsiwant)
if (all(c('xmn', 'xmx', 'ymn', 'ymx') %in% ls()))
  validate_footprint_extent(xmn, xmx, ymn, ymx)