# STILT Dependency Loader
# For documentation, see https://github.com/uataq/stilt
# Ben Fasoli
# add Dien's OCO-2/XSTILT source codes, 05/23/2018, DW

####
if (!'stilt_wd' %in% ls())
  stilt_wd <- getwd()

### Source r/src R scripts for Ben's STILT-R v2
rsc <- dir(file.path(stilt_wd, 'r', 'src'), pattern = '.*\\.r$',
  full.names = T)

### Source r/src R scripts for Dien's X-STILT, 05/23/2018, DW
rsc <- c(rsc, dir(file.path(stilt_wd, 'r', 'src', 'oco2-xstilt'),
  pattern = '.*\\.r$', full.names = T, recursive = T))

invisible(lapply(rsc, source))

### Load external libraries
if (!'lib.loc' %in% ls())
  lib.loc <- NULL
load_libs('dplyr', 'ncdf4', 'parallel', 'raster', 'readr', 'rslurm', 'uataq',
          'ggplot2', 'Hmisc', 'ggmap', 'geosphere', 'reshape2', 'PBSmapping',
          'sp', 'MASS', 'stringr',
          lib.loc = lib.loc)
# add few other packages for OCO-2/XSTILT, DW, 05/23/2018

### Load permute fortran dll for footprint matrix permutation
permute_exe <- file.path(stilt_wd, 'r/src/permute.so')
if (!file.exists(permute_exe))
  stop('calc_footprint(): failed to find permute.so in r/src/')
dyn.load(permute_exe)
