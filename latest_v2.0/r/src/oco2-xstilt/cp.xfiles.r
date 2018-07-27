
# script to copy files
# v for files of footprint, traj, or weighted traj
cp.xfiles <- function(workdir, v, nhrs.back, dpar, site, cpTF = F){

  from.pat <- c('_foot.nc', '_traj.rds', '_wgttraj.rds')
  from.var <- c('footprints', 'particles', 'by-id/*_X')

  from.path <- file.path(workdir, 'out', from.var[v])
  from.file <- list.files(path = from.path, pattern = from.pat[v])

  to.var   <- c('Footprints', 'Particles', 'Particles')
  to.path <- file.path(
    '/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2',
    to.var[v], site)

  if (cpTF) {
    cat(paste('cp.xfiles(): Copy from', from.path, '\n paste to', to.path, '\n\n'))
    system(paste0('cp ', from.path, '/*', from.pat[v], ' ', to.path))
  }

  #to.file <- list.files(path = to.path, pattern = from.pat[v])
  return(to.path)
}
