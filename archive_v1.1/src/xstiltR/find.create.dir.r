# subroutine to search a directory
# if not exist, create one
# Dien Wu, 05/02/2018

# workdir: original directory
# path and dir: searching/creating a folder under a given path
find.create.dir <- function(path, dir, workdir){

  if(length(list.files(path = path, pattern = dir)) == 0){
    setwd(path)   # move to 'path'
    system(paste0('mkdir ', dir))
    cat(paste(dir, 'not found; just created...\n'))

    # remember to go back to original working directory
    setwd(workdir)

  }else{
    cat(paste(dir,'existed...\n'))

  } # end if

  return(file.path(path, dir))   # return the whole path
}
