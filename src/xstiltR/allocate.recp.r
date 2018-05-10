# subroutine to allocate and assign receptors, given multiple nummodel
# write into control file for X-STILT
# DW, 05/09/2018

allocate.recp <- function(namelist){

  ## before writing, delete previous run_XSTILT and each trajlist* and CONTROL*
  if(length(list.files(namelist$workdir, "run_XSTILT"))>0){
    file.remove(list.files(namelist$workdir, "run_XSTILT"))
  }

  if(length(list.files(namelist$workdir, "XSTILT_CONTROL"))>0){
    file.remove(list.files(namelist$workdir, "XSTILT_CONTROL"))
  }

  if(length(list.files(namelist$workdir, "trajlist_"))>0){
    file.remove(list.files(namelist$workdir, "trajlist_"))
  }

  # create new names
  namelist$RUN <- paste("run_XSTILT_", namelist$timestr, sep="")
  namelist$CONTL <- paste("XSTILT_CONTROL_", namelist$timestr, "_copy", sep="")

  # store namelist
  str <- "backward"; if(namelist$forwardTF)str <- "forward"
  namelist$TRAJLIST <- paste("trajlist_", namelist$site, "_", namelist$timestr,
                       "_", str, sep="")
  trajlist <- file.path(namelist$workdir, namelist$TRAJLIST)  # link to path

  # how many copies?
  ncp <- length(namelist$nummodel)
  all.copy <- namelist$nummodel

  if(ncp > 1){

     # break receptors into different groups
     common.list <- namelist[-grep("recp.", names(namelist))]
     recp.list <- namelist[grep("recp.", names(namelist))]
     recp.list <- recp.list[names(recp.list)!= "recp.index"]

     # write CONTROL files for X-STILT, not STILT
     nrecp <- length(recp.list[[1]])
     sep.recp <- round(nrecp/ncp)
     group.index <- c(seq(1, nrecp, sep.recp), nrecp+1)

     for(n in 1:ncp){

       list.index <- group.index[n]:(group.index[n+1]-1)
       group.namelist <- sapply(recp.list, "[", list.index)

       # merge selected receptor info back to common.list
       merge.list <- common.list
       merge.list$recp.lat <- group.namelist[, "recp.lat"]
       merge.list$recp.lon <- group.namelist[, "recp.lon"]
       merge.list$recp.yr <- group.namelist[, "recp.yr"]
       merge.list$recp.mon <- group.namelist[, "recp.mon"]
       merge.list$recp.day <- group.namelist[, "recp.day"]
       merge.list$recp.hr <- group.namelist[, "recp.hr"]
       merge.list$nummodel <- all.copy[n]  # select copy number

       # output namelist into rds file
       group.filenm <- paste(trajlist, "_cp", n, ".rds", sep="")
       saveRDS(object=merge.list, file=group.filenm)

       # write startup r code as well
       cont <- paste("homedir='", namelist$homedir, "'\nworkdir='",
                     namelist$workdir,"'\nsetwd(workdir)\n",
                     "filenm='", group.filenm,"'\n",
                     "source(file.path(workdir, 'src', 'sourceall.r'))\n",
                     "namelist=readRDS(filenm)\n",
                     "run.backward.trajec(namelist=namelist)", sep="")
       cont.filenm <- file.path(namelist$workdir,
                              paste(namelist$CONTL, all.copy[n], ".r", sep=""))
       write(cont, cont.filenm)
     } # end for n

  }else{

    # if using onely 1 copy
    one.filenm <- paste(trajlist, "_cp", n, ".rds", sep="")
    saveRDS(object=namelist, file=one.filenm)

    # write startup r code as well
    cont <- paste("homedir='", namelist$homedir, "'\nworkdir='",
                  namelist$workdir,"'\nsetwd(workdir)",
                  "filenm='", one.filenm,"'\n",
                  "source(file.path(workdir, 'src', 'sourceall.r'))\n",
                  "namelist=readRDS(filenm)\n",
                  "run.backward.trajec(namelist=namelist)", sep="")
    cont.filenm <- file.path(namelist$workdir, paste(namelist$CONTL, all.copy,
                                                     ".r", sep=""))
    write(cont, cont.filenm)
  }  # end if one or multiple nodes

  cat(paste("allocate.recp(): Done with X-STILT start-up r code ...\n"))

  #### next, write the parallel run file as well
  parfile <- file.path(namelist$workdir, namelist$RUN)

  write("#!/bin/bash\n", parfile)
  write('echo "user is $USER"\n', parfile, append=T)

  # write lines for loading CHPC modules
  mods <- c("R/3.4.3", "ncl", "netcdf-f/4.4.4", "pgi/17.4", "mvapich2/2.2.p",
             "gcc/4.8.5")
  module <- paste("module load", mods)
  write(module, file=parfile, sep="\n", append=T)

  # go to working directory
  write(paste('\ncd ', namelist$workdir, '\n', sep=""), parfile, append=T)

  # point to different CONTROL R scripts
  nodes <- seq(0, ncp-1, 1)
  r.scripts <- paste("export script_", nodes, "=",
                     namelist$CONTL, all.copy, ".r", sep="")
  write(r.scripts, parfile, append=T)

  # print some info
  echo <- paste('\necho "Calculation (Core:$1) started at `date`"\necho $1',
                '\necho "You are in `pwd`"\n', sep="")
  write(echo, parfile, append=T)

  # excutable R codes
  exec <- paste('\nRscript $(eval echo \\$script_${1})\n',
                'echo "Calculation (Core:$1 ) stopped at `date`"\n', sep="")
  write(exec, parfile, append=T)
  cat(paste("allocate.recp(): Done with run_X-STILT...\n"))


  # finally, turn on permission of all files
  system(paste("chmod u+x ", parfile, sep=""))
  cont.files <- list.files(namelist$workdir, "XSTILT_CONTROL")
  Sys.chmod(cont.files, mode = "744")  # x for me

  return(namelist)
}  # end subroutine
