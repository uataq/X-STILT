# subroutine to allocate and assign receptors, given multiple nummodel
# write into control file for X-STILT
# DW, 05/09/2018

allocate.recp <- function(namelist, run.name="run_XSTILT"){

  # store namelist to output directory
  filenm <- paste("trajlist_", namelist$site, "_", namelist$timestr, sep="")
  if(namelist$forwardTF == F)filenm <- paste(filenm, "_backward.txt", sep="")
  if(namelist$forwardTF)filenm <- paste(filenm, "_forward.txt", sep="")

  filenm <- file.path(namelist$workdir, filenm)  # link to path

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
     group.index <- seq(1, nrecp, sep.recp)

     for(n in 1:ncp){

       #cat(paste("working on allocating copy #", all.copy[n],"...\n"))

       list.index <- group.index[n]:(group.index[n+1]-1)
       group.namelist <- sapply(recp.list, "[", list.index)

       merge.list <- common.list
       merge.list$recp.lat <- group.namelist[, "recp.lat"]
       merge.list$recp.lon <- group.namelist[, "recp.lon"]
       merge.list$recp.yr <- group.namelist[, "recp.yr"]
       merge.list$recp.mon <- group.namelist[, "recp.mon"]
       merge.list$recp.day <- group.namelist[, "recp.day"]
       merge.list$recp.hr <- group.namelist[, "recp.hr"]
       merge.list$nummodel <- all.copy[n]  # select copy number

       # output namelist into rds file
       group.filenm <- gsub(".txt", paste("_cp", n, ".rds", sep=""), filenm)
       saveRDS(object=merge.list, file=group.filenm)

       # write startup r code as well
       cont <- paste("homedir='", namelist$homedir, "'\nworkdir='",
                     namelist$workdir,"'\nsetwd(workdir)\n",
                     "filenm='", group.filenm,"'\n",
                     "source(file.path(workdir, 'src', 'sourceall.r'))\n",
                     "namelist=readRDS(filenm)\n",
                     "run.backward.trajec(namelist=namelist)", sep="")
       cont.filenm <- file.path(namelist$workdir,
                                paste("XSTILT_CONTROL_copy", all.copy[n],
                                      ".r", sep=""))
       write(cont, cont.filenm)
       #cat(paste("allocate.recp(): Done writing X-STILT start-up r code ...\n"))
     } # end for n

  }else{

    # if using onely 1 copy
    one.filenm <- gsub(".txt", paste("_cp", n, ".rds", sep=""), filenm)
    saveRDS(object=namelist, file=one.filenm)

    # write startup r code as well
    cont <- paste("homedir='", namelist$homedir, "'\nworkdir='",
                  namelist$workdir,"'\nsetwd(workdir)",
                  "filenm='", one.filenm,"'\n",
                  "source(file.path(workdir, 'src', 'sourceall.r'))\n",
                  "namelist=readRDS(filenm)\n",
                  "run.backward.trajec(namelist=namelist)", sep="")
    cont.filenm <- file.path(namelist$workdir, paste("XSTILT_CONTROL_copy",
                                                    all.copy, ".r", sep=""))
    write(cont, cont.filenm)
    #cat(paste("allocate.recp(): Done writing X-STILT start-up r code ...\n"))
  }  # end if
  cat(paste("allocate.recp(): Done with X-STILT start-up r code ...\n"))

  #### next, write the parallel txt file as well
  parfile <- file.path(namelist$workdir, run.name)

  write("#!/bin/bash\n", parfile)
  write('echo "user is $USER"\n', parfile, append=T)

  mods <- c("R/3.4.3", "ncl", "netcdf-f/4.4.4", "pgi/17.4", "mvapich2/2.2.p",
             "gcc/4.8.5")
  module <- paste("module load", mods)
  write(module, file=parfile, sep="\n", append=T)

  write(paste('\ncd ', namelist$workdir, '\n', sep=""), parfile, append=T)

  # point to different CONTROL R scripts
  nodes <- seq(0, ncp-1, 1)
  export.r <- paste("export script_", nodes, "=XSTILT_CONTROL_copy",
                    all.copy, ".r", sep="")
  write(export.r, parfile, append=T)

  # excutable codes
  #Rscript $(eval echo \$script_${1})

  echo <- paste('\necho "Calculation (Core:$1) started at `date`"\necho $1',
                '\necho "You are in `pwd`"\n', sep="")

  write(echo, parfile, append=T)

  exec <- paste('\nRscript $(eval echo \\$script_${1})\n',
                'echo "Calculation (Core:$1 ) stopped at `date`"\n', sep="")
  write(exec, parfile, append=T)
 cat(paste("allocate.recp(): Done with run_X-STILT...\n"))

}  # end subroutine
