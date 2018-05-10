# subroutine for running X-STILT trajec
# Written by DW, update on 05/19/2017 (add GDAS-driven STILT)
# ADD fixed agl run, 07/07/2017, DW
# clear up, 04/30/2017, DW

# read variables from output of "create_namelist_trajec.r"
run.backward.trajec <- function(namelist){

  # ------------------- CREATE outname for trajec -------------------- #
  if(namelist$columnTF){

    # If run with multiple locations, outname needs to be specified first
    # Set traj name convention, ONLY for uneven AGLs,
    # ONLY for if yr, mon, day, hr are the same
    #"2014x12x29x10x24.9836Nx47.0343Ex00000-03000by00100&03500-06000by00500x07400P"
    outname <- paste(namelist$filestr,"x", formatC(namelist$recp.hr[1],width=2,flag=0),"x",
                     formatC(namelist$recp.lat, width=5,format='f',digits=4,flag=0), "Nx",
                     formatC(namelist$recp.lon, width=5,format='f',digits=4,flag=0), "Ex",
                     formatC(namelist$minagl, width=5, flag=0), "-",
                     formatC(namelist$zi, width=5,flag=0), "by",
                     formatC(namelist$dh[1], width=5, flag=0), "&",
                     formatC(namelist$zi + namelist$dh[2], width=5,flag=0), "-",
                     formatC(max(namelist$agl), width=5, flag=0), "by",
                     formatC(namelist$dh[2], width=5, flag=0), "x",
                     formatC(namelist$npar, width=5, flag=0), "P", sep="")

  }else{

    # if fixed agl, no need to assign outname
    outname <- paste(namelist$filestr,"x", formatC(namelist$recp.hr[1],width=2,flag=0),"x",
                     formatC(namelist$recp.lat, width=5,format='f',digits=4,flag=0), "Nx",
                     formatC(namelist$recp.lon, width=5,format='f',digits=4,flag=0), "Ex",
                     formatC(namelist$agl, width=5, flag=0),"x",
                     formatC(namelist$npar, width=5, flag=0),"P",sep="")

  }  # end if column or not


  # -------------------  START running STILT trajec ------------------- #
  ### looping over selected & sorted soundings
  metlib <- paste(namelist$metpath,"/",sep="")

  for (t in 1:length(namelist$recp.lat)){

    cat("#------ Generating Traj", signif(t/length(namelist$recp.lat)*100), "% ------#\n")

    # GENERATE normal STILT trajectories, using "Trajec.r"
    # CALL trajec, return a list of input settings & store trajec in RData file
    # fixed receptors, and can loop over multiple fixed AGLs
    if(namelist$columnTF == FALSE){

      for(h in 1 : length(namelist$agl)){
       ident <- outname[h]
       info <- Trajec(yr=namelist$recp.yr[t], mon=namelist$recp.mon[t],
                      day=namelist$recp.day[t], hr=namelist$recp.hr[t],
                      lat=namelist$recp.lat[t], lon=namelist$recp.lon[t],
                      conv=namelist$convect, numpar=namelist$npar,
                      agl=namelist$agl[h], outname=ident, metd=namelist$metd,
                      mgmin=namelist$mgmin, nhrs=namelist$nhrs,
                      delt=namelist$delt, doublefiles=T, rundir=namelist$rundir,
                      nummodel=namelist$nummodel, veght=namelist$veght,
                      varsout=namelist$varstrajec, overwrite=namelist$overwrite,
                      metlib=metlib, metfile=namelist$metfile,
                      outpath=namelist$orig.trajpath)
      } # end loop h

    }else{   # column receptors

      # if multiple locations,
      # lat, lon, & agl have to be vectors with the same length!
      # multiple heights at the same lat, lon coordinate
      tmp.lat <- rep(namelist$recp.lat[t], length(namelist$agl))
      tmp.lon <- rep(namelist$recp.lon[t], length(namelist$agl))
      ident   <- outname[t]  # trajname

      # run trajec()
      info <- Trajec(yr=namelist$recp.yr[t], mon=namelist$recp.mon[t],
                     day=namelist$recp.day[t], hr=namelist$recp.hr[t],
                     lat=tmp.lat, lon=tmp.lon, agl=namelist$agl, outname=ident,
                     nhrs=namelist$nhrs, numpar=namelist$npar,veght=namelist$veght,
                     metlib=metlib, metfile=namelist$metfile, doublefiles=T,
                     delt=namelist$delt, metd=namelist$metd, mgmin=namelist$mgmin,
                     conv=namelist$convect, nummodel=namelist$nummodel,
                     overwrite=namelist$overwrite, varsout=namelist$varstrajec,
                     outpath=namelist$outpath, rundir=namelist$rundir)
    }  # end if columnTF

  }	# end t

} # end of subroutine
