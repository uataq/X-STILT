#***************************************************************************************************
# Function that loops over all starting times
#***************************************************************************************************

Trajecmod <- function(partarg=NULL, totpartarg=NULL, nodeoffset=NULL) {

#---------------------------------------------------------------------------------------------------
# Calls 'Trajec' for each starting time
# arguments assigned from call to setStiltparam.r
#
# $Id: Trajecmod.r,v 1.32 2013/12/10 17:58:31 trn Exp $
#---------------------------------------------------------------------------------------------------


# need to assign parameters; also save parameter setting in archive file with date in name
source("setStiltparam.r")
savename <- gsub(" ", ".", date())
savename <- substring(savename,4)
runs.done.dir <- NULL
if (file.exists('./Runs.done')) runs.done.dir <- './Runs.done/'
if (is.null(runs.done.dir) && file.exists(paste(sourcepath,'Runs.done',sep='')))
  runs.done.dir <- paste(sourcepath,'/Runs.done/',sep='')
if (is.null(runs.done.dir) &&
    substring(path, nchar(path)-nchar("Runs.done"), nchar(path)) == "Runs.done/")
  runs.done.dir <- sourcepath
if (!is.null(runs.done.dir)) {
   file.copy("setStiltparam.r",
             paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""),
             overwrite=T)
   cat("Saving copy of setStiltparam.r in ",
       paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""),
       "\n",sep="")
} else {
   cat("No Runs.done directory; parameter information not saved\n")
}



totpart <- 1
if (!is.null(totpartarg)) {
  cat('resetting totpart=', totpart, ' to totpartarg=', totpartarg, '\n', sep='')
  totpart <- totpartarg
}
part <- 1
if (!is.null(partarg)) {
  cat('Using totpart=', totpart, ' resetting part=', part, ' to partarg=', partarg, '\n', sep='')
  part <- partarg
}
if (!is.null(nodeoffset)) {
  nummodel <- part+nodeoffset
  cat('Using nodeoffset= ', nodeoffset, ' ,results in nummmodel= ', nummodel, '\n', sep='')
} else {
  nummodel <- part
}



# get Starting info
if (!existsr(Timesname, path)) stop(paste("cannot find object ", Timesname, " in directory ", path, sep=""))
StartInfo <- getr(paste(Timesname, sep=""), path) # object containing fractional julian day, lat, lon,
                                                  # agl, and possibly nhrsreceptor for starting position and time
# SELECTION OF A FEW Receptors for testing!
if (Times.startrow > 0) StartInfo <- StartInfo[Times.startrow:Times.endrow,, drop=FALSE] # can be just one (Times.startrow=Times.endrow)
nreceptors <- dim(StartInfo)[1]

# divide job into "totpart" parts to speed up
if (dim(StartInfo)[1] < totpart) {
  cat ('Warning: resetting totpart=', totpart, ' to dim(StartInfo)[1]=', dim(StartInfo)[1], '\n', sep='')
  totpart <- dim(StartInfo)[1]
}
if (part > totpart) {
  stop.message <- paste('Specified part=', part, ' > totpart=', totpart, ', stopping\n')
  cat(stop.message)
  stop(stop.message)
}
# round up the number of receptors per job:
start.rows <- c(1 + (0:(totpart-1))*ceiling(nreceptors/totpart))
# keep only those starting rows <= nreceptors:
start.rows <- start.rows[start.rows <= nreceptors]
# now add start.rows for final job:
start.rows <- c(start.rows,nreceptors+1)
StartInfo <- StartInfo[start.rows[part]:(start.rows[part+1]-1),, drop=FALSE]
dimnames(StartInfo)[[2]] <- toupper(dimnames(StartInfo)[[2]])

if (biomassburnTF) {
   biomassburnoutmat <- matrix(nrow=dim(StartInfo)[[1]], ncol=2)
   dimnames(biomassburnoutmat) <- list(NULL, c("ident", "CO"))
}

# OVERWRITE WARNING
if(existsr(paste("stiltresult",part,sep=""),path=path)) {
   warning("You are attempting to overwrite an existing stiltresult object")
   warning("Notice: If you have changed parameters and Trajecmod fails, first try to move or remove the existing stiltresult object")
}

nrows <- length(StartInfo[,1]) # 1 row for each trajectory
rownum <- 1
firsttraj <- T
firstflux <- T
l.remove.Trajfile <- FALSE
if (exists('remove.Trajfile')) l.remove.Trajfile <- remove.Trajfile
l.ziscale <- NULL
if (exists('ziscale')) l.ziscale <- ziscale
l.zsg.name <- NULL
if (exists('zsg.name')) l.zsg.name <- zsg.name
l.create.X0 <- FALSE
if (exists('create.X0')) l.create.X0 <- create.X0
l.use.multi <- TRUE
if (exists('use.multi')) l.use.multi <- use.multi
l.hymodelc.exe <- NULL
if (exists('hymodelc.exe')) l.hymodelc.exe <- hymodelc.exe
l.write.r <- TRUE
if (exists('write.r')) l.write.r <- write.r
l.write.nc <- FALSE
if (exists('write.nc')) l.write.nc <- write.nc
# set options for pos2id encoding (defaults are for old-style)
l.encode.minutes <- FALSE
if (exists('encode.minutes')) l.encode.minutes <- encode.minutes
l.digits.latlon <- 2
if (exists('digits.latlon')) l.digits.latlon <- digits.latlon
l.digits.agl <- 0
if (exists('digits.agl')) l.digits.agl <- digits.agl
# newer trajecmulti options:
l.MaxHymodelcDtmin <- 18*60 # no more than 18 hours discrepancy in receptor start times per hymodelc run
if (exists('MaxHymodelcDtmin')) l.MaxHymodelcDtmin <- MaxHymodelcDtmin
l.partperhymodelc <- 1
l.use.minutes <- FALSE
l.dxyp <- 0.
l.dzp <- 0.
if (l.use.multi) {
  l.setup.list <- list()
  if (exists('setup.list')) l.setup.list <- setup.list
  if (exists('use.minutes')) l.use.minutes <- use.minutes
  if (exists('partperhymodelc')) l.partperhymodelc <- partperhymodelc
  if (exists('dxyp')) l.dxyp <- dxyp
  if (exists('dzp')) l.dzp <- dzp
}
# Clean out old hymodelc.out file from previous run, if any:
if (file.exists(paste(rundir,"Copy",nummodel,"/","hymodelc.out",sep="")))
  file.remove(paste(rundir,"Copy",nummodel,"/","hymodelc.out",sep=""))
j1 <- 1
while (j1 <= nrows) {
  j2 <- min(nrows,j1+l.partperhymodelc-1)
  if (j2-j1 > 0) {
    # Only lump together receptors that are within MaxHymodelcDtmin minutes of each other
    # and that have the same value of qhrs:
    timediffs <- abs(StartInfo[j1,1]-StartInfo[(j1+1):j2,1])*24*60
    largediffs <- timediffs > l.MaxHymodelcDtmin
    if ('QHRSRECEPTOR' %in% colnames(StartInfo))
      largediffs <- largediffs | (StartInfo[(j1+1):j2,'QHRSRECEPTOR'] !=
                                  StartInfo[j1,'QHRSRECEPTOR'])
    if (any(largediffs))
      j2 <- j1 + min((1:length(largediffs))[largediffs]) - 1
  }

  lat <- StartInfo[j1:j2, "LAT"]; lon <- StartInfo[j1:j2, "LON"]; agl <- StartInfo[j1:j2, "AGL"]
  identname <- pos2id(StartInfo[j1:j2,1], lat, lon, agl, encode.minutes=l.encode.minutes,
                      digits.latlon=l.digits.latlon,digits.agl=l.digits.agl)
  for (xident in identname) cat("Trajecmod(): ", xident, " running at ", date(), "\n", sep="")
  if ('NHRSRECEPTOR' %in% colnames(StartInfo)) {
    this.nhrs <- StartInfo[j1:j2,'NHRSRECEPTOR']   # use parameter as set in receptor list
  }else{
    this.nhrs <- rep(nhrs,j2-j1+1)  # use parameter as set in setStiltparam.r
  }
  if ('QHRSRECEPTOR' %in% colnames(StartInfo)) {
    qhrs <- StartInfo[j1:j2,'QHRSRECEPTOR']   # use parameter as set in receptor list
  }else{
    qhrs <- rep(1./3600.,j2-j1+1)  # use default: 1 sec (always < dt) for instantaneous receptors
  }
  if ('DXYPRECEPTOR' %in% colnames(StartInfo)) {
    this.dxyp <- StartInfo[j1:j2,'DXYPRECEPTOR']   # use parameter as set in receptor list
  }else{
    this.dxyp <- rep(l.dxyp,j2-j1+1)  # use parameter as set in setStiltparam.r
  }
  if ('DZPRECEPTOR' %in% colnames(StartInfo)) {
    this.dzp <- StartInfo[j1:j2,'DZPRECEPTOR']   # use parameter as set in receptor list
  }else{
    this.dzp <- rep(l.dzp,j2-j1+1)  # use parameter as set in setStiltparam.r
  }
  dat <- month.day.year(floor(StartInfo[j1:j2,1])) # from julian to mmddyy
  yr4 <- dat$year # 4 digit year
  yr <- yr4%%100 # 2 digit year (or 1 digit...)
  mon <- dat$month
  day <- dat$day
  hr.float <- (StartInfo[j1:j2,1]-floor(StartInfo[j1:j2,1]))*24
  hr <- round(hr.float)
  mn <- rep(0,j2-j1+1)
  if (l.use.minutes) {
    # NOTE: mn here is now the actual minutes of the receptor time:
    hr <- floor(hr.float)
    mn <- round(60.*(hr.float-hr))
  }
  if (l.use.multi) {
    info <- Trajecmulti(yr=yr, mon=mon, day=day, hr=hr, mn=mn, lat=lat, lon=lon, agl=agl,
                        outname=identname,nhrs=this.nhrs,emisshrs=qhrs,dxyp=this.dxyp,dzp=this.dzp,
                        numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                        conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                        nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                        create.X0=l.create.X0,setup.list=l.setup.list,hymodelc.exe=l.hymodelc.exe,
                        write.r=l.write.r,write.nc=l.write.nc,
                        siguverr=siguverr,TLuverr=TLuverr,zcoruverr=zcoruverr,horcoruverr=horcoruverr,
                        sigzierr=sigzierr,TLzierr=TLzierr,horcorzierr=horcorzierr)
  } else {
    info <- Trajec(yr=yr, mon=mon, day=day, hr=hr, lat=lat, lon=lon, agl=agl, nhrs=nhrs,
                   maxdist=stepsize, numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                   conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                   nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                   create.X0=l.create.X0)
  }
  if (is.null(dim(info))) {
    infonames <- names(info)
  } else {
    infonames <- dimnames(info)[[2]]
  }
  if (firsttraj) { # set up array for run info
    run.info <- matrix(NA, nrow=nrows, ncol=length(infonames))
    dimnames(run.info) <- list(NULL, infonames)
    firsttraj <- F
  } else {
    havenames <- dimnames(run.info)[[2]]
    for (nm in infonames) {
      ndx <- which(havenames == nm)[1] # check if there are new column names
      if (is.na(ndx)) { # new column name, need to add column
        run.info <- cbind(run.info, rep(NA, dim(run.info)[1]))
        dimnames(run.info)[[2]][dim(run.info)[2]] <- nm
      }
    }
  }
  run.info[j1:j2, infonames] <- info

  #########################################################################
  ##### map trajectories to flux grids and vegetation maps ################
  ##### calculate mixing ratios at receptor points, save in result ########
  if (fluxTF) {
     if (j2>j1) stop ('non-unity partperhymodelc not supported for fluxTF=TRUE')
     print(paste("Trajecmod(): rownumber j:", j1))


     traj <- Trajecvprm(ident=identname, pathname=path, tracers=fluxtracers, coarse=aggregation,
                dmassTF=T, nhrs=nhrs, vegpath=vegpath, evilswipath=evilswipath,
                vprmconstantspath=vprmconstantspath, vprmconstantsname=vprmconstantsname, nldaspath=nldaspath,
                nldasrad=usenldasrad, nldastemp=usenldastemp, pre2004=pre2004,
                keepevimaps=keepevimaps, detailsTF=detailsTF, bios=fluxmod, landcov=landcov,
                numpix.x=numpix.x, numpix.y=numpix.y, lon.ll=lon.ll, lat.ll=lat.ll,
                lon.res=lon.res, lat.res=lat.res)


     # 'traj' is a vector
     if (existsr(paste("stiltresult", part, sep=""), path=path)) {
        result <- getr(paste("stiltresult", part, sep=""), path=path)
        if (dim(result)[1] != nrows) {
           if (firstflux) print("Trajecmod(): existing stiltresult has wrong dimension; creating new one.")
        } else {
           if (firstflux) print("Trajecmod(): found existing stiltresult, update rows in that.")
           firstflux <- FALSE
        }
     }
     if (firstflux) { # at beginning create result object
        ncols <- length(traj) # all from Trajec(), + 3 from StartInfo (agl, lat, lon)
        result <- matrix(NA, nrow=nrows, ncol=ncols)
        firstflux <- F
     }
     result[rownum, ] <- traj
     dimnames(result) <- list(NULL, c(names(traj)))
     dimnames(result) <- list(NULL, dimnames(result)[[2]])
     # write the object into default database; object names are, e.g., "Crystal.1"
     assignr(paste("stiltresult", part, sep=""), result, path=path)
  }
  rownum <- rownum+1

  ##### calculate footprint, assign in object ########
  if (footprintTF) {
     if (j2>j1) stop ('non-unity partperhymodelc not supported for footprintTF=TRUE')
     print(paste("Trajecmod(): ", identname, " running footprint at ", unix("date"), sep=""))
     print(paste("Trajecmod(): memory in use:", memory.size()[1]))
     foot <- Trajecfoot(identname, pathname=path, foottimes=foottimes, zlim=c(zbot, ztop),
                        fluxweighting=NULL, coarse=1, vegpath=vegpath,
                        numpix.x=numpix.x, numpix.y=numpix.y,
                        lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res)
     assignr(paste("foot", identname, sep=""), foot, path)
     print(paste("Trajecmod(): foot", identname, " assigned", sep=""))
  } # if (exists(foottimes))

  ##### plot footprint ########
  if (footplotTF) { # plot footprints
     if (j2>j1) stop ('non-unity partperhymodelc not supported for footplotTF=TRUE')
    foot <- getr(paste("foot", identname, sep=""), path)
    footplot(foot,identname,lon.ll,lat.ll,lon.res,lat.res)
    for (foottimespos in 1:(length(foottimes)-1)) {
    ############# NOT READY YET ###############
    }
  }

  # Specify the function parameters
  if (biomassburnTF) {
     if (j2>j1) stop ('non-unity partperhymodelc not supported for biomassburnTF=TRUE')
     biomassburnoutmat[j1, ] <- biomassburn(timesname=StartInfo, burnpath=burnpath, endpath=path, pathname=path, nhrs=nhrs, timesrow=j1)
     print(paste("Biomassburning influence calculated to be ", biomassburnoutmat[j1,2], " ppbv. Inserted into fireinfluence matrix row ", j1, sep=""))
  }

  if(l.remove.Trajfile)unix(paste("rm -f ",paste(path,".RData",identname,sep=""),sep=""))

  j1 <- j2+1 #increment row counter
}                                                           # for (j in 1:nrows)


# Wrap up all of the CO biomassburning calculations
if (biomassburnTF)
  write.table(biomassburnoutmat, file=paste(path, "fireinfluencex", nhrs, "hr_", part, ".txt", sep=""), row.names=F)

##### save mixing ratios at receptor points in file, e.g. stiltresult1.csv for part=1 ########
if (fluxTF) {
   dimnames(result) <- list(NULL, dimnames(result)[[2]])
   # write the object into default database; object names are, e.g., "Crystal.1"
   assignr(paste("stiltresult", part, sep=""), result, path=path)
   print(paste("stiltresult", part, " assigned in ", path, sep=""))
   write.table(result, file=paste(path, "stiltresult", part, ".csv", sep=""), na="", row.names=F)
}

# If evi and lswi maps from vprm calculations is saved to the global environment; it should be removed here

rm(list=objects(pattern="GlobalEvi"), envir=globalenv())
rm(list=objects(pattern="GlobalLswi"), envir=globalenv())
gc(verbose=F)

return(run.info)

}
