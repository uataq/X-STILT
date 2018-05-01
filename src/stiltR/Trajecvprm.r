#***************************************************************************************************
# calculate tracer concentrations in trajectories by mapping trajectories onto surface fluxes and
# adding CO2 boundary conditions
#***************************************************************************************************

Trajecvprm <- function(ident, pathname="", tracers=c("CO", "CO2"), coarse=1, dmassTF=T,
                       nhrs=NULL,
                       vegpath="/home/dmatross/ModelLab/VPRM/DevanVegRevise/",
                       evilswipath="/deas/group/cobra/DMMvprmTest/EVILSWI2004FULL/",
                       vprmconstantspath="/Net/Groups/BSY/BSY_3/cgerbig/RData/CarboEurope/",vprmconstantsname="vprmConstants.optCE",
                       nldaspath="/deas/group/cobra/DMMvprmTest/Radiation/NLDAS/",
                       nldasrad=FALSE, nldastemp=FALSE, pre2004=FALSE, keepevimaps=FALSE,
                       detailsTF=FALSE, linveg=FALSE,
                       numpix.x=376, numpix.y=324, lon.ll=-145, lat.ll=11,
                       lon.res=1/4, lat.res=1/6, bios="VPRM", landcov="DVN") {

# --------------------------------------------------------------------------------------------------
# Interface
# =========
#
# ident         character value specifying the trajectory ensemble to look at
# pathname      path where object with particle locations is saved
# tracers       vector of names for which mixing ratios are wanted; any subset of
#               c("co", "co2", "ch4", "h2", "n2o")
# coarse        degrade resolution (for aggregation error): 0: only 20 km resolution;
#               1-16: dynamic resolution, but limited highest resolution
#               coarse:         (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#               coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) # factors by which grids have been made coarser
#               coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)# e.g., 4 means grid is 4 times coarser
# dmassTF       if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
# vegpath       path under which vegetation and flux grids are stored
# evilswipath   path under which evi and lswi grids are stored
# vprmconstantspath     path under which vprm constants are stored (e.g. "/Net/Groups/BSY/BSY_3/cgerbig/RData/CarboEurope/")
# vprmconstantsname     name under which vprm constants are stored (e.g. "vprmconstants.12" or "vprmconstants.optCE" for EU) 
# nldaspath     path where gridded NLDAS temperature and radiation are stored
# nldasrad, nldastemp   use nldas inventory radiation and temperature instead of assimilated meteorology
# pre2004       is the year prior to 2004 or 2004+
# keepevimaps   this assigns evi and lswi maps to the global environment -- useful if this function will be called repeatedly
#               BE CAREFUL when using this, as it will put a number of LARGE objects in your database. If dynamic memory
#               allocation problems are anticipated, this feature should NOT be used!!!
# detailsTF     if TRUE, for each particle the flux contribution is saved in a big object (same name as ident, with "result" attached)
# linveg        if TRUE, CO2 fluxes are represented as linear functions of temperature and radiation (otherwise GEE is non linear)
# numpix.x      number of pixels in x directions in grid
# numpix.y      number of pixels in y directions in grid
# lon.ll        lower left corner of grid (longitude of southwest corner of southwest corner gridcell)
# lat.ll        lower left corner of grid (latitude of southwest corner of southwest corner gridcell)
# lon.res       resolution in degrees longitude
# lat.res       resolution in degrees latitude
# bios          which biosphere to use?
#               GSB -- Gerbig stupid biosphere
#               VPRM -- Devans biosphere
# landcover     which landcover to use?
#               IGBP -- Original IGBP, from Christoph
#               GLCC -- Original VPRM, based on updated IGBP (GLC2.0)
#               DVN -- The updated Devan system
#               SYNMAP Martin Jung synmap product
#               SYNMAP.VPRM8 8 VPRM classes based on synmap product
#               SYNMAP.VPRM9 9 VPRM classes based on synmap product

#global variables (see setStiltparam.r): 
# inikind    trajectory initial values. Possible: "climat", "CT" (CarbonTracker), "TM3". Vector, specifies info for each desired tracer
# inifile    absolute file name. Vector, specifies info for each desired tracer
#
#
# ---------------------------------------------------------------------------------------------------
#
# Algorithm
# =========
#
#    Step 1: Initialize & Set up output
#    Step 2: Retrieve trajectory information
#    Step 3: Set-up Grid information
#    Step 4: Accounting--mass correction and organizing particle info into only what is relevant to
#            our ground
#    Step 5: OH effects on CO
#    Step 6: Further reduce particle information to surface influence
#    Step 7: Define Emission grid wrt dynamic gridding
#    Step 8: Setup vegetation grid loop
#    Step 10: Adjust Fossil Fuel Emissions
#    Step 11: Reclassify vegetation, calculate a priori vegetation flux
#    Step 12: add CO2 boundary values
#
#
# ---------------------------------------------------------------------------------------------------
#
# History
# =======
#
# This is an attempt a rebuild Trajecflux to include the VPRM (Vegetation and Photosynthesis Model)
# into the STILT code. It represents a major departure from the GSB (Gerbig's Stupid Biosphere) originally written into
# the Trajecflux code.
#
# 1) We will now use the Xiao VPM model (Xiao et al. 2004) to get an a priori flux estimate.
#
# NEE = GEE + R
# GEE ~= GPP
# GPP = eta * FAPAR.pav * PAR
#      eta=eta.0*T.sc * W.sc * P.sc * PAR.sc, where eta.0 = 0.044 initially (according to Xiao)
#            T.sc = [(T-Tmin)(T-Tmax)]/[((T-Tmin)(T-Tmax)) - (T - Topt)^2]
#            W.sc = 1+LSWI/1+LSWImax
#           P.sc = 1+LSWI/2 {Bud Burst period} -OR- P.sc = 1 {rest of year}
#            PAR.sc = 1/(1 + PAR/PAR.0)
#      FAPAR.pav = a*EVI where a~=1 initially 
#
#     this translates to:
# GPP = lambda.par * T.sc * P.sc * W.sc * PAR.sc * EVI * PAR
#
# R = alpha*T + beta
#   
# Given an increase in the number of vegetation classes {up to 11} there will be an individual lambda parameter,
# where lambda is a scaling factor. EVI will vary over time, and we will assume a functional form with fitted 
# parameters, based on Julian time. The maximum values and minimum values will be derived from the table
# which comes from that functional form.
#
#
# 2) New capability has been added to make it possible to use Devan's vegetation classes (moniker= "DVN").
# With this, there are 11 vegetation classes, 4 evergreen and 7 others. The new vegetation scheme is as follows:
#
#                 VPRM Class    STILT-VPRM class   
# Evergreen A                1A                1
# Evergreen B                1B                2
# Evergreen C                1C                3
# Evergreen D                1D                4
# Deciduous                   2                5
# Mixed forest                3                6
# Shrubland                   4                7
# Savanna                     5                8
# Cropland-Soy               6A                9
# Cropland-Corn              6B                9
# Grassland                   7               10
# Peat                        8               11
# Others                      9               12
#
# 3) Potential global variable capability for repeated runs, eliminating memory calls
#
# 4) Now include NLDAS temperature as well as NLDAS radiation. The NLDAS radiation fields and temperature
# fields remain our best option for inputs. The STILT-based temp0 and swrad can also be used, depending
# on met input, but NLDAS represents a reanalysis version.
#
# 5) Added 'pre2004' variable. For 2003 and earlier, our radiation conversion factors are different. This is mostly
# used in converting shortwave radiation to par.
#
#---------------------------------------------------------------------------------------------------
# $Id: Trajecvprm.r,v 1.26 2009/11/25 08:51:03 gerbig Exp $
#---------------------------------------------------------------------------------------------------


####################################################################################################
# Step 1: Initialize & Set up output
####################################################################################################

tracers <- tolower(tracers)                                 # use lower case
# get names for output

tracersini <- c(paste(tracers, "ini", sep=""))
if ("co"%in%tracers)
   tracersini <- c(tracersini, "coinio")                    # also want initial CO without OH...
names.output <- c("ident", "latstart", "lonstart", "aglstart", "btime", "late", "lone", "agle",
                  tracersini, "zi", "grdht", "nendinarea", paste("sd", tracersini, sep=""))
if ("co2"%in%tracers&bios!="") #add tracer release over water when have landuse available for biosphere
   names.output <- c(names.output, "inflwater", "co2ffm")
if ("co2"%in%tracers&bios=="")
   names.output <- c(names.output, "co2ffm")
if ("co"%in%tracers)
   names.output <- c(names.output, "coffm")

output.veg <- NULL
if ("co2"%in%tracers)  {
   if (bios != "GSB" & bios != "VPRM" & bios != "")
      stop("Trajecvprm: parameter bios incorrectly specified. Redefine as GSB or VPRM.")
   if (bios == "GSB" | bios == "VPRM"){
      out.type <- c("infl", "gee", "resp")
      cat(bios, "\n")
      if (bios == "GSB") {
         output.veg <- c("frst", "shrb", "crop", "wetl")
      }
      if (bios != "GSB"&(landcov == "GLCC"|landcov == "SYNMAP"|landcov == "SYNMAP.VPRM8"))
        output.veg <- c("evergreen", "decid", "mixfrst", "shrb", "savan", "crop", "grass", "peat") #"peat" is replaced by "others" in Jena VPRM preproc.
      if (bios != "GSB"&(landcov == "DVN"))
         output.veg <- c("evergreenA", "evergreenB", "evergreenC", "evergreenD", "decid", "mixfrst", "shrb", "savan", "crop", "grass", "peat")
      names.output <- c(names.output, as.vector(outer(out.type, output.veg, paste, sep="")))
   }
}
others<-tracers[tracers!="co2"&tracers!="co"&tracers!="cofire"]
if(length(others)>0)names.output <- c(names.output, paste(others,"ffm",sep=""))
if ("cofire"%in%tracers)
   names.output <- c(names.output, "cofire")



####################################################################################################
# Step 2: Retrieve trajectory information
####################################################################################################

# Check if object exists
if (existsr(ident,pathname)) {
   cat("Trajecvprm(): starting with ident=", ident, "\n")   #found object
   part <- getr(ident,pathname)                             #get it
} else {
   if (paste(pathname,".RData",ident,".gz",sep="") %in%
      dir(pathname,pattern=paste(".RData",ident,sep=""),all.files=TRUE)) {
      cat("Trajecvprm(): starting with .gz ident=",ident, "\n") #found .gz object
      part <- getr(ident,pathname,gz=TRUE)                  #get it
   } else {
    cat("object ", pathname,".RData", ident, " NOT FOUND\n", sep="")
    lastresult <- rep(NA,length(names.output)-1)
    lastresult <- c(ident,lastresult)
    names(lastresult) <- names.output
    return(lastresult)
   }
} #if exists or not


# get time and position information from name (ident)
pos <- id2pos(ident)
time <- month.day.year(floor(pos[1]))
yr4 <- time$year                                            # 4 digit year
yr <- yr4%%100                                              # 2 digit year (or 1 digit...)
mon <- time$month
day <- time$day
hr <- round((pos[1]-floor(pos[1]))*24)

# check if required particle information is available
rqdnames <- c("time", "lat", "lon", "agl", "zi", "index", "temp0", "foot")
if ("co2"%in%tracers) {
   rqdnames <- c(rqdnames, "swrad")
   # need pressure for hybrid coordinates
   if (inikind["co2"] == "CT" || inikind["co2"] == "TM3") rqdnames <- c(rqdnames, "pres")
}
if (dmassTF)
        rqdnames <- c(rqdnames, "dmass")
for (nm in rqdnames) {
   if (!(nm%in%dimnames(part)[[2]])) {
      cat("need column '", nm, "' for this run\n", sep="")
      lastresult <- rep(NA, length(names.output)-1)
      lastresult <- c(ident, lastresult)
      names(lastresult) <- names.output
      return(lastresult)
   }
}


####################################################################################################
# Step 3: Set-up Grid information
####################################################################################################

part[, "time"] <- (-part[, "time"]/60) # use time back from now on, and transform to hours
dimnames(part)[[2]][dimnames(part)[[2]] == "time"] <- "btime"
# get grid indices
# For horizontal grids (lower left corner of south-west gridcell: 11N,145W; resolution: 1/4 lon, 1/6 lat, 376 (x) times 324 (y))
gitx <- floor(1/lon.res*(part[, "lon"]-lon.ll)+1)
gity <- floor(1/lat.res*(part[, "lat"]-lat.ll)+1)
part <- cbind(part, gitx, gity)
dimnames(part) <- list(NULL, dimnames(part)[[2]])



####################################################################################################
# Step 4: Accounting--mass correction and organizing particle info into only what is relevant to our
# ground area, as defined in step 3.
####################################################################################################

if (dmassTF) {
  #remove particles with too strong dmass violation
  ind<-unique(part[part[,"dmass"]>1E3|part[,"dmass"]<1/1E3,"index"])
  if (length(ind) >= length(unique(part[, "index"]))/2){
    message("Trajecvprm(): ", length(ind), ' of ', length(unique(part[, "index"])), ' particles have mass defect; returning NA')
    lastresult <- rep(NA, length(names.output)-1)
    lastresult <- c(ident, lastresult)
    names(lastresult) <- names.output
    return(lastresult)
  }
  part<-part[!part[,"index"]%in%ind,]

  # get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
  # i.e. correction for average mass loss of particles, since they get attracted to areas of mass destruction
  mean.dmass <- tapply(part[, "dmass"], part[, "btime"], mean) # this gives for each btime a mean dmass
  # DMM
  # To account for situations where mean.dmass is zero (mass violation total), need to avoid division by zero to
  # avoid downstream problems.
  mean.dmass[which(mean.dmass == 0)] <- 0.00001

  # need to "merge" this with part; can't use array since not all paticles are left at very large btime
  nparleft <- rle(part[, "btime"])$length # number of times the same btime is repeated
  mean.dmass <- rep(mean.dmass, nparleft) # long vector of mean dmass
  # need to link this info to each particles dmass: normalize individual dmass by mean.dmass
  part[, "dmass"] <- part[, "dmass"]/mean.dmass             # Dan Matross gets problems with that
}  
# To allow short runs, cut off after nhrs
if (!is.null(nhrs))
  part <- part[abs(part[, "btime"])<=abs(nhrs), ]

# remove particles when they cross the longitude -145 West for the first time (that's where the climatology is valid)
inbgarea <- floor(part[, "gitx"])<1
part[inbgarea, c("gitx", "gity")] <- NA
dimnames(part) <- list(NULL, dimnames(part)[[2]])

# remove points when they enter the background area for the first time
sumx <- tapply(part[, "gitx"], part[, "index"], cumsum) # cumsum gives "NA" after first occurence of "NA"
ordert <- order(part[, "index"], part[, "btime"]) # order first by index, then by time
ordern <- order(part[ordert, "btime"], part[ordert, "index"]) # to recreate original order
sumx <- unlist(sumx)[ordern]
part <- part[!is.na(sumx), ]
dimnames(part) <- list(NULL, dimnames(part)[[2]])

# only keep points, when information is changed:
# 1. position and times at boundary of NGM grid
# 2. position and times when surface influences particles
# 3. position and times at start of trajectory

# 1. boundary
ordert <- order(part[, "index"], part[, "btime"]) # order first by index, then by time
ordern <- order(part[ordert, "btime"], part[ordert, "index"]) # to recreate original order
delbte <- c(diff(part[ordert, "btime"]), -1000)[ordern] # timestep will be negative at last obs. for each particle
selend <- delbte<0

# 2. surface influence over NGM area
inngm <- floor(part[, "gitx"])<=numpix.x&floor(part[, "gitx"])>=1&floor(part[, "gity"])<=numpix.y&floor(part[, "gity"])>=1
selinf <- part[, "foot"]>0&inngm

# 3. start of particle trajectory
delbte <- c(-1000, diff(part[ordert, "btime"]))[ordern] # timestep will be negative at first obs. for each particle
selfirst <- delbte<0



####################################################################################################
# Step 5: OH effects on CO
####################################################################################################

if ("co"%in%tracers) {
   #################### CO + OH losses########################################
   # first get OH at particle position, calculate COrel. loss and CH4 source
   # OH from SAS, parameterized as oh=oh0+oh1*p+oh2*p**2, with parameters ohi for each month and for 30 and 60 lat
   pmb <- 1013*exp(-part[, "agl"]/8000) # 8 km scale height
   ohm <- oh[oh[, "month"] == mon, ]
   oh60 <- ohm[ohm[, "lat"] == 60, "oh0"]+ohm[ohm[, "lat"] == 60, "oh1"]*pmb+ohm[ohm[, "lat"] == 60, "oh2"]*pmb**2 # OH at 60 north, same altitude
   oh30 <- ohm[ohm[, "lat"] == 30, "oh0"]+ohm[ohm[, "lat"] == 30, "oh1"]*pmb+ohm[ohm[, "lat"] == 30, "oh2"]*pmb**2 # OH at 30 north, same altitude
   oh60[oh60<0] <- 0; oh30[oh30<0] <- 0
   # local OH, interpolated
   ohl <- (oh60*(part[, "lat"]-30)+oh30*(60-part[, "lat"]))/(60-30)
   ohl[ohl<0] <- 0
   tair.k <- part[, "temp0"]-part[, "agl"]*6.5/1000 # average lapse rate
   # get CO loss
   kohCO <- 1.3E-13*(1+(0.6*pmb/1000))*300/tair.k # in 1/cm^3/sec
   tau <- 1/(kohCO*ohl)/3600/24; # in days
   # integral (k*OH*dt)
   # get dbtime in same format as kohCO and ohl
   delbt <- c(0, diff(part[ordert, "btime"]))[ordern][!selfirst]*60*60 # timestep in seconds
   kohdt <- kohCO[!selfirst]*ohl[!selfirst]*delbt
   # integral over dt: sum over timesteps for each particle individually
   Ikohdt <- tapply(kohdt, part[!selfirst, "index"], sum)
   CO.frac <- exp(-Ikohdt) # preliminary result: factors for each particles initial CO
   CO.fact <- part[, "btime"]*0 # initialize
   CO.fact[selend] <- CO.frac

   #################### CO from CH4 + OH ######################################
   kohCH4 <- 2.3E-12*exp(-1765/tair.k);
   kohCH4dt <- kohCH4[!selfirst]*ohl[!selfirst]*delbt
   IkohCH4dt <- tapply(kohCH4dt, part[!selfirst, "index"], sum)
   COfrCH4 <- 1.780*IkohCH4dt # CO from CH4 in ppm
   # not all will make it, account for some losses (about half the losses for COini)
   COfrCH4 <- COfrCH4*exp(-0.5*Ikohdt)
   CO.frCH4 <- part[, "btime"]*0 # initialize
   CO.frCH4[selend] <- COfrCH4
   part <- cbind(part, CO.fact, CO.frCH4) # add both to particle location object
}



####################################################################################################
# Step 6: Further reduce particle information to surface influence
####################################################################################################

# keep only particles where they matter: create flag for first, last, or when surface influence
# apply selection here (only first, last, or when surface influence)
part <- part[(selinf+selend+selfirst)>0, ]
# print("3"); print(dim(part)); print(sum(inngm)); print(sum(part[, "foot"]>0))
# move x and y position of final position to initialization area (gitx=1, gity= 1 to numpix.y), at least make sure they are not outside NGM
part[part[, "gitx"]>numpix.x, "gitx"] <- numpix.x
part[part[, "gity"]>numpix.y, "gity"] <- numpix.y
part[part[, "gitx"]<1, "gitx"] <- 1
part[part[, "gity"]<1, "gity"] <- 1

# get different resolutions for surface grids depending on range in x and y and on particle number for each timestep
# get selector for first and last row w/ a given btime
selfirst <- c(T, diff(part[, "btime"])>0)
selast <- c(diff(part[, "btime"])>0, T)

max.x <- part[order(part[, "btime"], part[, "gitx"]), ][selast>0, "gitx"]
min.x <- part[order(part[, "btime"], part[, "gitx"]), ][selfirst>0, "gitx"]
max.y <- part[order(part[, "btime"], part[, "gity"]), ][selast>0, "gity"]
min.y <- part[order(part[, "btime"], part[, "gity"]), ][selfirst>0, "gity"]
btime <- part[order(part[, "btime"], part[, "gity"]), ][selfirst>0, "btime"]
names(max.x) <- NULL; names(min.x) <- NULL; names(max.y) <- NULL; names(min.y) <- NULL

# now get information back in format for all timesteps and index-numbers
minmax.yx <- cbind(btime, max.x, min.x, max.y, min.y)
minmax.yx <- merge(part[, c("btime", "index")], minmax.yx, by="btime")
max.x <- minmax.yx[, "max.x"]
min.x <- minmax.yx[, "min.x"]
max.y <- minmax.yx[, "max.y"]
min.y <- minmax.yx[, "min.y"]
names(max.x) <- NULL; names(min.x) <- NULL; names(max.y) <- NULL; names(min.y) <- NULL



####################################################################################################
# Step 7: Define Emission grid wrt dynamic gridding
####################################################################################################

# Call 'getgrid' to get correct emission grid--necessary b/c emission grid is too large, so divided into several diff objects
# use getgridp.ssc function: don't allow the resolution to get finer at earlier backtime; use cummax(ran.x)

gridresult <- getgridp(min.x, max.x, min.y, max.y, numpix.x, numpix.y, coarse)
emissname <- paste(gridresult[, "xpart"], gridresult[, "ypart"], gridresult[, "gridname"], sep="")
# Extract appropriate emissions within each emission grid--do one grid at a time b/c reduces # of times grid has to be accessed
coarsex <- c(1,1,2,2,2,4,4,4,8,8,8,16,16,16,32,32)  # factors by which grids have been made coarser
coarsey <- c(1,2,1,2,4,2,4,8,4,8,16,8,16,32,16,32)  # e.g., '4' means grid is 4 times coarser



####################################################################################################
# Step 8: Setup vegetation grid loop
####################################################################################################

if (landcov == "IGBP") {
        veghead <- "veg."
} else if (landcov == "GLCC") {
        veghead <- "glcc."
} else if (landcov == "DVN") {
        veghead <- "devanveg."
} else if (landcov == "SYNMAP") {
        veghead <- "synmap."
} else if (landcov == "SYNMAP.VPRM8") {
        veghead <- "synvprm8."
} else {
        stop("Improperly specified Landcover format; Exiting now!")
        veghead <- "veg."
}
if(!"co2"%in%tracers)landcov <- "" # use dummy, so that columns vor veg. cov. don't get created
nBaseVeg <- length(output.veg)
if (landcov == "GLCC"|landcov == "IGBP") nBaseVeg <- 17
if (landcov == "SYNMAP") nBaseVeg <- 48
if (landcov == "SYNMAP.VPRM8") nBaseVeg <- 8
if (landcov == "DVN") nBaseVeg <- 12
if (bios == "") nBaseVeg <- 0
nReclss <- ifelse(bios == "GSB",5,8)
if (landcov == "DVN")
        nReclss <- 12
pth2o <- nBaseVeg            # water (assumed to be 'last' class)
ptco2ff <- nBaseVeg+1        # co2 from fossil fuels
ptcoff <- nBaseVeg+2         # co from fossil fuels
if(!ncdfTF["ch4"]){ #surface fluxes and modis indices not in netCDF format
  ptch4 <- (nBaseVeg+3)+0:11   # methane emission for Jan-Dec, monthly
  ptn2o <- (nBaseVeg+3)+12   # assume constant fluxes (they are not there yet)
  pth2 <- (nBaseVeg+3)+13   # assume constant fluxes (they are not there yet)
  ptcofire <- (nBaseVeg+3)+14   #fire emissions ncdf
} else {
  ptch4 <- (nBaseVeg+3)   # hourly methane emission in ncdf file
  ptn2o <- (nBaseVeg+4)   # hourly n2o emission in ncdf file
  pth2 <- (nBaseVeg+5)   # hourly h2 emission in ncdf file (don't exist yet)
  ptcofire <- (nBaseVeg+6)   # fire emissions ncdf
}

# loop over different surface grids
vegs <- NULL # initialize vector containing all flux grid numbers

if ("co2"%in%tracers&bios!="")
        vegs <- c(vegs,1:nBaseVeg, ptco2ff)
if ("co2"%in%tracers&bios=="")
        vegs <- c(vegs, ptco2ff)
if ("co"%in%tracers)
        vegs <- c(vegs, ptcoff)
if ("ch4"%in%tracers)
        vegs <- c(vegs, ptch4)
if ("n2o"%in%tracers)
        vegs <- c(vegs, ptn2o)
if ("cofire"%in%tracers)
        vegs <- c(vegs, ptcofire)
nveg <- length(vegs)

result <- cbind(part, matrix(ncol=nveg, nrow=length(part[, "btime"])))
dimnames(result) <- list(NULL, c(dimnames(part)[[2]], paste("v", vegs, sep="")))



# DMM lswi and evi
# This handles dynamic gridding
# Depends on formatting of EVI and LSWI to be in 376x324x72 array (or degraded to YxXx72)
if (("co2"%in%tracers & bios == "VPRM"&!ncdfTF["co2"])|nldasrad|nldastemp) {
   emssname.msbb <- gridresult[, "gridname"]

   if (bios == "VPRM") {
      eviset <- matrix(0, nrow=nrow(part), ncol=nReclss)
      lswiset <- matrix(0, nrow=nrow(part), ncol=nReclss)

      eviMaxVec <- matrix(1, nrow=nrow(part), ncol=nReclss)
      eviMinVec <- matrix(0, nrow=nrow(part), ncol=nReclss)
      lswiMaxVec <- matrix(1, nrow=nrow(part), ncol=nReclss)
      lswiMinVec <- matrix(-1, nrow=nrow(part), ncol=nReclss)
   }

   if (nldasrad) {
      nldasRadVec <- rep(NA, nrow(part))
   }

   if (nldastemp) {
      nldasTempVec <- rep(NA, nrow(part))
   }

   rTime <- month.day.year(floor(id2pos(ident)[1]))         # Month day year
   rDayFrac <- id2pos(ident)[1] - floor(id2pos(ident)[1])   # Fractional day--accounts for starting hour
   julDay <- rTime$day+switch(rTime$month,0,31,59,90,120,151,181,212,243,273,304,334)+
                           ifelse(rTime$year%%4 == 0&rTime$year != 2000&rTime$month>2,1,0) # leap year

   if (nldasrad|nldastemp) {
      nlddoy <- floor(julDay + rDayFrac - part[, "btime"]/24)
      nldstartday <- min(nlddoy, na.rm=T)
      nldendday <- max(nlddoy, na.rm=T)
      nldMaster <- array(NA, dim=c(numpix.y, numpix.x, length(nldstartday:nldendday)*24))
      nldTempMaster <- array(NA, dim=c(numpix.y, numpix.x, length(nldstartday:nldendday)*24))
      nlddayvec <- nldstartday:nldendday
      for (g1 in 1:length(nlddayvec)) {
         for (g2 in 0:23) {
            nldmon <- which.min(nlddayvec[g1]%/%(c(31,59,90,120,151,181,212,243,273,304,334,365)+
                            c(0, rep(ifelse(rTime$year%%4 == 0&rTime$year != 2000,1,0),11))))
            nldday <- nlddayvec[g1]%%((c(365,31,59,90,120,151,181,212,243,273,304,334)+
                            c(0,0, rep(ifelse(rTime$year%%4 == 0&rTime$year != 2000,1,0),10)))[nldmon])
            if (nldday == 0) {
                    nldmon <- nldmon-1
                    nldday <- c(31,28+ifelse(rTime$year%%4 == 0&rTime$year != 2000,1,0),
                            31,30,31,30,31,31,30,31,30,31)[nldmon]
            }
            nldyear <- rTime$year
            nldfilen <- paste("nldassw.", nldyear, ifelse(nchar(as.character(nldmon)) == 1, "0", ""), nldmon,
                            ifelse(nchar(as.character(nldday)) == 1, "0", ""), nldday,
                            ifelse(nchar(as.character(g2)) == 1, "0", ""), g2, sep="")

            nldTempfilen <- paste("nldastempK.", nldyear, ifelse(nchar(as.character(nldmon)) == 1, "0", ""), nldmon,
                            ifelse(nchar(as.character(nldday)) == 1, "0", ""), nldday,
                            ifelse(nchar(as.character(g2)) == 1, "0", ""), g2, sep="")
            if (nldasrad) {
                    nldhrmat <- getr(nldfilen, path=nldaspath)
                    nldMaster[, , (g1-1)*24+g2] <- nldhrmat
            }
            if (nldastemp) {
                    nldTemphrmat <- getr(nldTempfilen, path=nldaspath)
                    nldTempMaster[, , (g1-1)*24+g2] <- nldTemphrmat
            }# if (nldastemp)
         } # for g2
      } # for g1

   } # if nldas

   for (rs in unique(emssname.msbb)) {
      sel <- gridresult[, "gridname"] == rs
      doy <- floor(julDay+rDayFrac-part[sel, "btime"]/24)
      decday <- julDay+rDayFrac-part[sel, "btime"]/24
      gridname <- unique(gridresult[sel, "gridname"])   # gridname can be 1~16, representing diff. resolutions of grid
      x <- part[sel, "gitx"]
      y <- part[sel, "gity"]

      xbestres <- x
      ybestres <- y
      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridname]; shrink.y <- coarsey[gridname]
      # grids have NOT been divided
      x <- ceiling(x/shrink.x)
      y <- ceiling(y/shrink.y)

      if (nldasrad|nldastemp) {
         # NLDAS master is an array dynx, dyny, nday
         nldtmmark <- (floor(decday-nldstartday)*24+floor(((decday-nldstartday) - floor(decday-nldstartday))*24))+1
         nldtmmark[which(nldtmmark>dim(nldMaster)[3])] <- NA
         nldtmmark[which(nldtmmark<0)] <- NA
      }
      if (nldasrad) {
         nldasRadVec[sel] <- nldMaster[cbind(ybestres, xbestres, nldtmmark)]
      }

      if (nldastemp) {
         nldasTempVec[sel] <- nldTempMaster[cbind(ybestres, xbestres, nldtmmark)]
      }




      if (bios == "VPRM") {
      # Changed--to do by vegetation class
      # Now have to pick it up down below
      ## Need to add max and min values for evi and lswi from mapped values
          if (exists(paste("GlobalEviMaxMap.res", rs, sep=""), where=globalenv())) {
                 eviMaxMap <- get(paste("GlobalEviMaxMap.res", rs, sep=""), envir=globalenv())
         } else {
                 eviMaxMap <- getr(paste("eviMaxMap.res", rs, sep=""), path=evilswipath)
                 if (keepevimaps)
                         assign(paste("GlobalEviMaxMap.res", rs, sep=""), eviMaxMap, envir=globalenv())
         }

         if (exists(paste("GlobalEviMinMap.res", rs, sep=""), where=globalenv())) {
                 eviMinMap <- get(paste("GlobalEviMinMap.res", rs, sep=""), envir=globalenv())
         } else {
                 eviMinMap <- getr(paste("eviMinMap.res", rs, sep=""), path=evilswipath)
                 if (keepevimaps)
                         assign(paste("GlobalEviMinMap.res", rs, sep=""), eviMinMap, envir=globalenv())
         }

         if (exists(paste("GlobalLswiMaxMap.res", rs, sep=""), where=globalenv())) {
                 lswiMaxMap <- get(paste("GlobalLswiMaxMap.res", rs, sep=""), envir=globalenv())
         } else {
                 lswiMaxMap <- getr(paste("lswiMaxMap.res", rs, sep=""), path=evilswipath)
                 if (keepevimaps)
                         assign(paste("GlobalLswiMaxMap.res", rs, sep=""), lswiMaxMap, envir=globalenv())
         }

         if (exists(paste("GlobalLswiMinMap.res", rs, sep=""), where=globalenv())) {
                 lswiMinMap <- get(paste("GlobalLswiMinMap.res", rs, sep=""), envir=globalenv())
         } else {
                 lswiMinMap <- getr(paste("lswiMinMap.res", rs, sep=""), path=evilswipath)
                 if (keepevimaps)
                         assign(paste("GlobalLswiMinMap.res", rs, sep=""), lswiMinMap, envir=globalenv())
         }


      # UPDATE UPDATE--back to linear interpolation of smoothed, approx fitted 8-day curves. The 273 day matrix was
      # too unwieldy for bulk calculations

         setDys <- seq(1,361, by=8)
         lateMrk <- ((doy-1)%/%8 + 1)*8+1
         earlyMrk <- ((doy-1)%/%8)*8+1
         interpFrac <- (doy-earlyMrk)/8
         earlyMtch <- match(earlyMrk, setDys)
         lateMtch <- match(lateMrk, setDys)


         for (vt1 in 1:nReclss) {
            if (exists(paste("GlobalEviMaster.", vt1, ".res", rs, sep=""), where=globalenv())) {
                    eviMaster <- get(paste("GlobalEviMaster.", vt1, ".res", rs, sep=""), envir=globalenv())
            } else {
                    eviMaster <- getr(paste("eviMaster.", vt1, ".res", rs, sep=""), path=evilswipath)
                    if (keepevimaps)
                            assign(paste("GlobalEviMaster.", vt1, ".res", rs, sep=""), eviMaster, envir=globalenv())
            }

            if (exists(paste("GlobalLswiMaster.", vt1, ".res", rs, sep=""), where=globalenv())) {
                    lswiMaster <- get(paste("GlobalLswiMaster.", vt1, ".res", rs, sep=""), envir=globalenv())
            } else {
                    lswiMaster <- getr(paste("lswiMaster.", vt1, ".res", rs, sep=""), path=evilswipath)
                    if (keepevimaps)
                            assign(paste("GlobalLswiMaster.", vt1, ".res", rs, sep=""), lswiMaster, envir=globalenv())
            }

            vtvec <- rep(vt1, length(x))
            eviEarly <- eviMaster[cbind(y, x, earlyMtch)]
            eviLate <- eviMaster[cbind(y, x, lateMtch)]
            lswiEarly <- lswiMaster[cbind(y, x, earlyMtch)]
            lswiLate <- lswiMaster[cbind(y, x, lateMtch)]

            eviLate[which(lateMrk>313|lateMrk<41)] <- 0
            eviEarly[which(earlyMrk>313|earlyMrk<41)] <- 0
            lswiLate[which(lateMrk>313|lateMrk<41)] <- 0
            lswiEarly[which(earlyMrk>313|earlyMrk<41)] <- 0

            eviset[sel, vt1] <- (1-interpFrac)*eviEarly+interpFrac*eviLate
            lswiset[sel, vt1] <- (1-interpFrac)*lswiEarly+interpFrac*lswiLate

            # Recomplication--we're back to the maps

            eviMaxVec[sel, vt1] <- eviMaxMap[cbind(y, x, vtvec)]
            eviMinVec[sel, vt1] <- eviMinMap[cbind(y, x, vtvec)]
            lswiMaxVec[sel, vt1] <- lswiMaxMap[cbind(y, x, vtvec)]
            lswiMinVec[sel, vt1] <- lswiMinMap[cbind(y, x, vtvec)]

         } # for vt1
      } # of if (bios == "VPRM")
   } # for rs in emssnames
} # of if ((bios == "VPRM"&!ncdfTF)|nldasrad|nldastemp)        

# Note in the following that emissgrid is a misnomer, these are really vegetation fractions. Only emissions for
# fossil fuels
# Following is legacy code setupt. Loops over vegetation types at different resolutions. Should be rolled into the above loop, but hasn't been.
for (vegnum in vegs) {
   speci <- NULL
   if (vegnum <= nBaseVeg & ncdfTF["co2"]) next # read veg frac and/or modis indices from netcdf files
   if (vegnum == ptco2ff) speci <- "CO2"
   if (vegnum == ptcoff) speci <- "CO"
   if (any(vegnum == ptch4)) speci <- "CH4"
   if (vegnum == ptn2o) speci <- "N2O"
   if (vegnum == ptcofire) speci <- "cofire"
   if (!is.null(speci) && ncdfTF[tolower(speci)]) next # read surface fluxes from netcdf files

   for (rs in unique(gridresult[, "gridname"])) {
      # emissgrid <- getr(paste(veghead, vegnum, ".", name, sep=""), path=vegpath)  # select appropriate emission grid--note extension '.mat' for new format
      emissgrid <- getr(paste(veghead, vegnum, ".00", rs, sep=""), path=vegpath)
      if (vegnum == pth2o&!("co2"%in%tracers))
              emissgrid <- 1-emissgrid # use land mask for CO emission w/o prior
      sel <- gridresult[, "gridname"] == rs
      gridname <- unique(gridresult[sel, "gridname"])   # gridname can be 1~16, representing diff. resolutions of grid
      x <- part[sel, "gitx"]
      y <- part[sel, "gity"]

      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridname]; shrink.y <- coarsey[gridname]
      # NO DIVIDING
      # grids have NOT been divided
      x <- ceiling(x/shrink.x)
      y <- ceiling(y/shrink.y)

      yx <- cbind(y, x)
      # BUDGET
      # Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes,
      # in ppm/(micro-mol/m^2/s)
      EMCO <- emissgrid[yx]*part[sel, "foot"]
      # also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
      if (dmassTF) EMCO <- EMCO*part[sel, "dmass"]
      result[sel, paste("v", vegnum, sep="")] <- EMCO
   }# for rs
}# for different surface grid
for (vegnum in vegs[(nBaseVeg + 1):length(vegs)]) { # all fossil fuel emissions, do the ones with netCDF format
   speci <- NULL
   if (vegnum == ptco2ff) speci <- "CO2"
   if (vegnum == ptcoff) speci <- "CO"
   if (any(vegnum == ptch4)) speci <- "CH4"
   if (vegnum == ptn2o) speci <- "N2O"
   if (vegnum == ptcofire) speci <- "cofire"
   if (is.null(speci)) stop ('Failed to match speci based on vegnum')
   if (ncdfTF[tolower(speci)]) { # read surface fluxes from netcdf files
        # BUDGET
              # Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes,
              # in ppm/(micro-mol/m^2/s)
      # need time in YY MM DD hh mm, not local, but GMT!
      gmtime <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=rep(0, dim(result)[1])) # last column is weekday, sunday=0, monday=1 etc.
      fdate <- as.vector(rbind(gmtime[, "yr"], gmtime[, "mon"], gmtime[, "day"], gmtime[, "hr"], rep(0, dim(result)[1])))
      if (sum(is.na(fdate))>0) stop("NA in fdate!")
      # print(fdate)
      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridresult[, "gridname"]]; shrink.y <- coarsey[gridresult[, "gridname"]]
      x <- ceiling(part[, "gitx"]/shrink.x)
      y <- ceiling(part[, "gity"]/shrink.y)
      emiss <- NULL
      # print(pathname)
      # assignr("testncdf", list(gmtime=gmtime, fdate=fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, emiss=emiss), pathname)
      # stop("testing")
      # print("date 1")
      # print(date())
  
 #    if(speci!="cofire")emiss <- get.fossEU.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci)
  #    if(speci=="cofire")emiss <- get.fireBARCA.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci)
      print("DEBUG Now call get.Emisnetcdf")
      emiss <- get.Emis.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci,numpix_x=numpix.x,numpix_y=numpix.y)

      EMCO <- emiss*part[, "foot"]
      # also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
      if (dmassTF) EMCO <- EMCO*part[, "dmass"]
      result[, paste("v", vegnum, sep="")] <- EMCO
   } # of if (ncdf) {
} # for different surface grid

if (ncdfTF["co2"]&bios=="VPRM"&"co2"%in%tracers) { # VPRM with fortran call to directly access netcdf files
        # BUDGET
              # Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes,
              # in ppm/(micro-mol/m^2/s)
      # need time in YY MM DD hh mm, not local, but GMT!
      gmtime <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=rep(0, dim(result)[1])) # last column is weekday, sunday=0, monday=1 etc.
      fdate <- as.vector(rbind(gmtime[, "yr"], gmtime[, "mon"], gmtime[, "day"], gmtime[, "hr"], rep(0, dim(result)[1])))
      if (sum(is.na(fdate))>0) stop("NA in fdate!")
      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridresult[, "gridname"]]; shrink.y <- coarsey[gridresult[, "gridname"]]
      x <- ceiling(part[, "gitx"]/shrink.x)
      y <- ceiling(part[, "gity"]/shrink.y)
      if(bios=="VPRM"){
        vprmset <- NULL
        # get EVI parameters, vegetation specific
        vprmset  <- get.modis.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, dpath=evilswipath)
        eviset <- vprmset$EVI
        lswiset <- vprmset$LSWI
        eviMaxVec <- vprmset$EVI_amax
        eviMinVec <- vprmset$EVI_amin
        lswiMaxVec <- vprmset$LSWI_amax
        lswiMinVec <- vprmset$LSWI_amin
        veg.fra<-vprmset$VEG_FRA
      } else {
        veg.fra <- get.vegfrac.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, dpath=evilswipath)

      }
      infl.veg <- veg.fra*part[, "foot"] #influence from specific vegetation type
      # also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
      if (dmassTF) infl.veg <- infl.veg*part[, "dmass"]
      result[, paste("v", vegs[1:nBaseVeg], sep="")] <- infl.veg
} # of if (ncdf.fossTF) {

dimnames(result) <- list(NULL, dimnames(result)[[2]])

###################################


####################################################################################################
# Step 10: Adjust Fossil Fuel Emissions
####################################################################################################

   # time factor for CO emissions
   ltime <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=round(result[, "lon"]*24/360)) # last column is weekday, sunday=0, monday=1 etc.
           # hourly factors from NAPAP
   hfac <- c(0.272,0.231,0.214,0.226,0.322,0.705,1.22,1.39,1.35,1.40,1.49,1.54,1.58,1.59,1.65,1.74,1.69,1.43,1.05,0.789,0.688,0.592,0.483,0.370)
           # weekday-factors: sun, mon, tue, wed, thu, fri, sat, sun; napap-grid contains fluxes for summer-weekday, therefore weekday factor of 1
   dfac <- c(0.865,1,1,1,1,1,0.875)

   if ("co"%in%tracers&!ncdfTF["co"]) { #not for netCDF files with hourly emissions
           emfacco <- hfac[ltime[, "hr"]+1]*dfac[ltime[, "weekd"]+1]
   # DMM          result[, "v19"] <- result[, "v19"]*emfacco
           result[, paste("v", ptcoff, sep="")] <- result[, paste("v", ptcoff, sep="")]*emfacco
   #          if (!("co2"%in%tracers)) result[, "v17"] <- result[, "v17"]*emfacco # land mask weighted w/ CO emissions
   }

   # time factor for CO2 emissions
   if ("co2"%in%tracers&!ncdfTF["co2"]) { #not for netCDF files with hourly emissions
             # use hourly factors from CO, but with reduced amplitude (0.4 reduction)
             hfac <- 1+0.4*(hfac-1)
             # weekday-factors: use factors from CO, but w/ reduced amplitude, and with mean == 1 (CO2 emissions are for average day, not for weekday)
             dfac <- dfac/mean(dfac)
             dfac <- 1+0.4*(dfac-1)
             emfacco2 <- hfac[ltime[, "hr"]+1]*dfac[ltime[, "weekd"]+1]
   # DMM                  result[, "v18"] <- result[, "v18"]*emfacco2
             result[, paste("v", ptco2ff, sep="")] <- result[, paste("v", ptco2ff, sep="")]*emfacco2
   }



####################################################################################################
# Step 11: Reclassify vegetation, calculate a priori vegetation flux
####################################################################################################

if (landcov == "IGBP"&bios == "GSB") {
   nReclss <- 5
   reclss <- c(1,1,1,1,1,2,2,2,2,2,4,3,0,3,5,2,5)

   # GSB is based on grouping IGBP classes
   # -------------------------------------
   # vgroup1: Forrests: IGBP 1, 5, 4, 2, 3 (not existent)
   # vgroup2: shrublands etc.:7,10,8,16,6,9
   # vgroup3: Crops etc.: 14,12
   # vgroup4: wetland: 11
   # vgroup5: Water IGBP 17 (water), 15 (Snow and Ice)
   # Rest:13 (urban) -- are assumed to have zero vegetative flux influence
}

if (landcov == "GLCC"&bios != "GSB") {
   nReclss <- 8
   reclss <- c(1,1,2,2,3,4,4,5,5,7,7,6,8,7,8,7,8)
   # CHECK THIS
   # NEEDS WORK TO DEFINE VEGETATION GROUPS
   ## Value        Devan Category
   # 1        Evergreen Forest
   # 2        Deciduous Forest
   # 3        Mixed Forest
   # 4        Shrublands
   # 5        Savannas
   # 6        Croplands
   # 7        Grasslands
   # 8        Water, Snow/ice, other
   # Rest(urban, barren) assumed to have zero vegetative flux influence
}

if (landcov == "SYNMAP.VPRM8"&bios == "GSB") {
   nReclss <- 5
   reclss <- c(1,1,1,2,2,3,2,5)
}
if (landcov == "SYNMAP.VPRM8"&bios == "VPRM") {
   nReclss <- 9
   reclss <- 1:9
}
# Devan vegetation classes, updated
#                 VPRM Class    STILT-VPRM class
# Evergreen A                1A                1
# Evergreen B                1B                2
# Evergreen C                1C                3
# Evergreen D                1D                4
# Deciduous                2                5
# Mixed forest                3                6
# Shrubland                4                7
# Savanna                5                8
# Cropland-Soy                6A                9
# Cropland-Maize                6B                9
# Grassland                7                10
# Peatlands                8                11
# Others                        9                12
#

if (landcov == "DVN") {
   nReclss <- 12
   reclss <- 1:12
}

result.ready <- result[, substring(dimnames(result)[[2]],1,1) != "v"]

if ("co2"%in%tracers&bios!="") {
   result.add <- cbind(result[, paste("v", pth2o, sep="")], result[, paste("v", ptco2ff, sep="")])
   dimnames(result.add) <- list(NULL, paste("v", nBaseVeg+0:1, sep=""))
   result.ready <- cbind(result.ready, result.add)
}
if ("co2"%in%tracers&bios=="") {
   result.add <- matrix(result[, paste("v", ptco2ff, sep="")], ncol=1, dimnames=list(NULL, paste("v", nBaseVeg+2, sep="")))
   result.ready <- cbind(result.ready, result.add)
}
if ("co"%in%tracers) {
   result.add <- matrix(result[, paste("v", ptcoff, sep="")], ncol=1, dimnames=list(NULL, paste("v", nBaseVeg+2, sep="")))
   result.ready <- cbind(result.ready, result.add)
}

if ("co2"%in%tracers&bios!="") {
   # Conversion factor of SW radiation (W/m2) to par for use in VPRM. Based on accumulated fit of tower data
   # In NLDAS, EDAS gets ingested to supplement GOES at low sun angles. EDAS is an overestimation,
   if (pre2004) {
      swtoparnldas <- 1.60
      swtoparedas <- 1.5 # placeholder
   } else {
      swtoparnldas <- 1.89
      swtoparedas <- 1.61
   }
        
   # Get light and temperature
   if (nldasrad) {
      swrad <- nldasRadVec
      swradpar <- swtoparnldas*swrad
      radsubpts <- which(is.na(swrad))
      swrad[radsubpts] <- result[radsubpts, "swrad"]
      swradpar[radsubpts] <- swtoparedas*result[radsubpts, "swrad"]


      # Note: We fill in the missing points in radiation with modeled values
   } else {
      swrad <- result[, "swrad"]
      swradpar <- swtoparedas*swrad
   }

   if (nldastemp) {
      tempAir <- nldasTempVec-273.15
      tempAir[which(tempAir<(-100)|tempAir>(100)|is.nan(tempAir))] <- NA
      tempAir[which(is.na(tempAir))] <- (result[which(is.na(tempAir)), "temp0"]-273.15)
   } else {
      tempAir <- result[, "temp0"]-273.15
      tempAir[which(tempAir<(-100)|tempAir>(100)|is.nan(tempAir))] <- NA
   }
                
   if (bios == "GSB") {
#! GSB functionality has only a few land cover types
#                if (landcov == "IGBP") {
      if (!linveg) {
         dresp.dT <- c(dlambda.veg[, "drdt"], dlambda.veg[2, "drdt"])
         a3 <- c(dlambda.veg[, "a3"], dlambda.veg[2, "a3"])
         a4 <- c(dlambda.veg[, "a4"], dlambda.veg[2, "a4"])
         # get parameterized GEE and R
      } else { # fluxes linear in temp and radiation
         dresp.dT <- c(dlambda.simp.veg[, "drdt"], dlambda.simp.veg[2, "drdt"])
         a3 <- c(dlambda.simp.veg[, "a3"], dlambda.simp.veg[2, "a3"])
      }# of if not linear fluxes as fct of temp and radiation
#                } # if IGBP

      # calculate t and par dependent fluxes
      for (k in 1:(nReclss-1)) {
         # Assume last reclass category is water
         # Determine total influence for all vegetation from this uberclass
         # v.tot <- apply(result[, paste("v", as.character(which(reclss == k)), sep="")],1, sum)
         selCols <- paste("v", as.character(which(reclss == k)), sep="")
         if (length(selCols)>1)
                 v.tot <- apply(result[, selCols],1, sum)
         else if (nchar(selCols)>1)
                 v.tot <- result[, selCols]
         else # not all classes are present in all areas
#                 stop("Sel Cols missing--Veg types off")
                 v.tot <- result[,1]*NA
         # ENDFIX
         if (!linveg)
                 v.gee <- a3[k]*swrad/(swrad+a4[k])*v.tot
         if (linveg)
                 v.gee <- a3[k]*swrad*v.tot

         v.resp <- dresp.dT[k]*tempAir*v.tot

         cnameAdd <- paste("v", (nBaseVeg+3)+(k-1)*3+0:2, sep="")
         gsbMat <- cbind(v.tot, v.gee, v.resp)
         dimnames(gsbMat) <- list(NULL, cnameAdd)

         result.ready <- cbind(result.ready, gsbMat)
      }
   }

   if (bios == "VPRM") {
      # read in values for all classes
      if (landcov != "DVN") {
              vprmConstants <- getr(vprmconstantsname, path=vprmconstantspath)
      } else {

              vprmConstants <- getr("publishedvprmconstants.12", path=vprmconstantspath)
      }


      for (k in 1:(nReclss-1)) {

         # Determine total influence for all vegetation in this uberclass
         selCols <- paste("v", as.character(which(reclss == k)), sep="")
         if (length(selCols)>1)
                 v.tot <- apply(result[, selCols],1, sum)
         else if (length(selCols)>0)
                 v.tot <- result[, selCols]
         else
                 stop("Sel Cols missing--Veg types off")
         # get parameters
         if(names(vprmConstants)[1]=="lambdaGPP"){ #devan style 
            lambdaGPP <- vprmConstants$lambdaGPP[k]
            alphaResp <- vprmConstants$alphaResp[k]
            betaResp <- vprmConstants$betaResp[k]
            parZero <- vprmConstants$parZero[k]
            radScalar <- swradpar/(1 + (swradpar/parZero)) # Include radiation factor
         } else { #Jena style, differentiate between PAR and SW, and use SW
            lambdaGPP <- vprmConstants$lambdaGPP.sw[k]
            alphaResp <- vprmConstants$alphaResp[k]
            betaResp <- vprmConstants$intResp[k]
            parZero <- vprmConstants$swradZero[k]
            radScalar <- swrad/(1 + (swrad/parZero)) # Include radiation scalar
         }
         # get EVI parameters, vegetation specific
         evi <- eviset[, k]
         lswi <- lswiset[, k]

         # get temperature
         # temperature parameters, assume calculations in degrees C
         tempOpt <- vprmConstants$tempOpt[k]
         tempMax <- vprmConstants$tempMax[k]
         tempMin <- vprmConstants$tempMin[k]
         if("tempRLow"%in%names(vprmConstants)) #min temp with resp 
            tempRLow <- vprmConstants$tempRLow[k]
         else
            tempRLow <- 0*vprmConstants$tempMin[k] #set to zero

         tempScalar <- ((tempAir-tempMin)*(tempAir-tempMax))/
                     (((tempAir-tempMin)*(tempAir-tempMax))-(tempAir-tempOpt)^2)
         tempScalar[which(tempAir>tempMax | tempAir < tempMin)] <- 0

         eviMax <- eviMaxVec[, k]
         eviMin <- eviMinVec[, k]
         lswiMax <- lswiMaxVec[, k]
         lswiMin <- lswiMinVec[, k]

         # Need to eliminate ALL noise in evi and lswi--when greater than max values, set to max value, same for min (on evi)
         # with extremely low values.

#remove later
#print(paste("k:",k))
#print(paste("Trajecvprm: checking if evi threshold required",sum(evi>eviMax),sum(evi<eviMin)))
#if(sum(evi>eviMax)>0)print(paste("locations evi>max:",which(evi>eviMax, arr.ind = TRUE)))
#print(paste("Trajecvprm: checking if lswi threshold required",sum(lswi>lswiMax),sum(lswi<lswiMin)))
#end remove later

         evi[which(evi>eviMax)] <- eviMax[which(evi>eviMax)]
         evi[which(evi<eviMin)] <- eviMin[which(evi<eviMin)]
         lswi[which(lswi>lswiMax)] <- lswiMax[which(lswi>lswiMax)]
         lswi[which(lswi<lswiMin)] <- lswiMin[which(lswi<lswiMin)]

         # modification for so-called "xeric systems", comprising shrublands and grasslands
         # these have different dependencies on ground water.

         if ((landcov == "DVN"&is.element(k, c(7,10)))|(landcov != "DVN"&is.element(k, c(4,7)))) {
                 wScalar <- (lswi-lswiMin)/(lswiMax-lswiMin)
         } else {
                 wScalar <- (1+lswi)/(1+lswiMax)
         }


         # Case 1: Evergreens--pScalar always = 1
         # Case 2: Grasslands and Savannas -- pScalar never equals 1
         # Case 3: Others -- pScalar = 1 for growing season only as determined from evi threshold
         #        evithreshold = eviMin  + (0.55*range); range = eviMax - eviMin

         pScalar <- (1+lswi)/2
         if ((landcov == "DVN"&is.element(k,1:4))|(landcov != "DVN"&k == 1)) { # if evergreen
                 phenologyselect <- 1:length(evi)
         }
         if ((landcov == "DVN"&is.element(k, c(5,6,7,9,12)))|(landcov != "DVN"&is.element(k, c(2,3,4,6,8)))) { # if decid, mixed, shrub, crop, or other
                 threshmark <- 0.55
                 evithresh <- eviMin+(threshmark*(eviMax-eviMin))
                 phenologyselect <- which(evi>evithresh)
         }
         # by default, grasslands and savannas and peatlands never have pScale=1

         pScalar[phenologyselect] <- 1

         # Error Check evi, lswi, tempScalar
         evi[which(is.na(evi))] <- 0   #checked for Jena implementation: na check not really required
         lswi[which(is.na(lswi))] <- 0 #checked for Jena implementation: na check not really required
         tempScalar[which(is.na(tempScalar))] <- 0  #checked for Jena implementation: na check not really required
         tempScalar[which(tempScalar<0)] <- 0
         # Determine GPP
         # NOTE UNITS--VPRM outputs GPP and Respiration in umol/m2/s (conveniently, what is needed here); when multiplied by
         #                influence (ppm/(umol/m2/s)) get ppm
         xiaoGPP <- lambdaGPP*tempScalar*wScalar*pScalar*evi*radScalar*v.tot*-1

         # want vegetative uptake to be negative with respect to the atmosphere, so multiply by negative one
         # for symmetry
         xiaoGPP[which(is.na(xiaoGPP))] <- 0

         v.gee <- xiaoGPP

         tempForR <- tempAir
         tempForR[which(tempForR<tempRLow)] <- tempRLow

         devanResp <- alphaResp*tempForR*v.tot+betaResp*v.tot
         devanResp[which(is.na(devanResp))] <- 0
         v.resp <- devanResp


         # Output
         cnameAdd <- paste("v", (nBaseVeg+3)+(k-1)*3+0:2, sep="")
         msbMat <- cbind(v.tot, v.gee, v.resp)
         dimnames(msbMat) <- list(NULL, cnameAdd)

         result.ready <- cbind(result.ready, msbMat)
      } # for k
   } # of if (bios == "VPRM")
}# if co2

cnameAdd <- paste("v", nBaseVeg+2+ifelse("co"%in%tracers&!("co2"%in%tracers),0, length(output.veg)*3+1), sep="")
if ("ch4"%in%tracers&!ncdfTF["ch4"]) { # also put ch4 fluxes....
   # CH4: need month selector
   monthlong <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60)[, "mon"]
     v.ch4 <- rep(0, nrow(result)) # initialize
     for (monch4 in 1:12) { # loop over 12 months
               v.ch4[monthlong == monch4] <- result[monthlong == monch4, paste("v", ptch4[monch4], sep="")]
     }
   v.ch4 <- matrix(v.ch4, ncol=1, dimnames=list(NULL, cnameAdd))

   result.ready <- cbind(result.ready, v.ch4)
}
if ("ch4"%in%tracers&ncdfTF["ch4"]) { # also put ch4 fluxes....
   result.add <- matrix(result[, paste("v", ptch4, sep="")], ncol=1, dimnames=list(NULL, cnameAdd))
   result.ready <- cbind(result.ready, result.add)
}

if ("n2o"%in%tracers&ncdfTF["n2o"]) { # also put ch4 fluxes....
   cnameAdd <- paste("v", as.numeric(substring(cnameAdd,2))+1, sep="")
   result.add <- matrix(result[, paste("v", ptn2o, sep="")], ncol=1, dimnames=list(NULL, cnameAdd))
   result.ready <- cbind(result.ready, result.add)
}

if ("h2"%in%tracers&!ncdfTF["h2"]) { # H2: use soil moisture parameterization in presence of vegetation
   # For volumetric soil moisture 0- 0.15:
   # -Flux=80*soil moisture
   # For soil moisture 0.15-0.50
   # -Flux=-34.3*soil moisture + 17.14
   # The units are nmol m-2 s-1.  You get max flux of 12 at 0.15 linearly
   # decreasing to 0 at soil moisture values of 0 and 0.5.

   idveg <- which(is.element(reclss,1:(nReclss-1))) # Assume vegetation in all categories except the 'final' landclass

   v.h2 <- apply(result[, paste("v", as.character(idveg), sep='')],1, sum) # add all influences from this vegetation class
   if (!("solw"%in%dimnames(result)[[2]]))
          stop("need soil moisture for H2")
   vsm <- result[, "solw"]/1000
   v.h2 <- v.h2*( (vsm<0.15)*(80*vsm) + (vsm>=0.15&vsm<0.5)*(17.14-34.3*vsm))/1000 # from nmol to mumol

   cnameAdd <- paste("v", as.numeric(substring(cnameAdd,2))+1, sep="")
   v.h2 <- matrix(v.h2, ncol=1, dimnames=list(NULL, cnameAdd))

   result.ready <- cbind(result.ready, v.h2)
}
if ("h2"%in%tracers&ncdfTF["h2"]) { # also put h2 fluxes....
   cnameAdd <- paste("v", as.numeric(substring(cnameAdd,2))+1, sep="")
   result.add <- matrix(result[, paste("v", pth2, sep="")], ncol=1, dimnames=list(NULL, cnameAdd))
   result.ready <- cbind(result.ready, result.add)
}
if ("cofire"%in%tracers&ncdfTF["cofire"]) { # also put ch4 fluxes....
   cnameAdd <- paste("v", as.numeric(substring(cnameAdd,2))+1, sep="")
   result.add <- matrix(result[, paste("v", ptcofire, sep="")], ncol=1, dimnames=list(NULL, cnameAdd))
   result.ready <- cbind(result.ready, result.add)
}

# rename all flux inputs, start with v1
dimnames(result.ready)[[2]][substring(dimnames(result.ready)[[2]],1,1) == "v"] <- paste("v",
                1:length(dimnames(result.ready)[[2]][substring(dimnames(result.ready)[[2]],1,1) == "v"]), sep="")
# this means, v1 is water influence, v2 and v3 are CO2 and CO ff emission inputs, then influence, gee and R due to different vegetations, 
# then (if req.) CH4, then H2, then biomass burning 
# For reinsertion
result <- result.ready

        # get selector for last observation for each particle
ordert <- order(result[, "index"], result[, "btime"]) # order first by index, then by time
ordern <- order(result[ordert, "btime"], result[ordert, "index"]) # to recreate original order
delbte <- c(diff(result[ordert, "btime"]), -1000)[ordern] # timestep will be negative at last obs. for each particle
selend <- delbte<0

        # get selector for first observation for each particle
delbte <- c(-1000, diff(part[ordert, "btime"]))[ordern] # timestep will be negative at first obs. for each particle
selfirst <- delbte<0

        # get initial conditions for CO and CO2, read with read <- co2 <- bg <- aug00.ssc and read <- co <- bg <- aug00.ssc
        # co2.ini and co.ini objects are 3d arrays, agl*lat*sasdate; sasdate is 0 at 1/1/1960 (i.e. elapsed days since 1/1/1960, same as julian(m, d, y))
aglg <- round(result[selend, "agl"]/500) # get agl in 500 m intervals
latg <- round((result[selend, "lat"]-10)/2.5) # get lat in 2.5 deg intervals, starting at 10 deg

for(tr in tracers){ #default: no boundary condition specified
   cini <- result[,1]*0                                     # initialize vector with the right length
   result <- cbind(result, cini)
}
dimnames(result)[[2]][(dim(result)[2]-length(tracers)+1):dim(result)[2]] <- paste(tracers, "ini", sep="") # give advected boundary mixing ratio names
tracers.clim<-tracers[inikind[tracers]=="climat"]
if(length(tracers.clim)>0){ #are there climatological boundary files to be used
        # get 1st boundary field to set things up
   if (!existsr(paste(tracers.clim[1], ".ini", sep=""), pathname)){
           print("need to use read.bg() to get boundary condition")
           stop(paste("need ", tracers.clim[1], " boundary condition", sep=""))
   }
   ini <- getr(paste(tracers.clim[1], ".ini", sep=""), pathname)

   startday <- as.numeric(dimnames(ini)[[3]][1]); delday <- as.numeric(dimnames(ini)[[3]][2])-startday

   sasdate <- julian(mon, day, yr4)-result[selend, "btime"]/24-startday # days elapsed since beginning of .ini file
   pointer <- cbind(aglg+1, latg+1, round(sasdate/delday)+1) # array indices must start with 1
           # use constant ini field when no initial data available
   if (any(pointer[,3]>dim(ini)[3]))
           cat("Trajecvprm(): extrapolating ", tracers.clim[1], ".ini, need later times\n", sep="")
   if (any(pointer[,3]<1))
           cat("Trajecvprm(): extrapolating ", tracers.clim[1], ".ini, need earlier times\n", sep="")
   pointer[pointer[,3]>dim(ini)[3],3] <- dim(ini)[3]
   pointer[pointer[,3]<1,3] <- 1

           # limit to lat range of initial field
           # limit pointer to valid values for latg and aglg
   pointer[pointer[,1]>dim(ini)[1],1] <- dim(ini)[1]
   pointer[pointer[,2]>dim(ini)[2],2] <- dim(ini)[2]
   pointer[pointer[,1]<1,1] <- 1
   pointer[pointer[,2]<1,2] <- 1

   for (tr in tracers.clim[1:length(tracers.clim)]) {
      if (tr!=tracers.clim[1]) {
         if (!existsr(paste(tr, ".ini", sep=""), pathname)) {
            stop(paste("warning: need ", tr, " boundary condition", sep=""))
         }
         ini <- getr(paste(tr, ".ini", sep=""), pathname)
      }
      cini <- result[,1]*0                                     # initialize vector with the right length
      cini[selend] <- ini[pointer]
      result[,paste(tr, "ini", sep="")] <- cini
   }
} #if(length(tracers.clim)>0){ #are there climatological boundary files to be used


####################################################################################################
# Step 12: add CO2 boundary values
####################################################################################################

if ("co2"%in%tracers & inikind["co2"] == "CT") {
   cat("Trajecvprm: using CarbonTracker initial values.\n")
   result <- get.CarbonTracker.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                              result=result, result.sel=selend)
} else if ("co2"%in%tracers & inikind["co2"] == "TM3") {
   cat("Trajecvprm: using TM3 initial values.\n")
   ftype<-substring(inifile["co2"],nchar(inifile["co2"])-1,nchar(inifile["co2"]))
   if(ftype==".b")result <- get.TM3.bin(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
   if(ftype=="nc")result <- get.TM3.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
} else if ("co2"%in%tracers & inikind["co2"] == "LMDZ") {
   cat("Trajecvprm: using LMDZ initial values.\n")
   result <- get.LMDZ.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
}


if ("co"%in%tracers & inikind["co"] == "climat") {                                      # CO: need also advected boundary value with no chemistry
   coinio <- result[, "coini"]                              # no chemistry yet
   result[selend, "coini"] <- result[selend, "CO.fact"]*coinio[selend]+result[selend, "CO.frCH4"]
   result <- cbind(result, coinio)
} else if ("co"%in%tracers & inikind["co"] == "GEMS") {
   cat("Trajecvprm: using GEMS CO initial values.\n")
   result <- get.GEMS_CO.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co"],
                    result=result, result.sel=selend,spec="coini")
}
if ("cofire"%in%tracers & inikind["cofire"] == "GEMS") {                                      # CO: need also advected boundary value with no chemistry
   cat("Trajecvprm: using GEMS CO initial values for cofire.\n")
   result <- get.GEMS_CO.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["cofire"],
                    result=result, result.sel=selend,spec="cofireini")
}


if ("ch4"%in%tracers & inikind["ch4"] == "TM3") {
   cat("Trajecvprm: using TM3 initial values for CH4.\n")
   result <- get.TM3CH4.bin(yr4=yr4, mon=mon, day=day, hr=hr, ch4inifile=inifile["ch4"],
                    result=result, result.sel=selend)
}

dimnames(result) <- list(NULL, dimnames(result)[[2]])

# add all influences
v <- result[, substring(dimnames(result)[[2]],1,1) == "v"]
sv <- v*0 # initialize sum of veg. influence
if(is.matrix(v)){ #multiple tracers
   dimnames(sv) <- list(NULL, paste("s", dimnames(sv)[[2]], sep=""))
   for (veg in 1:length(dimnames(sv)[[2]])) {
      cumsamp <- tapply(v[, paste("v", veg, sep="")], result[, "index"], cumsum)
      cumsamp <- unlist(cumsamp)[ordern]
      sv[, paste("sv", veg, sep="")] <- cumsamp
   }# for veg in different surface flux variables
   result <- cbind(result, sv)
} else { #only single tracer
   cumsamp <- tapply(v, result[, "index"], cumsum)
   cumsamp <- unlist(cumsamp)[ordern]
   sv1 <- cumsamp
   result <- cbind(result, sv1)
}
#################################
dimnames(result) <- list(NULL, dimnames(result)[[2]])

# drop all vegetation fluxes, keep only the accumulated fluxes
resulto<-result
result <- result[, substring(dimnames(result)[[2]],1,1) != "v"]
dimnames(result) <- list(NULL, dimnames(result)[[2]])


# keep only information which is required: last entry position and accumulated fluxes for coarse vegetation classes, and CO and CO2 fossil fuel emissions
lastresult <- result[selend>0, ]
dimnames(lastresult) <- list(NULL, dimnames(lastresult)[[2]])

# keep also certain parameters at beginning of particle trajectory
firstresult <- result[selfirst>0, ]

# give each particle the same weight throughout the run,
# discard a complete particle run when particle stops before max. btime and before it reaches boundary: done with object "part" higher up
# get stdeviation for co2ini, coini
if (length(tracersini) > 1) {
   sdini <- apply(lastresult[, tracersini], 2, stdev, na.rm=T)
} else {
   sdini <- stdev(lastresult[, tracersini])
}

# get number of particles ending insight area before btime=btimemax; allow for ~2.5 degree band
# at border and 1 day time
endinarea <- (lastresult[, "gitx"]>10&lastresult[, "gitx"]<(numpix.x-10))&(lastresult[, "gity"]>10&lastresult[, "gity"]<(numpix.y-10))&lastresult[, "btime"]<nhrs-0.05*nhrs
endinarea <- sum(endinarea)
# get averages
lastresult <- apply(lastresult,2, mean)
firstresult <- apply(firstresult,2, mean)

# keep only ident, btime, lat, lon, co2ini, coini, mixed layer height and ground height, at start,
# at ending in area, sdev. co and co2 ini, influence from water, and veg. influences
output <- c(pos, lastresult[c("btime", "lat", "lon", "agl", tracersini)], firstresult[c("zi", "grdht")], endinarea,
            sdini, lastresult[substring(names(lastresult),1,2) == "sv"])
names(output) <- names.output

if (detailsTF) {
  resulthr<-get.mean.traj(resulto,hr=2) #only for looking at data...
  selnam <- substring(dimnames(resulthr)[[2]],1,1) == "v"
  dimnames(resulthr)[[2]][selnam] <- names.output[(length(names.output)-sum(selnam)+1):length(names.output)]
  assignr(paste(ident,"resulthr", sep=""), resulthr,pathname, printTF=TRUE)

  selnam <- substring(dimnames(result)[[2]],1,2) == "sv"
  dimnames(result)[[2]][selnam] <- names.output[(length(names.output)-sum(selnam)+1):length(names.output)]
  if(ncdfTF["co2"]&bios=="VPRM"){
    vprmstuff<-cbind(eviset,eviMaxVec,eviMinVec)
    dimnames(vprmstuff)[[2]]<-c(paste(output.veg,".evi",sep=""),paste(output.veg,".evimax",sep=""),paste(output.veg,".evimin",sep=""))
    result<-cbind(result,vprmstuff)
  }
  assignr(paste(ident, "result", sep=""), result, pathname, printTF=TRUE)
}

gc(verbose=F)
return(output)

}
