# plot lat-series
# by DW

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- "/uufs/chpc.utah.edu/common/home"
workdir <- file.path(homedir, "lin-group1/wde/github/XSTILT")
setwd(workdir) # move to working directory

 # source all functions
source(file.path(workdir, "src/sourceall.r"))
library(ncdf4); library(ggplot2)
library(reshape)

index <- 4
site  <- c("Riyadh", "Cairo", "PRD", "Jerusalem")[index]
track.timestr <- "2015061510"
filestr <- paste(substr(track.timestr,1,4), "x", substr(track.timestr,5,6), "x",
                 substr(track.timestr,7,8), sep="")
met <- c('1km', 'gdas0p5')[2]

# spatial domains placing receptors and city center, help select OCO-2 data ***
# in form of "lat.lon <- c(minlat, maxlat, minlon, maxlon, city.lat, city.lon)"
# oco2.hr for overpass hour in UTC
if(site == "Riyadh"){lat.lon <- c(23, 26, 45, 50, 24.71, 46.75); oco2.hr <- 10}
if(site == "Cairo"){lat.lon <- c(29, 31, 30, 32, NA, NA); oco2.hr <- NA}
if(site == "PRD"){lat.lon <- c(22, 23, 112, 115, NA, NA); oco2.hr <- NA}
if(site == "Jerusalem"){lat.lon <- c(31, 33, 35, 36, 31.78, 35.22); oco2.hr<-10}

odiac.vname <- "2017"
txtfile <- paste("STILT_OCO2_XCO2_", site, "_", track.timestr, "_odiac",
                  odiac.vname,"_",met,"v5.txt", sep="")

outdir <- file.path(workdir, 'output', txtfile)
sim <- read.table(outdir, header=T, sep=",")

# path for input data, OCO-2 Lite file
oco2.version <- c("b7rb","b8r")[2]        # OCO-2 version
ocostr   <- paste("OCO-2/OCO2_lite_", oco2.version, sep="")
ocopath  <- file.path(homedir, "lin-group1/wde/STILT_input", ocostr)

# grab OCO-2 data
obs <- grab.oco2(ocopath=ocopath, timestr=track.timestr, lat.lon=lat.lon)

t1 <- ggplot() + geom_point(data=obs, aes(x=lat, y=xco2))

p1 <- ggplot.sim.obs(site=site, timestr=track.timestr, met=met, lat.lon=lat.lon,
                     oco2.version=oco2.version, sim.gdas=sim, all.obs=obs,
                     plotTF=F, rmTF=F, rawTF=F, covTF=F, obs.uncertTF=F)
