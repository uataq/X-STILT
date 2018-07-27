# plot footprint using ggplot2
# Dien Wu

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

# plot which variables
p <- 1
plot.var <- c("foot_anthro_", "foot_bio_", "foot_ocean_", "intfoot")[p]
nc.var <- c("foot_anthro", "foot_bio", "foot_ocean", "footprint")[p]
plot.str <- c("XCO2.anthro", "XCO2.bio", "XCO2.ocean", "Footprint")[p]
outdir <- c('akpw_xco2', 'akpw_xco2', 'akpw_xco2', 'akpw_intfoot')[p]
limits <- list(c(1E-6, 1E-4, 0.05, -3), c(1E-6, 1E-4, 0.01, -2),
               c(1E-6, 1E-4, 0.01, -2), c(1E-8, 1E-6, 8E-4, -3))[[p]]
signal <- limits[1]
plot.min <- limits[2]
plot.max <- limits[3]
mp <- limits[4]  # midpoint in log

# paths and files
outdir <- file.path(workdir, 'output', 'NetCDF', outdir)
plotdir <- find.create.dir(path=file.path(workdir,'output'), dir='Plots',
                           workdir=workdir)

file <- list.files(outdir, paste(plot.var, filestr, sep=""))
ident <- substr(file, nchar(plot.var)+1, nchar(file)-3)
info <- ident.to.info(ident, aglTF=F)[[1]]

# select regions data
minlat <- 28
maxlat <- 34
minlon <- 32
maxlon <- 40
nhrs <- -72
center.lat <- 31
center.lon <- 36

#### plotting variables
font.size = rel(1.3); alpha = 0.8; zoom <- 7
# load google maps
mm <- ggplot.map(map="ggmap", center.lat=center.lat,
                 center.lon=center.lon, zoom=zoom)

m1 <- mm[[1]] + theme_bw()+ coord_equal(1.1)
shift.lat <- mm[[2]]; shift.lon <- mm[[3]]

pp<-list()

# loop over all exceeding events--
for (f in 1:length(file)){

  tmp.file <- file[f]

  cat("Working on time", info$timestr[f], "recp", info$recp.lat[f],"N...\n")
  dat <- nc_open(file.path(outdir, tmp.file))

  # read data and lat/lon
  tmp.dat <- ncvar_get(dat, nc.var)  # [lon,lat]
  tmp.lat  <- ncvar_get(dat, 'Lat')   # cell center lat
  tmp.lon  <- ncvar_get(dat, 'Lon')   # cell center lon
  dimnames(tmp.dat) <- list(tmp.lat, tmp.lon)
  print(max(tmp.dat))
  nc_close(dat)

  melt.dat <- melt(tmp.dat)
  colnames(melt.dat) <- list('lat', 'lon', 'value')
  sel.dat <- melt.dat[melt.dat$lat > minlat & melt.dat$lat < maxlat &
                      melt.dat$lon > minlon & melt.dat$lon < maxlon &
                      melt.dat$value > signal,]
  print(max(sel.dat$value))

  sel.loc <- data.frame(lon=info$recp.lon[f]+shift.lon,
                        lat=info$recp.lat[f]+shift.lat, fac=info$recp.lat[f])

  # plotting...
  title <- paste("Maps for receptor", info$recp.lat[f], "N over", site,
                 "on", info$timestr[f])
  p1 <- m1 + labs(x="LONGITUDE [E]",y="LATITUDE [N]", title=title)
  p2 <- p1 + geom_raster(data=sel.dat,alpha=alpha,
                         mapping=aes(x=lon+shift.lon, y=lat+shift.lat, fill=value))

  # #ffffff • #d0e1f9 • #4d648d • #283655 • #1e1f26
  #col<-data.frame(low="#d0e1f9",mid="#4d648d",high="#283655", stringsAsFactors = FALSE)
  #col<-data.frame(low="#bdeaee",mid="#76b4bd",high="#58668b", stringsAsFactors = FALSE)
  #col<-data.frame(low="#ffaf7b",mid="#d76d77",high="#3a1c71", stringsAsFactors = FALSE)
  col <- data.frame(low="yellow", mid="orange", high="red",
                    stringsAsFactors = FALSE)

  p2 <- p2 + scale_fill_gradient(limits=c(plot.min, plot.max), name=plot.str,
                                  low=col$low, high=col$high, trans="log10")

  p3 <- p2 + theme(legend.position="bottom", legend.key.width=unit(3, "cm"),
                   legend.key.height=unit(0.5, "cm"),
                   legend.text=element_text(size=font.size),
                   legend.key = element_blank(),
                   axis.title.y=element_text(size=font.size,angle = 90),
                   axis.title.x=element_text(size=font.size,angle=0),
                   axis.text=element_text(size=font.size),
                   axis.ticks=element_line(size=font.size),
                   title=element_text(size=font.size))

  p3 <- p3 + geom_text(data=sel.loc,label="receptor", size=6,
                       mapping=aes(x=lon+0.5, y=lat), colour="purple")

  p3 <- p3 + geom_point(data=sel.loc, size=2, pch=21, colour="purple",
                        mapping=aes(x=lon, y=lat))

  pp[[f]] <- p3
  filenm <- paste(nc.var, "_", info$recp.lat[f], "N_", site,"_",
                  info$timestr[f],"_", met,".png",sep="")
  ggsave(p3, file=file.path(plotdir, filenm), width=12, height=12)
}

#library(ggpubr)
#mm <- ggarrange(plotlist=pp, labels = paste(seq(1,f,1),")",sep=""), ncol=4, nrow=2)
#filenm <- paste("Int_Foot_Yuma_", met, "_", year, "_", abs(nhrs)/24, "days.png", sep="")
#ggsave(mm, file=filenm, width=40, height=20)

# end of script

if(F){

  ##### read trajec
  traj.files <- list.files(pattern='traj', path=traj.path)
  traj.files <- traj.files[substr(traj.files, 1, 4)==as.character(year)]
  ident <- matrix(unlist(strsplit(traj.files,'_')), ncol = 5, byrow=T)
  ident <- data.frame(ident, stringsAsFactors = FALSE)
  colnames(ident) <- list('timestr', 'lon', 'lat', 'agl', 'format')

  # loop over all files
  each.hr <- 24   # plot every ? hours
  for(t in 1:length(traj.files)){

    # get trajectory info
    traj.info <- readRDS(file=file.path(traj.path, traj.files[t]))
    trajdat <- traj.info$particle
    recp.info <- traj.info$receptor
    print(range(trajdat$time))

    # get met file types
    met <- "EDAS40km"
    if(ident$timestr[t] > as.character("2015061500"))met <- "HRRR"

    # plot trajec
    zoom <- 5
    mm <- ggplot.map(map="ggmap", zoom=zoom,
                     center.lat=recp.info$lati, center.lon=recp.info$long)
    m1 <- mm[[1]]
    sel.trajdat <- trajdat[trajdat$time %% (each.hr*60) ==0, ]

    p1 <- m1+geom_point(data=sel.trajdat, aes(x=long, y=lati, colour=as.factor(time)), pch=21, size=1)
    p1 <- p1+scale_colour_manual(values=ggcol.def(n=length(unique(sel.trajdat$time))),
                                 breaks=unique(sel.trajdat$time), name="hrs back",
                                 labels=unique(sel.trajdat$time)/60)
    ggsave(p1, file=paste("Traj_Yuma_",ident$timestr[t],"_",met,".png",sep=""), width=12, height=12)
  }

}
