### subroutine to draw STILT particle locations on 2D or 3D plot
# by Dien Wu, 06/21/2017
# convert kernel density to percentile, DW, 11/15/2017
# return latitude range, DW, 11/16/2017

ggplot.forward.trajec <- function(ident, trajpath, site, timestr, ocopath, zoom,
                                   lat.lon){

  # need joinPolys() and point.in.polygon() functions
  library(ggplot2); library(ncdf4); library(MASS)
  library(PBSmapping); library(sp)

  font.size <- rel(1.2)

  # call grab.oco2() to read in observations
  obs.all <- grab.oco2(ocopath=ocopath, timestr=timestr, lat.lon=lat.lon)

  # read in forward trajec
  trajdat <- NULL
  for(f in 1:length(ident)){
    tmp.trajdat <- getr(xname=ident[f], path=trajpath)
    trajdat <- rbind(trajdat, tmp.trajdat)
  }
  trajdat <- data.frame(trajdat)

  # before subsetting trajec, grab box receptor
  recp.trajdat <- trajdat[trajdat$time == min(trajdat$time), ]
  llon <- c(min(recp.trajdat$lon), max(recp.trajdat$lon))
  llat <- c(min(recp.trajdat$lat), max(recp.trajdat$lat))
  box.recp <- data.frame(lon=c(llon[1], llon[2], llon[2], llon[1]),
                         lat=c(llat[1], llat[1], llat[2], llat[2]))

  # only allow for an hour duration
  info <- ident.to.info(ident=ident)[[1]]
  min.recp.mm <- min(info$recp.hour)*60

  # compute overpass time
  hhmm <- paste(formatC(obs.all$hr, width=2, flag=0),
                formatC(obs.all$min, width=2, flag=0), sep="")

  min.xtime <- min(hhmm, na.rm=T)
  max.xtime <- max(hhmm, na.rm=T)
  dmin <- as.numeric(substr(min.xtime,1, 2)) * 60 +
          as.numeric(substr(min.xtime, 3, 4)) - min.recp.mm

  min.xtime <- dmin - 2         # allow for 2 mins
  max.xtime <- dmin + 2
  trajdat <- trajdat[trajdat$time <= max.xtime & trajdat$time >= min.xtime,]

  # load map, from https://susanejohnston.wordpress.com/2012/07/03/creating-a-
  # large-scale-map-using-ggplot2-a-step-by-step-guide/
  # 5th and 6th elements of lat.lon are city center
  gg <- ggplot.map(center.lat=lat.lon[5], center.lon=lat.lon[6], map="ggmap", zoom=zoom)
  m1 <- gg[[1]]; shift.lat=gg[[2]]; shift.lon=gg[[3]]

  # draw receptor box
  m1 <- m1+geom_polygon(data=box.recp,aes(x=lon,y=lat),linetype=3,fill="gray30",alpha=0.5)

  # calculate 2D kernel density
  cat("\n\nggplot.forward.trajec(): calculating kernel density...\n")
  dens <- kde2d(trajdat$lon, trajdat$lat, h=c(0.1,0.1), n=100)
  densf <- data.frame(expand.grid(lon=dens$x, lat=dens$y), prob=as.vector(dens$z))
  densf$norm.prob <- densf$prob/max(densf$prob)
  lab.norm <- c(0.05, seq(0, 1, 0.1)[-1])

  p1 <- m1+geom_contour(data=data.frame(densf),aes(x=lon,y=lat,z=norm.prob,colour=..level..),
                        breaks=lab.norm, size=1.3)
  p1 <- p1+scale_colour_gradient(name="Normalized\nKernel\nDensity",
                                 low="lightblue", high="purple", breaks=lab.norm,
                                 labels=lab.norm, limits=c(0,max(lab.norm)))

  get.info <- ggplot_build(p1)$data[[5]]
  uni.level <- unique(get.info$level)

  # assign level to particles
  assign.level <- rep(0,nrow(trajdat))
  trajdat$assign.level <- assign.level

  for(u in 1:length(uni.level)){
      tmp.bound <- get.info[get.info$level == uni.level[u], ]
      tmp.index <- point.in.polygon(trajdat$lon, trajdat$lat, tmp.bound$x, tmp.bound$y)
      trajdat[tmp.index > 0, "assign.level"] <- uni.level[u]
  }

  p1 <- p1+geom_point(data=trajdat, aes(x=lon,y=lat,colour=assign.level),size=0.4,alpha=0.5)

  # add observed soundings
  max.y <-406; min.y <- 394
  if(timestr >= "20160101"){max.y <- 410; min.y <- 400}

  p2 <- p1+geom_point(data=obs.all,aes(x=lon,y=lat,fill=xco2),shape=21,colour="gray90")+
           scale_fill_gradientn(colours=default.col(),name="XCO2",limits=c(min.y,max.y),
                                breaks=seq(min.y,max.y,2),labels=seq(min.y,max.y,2))

  # compute outmost boundary
  wide.bound <- get.info[get.info$level==min(get.info$level),]
  bound.traj <- data.frame(PID=rep(1,nrow(wide.bound)),POS=1:nrow(wide.bound),
                           X=wide.bound[,"x"],Y=wide.bound[,"y"])

  td<-0.05         # threshold for cutting polygons
  # if diff is too large, meaning there are more than one polygon
  if(length(which(abs(diff(bound.traj$Y)) >td)) >0){
    bound.traj <- bound.traj[1:which(abs(diff(bound.traj$Y))>td),]
  }

  if(length(which(abs(diff(bound.traj$X)) >td)) >0){
    bound.traj <- bound.traj[1:which(abs(diff(bound.traj$X))>td),]
  }

  p3 <- p2+geom_polygon(data=bound.traj,aes(x=X,y=Y),colour="gray10",
                        linetype=1,fill=NA,size=0.9,alpha=0.5)

  # further select OCO2 soundings  over OCO-2 track
  SEL <- obs.all$lat < max(bound.traj$Y) & obs.all$lat > min(bound.traj$Y)
  SEL <- SEL & obs.all$lon> min(bound.traj$X) & obs.all$lon< max(bound.traj$X)
  sel.obs <- obs.all[SEL,]

  if(nrow(sel.obs)==0){
    cat("ggplot.forward.trajec(): no intersection with obs...\n")
    return(p3)

  }else{

    hpts <- chull(x=sel.obs$lon, y=sel.obs$lat)
    hpts <- c(hpts,hpts[1])
    narrow.bound <- data.frame(sel.obs[hpts,])
    bound.obs <- data.frame(PID=rep(2,nrow(narrow.bound)), POS=1:nrow(narrow.bound),
                            X=narrow.bound$lon, Y=narrow.bound$lat)

    # find overlapping region and find all points in overlapping polygon
    joint.bound <- joinPolys(bound.traj,bound.obs) # in PBS Mapping package
    pol.bound   <- joint.bound

    pol.index <- point.in.polygon(obs.all$lon, obs.all$lat, pol.bound$X, pol.bound$Y)
    pol.obs <- obs.all[pol.index >0,]

    p3 <- p3+geom_polygon(data=pol.bound, aes(x=X, y=Y), colour="gray50", fill="gray40", alpha=0.5)
    p4 <- p3+geom_point(data=pol.obs, aes(x=lon, y=lat, fill=xco2), shape=21, colour="gray50")

    p5 <- p4+annotate("text", x=unique(info$recp.lon)+shift.lon, y=unique(info$recp.lat)+shift.lat-0.05,
                      label=site, fontface=1, size=6)+annotate("point", x=info$recp.lon+shift.lon,
                      y=info$recp.lat+shift.lat, size=2, shape=17)}

    p5 <- p5+theme(legend.position="right",legend.key.width=unit(0.5, "cm"),
                   legend.key.height=unit(1.5, "cm"),
                   legend.text=element_text(size=font.size),
                   legend.key = element_blank(),
                   axis.title.y=element_text(size=font.size,angle = 90),
                   axis.title.x=element_text(size=font.size,angle=0),
                   axis.text=element_text(size=font.size),
                   axis.ticks=element_line(size=font.size),
                   title=element_text(size=font.size))

    # compute mean background for both northern and southern
    if(nrow(pol.obs) > 0){

      scr.obs <- obs.all[obs.all$qf==0, ]
      north.bg.ind <- scr.obs$lat > max(pol.obs$lat) & scr.obs$lat < lat.lon[2]
      south.bg.ind <- scr.obs$lat < min(pol.obs$lat) & scr.obs$lat > lat.lon[1]

      north.bg <- mean(scr.obs[north.bg.ind, "xco2"])
      south.bg <- mean(scr.obs[south.bg.ind, "xco2"])

      mean.bg <- c(north.bg + south.bg)/2
      sd.bg <- sd(c(north.bg, south.bg))

      cat(paste("North:", signif(north.bg,5),"ppm; South:",
                          signif(south.bg,5),"ppm; Final:",
                          signif(mean.bg,5),"ppm..\n"))

      min.lat <- min(pol.obs$lat)
      max.lat <- max(pol.obs$lat)
      cat(paste("Enhanced lat range:", signif(min.lat,4), signif(max.lat,4),"N\n"))

    }else{  # if no intersection
      north.bg <- NA; south.bg <-NA; mean.bg <- NA
      cat("ggplot.forward.trajec(): No intersection with OCO-2 track...\n")
      min.lat<-NA;max.lat<-NA
    }

    title <- paste("Forward-time urban plume for overpass on",timestr)
    p5 <- p5+labs(x="LONGITUDE",y="LATITUDE",title=title)

    picname <- paste("urban_plume_forward_",site,"_",timestr,".png",sep="")
    picname <- file.path(trajpath, picname)
    ggsave(p5, filename=picname, width=12, height=12)

    # also, return the max min latitude ranges for polluted range
    bg.lat <- data.frame(north.bg, south.bg, mean.bg, min.lat, max.lat)
    return(list(p5, bg.lat, bound.traj))

  }

 # end of subroutine
