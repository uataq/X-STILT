# subroutine for running X-STILT trajec
# Written by DW, update on 05/19/2017 (add GDAS-driven STILT)
# ADD fixed agl run, 07/07/2017, DW
# clear up, 04/30/2017, DW

# read variables from output of "create_namelist_trajec.r"
get.more.namelist <- function(namelist, plotTF){

  # ------------------- Step 1. READ IN OCO-2 LITE FILES --------------------- #
  library(ncdf4)

  YYYYMMDD <- substr(namelist$timestr, 1, 8)
  YYMMDD <- substr(namelist$timestr, 3,8)
  ocofile <- list.files(pattern = YYMMDD, path = namelist$ocopath)
  ocodat <- nc_open(file.path(ocopath, ocofile))

  # grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  # may need to reverse later
  ocolevel <- ncvar_get(ocodat, "levels")
  ocolat <- ncvar_get(ocodat, "latitude")
  ocolon <- ncvar_get(ocodat, "longitude")
  ocofoot <- ncvar_get(ocodat, "Sounding/footprint")

  xco2 <- ncvar_get(ocodat, "xco2")
  xco2.uncert <- ncvar_get(ocodat, "xco2_uncertainty")  # posterior uncertainty

  # grabbing time for STILT receptors
  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id <- as.character(ncvar_get(ocodat,"sounding_id"))
  YYYY <- substr(id,1,4)
  MM <- substr(id,5,6)
  DD <- substr(id,7,8)
  HH <- substr(id,9,10)
  mm <- substr(id,11,12)
  ss <- substr(id,13,14)
  ms <- substr(id,15,15)
  f <- substr(id,16,16)
  YYYYMMDDHHmmss <- substr(id,1,14)

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition
  OM <- ncvar_get(ocodat, "Sounding/operation_mode")
  LF <- ncvar_get(ocodat, "Sounding/land_fraction")   # > 80%: land, < 20%: sea
  QF <- ncvar_get(ocodat, "xco2_quality_flag")
  WL <- ncvar_get(ocodat, "warn_level")

  # operation modes:
  mode <- rep(NULL, length(OM))
  mode[LF> 80 & OM==0] <- "Land_Nadir"
  mode[LF< 20 & OM==0] <- "Sea_Nadir"
  mode[LF> 80 & OM==1] <- "Land_Glint"
  mode[LF< 20 & OM==1] <- "Sea_Glint"
  mode[LF> 80 & OM==2] <- "Land_Target"
  mode[LF< 20 & OM==2] <- "Sea_Target"
  mode[LF> 80 & OM==3] <- "Land_Transition"
  mode[LF< 20 & OM==3] <- "Sea_Transition"

  ## SELECT REGION
  lat.lon <- namelist$lat.lon   # minlat, maxlat, minlon, maxlon
  region.index <- which(ocolat >= lat.lon[1] & ocolat < lat.lon[2] &
                        ocolon >= lat.lon[3] & ocolon < lat.lon[4])

  sel.lat <- ocolat[region.index]
  sel.lon <- ocolon[region.index]
  sel.xco2 <- xco2[region.index]
  sel.qf <- QF[region.index]
  sel.wl <- WL[region.index]
  sel.id <- id[region.index]

  sel.foot <- ocofoot[region.index]  # OCO-2 footprint
  sel.time <- substr(id,1,10)[region.index]
  sel.YYYY <- YYYY[region.index]
  sel.MM <- MM[region.index]
  sel.DD <- DD[region.index]
  sel.HH <- HH[region.index]
  sel.mode <- mode[region.index]
  cat("Operational Modes:", unique(sel.mode),"\n\n")

  # whether plotting XCO2 from OCO-2
  if(plotTF){
    library(ggmap); library(ggplot2)
    xco2.obs <- data.frame(lat=as.numeric(sel.lat), lon=as.numeric(sel.lon),
                           xco2=as.numeric(sel.xco2))
    # lat/lon for city
    obs.lat <- c(24.63-0.1,30.05+0.2,22.40+1)[1]
    obs.lon <- c(46.72+0.1,31.23-0.2,114.11-0.4)[1]
    zoom <- 8; alpha <- 1; font.size <- rel(1.2)
    col <- c("black", '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
            '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')

    # plot google map
    sitemap <- get_map(location=c(lon=obs.lon,lat=obs.lat), zoom=zoom,
                       maptype='roadmap')
    m1 <- ggmap(sitemap)+theme_bw()

    c1 <- m1+geom_point(data=xco2.obs,aes(x=lon,y=lat,colour=xco2))
    c1 <- c1+labs(x="Longitude",y="Latitude")
    c1 <- c1+labs(title=paste("OCO-2 XCO2 [ppm] for",site,"on",YYYYMMDD))
    c1 <- c1+scale_colour_gradientn(name="OCO-2 XCO2 [ppm]",colours=col,
                                   breaks=seq(380,420,2), labels=seq(380,420,2))

    c2 <- c1+theme(legend.position="bottom",
                   legend.text=element_text(size=font.size),
                   legend.key=element_blank(),
                   legend.key.height = unit(0.5, "cm"),
                   legend.key.width = unit(3, "cm"),
                   axis.title.y=element_text(size=font.size,angle=90),
                   axis.title.x=element_text(size=font.size,angle=0),
                   axis.text=element_text(size=font.size),
                   axis.ticks=element_line(size=font.size),
                   title=element_text(size=font.size))

    ggsave(c2,filename=paste("ggmap_xco2_",site,"_",YYYYMMDD,".png",sep=""),
              width=11,height=12)
  }  # end if plotTF

  # ------------------- Step 2. SET UP the STILT receptors ------------------- #
  # e.g., Time, Lat, Lon, and Height based on observation
  yr <- as.numeric(sel.YYYY)-2000
  mon <- as.numeric(sel.MM)
  day <- as.numeric(sel.DD)
  hr <- as.numeric(sel.HH)

  ### for most cases, no need to modify STILT variables below
  # T for convection (RAMS winds: grell convection scheme,EDAS and FNL: simple
  # redistribution within vertical range with CAPE>0)
  namelist$convect <- F

  # Enforces Courant criterion. When >0: maximum horizontal distances travelled
  # by a particle in a single timestep, in km.
  namelist$stepsize <- 0

  # Set to 2000 to keep particles from dropping out of high res WRF runs
  namelist$mgmin <- 2000
  namelist$veght <- 0.5         # surface layer, half of PBL
  namelist$metd <- c("fnl","awrf") # meteorological file types
  namelist$varstrajec <- c("time","index","lat","lon","agl","grdht", "foot",
                           "sampt","dmass","zi","pres")

  # update release AGLs
  recp.info <- list()
  if(namelist$columnTF){
    namelist$agl <- namelist$agl
    namelist$npar <- namelist$npar   # total number of particles
  }

  ### 'lat','lon',&'agl' can be a VECTOR of the same length
  # filter by quality flag, more receptors when XCO2 is high
  if(namelist$oco2.version == "b8r")sort.lat <- sort(sel.lat[sel.wl == 0])
  if(namelist$oco2.version == "b7rb")sort.lat <- sort(sel.lat[sel.qf == 0])
  recp.index <- namelist$recp.index
  recp.lat <- sort.lat[findInterval(recp.index, sort.lat)]
  match.index <- unique(match(recp.lat, sel.lat))

  # change order for longitude, and time
  namelist$recp.lat <- signif(sel.lat[match.index], 6)
  namelist$recp.lon <- signif(sel.lon[match.index], 7)
  namelist$recp.yr  <- yr [match.index]
  namelist$recp.mon <- mon[match.index]
  namelist$recp.day <- day[match.index]
  namelist$recp.hr  <- hr [match.index]

  # return namelist
  return(namelist)

} # end of subroutine
