# Script to make footprint plots with OCO-2 observations, DW, 05/12/2020

homedir <- '/uufs/chpc.utah.edu/common/home'
xstilt_wd <- file.path(homedir, 'lin-group7/wde/X-STILT') # current dir
register_google(key = '')

# choose the site
site <- 'Seoul'
input.path <- file.path(homedir, 'lin-group7/wde/input_data')
oco2.ver  <- c('b7rb', 'b8r', 'b9r')[3]           # OCO-2 version
oco2.path <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))

# footprint path
store.path  <- file.path(homedir, 'lin-group7/wde/output', site)
out.path    <- list.files(store.path, 'out_20', full.names = T)
all.timestr <- strsplit.to.df(basename(out.path))$V2

# **** choose a time stamp and load all footprint 
tt <- 3
timestr   <- all.timestr[tt]
byid.path <- file.path(store.path, paste0('out_', timestr, '_', site), 'by-id')
foot.fns  <- list.files(byid.path, 'X_foot.nc', full.names = T, recursive = T)
recp.info <- strsplit.to.df(basename(foot.fns))

# for example we only plot one footprint
# alternatively, you can merge multiple footprints and plot together
# **** choose which footprint file to plot
indx <- 35  
foot.rt  <- raster(foot.fns[indx])
recp.lat <- recp.info$V3[indx]
recp.lon <- recp.info$V2[indx]

# convert footprint from rasterLayer to data frame and rename it
foot <- as.data.frame(foot.rt, xy = T); colnames(foot) <- c('lon', 'lat', 'foot')
picname <- file.path(store.path, paste0('footprint_', site, '_', timestr, 
                                        '_', recp.lat, '_', recp.lon, '.png'))

# load map and plot footprint
mm <- ggplot.map(map = 'ggmap', maptype = 'roadmap', center.lat = recp.lat
                 center.lon = recp.lon, zoom = 8)

# ggmap.xfoot.obs will automatically store figures in store.path 
title <- paste('Spatial footprint for', site, 'on', timestr)

f1 <- ggmap.xfoot.obs(mm, site, oco2.ver, oco2.path, timestr, recp.lon, recp.lat, 
                      foot, min.foot.sig = 1E-8, max.foot.sig = 1E-2, qfTF = T, 
                      title, picname, storeTF = T, leg.pos = 'bottom')

