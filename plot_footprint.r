# Script to make footprint plots with OCO-2 observations, DW, 05/12/2020

homedir = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'models/X-STILT') # current dir
source(file.path(xstilt_wd, 'r/dependencies.r'))

api.key <- readLines('insert_ggAPI.csv')
if (api.key == '') stop('Missing googleAPI, insert your API in insert_ggAPI.csv\n')
register_google(key = api.key)


# choose the site
site = 'Paris'
input.path  = '/central/groups/POW'
store.path  = file.path(homedir, 'postdoc_proj/XSTILT_output', site)
out.path    = list.files(store.path, 'out_20', full.names = T)
all.timestr = strsplit.to.df(basename(out.path))$V2


# **** choose a time stamp and load all footprint 
tt = 1
timestr   = all.timestr[tt]
byid.path = file.path(out.path, 'by-id')
foot.fns  = list.files(byid.path, 'X_foot.nc', full.names = T, recursive = T)
recp.info = strsplit.to.df(basename(foot.fns))


## whether this is an ideal simulations
idealTF = T 
if (!idealTF) {
    oco.sensor  = c('OCO-2', 'OCO-3')[2]
    oco.ver     = c('V7rb', 'V8r', 'V9r', 'VEarlyR')[4]   # retrieval algo version
    oco.path    = file.path(input.path, oco.sensor, paste0('L2_Lite_FP_', oco.ver))

} else oco.sensor = oco.ver = oco.path = NULL   # set oco.* to NULL


## plot mean footprint from all receptors or footprint from one single receptor
meanTF = TRUE 
indx   = 35       # if plotting foot from one receptor, give it an index


if (meanTF) {   # 1) If you wanna plot spatial mean footprint -----------------

    # load all footprints
    foot.stk = stack(foot.fns); names(foot.stk) = basename(foot.fns)

    # take spatial average, this could take a while given large # of receptors
    foot.rt = mean(foot.stk)  
    picname = file.path(store.path, paste0('footprint_', site, '_', timestr, '_mean.png'))
    recp.lat = recp.lon = NULL 

    # load map
    mm = ggplot.map(map = 'ggmap', maptype = 'roadmap', zoom = 7, 
                    center.lat = mean(as.numeric(recp.info$V3)),
                    center.lon = mean(as.numeric(recp.info$V2)))

} else {    # 2) OR if plotting footprint from one receptor -------------------

    # load footprint from one receptor indicated by `indx`
    foot.rt  = raster(foot.fns[indx])
    recp.lat = as.numeric(recp.info$V3[indx])
    recp.lon = as.numeric(recp.info$V2[indx])
    picname = file.path(store.path, paste0('footprint_', site, '_', timestr, '_', 
                                            recp.lat, '_', recp.lon, '.png'))
    
    mm = ggplot.map(map = 'ggmap', maptype = 'roadmap', zoom = 7, 
                    center.lat = recp.lat, center.lon = recp.lon)

}   # end if mean TF


# plot footprint -------------------------------------------------
# convert footprint from rasterLayer to data frame and rename it
foot = as.data.frame(foot.rt, xy = T); colnames(foot) = c('lon', 'lat', 'foot')


# ggmap.xfoot.obs will automatically store figures in store.path 
title = paste('Spatial footprint for', site, 'on', timestr)
f1 = ggmap.xfoot.obs(mm, site, oco.ver, oco.path, timestr, recp.lon, recp.lat, 
                     foot, min.foot.sig = 1E-8, max.foot.sig = 1E-2, qfTF = T, 
                     title, picname, storeTF = T, width = 7, height = 7, 
                     leg.pos = 'bottom', scale.coord = 1.5)

