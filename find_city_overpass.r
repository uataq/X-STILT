# script to search for cities, Dien Wu, 07/03/2018
# update on 10/12/2018, DW

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT')
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

# insert your API for the use of ggplot and ggmap
api.key <- ''
register_google(api.key)

# OCO-2 version, path
oco2.ver    <- c('b7rb', 'b8r', 'b9r')[3]  # OCO-2 version
input.path  <- file.path(homedir, 'lin-group5/wde/input_data')
oco2.path   <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path    <- file.path(input.path, 'OCO-2/L2/OCO2_lite_SIF_b8r')
txtpath     <- file.path(input.path, 'OCO-2/overpass_city')
anthro.path <- file.path(input.path, 'anthromes/v2')
pop.path    <- file.path(input.path, 'GPWv4')


#---------------------------------------------------------------------------- #
# one can add more urban regions here
site.list <- c(   
  # middle east and SE Asia (except for China and Japan)
  'Riyadh',  'Medina',  'Mecca',     'Cairo',     'Jerusalem', 
  'Jeddah',  'Karachi', 'Tehran',    'Baghdad',   'Ankara',    
  'Istanbul, Turkey', 'Mumbai',  'Delhi',   'Bangalore', 'Hyderabad', 
  'Ahmedabad', 'Manila', 'Jakarta', 'Bangkok', 'HoChiMinh', 
  'Singapore', 'KualaLumpur',

  # China and Japan
  'Nanjing',   'Suzhou',   'Shanghai', 'Beijing', 'Tianjin', 'Xian',  'Lanzhou', 
  'Zhengzhou', 'HongKong', 'Seoul',    'Busan',   'Nagoya',  'Tokyo', 'Osaka',

  # Africa and Europe
  'Lagos',  'Luanda', 'Kinshasa', 'CapeTown',  'Johannesburg', 'Moscow',
  'Paris',  'London', 'Madrid',   'Rome',      'Berlin',       'Milan',
  'Athens', 'StPetersburg, Russia', 'Barcelona, Spain',

  # North America
  'Indianapolis', 'Chicago',  'Phoenix', 'Denver',     'SaltLakeCity',
  'LasVegas',     'Houston',  'Dallas',  'LosAngeles', 'Albuquerque',
  'NY',           'DC',       'Miami',   'Atlanta',    'Seattle',
  'Toronto',      'Montreal', 'Vancouver',

  # central & south America and Australia
  'SaoPaulo', 'Lima', 'Perth', 'Sydney', 'Brisbane', 'Melbourne', 'Adelaide',

  # updates on 10/12/2018
  'Guadalajara', 'Bogota',      'RiodeJaneiro', 'Santiago',  'MexicoCity', 
  'Caracas',     'BuenosAires', 'Salvador',     'Brasilia',  'Fortaleza', 
  'Venice',      'Doha',        'Chengdu',      'Hangzhou',  'Wuhan', 
  'Shenyang',    'Chongqing',   'Calgary',      'Edmonton',  'Ottawa',   
  'Quebec',      'Dhaka',       'Kolkata',      'Khartoum',  'Kiev', 
  'Vienna',      'Bucharest',   'Darwin, Australia',  'Pretoria', 'Buffalo'
)

# define urban box and data range, default is 1x1 deg 
dlon.urban <- 0.5; dlat.urban <- 0.5; date.range <- c('20140901', '20181231')

# txtfile for all cities and tracks
fn <- paste0('city_track_veg_', oco2.ver, '_v', Sys.Date(),'.txt')


# ---------------------------------------------------------------------------
info <- NULL 
plotTF <- T  # whether plot XCO2 and SIF

for (s in 1 : length(site.list)) {

    site <- site.list[s]
    cat(paste('\n# ----- working on overpass on', site, '----- #\n'))

    # call get.site.info() to get lon.lat and OCO2 overpasses info
    # PLEASE add lat lon info in 'get.site.track'
    lon.lat <- get.lon.lat(site, dlon = 2, dlat = 2)
    print(lon.lat)

    # if found txtfile, no need to search or plot again, unless one wants to
    txtfile  <- paste0('oco2_overpass_', site, '_', oco2.ver, '.txt')
    searchTF <- !file.exists(file.path(txtpath, txtfile))

    # get OCO-2 overpasses
    site.info <- get.site.track(site, oco2.ver, oco2.path, searchTF,
                                date.range, thred.count.per.deg = 100, 
                                lon.lat, urbanTF = T, dlon.urban, dlat.urban, 
                                thred.count.per.deg.urban = 50, txtpath)

    # remove overpasses with target modes
    oco2.track  <- site.info #%>% filter(tot.count < 2000)
    ntrack      <- nrow(oco2.track)
    cat(paste('Found total', ntrack, 'overpasses\n'))

    # one can further subset 'oco2.track' based on quality flag
    # add QF filtering (at least 50 screened soundings)
    screen.oco2.track <- oco2.track %>% filter(qf.urban.count > 50)
    ss.oco2.track     <- oco2.track %>% filter(qf.urban.count > 100)
    all.timestr       <- screen.oco2.track$timestr

    mm <- NULL
    if (plotTF & nrow(screen.oco2.track) > 0) {
        # create plotting directory
        plotdir <- file.path(homedir, 'lin-group5/wde/github/stilt', 
                            gsub(' ', '', lon.lat$regid), site, 'plot', oco2.ver)
        dir.create(plotdir, showWarnings = F, recursive = T)
        mm <- ggplot.map(map = 'ggmap', zoom = 8, center.lon = lon.lat$citylon, 
                         center.lat = lon.lat$citylat, api.key = api.key)

        # plot SIF and XCO2
        for (t in 1 : nrow(screen.oco2.track)) {
            timestr <- all.timestr[t]
            ggmap.obs.xco2(site, timestr, oco2.ver, oco2.path, lon.lat, workdir,
                           plotdir, zoom = 8, qfTF = T, box.dlat = dlat.urban, 
                           box.dlon = dlon.urban)
            ggmap.obs.sif (site, timestr, sif.path, lon.lat, workdir, plotdir, 
                           zoom = 8, box.dlat = dlat.urban, box.dlon = dlon.urban)
        }  # end for t
    }  # end if plotTF

    ### -------------------- grab anthromes based on lon.lat ----------------- #
    atm.info <- ggmap.anthromes(lon.lat, anthro.path, mm, site, picpath = plotdir) 

    # define bio as 'cropland' and 'woodland' and estimate vegetation cover %
    new.atm <- atm.info %>% mutate(bio.flag = grepl('woodlands', type ) | 
                                              grepl('croplands', type))
    veg.frac <- length(which(new.atm$bio.flag)) / nrow(new.atm) * 100

    ### -------------------- grab and plot population data ------------------- #
    pop.info <- ggmap.gpw(lon.lat, pop.path, mm, site, picpath = plotdir)

    # select pop density over urban region
    urban.pop <- pop.info %>% filter(lon >= (lon.lat$citylon - dlon.urban),  
                                     lon <= (lon.lat$citylon + dlon.urban), 
                                     lat >= (lon.lat$citylat - dlat.urban),  
                                     lat <= (lon.lat$citylat + dlat.urban)) %>% 
                              na.omit()

    # store all info
    tmp.info <- data.frame(lon.lat, count.all = ntrack, 
                           count.qf50 = nrow(screen.oco2.track), 
                           count.qf100 = nrow(ss.oco2.track), 
                           veg.frac = veg.frac, tot.pop = sum(urban.pop$pop), 
                           stringsAsFactors = F)
    info <- rbind(info, tmp.info)
}  # end for s

write.table(info, file = fn, sep = ',', row.names = F, quote = F)

# end of script
