# subroutine to plot raob stations on map and wind error on timeseries
# by Dien Wu, 12/03/2018

# selTF for selecting certain raob stations, DW
# add the calculation of mean u, v winds, DW, 01/17/2019 

ggplot.raob.uverr <- function(timestr, site, lon.lat, map = NULL, raob.path, 
                              err.path = NULL, err.file = NULL, nhrs = -72, 
                              selTF = F, font.size = rel(1.0), plotTF = T) {
    
    # 1. plot raob stations on a map 
    if (plotTF) {

        # grab RAOB
        raob <- grab.raob(raob.path, timestr, workdir = err.path, nhrs = nhrs, 
                          format = 'fsl', overwrite = F) %>% 
                mutate(loc = paste0(lon,' ', lat))

        if (selTF) raob <- raob %>% filter(lat >= lon.lat$minlat, 
                                           lat <= lon.lat$maxlat, 
                                           lon >= lon.lat$minlon, 
                                           lon <= lon.lat$maxlon)

        loc.count <- as.data.frame(table(raob$loc), stringsAsFactors = F)
        colnames(loc.count) <- c('loc', 'count')
        loc.str = matrix(unlist(strsplit(loc.count$loc, ' ')), ncol = 2, byrow = T)

        norm.count <- 6 * 16 # 3 days, 16 vertical levels
        loc.df <- data.frame(lon = as.numeric(loc.str[,1]), 
                             lat = as.numeric(loc.str[,2]), 
                             count = as.numeric(loc.count$count)) %>% 
                  mutate(frac = count/norm.count)

        if (is.null(map)) 
            map <- ggplot.map(map = 'ggmap', center.lat = lon.lat$citylat, 
                              center.lon = lon.lat$citylon, zoom = 5)[[1]]
                         
        r1 <- map + #geom_path(data = cir, aes(x, y), linetype = 2, colour = 'gray50') +
            geom_point(data = lon.lat, aes(citylon, citylat), size = 5,
                       fill = NA, colour = 'red', shape = 21) +
            geom_point(data = loc.df, aes(lon, lat, colour = frac * 100), 
                       size = 2, shape = 17) + 
            scale_colour_gradient(low = 'yellow', high = 'red',
                                  name = 'RAOB Available\nData Fraction [%]') + 
            labs(x = 'LONGITUDE', y = 'LATITUDE', 
                 title = paste('Map of RAOB stations for', site, 'on', timestr)) + 
            theme(legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
                  legend.key.height = unit(0.5, 'cm'),
                  legend.text = element_text(size = font.size),
                  legend.key = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.title.y = element_text(size = font.size, angle = 90),
                  axis.title.x = element_text(size = font.size, angle = 0),
                  axis.text = element_text(size = font.size),
                  axis.ticks = element_line(size = font.size),
                  title = element_text(size = font.size))

        png <- file.path(err.path, paste0('map_', site, '_', timestr, '.png'))
        ggsave(r1, filename = png, width = 6, height = 6)
    } # end of plotTF

    # 2. plot time series of uverr
    if (!is.null(err.file)) {

        date1 <- as.POSIXct(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
        date2 <- date1 + nhrs * 60 * 60
        min.date <- min(date1, date2)
        max.date <- max(date1, date2)

        err.dat <- read.table(file.path(err.path, err.file), header = T, sep = ',') 
        err.dat <- err.dat %>% 
                   mutate(date = as.POSIXct(as.character(timestr), 
                                            format = '%Y%m%d%H', tz = 'UTC'), 
                          agl.str = ifelse(hgt <= 3000, 'a) 0-3km', 
                                    ifelse(hgt <= 6000, 'b) 3-6km', 
                                    ifelse(hgt <= 10000, 'c) 6-10km', 'd) > 10km'))), 
                          dist = as.numeric(rdist.earth(
                                            x1 = lon.lat[c('citylon', 'citylat')], 
                                            x2 = err.dat[c('lon', 'lat')], 
                                            miles = F))) %>% 

                   # remove wind error outlier, abs() > 40 m/s
                   filter(abs(u.err) <= 40, abs(v.err) <= 40, dist <= 2000) %>% 

                   # select wind errors based on 3-day-time
                   filter(date >= min.date, date <= max.date) %>% na.omit()

        if (selTF) err.dat <- err.dat %>% filter(lat >= lon.lat$minlat, 
                                                 lat <= lon.lat$maxlat, 
                                                 lon >= lon.lat$minlon, 
                                                 lon <= lon.lat$maxlon)
                                  
        # calculate a distance weighted RMSE
        err.stat <- err.dat %>% group_by(agl.str) %>% 
                    dplyr::summarize(mean.uv.met = mean(abs(c(u.met, v.met))), 
                                     mean.u.met = mean(abs(u.met)), 
                                     mean.v.met = mean(abs(v.met)), 

                                     mean.uv.raob = mean(abs(c(u.raob, v.raob))), 
                                     mean.u.raob = mean(abs(u.raob)), 
                                     mean.v.raob = mean(abs(v.raob)), 

                                     siguverr = sqrt(mean(c(u.err^2, v.err^2))), 
                                     siguerr = sqrt(mean(u.err^2)), 
                                     sigverr = sqrt(mean(v.err^2)), 

                                     mbuverr = mean(c(u.err, v.err)), 
                                     mbuerr  = mean(u.err), 
                                     mbverr  = mean(v.err)) %>% 
                    mutate(timestr = timestr, site = site)

        # start plotting
        if (plotTF) {

            title <- paste0('Time series of u,v errors before the overpass time ', 
                            timestr, ' for ', site, '\nRMSE = ', 
                            signif(err.stat$siguverr[1], 3), ' m/s (0-3 km); ', 
                            signif(err.stat$siguverr[2], 3), ' m/s (3-6 km); ', 
                            signif(err.stat$siguverr[3], 3), ' m/s (6-10 km); ', 
                            signif(err.stat$siguverr[4], 3), ' m/s (> 10 km)')
            brk <- sort(c(signif(min(err.dat$dist), 2), 500, 1000, 1500, 2000, 
                            signif(max(err.dat$dist), 4)))
                        
            e1 <- ggplot(data = err.dat) + theme_bw() + 
                  facet_wrap(~agl.str, nrow = 1) +
                  geom_hline(yintercept = 0, linetype = 2) + #ylim(c(-15, 15)) +
                  geom_point(aes(date, u.err, colour = dist), size = 0.5, 
                             position = position_jitterdodge()) + 
                  #geom_boxplot(aes(date, u.err, group = date), notchwidth = 1) + 
                  scale_colour_gradient(low = 'purple', high = 'lightblue', 
                                        limits = c(min(err.dat$dist) - 5, 
                                                   max(err.dat$dist) + 5), 
                                        breaks = brk, labels = brk) + 
                  labs(x = 'DATE (days-back)', y = 'Wind errors [m/s]', 
                       colour = 'Distance between RAOB stations\nand city center [km]', 
                       title = title) + 
                  theme(legend.key.height = unit(0.5, 'cm'),
                        legend.key.width = unit(1, 'cm'), 
                        legend.text = element_text(angle = 45, hjust = 1))
                
            e2 <- ggplot(data = err.dat) + theme_bw() + 
                  facet_wrap(~agl.str, nrow = 1) +
                  geom_hline(yintercept = 0, linetype = 2) + #ylim(c(-15, 15)) +
                  geom_point(aes(date, v.err, colour = dist), size = 0.5, 
                             position = position_jitterdodge()) + 
                  #geom_boxplot(aes(date, v.err, group = date), notchwidth = 1) +
                  labs(x = 'DATE (3-day-back)', y = 'Wind errors [m/s]', 
                       colour = 'Distance between RAOB stations\nand city center [km]') +
                  scale_colour_gradient(low = 'purple', high = 'lightblue', 
                                        limits = c(min(err.dat$dist) - 5, 
                                                   max(err.dat$dist) + 5), 
                                        breaks = brk, labels = brk) + 
                  theme(legend.key.height = unit(0.5, 'cm'),
                        legend.key.width = unit(1, 'cm'))

            ee <- ggarrange(e1, e2, labels = c('u', 'v'), nrow = 2, 
                            common.legend = T, heights = c(1.2, 1.0))
            ree <- ggarrange(r1, ee, ncol = 2, widths = c(0.5, 1))
            err.png <- file.path(err.path, paste0('map_uverr_', site, '_', timestr, '.png'))
            ggsave(ree, filename = err.png, width = 16, height = 6)
        }  # end of plotTF

        return(err.stat)
    } 

} # end of subroutine

