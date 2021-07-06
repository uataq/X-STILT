# script to plot urban plume, obs and latitude series, DW, 03/04/2021

plot.urban.plume = function(site, site_lon, site_lat, sensor, sensor_gas, 
                            recp_box, recp_info, sel_traj, densf, obs_df, 
                            plm_df, intersectTF, bg_df, bg_side = 'south', 
                            bg.dlat, bg.dlon, map, td, font.size = rel(0.9), 
                            pp_fn = NULL ) {
    
    # -------------------------------------------------------------------------
    uni_times = unique(substr(obs_df$time_utc, 1, 10))
    uni_times = uni_times[!is.na(uni_times)]
    if (length(uni_times) > 1) {
      time_string = paste('during', min(uni_times), '-', max(uni_times))
    } else time_string = paste('on', uni_times)

    if (length(unique(recp_info$recp.lat)) > 1) 
      cat('***Found trajec released from multiple locations...\n')
      
    # -------------------------------------------------------------------------    
    # plot trajectories within overpass duration and observations
    lab.kde = c(td, seq(0.1, 1, 0.1))
    title = paste('Urban plume with', sensor, sensor_gas, 'data over', site, time_string)
    p1 = map + labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) + 
         geom_point(data = sel_traj, aes(lon, lat, colour = dens.level),
                    size = 0.1, alpha = 0.2) + 
         
         # plot normalized 2D kernel density (between 0 and 1) as contours
         geom_contour(data = densf, aes(lon, lat, z = norm.prob, colour = ..level..), 
                      breaks = lab.kde, size = 1.3) +
         scale_colour_gradient(name = 'Normalized\nKernel\nDensity',
                               low = 'lightblue', high = 'purple', 
                               breaks = lab.kde, labels = lab.kde,
                               limits = c(0, max(lab.kde))) + 
        
         # plot receptor box in green 
         geom_polygon(data = recp_box, aes(lon, lat), color = 'darkgreen', 
                      size = 1, linetype = 3, fill = NA) 

    # -------------------------------------------------------------------------
    # plot observations as polygons 
    # and the outmost curve of interpolated urban plume with black line
    if (grepl('OCO', sensor)) { alphas = c(0.5, 0.9); unit = ' [ppm]' }
    if (grepl('TROPOMI', sensor)) { alphas = c(0.4, 0.8); unit = ' [ppb]' }
    p2 = p1 + geom_polygon(data = obs_df, aes(lons, lats, fill = val, 
                                              group = polygon, alpha = plmTF), 
                           colour = 'gray60', size = 0.2) +
              scale_alpha_manual(name = NULL, values = alphas) + guides(alpha = F) + 
              scale_fill_gradientn(colours = def.col(), 
                                   name = paste0('X', sensor_gas, unit)) + 
              geom_polygon(data = plm_df, aes(X, Y), colour = 'gray10', 
                           fill = NA, size = 0.9, alpha = 0.5) 
    
    if (!intersectTF) return(p2) 

    # ------------------------------------------------------------------------- 
    # if there is an intersection between the satellite soundings and urban plume
    # plot overlap polygon and polluted obs (only screened data), DW, 08/20/2018
    p3 = p2 + annotate('text', x = unique(recp_info$recp.lon), label = site, 
                               y = unique(recp_info$recp.lat) - 0.05, size = 5) +
              annotate('point', x = recp_info$recp.lon, size = 2, 
                                y = recp_info$recp.lat, shape = 17) + 
              theme(legend.position = 'right', legend.key = element_blank(),
                    legend.key.width = unit(0.5, 'cm'),
                    legend.key.height = unit(1.2, 'cm'),
                    legend.text = element_text(size = font.size),
                    axis.title.y = element_text(size = font.size, angle = 90),
                    axis.title.x = element_text(size = font.size, angle = 0),
                    axis.text = element_text(size = font.size),
                    axis.ticks = element_line(size = font.size),
                    title = element_text(size = font.size))
    
    # if there is background, draw the background regions and display background values 
    if (!is.null(bg_df)) {
      
      # obs within and outside the urban plume as indicated by the convex hull
      blist = select.obs.side(obs_df, sensor, bg.side = bg_side, bg.dlat, bg.dlon)
      obs_bg = blist$obs_bg
      p4 = p3 + geom_polygon(data = obs_bg, aes(lons, lats, fill = val, group = polygon), 
                             color = 'gray20', alpha = 0.3, size = 0.2) 
      
      # ------------------------------------------------------------------------- 
      # plot enhanced concentrations, DW, 05/25/2021
      # ------------------------------------------------------------------------- 
      # get background stats for the right side
      bg_df = bg_df[bg_df$bg.side == bg_side, ]
      bg_na = mean(bg_df$bg.median) # bg for places outside plume and over side other than bg_side

      # merge obs and bg stats, gap fill bg, and calc ffgas; now as 'obs_rev'
      obs_rev = blist$obs_df %>% left_join(bg_df, by = c('swath', 'bin'))
      obs_rev[is.na(obs_rev$bg.median), 'bg.median'] = bg_na
      obs_rev = obs_rev %>% mutate(ff = val - bg.median)

      # sanity check, check mean ff in background region; should be close to 0
      #ply_bg = unique(obs_bg$polygon)
      #chk_bg = obs_rev %>% filter(polygon %in% ply_bg); mean(chk_bg$ff) 
      
      # also plot enhancements 
      title = paste0('Urban plume with ', sensor, ' \u0394', sensor_gas, 
                     ' over ', site, ' ', time_string, ' (', bg_side, ')')
      if (sensor_gas == 'CO' ) { minz = -10; dz = 10; ts = 1.7; signum = 2; bg.col = 'gray20'}
      if (sensor_gas == 'NO2') { minz = -0.2; dz = 0.2; ts = 1.3; signum = 1; bg.col = 'gray30'}
      if (sensor_gas == 'CO2') { minz = -1; dz = 1; bg.col = 'gray40'}

      e1 = map + labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) +
           geom_polygon(data = obs_rev, aes(lons, lats, fill = ff, 
                                           group = polygon, alpha = plmTF), 
                        colour = 'gray60', size = 0.2) +
           geom_polygon(data = plm_df, aes(X, Y), colour = 'gray10', 
                        fill = NA, size = 0.9, alpha = 0.5) +
           geom_polygon(data = obs_bg, aes(lons, lats, group = polygon), 
                        fill = NA, color = bg.col, size = 0.3) +
           scale_alpha_manual(name = NULL, values = alphas) + guides(alpha = F) + 
           scale_fill_gradientn(colours = def.col(), 
                                limits = c(minz, max(obs_rev$ff)), 
                                breaks = seq(minz, max(obs_rev$ff), dz), 
                                labels = seq(minz, max(obs_rev$ff), dz), 
                                name = paste0('\u0394X', sensor_gas, unit)) + 
           theme(legend.position = 'right', legend.key = element_blank(),
                  legend.key.width = unit(0.5, 'cm'),
                  legend.key.height = unit(1., 'cm'),
                  legend.text = element_text(size = font.size),
                  axis.title.y = element_text(size = font.size, angle = 90),
                  axis.title.x = element_text(size = font.size, angle = 0),
                  axis.text = element_text(size = font.size),
                  axis.ticks = element_line(size = font.size),
                  title = element_text(size = font.size))
    
      # add power plant locations
      if (!is.null(pp_fn)) {

        pp_df = read.csv(pp_fn) %>% mutate_if(is.factor, as.character) %>% 
                filter(longitude >= min(map$data$lon), longitude <= max(map$data$lon), 
                       latitude  >= min(map$data$lat), latitude  <= max(map$data$lat))

        if (nrow(pp_df) > 0) {
          num_fuel = length(unique(pp_df$primary_fuel))
          vals = c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12)[1:num_fuel]
          e1 = e1 + geom_point(data = pp_df %>% filter(capacity_mw < 1000), 
                              aes(longitude, latitude, shape = primary_fuel), 
                              color = 'gray90', size = 2) + 
                    geom_point(data = pp_df %>% filter(capacity_mw > 1000), 
                              aes(longitude, latitude, shape = primary_fuel), 
                              color = 'gold', size = 2) + 
                    scale_shape_manual(name = 'Fuel', values = vals) + 
                    theme(legend.key = element_rect(fill = 'gray10')) + 
                    geom_label(data = pp_df %>% filter(capacity_mw >= 3200, 
                                                       primary_fuel == 'Coal'), 
                               aes(longitude, latitude + 0.05, label = capacity_mw), 
                               color = 'black', fill = alpha(c('white'), 0.4), size = 1.5)
        } else cat('no power stations found...\n')
      } # end if pp_fn

      if (sensor_gas != 'CO2') 
        e1 = e1 + geom_text(data = obs_rev, aes(lon, lat, label = signif(ff, signum)), size = ts)
      
      pp = list(total = p4, delta = e1)
    } else pp = p3

    return(pp)
}
