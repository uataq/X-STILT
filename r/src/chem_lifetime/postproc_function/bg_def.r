
#' @param bg_quad 'left', 'right', 'top', 'bottom'
#' use a broader spatial region for calculating background 

bg_def = function(pr_edgar, pr_epa = NULL, bg_quad, lon_lat, 
                  xbuf = 0.2, ybuf = 0.2, obs_tno2_df = NULL,
                  xmn = NULL, xmx = NULL, ymn = NULL, ymx = NULL) {

     if (is.null(xmn) | is.null(xmx) | is.null(ymn) | is.null(ymx)) {

          if (bg_quad == 'right') {
               xmn = lon_lat$site_lon + xbuf
               xmx = lon_lat$maxlon
               ymn = lon_lat$minlat
               ymx = lon_lat$maxlat

          } else if (bg_quad == 'left') {
               xmn = lon_lat$minlon
               xmx = lon_lat$site_lon - xbuf
               ymn = lon_lat$minlat
               ymx = lon_lat$maxlat

          } else if (bg_quad == 'top') {

               xmn = lon_lat$minlon
               xmx = lon_lat$maxlon 
               ymn = lon_lat$site_lat + ybuf
               ymx = lon_lat$maxlat

          } else if (bg_quad == 'bottom') {
               xmn = lon_lat$minlon
               xmx = lon_lat$maxlon 
               ymn = lon_lat$minlat
               ymx = lon_lat$site_lat - ybuf

          } else if (bg_quad == 'bottom-left') {
               xmn = lon_lat$minlon
               xmx = lon_lat$site_lon - xbuf 
               ymn = lon_lat$minlat
               ymx = lon_lat$site_lat - ybuf

          } else if (bg_quad == 'top-left') {
               xmn = lon_lat$minlon
               xmx = lon_lat$site_lon - xbuf 
               ymn = lon_lat$site_lat + ybuf
               ymx = lon_lat$maxlat

          } else if (bg_quad == 'top-right') {
               xmn = lon_lat$site_lon + xbuf
               xmx = lon_lat$maxlon 
               ymn = lon_lat$site_lat + ybuf
               ymx = lon_lat$maxlat

          } else if (bg_quad == 'bottom-right') {
               xmn = lon_lat$site_lon + xbuf
               xmx = lon_lat$maxlon
               ymn = lon_lat$minlat
               ymx = lon_lat$site_lat - ybuf
          }
     } 

     if (is.null(xmn) | is.null(ymn) | is.null(xmx) | is.null(ymx)) 
          stop('WRONG background side...\n')

     tno2_bg = NA
     if (!is.null(obs_tno2_df)) {
          tno2_bg_df = obs_tno2_df %>% filter(lats >= ymn, lats <= ymx, 
                                              lons >= xmn, lons <= xmx)
          tno2_bg = median(tno2_bg_df$tno2)
     }

     if (!is.null(pr_edgar)) {
          tno2_bg_df = pr_edgar %>% filter(lat >= ymn, lat <= ymx, 
                                           lon >= xmn, lon <= xmx)
          if (is.na(tno2_bg)) tno2_bg = median(tno2_bg_df$obs_tno2)

          # overwrite FF signal 
          pr_edgar = pr_edgar %>% 
          mutate(sim_tno2_ff = sim_tno2_mix - median(tno2_bg_df$sim_tno2_mix),
                 obs_tno2_ff = obs_tno2 - tno2_bg)
     } 
     
     if (!is.null(pr_epa)) {
          tno2_bg_df = pr_epa %>% filter(lat >= ymn, lat <= ymx, 
                                         lon >= xmn, lon <= xmx)
          if (is.na(tno2_bg)) tno2_bg = median(tno2_bg_df$obs_tno2)
          
          pr_epa = pr_epa %>% 
          mutate(sim_tno2_ff = sim_tno2_mix - median(tno2_bg_df$sim_tno2_mix),
                 obs_tno2_ff = obs_tno2 - tno2_bg)
     } 

     return(list(pr_edgar = pr_edgar, pr_epa = pr_epa, xmn = xmn, xmx = xmx, 
                 ymn = ymn, ymx = ymx))
}

