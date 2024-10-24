# 
plot.bg = function(site, site_lon, site_lat, sensor, sensor_gas, recp_box, 
                   recp_info, sel_traj, densf, obs_df, plm_df, intersectTF, 
                   bg_df, bg_side = NA, bg_deg, bin_deg, map, td, picname, 
                   font.size, pp_fn = NULL) {
    
    # plot map first 
    uni_sides = unique(bg_df$bg.side)
    print(uni_sides)
    if (!is.na(bg_side)) uni_sides = bg_side
    width = 9; height = 9

    for ( bg_side in uni_sides ) {
        pp = plot.urban.plume(site, site_lon, site_lat, sensor, sensor_gas, 
                              recp_box, recp_info, sel_traj, densf, obs_df, 
                              plm_df, intersectTF, bg_df, bg_side, bg_deg, 
                              bin_deg, map, td, font.size, pp_fn) 

        if ( 'list' %in% class(pp) ) {
            e1 = pp$delta 
            picname_delta = gsub('forward_plume', 'forward_plume_delta', picname)
            picname_delta = gsub('.png', paste0('_', bg_side, '.png'), picname_delta)
            
            ggsave(e1, filename = picname_delta, width = width, height = height)
            p1 = pp$total
        } else p1 = pp

        ggsave(p1, filename = picname, width = width, height = height)
    }   # end for

}   # end of function



if (F) {

    # if there is an intersection, plot latitude series
    l1 = plot.bg.3d(site, timestr, obs_df, bg_df)

    # merge map of forward plume and latitude series, DW, 10/30/2018
    pl = ggarrange(plotlist = list(p1, l1), heights = c(2, 1), nrow = 2, labels = c('a)', 'b)'))
  
}
