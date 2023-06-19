
load_eco2 = function(timestr, invent = c('edgar', 'epa', 'odiac'), 
                     eco2_edgar_fn, eco2_odiac_fn = NA, epa_fn = NA, 
                     epa_name = NA, epa_tz = NA, xmin, xmax, ymin, ymax) {
    
    # loading EDGAR emission ------------------------------------------
    if (grepl('.tif', eco2_edgar_fn)) {  
        
        # if using pre-processed tif file for EDGAR, 
        # already converted emiss unit in umol m-2 s-1
        edgar_rt = raster(eco2_edgar_fn)
        crs(edgar_rt) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
    
    } else {
        
        dat = nc_open(eco2_edgar_fn)

        # center lat/lon, create our own to avoid floating numbers
        lon = seq(0.05, 359.95, 0.1)
        lat = seq(-89.95, 89.95, 0.1)
        #lon = ncvar_get(dat, 'lon') - ncvar_get(dat, 'lon')[1]      
        #lat = ncvar_get(dat, 'lat') - ncvar_get(dat, 'lon')[1]    

        # standard_name: tendency_of_atmosphere_mass_content_of_carbon_dioxide_due_to_emission
        emiss_mtrx = ncvar_get(dat, 'emi_co2') * 1E3 / 44 * 1E6
        dimnames(emiss_mtrx) = list(lon, lat)
        nc_close(dat)

        # convert to raster and crop emissions
        edgar_rt = reshape2::melt(emiss_mtrx) %>% 
                   rename(lon = Var1, lat = Var2, eco2 = value) %>% 
                   mutate(lon = ifelse(lon >= 180, lon - 360, lon)) %>% 
                   rasterFromXYZ(., crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
    }   # end if data format

    # if using ODIAC, further load ODIAC --------------------------------
    if (invent == 'odiac') {
        
        # ODIAC unit in tonne carbon/cell (monthly total)
        # select EDGAR emissions and convert unit to umol m-2 s-1
        odiac_rt = raster(eco2_odiac_fn)
        crs(odiac_rt) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        odiac_rt = crop(odiac_rt, extent(xmin, xmax, ymin, ymax)) 

        area_rt = raster::area(odiac_rt) * 1E6    # convert km2 to m2
        mod = Hmisc::monthDays(as.Date(paste0(substr(timestr, 1, 4), '-',
                                              substr(timestr, 5, 6), '-', 
                                              substr(timestr, 7, 8))))

        # unit conversion
        odiac_rt = odiac_rt * 1E6 / 12 * 1E6 / mod / 24 / 3600 / area_rt
    } else odiac_rt = NULL  # end if
    

    # if scaled by EPA -------
    # ad hoc fix for Hunter Power plant, simply shift the EDGAR emission grid eastward by 0.1degree, DW, 03/27/2022
    if (epa_name == 'Hunter' & invent == 'epa') {
        edgar_rt = crop(edgar_rt, extent(xmin - 0.1, xmax + 0.1, ymin, ymax))  
        extent(edgar_rt)[1:2] = extent(edgar_rt)[1:2] + 0.1

    } else edgar_rt = crop(edgar_rt, extent(xmin, xmax, ymin, ymax))  

    if (!is.na(epa_fn) & invent == 'epa') {
        epa_list = correct_emiss_epav2(edgar_rt, epa_name, epa_tz, epa_fn,'CO2')
        epa_df = epa_list$epa_df 
        edgar_pp = epa_list$emiss_pp 

    } else epa_df = edgar_pp = NULL
    

    # returning ------
    emiss_list = list(edgar_rt = edgar_rt, odiac_rt = odiac_rt, 
                      epa_df = epa_df, edgar_pp = edgar_pp)

    return(emiss_list)
}





if (F) {
    
    # always load EDGAR emission even using EPA scaled or ODIAC based
    # convert kg-CO2 m-2 s-1 to umol m-2 s-1
    edgar_rt = raster(eco2_edgar_fn) * 1E3 / 44 * 1E6

    # fix (0, 360) to (-180, 180) for EDGAR's longitude
    if (extent(edgar_rt)[2] > 181) {
        east = crop(edgar_rt, extent(0, 180, -90, 90))
        west = crop(edgar_rt, extent(180, 360, -90, 90))

        # then change extent of west to negative long
        extent(west) = c(-180, 0, -90, 90)
        edgar_rt = merge(west, east)
    }

    # select EDGAR emissions
    edgar_rt = crop(edgar_rt, extent(xmin, xmax, ymin, ymax))  
    
    # select EDGAR emissions
    m1 = ggplot.map(map = 'ggmap', maptype = 'hybrid', zoom = 10, 
                    center.lon = epa_loc$lon, center.lat = epa_loc$lat)[[1]] + 
         coord_cartesian() + 
         geom_raster(data = edgar_rt %>% as.data.frame(xy = T), 
                    aes(x, y, fill = eco2), alpha = 0.5) + 
         geom_point(data = epa_loc, aes(lon, lat), color = 'white')
    ggsave(m1, filename = 'debug_emiss.png')

}