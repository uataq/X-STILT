
# subroutine to load TCR or EDGAR emissions 
# add an option using ODIAC, but with CO2-NO emission ratio from EDGAR 

load_eno = function(timestr, invent = c('tcr', 'edgar', 'epa')[3], 
                    emiss_fn, epa_fn = NA, epa_name = NA, epa_tz = NA, 
                    xmin, xmax, ymin, ymax, dmn = NULL, dmx = NULL) {

    if (invent == 'tcr') {
        #dmn = paste0(substr(min(run_dates), 1, 7), '-01')
        #dmx = paste0(substr(max(run_dates), 1, 7), '-31')
        tcr_nox = grab_tcr(emiss_fn, xmin - 1, xmax + 1, ymin - 1, ymax + 1, 
                           tmn = dmn, tmx = dmx, tcr.varname = 'nox', 
                           levelTF = F, monTF = T, tformat = '%Y-%m-%d') %>% 
                   mutate(time = as.character(time), date = as.Date(time)) %>% 
                   arrange(desc(lat))

        tcr_rt = raster(xmin = min(tcr_nox$lon), 
                        xmax = max(tcr_nox$lon) + 1.125, 
                        ymin = min(tcr_nox$lat), 
                        ymax = max(tcr_nox$lat) + 1.122, 
                        nrow = length(unique(tcr_nox$lat)), 
                        ncol = length(unique(tcr_nox$lon)), vals = tcr_nox$nox)
        emiss_rt = crop(tcr_rt, extent(xmin, xmax, ymin, ymax))
    
    } else if (invent == 'edgar' | invent == 'epa' | invent == 'odiac') {
        
        if (grepl('.tif', emiss_fn)) {  
            
            # if using pre-processed tif file, already converted emiss unit in umol m-2 s-1, DW, 11/17/2022
            emiss_rt = raster(emiss_fn)
            crs(emiss_rt) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        } else {

            # convert kg-NOx m-2 s-1 to umol-NOx(or NO) m-2 s-1, DW, 03/27/2022
            # correct for molecular weights of NOx in EDGAR
            dat = nc_open(emiss_fn)

            # center lat/lon, create our own to avoid floating numbers
            lon = seq(0.05, 359.95, 0.1)
            lat = seq(-89.95, 89.95, 0.1)
            #lon = ncvar_get(dat, 'lon') - ncvar_get(dat, 'lon')[1]      
            #lat = ncvar_get(dat, 'lat') - ncvar_get(dat, 'lon')[1]    

            # standard_name: tendency_of_atmosphere_mass_content_of_nox_expressed_as_nitrogen_dioxide_due_to_emission
            emiss_mtrx = ncvar_get(dat, 'emi_nox') * 1E3 / 46 * 1E6
            dimnames(emiss_mtrx) = list(lon, lat)
            nc_close(dat)

            emiss_rt = reshape2::melt(emiss_mtrx) %>% 
                       rename(lon = Var1, lat = Var2, enox = value) %>% 
                       mutate(lon = ifelse(lon >= 180, lon - 360, lon)) %>% 
                       rasterFromXYZ(., crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

        }   # end if data format
     
        # ad hoc fix for Hunter Power plant, simply shift the EDGAR emission grid eastward by 0.1degree, DW, 03/27/2022
        if (epa_name == 'Hunter' & invent == 'epa') {
            emiss_rt = crop(emiss_rt, extent(xmin - 0.1, xmax + 0.1, ymin,ymax))
            extent(emiss_rt)[1:2] = extent(emiss_rt)[1:2] + 0.1

        } else emiss_rt = crop(emiss_rt, extent(xmin, xmax, ymin, ymax))  

    } else stop("load_eno(): can only try among 'tcr', 'edgar', 'epa, check your @param invent...\n")


    # hourly EPA data at the facility level -----------------------
    if (!is.na(epa_fn) & invent == 'epa') {
        epa_list = correct_emiss_epav2(emiss_rt, epa_name, epa_tz, epa_fn,
                                       epa_species = 'NOx')
        epa_df = epa_list$epa_df 
        emiss_pp = epa_list$emiss_pp 
        
    } else epa_df = emiss_pp = NULL

    emiss_list = list(emiss_rt = emiss_rt, epa_df = epa_df, emiss_pp = emiss_pp)
    return(emiss_list)
}  




if (F) {
    
    # used to use 30 g/mol for NOx
    emiss_rt = raster(emiss_fn) * 1E3 / 46 * 1E6
    # fix (0, 360) to (-180, 180) for EDGAR's longitude
    if (extent(emiss_rt)[2] > 181) {
        east = crop(emiss_rt, extent(0, 180, -90, 90))
        west = crop(emiss_rt, extent(180, 360, -90, 90))

        # then change extent of west to negative long
        extent(west) = c(-180, 0, -90, 90)
        emiss_rt = merge(west, east)
    }

    # select EDGAR emissions
    m1 = ggplot.map(map = 'ggmap', maptype = 'hybrid', zoom = 10, 
                    center.lon = epa_loc$lon, center.lat = epa_loc$lat)[[1]] + 
         coord_cartesian() + 
         geom_raster(data = emiss_df %>% filter(lon >= -111.5, lon <= -110.5, lat >= 38.8, lat <= 39.6), 
                    aes(lon + 0.1, lat, fill = enox), alpha = 0.5) + 
         geom_point(data = epa_loc, aes(lon, lat), color = 'white')
    ggsave(m1, filename = 'debug_emiss.png')

}