
# subroutine to load EDGAR emissions
# always grab EDGAR emissions even using ODIAC option (since we need to extract the emission ratio)
load_eco = function(invent = 'EDGAR', emiss_fn, xmin, xmax, ymin, ymax) {
    
    # loading EDGAR emission ------------------------------------------
    if (grepl('.tif', emiss_fn)) {  
        
        # if using pre-processed tif file, emiss unit already in umol m-2 s-1
        emiss_rt = raster(emiss_fn)
        crs(emiss_rt) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

    } else {
        
        # convert kg-CO m-2 s-1 to umol-CO(or NO) m-2 s-1, DW, 03/27/2022
        # correct for molecular weights of CO in EDGAR
        dat = nc_open(emiss_fn)

        # center lat/lon, create our own to avoid floating numbers
        lon = seq(0.05, 359.95, 0.1)
        lat = seq(-89.95, 89.95, 0.1)
        #lon = ncvar_get(dat, 'lon') - ncvar_get(dat, 'lon')[1]      
        #lat = ncvar_get(dat, 'lat') - ncvar_get(dat, 'lon')[1]    

        # standard_name: tendency_of_atmosphere_mass_content_of_carbon_monoxide_due_to_emission
        emiss_mtrx = ncvar_get(dat, 'emi_co') * 1E3 / 28 * 1E6
        dimnames(emiss_mtrx) = list(lon, lat)
        nc_close(dat)

        # convert to raster and crop emissions
        emiss_rt = reshape2::melt(emiss_mtrx) %>% 
                   rename(lon = Var1, lat = Var2, eco = value) %>% 
                   mutate(lon = ifelse(lon >= 180, lon - 360, lon)) %>% 
                   rasterFromXYZ(., crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
    }   # end if
    
    emiss_rt = crop(emiss_rt, extent(xmin, xmax, ymin, ymax))  
    return(emiss_rt)
}  



if (F) {
    # convert kg-CO m-2 s-1 to umol m-2 s-1
    emiss_rt = raster(emiss_fn) * 1E3 / 28 * 1E6

    # fix (0, 360) to (-180, 180) for EDGAR's longitude
    if (extent(emiss_rt)[2] > 181) {
        east = crop(emiss_rt, extent(0, 180, -90, 90))
        west = crop(emiss_rt, extent(180, 360, -90, 90))

        # then change extent of west to negative long
        extent(west) = c(-180, 0, -90, 90)
        emiss_rt = merge(west, east)
    }

    # select EDGAR emissions
    emiss_cp = crop(emiss_rt, extent(xmin, xmax, ymin, ymax))  

}