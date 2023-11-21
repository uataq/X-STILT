
load_eco2_odiac = function(odiac_fn, timestr, xmin, xmax, ymin, ymax) {

    library(raster)

    # ODIAC unit in tonne carbon/cell (monthly total)
    # select EDGAR emissions and convert unit to umol m-2 s-1
    odiac_rt = raster(odiac_fn)
    crs(odiac_rt) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

    if ( !is.na(xmin) & !is.na(xmax) & !is.na(ymin) & !is.na(ymax) )
        odiac_rt = crop(odiac_rt, extent(xmin, xmax, ymin, ymax)) 

    area_rt = raster::area(odiac_rt) * 1E6    # convert km2 to m2
    mod = Hmisc::monthDays(as.Date(paste0(substr(timestr, 1, 4), '-',
                                          substr(timestr, 5, 6), '-', 
                                          substr(timestr, 7, 8))))

    # unit conversion
    odiac_rt = odiac_rt * 1E6 / 12 * 1E6 / mod / 24 / 3600 / area_rt
    return(odiac_rt)
}