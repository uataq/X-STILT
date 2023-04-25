#' subroutine to grab TROPOMI data, given a spatial domain and a date
#' @author: Dien Wu

### ------------------------------- 
# tropomi module used in main script, config_tropomi() being called by config_xstilt()
config_tropomi = function(timestr, tropomi.speci, tropomi.path, lon.lat) {

    # get TROPOMI info, if tropomi overpass hour differs from OCO, DW, 10/26/2020
    if (!NA %in% tropomi.speci) {
        
        oco.hr = substr(timestr, 9, 10)
        tropomi.fn = NULL 
        for (ss in 1 : length(tropomi.speci)) {

            cat(paste('- Obtaining TROPOMI', tropomi.speci[ss], 'info...\n'))
            tmp.path = tropomi.path[grep(tropomi.speci[ss], tropomi.path)]
            tropomi.info = find.tropomi(tropomi.path = tmp.path, 
                                        timestr = substr(timestr, 1, 8), 
                                        lon.lat = lon.lat)
            
            if (ss == length(tropomi.speci)) {

                # if TRUE, use TROPOMI overpass hour instead of OCO
                # while still keep the same lat/lon soundings of OCO as receptor locations
                # if FALSE, combine trajec simulation as one run, given their same overpass hour
                tropomi.hr = substr(tropomi.info[, c(1, 2)], 10, 11)
                cat(paste('\nOCO overpass hour:', oco.hr, 
                          'UTC; TROPOMI overpass hour:', tropomi.hr[1], 
                          'to', tropomi.hr[2], 'UTC\n'))

                tropomiTF = !as.numeric(oco.hr) %in% seq(as.numeric(tropomi.hr[1]), 
                                                         as.numeric(tropomi.hr[2]), 1)
            }   # end if ss

            tropomi.fn = c(tropomi.fn, tropomi.info$fn)
        }   # end for ss
    } # end if NA 

    return(list(tropomiTF = tropomiTF, tropomi.fn = tropomi.fn))
}
