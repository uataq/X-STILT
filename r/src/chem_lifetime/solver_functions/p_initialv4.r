
# initialize NOx from last or initial time stamp, DW, 2021/10/25

# remove tracking [NOx] before and after chem, emiss...etc. 

# only pass in particle info for the current (t) and the last time stamp (t - 1)
p_initialv4 = function(p_curr, p_last = NULL, mx_res = NA, uni_sza, 
                       bg_nox, a_rto) {
    
    library(dplyr)
    if (is.null(p_last)) p_last = p_init_endpointv4(p_curr, bg_nox, a_rto = 1)
    indx_curr = sort(unique(p_curr$indx))
    indx_last = sort(unique(p_last$indx))
    indx_miss = indx_curr[!indx_curr %in% indx_last]
    
    # ------------------------------------------------------------------------ #
    # work on particles that showed up at time = t, but not at time = t - 1
    # need to initialize those missing particles an IC/BC 
    # (treat like the first time stamp, i.e., either grab IC from TM5, or assign as a const), DW, 01/24/2022
    if ( length(indx_miss) > 0 ) {  
        # particle indices not shown up at the last timestemp, 
        # need to initialize those those
        p_curr_miss = p_curr %>% filter(indx %in% indx_miss)
        p_miss = p_init_endpointv4(p_curr_miss, bg_nox, a_rto)

        # and put the initial conditions for the missing particles back
        p_last_sub = p_last[colnames(p_last) %in% colnames(p_miss)]
        p_last = rbind(p_last_sub, p_miss) %>% arrange(indx)
    }   # end if for paricles that showed up for the first time 

    # create grid for mixing, use lower left as coordinate
    if (is.na(mx_res)) stop('p_initialv4(): no resolution for mixing...\n')
    grd_lon = seq(floor(min(p_curr$long)), ceiling(max(p_curr$long)), 
                  by = mx_res / 111)
    grd_lat = seq(floor(min(p_curr$lati)), ceiling(max(p_curr$lati)), 
                  by = mx_res / 111)
    
    # ------------------------------------------------------------------------ #
    # only subset columns that contain [NOx], [NO2], [CO2], [CO], and [NO2-NOx ratio] and mlht (rename as mlht_last for the entrainment module)
    p_ll = p_last %>% dplyr::select(indx, mlht_last = mlht, 
                                    starts_with('p_'), -contains('_OX'))
    
    # merge particle info from last timestep as initial condition for the current timestep, also locate the corresponding lon, lat, SZA
    p_curr = p_curr %>% left_join(p_ll, by = 'indx') %>%
             mutate(gsza = uni_sza[findInterval(psza, uni_sza)], 
                    glon = grd_lon[findInterval(long, grd_lon)], 
                    glat = grd_lat[findInterval(lati, grd_lat)])   
    
    # check if p_curr$p_* has NA values, if so, debug why??
    if (NA %in% p_curr$p_nox_nomix | NA %in% p_curr$p_nox_mix) 
        cat(paste('NA deteched for time = ', t, '\n'))
    
    return(p_curr)
}




# ------------------------------------------------------------------------ #
# get boundary conditions for trajec endpoint
p_init_endpointv4 = function(curr_pex, bg_nox, a_rto) {

    curr_mlht = mean((curr_pex %>% as.data.frame())$mlht)  # in meter

    # initialize NOx if constant bg_nox is prescribed, 
    # b and a for before and after chemical changes
    if (!is.na(bg_nox)) {
        p_last = data.frame(indx = unique(curr_pex$indx),    
                            p_nox_nochem = bg_nox,  # as passive tracer
                            p_nox_nomix = bg_nox,   # NOx without mixing
                            p_nox_mix = bg_nox,     # NOx with parcel mixing
                            mlht = curr_mlht)

    } else   
    
        # initialize particle-level NOx from 1x1 degree TM5 prior profile
        # assuming NO2-NOx ratio being 1 at trajec endpoint (near midnight)
        p_last = curr_pex %>% dplyr::select(indx, p_nox_nochem = no2_ic, 
                                            p_nox_nomix = no2_ic, 
                                            p_nox_mix = no2_ic)

    # initialize a few more variables
    p_last = p_last %>% mutate(p_co2_nomix = 0, p_co2_mix = 0, 
                               p_co_nomix = 0, p_co_mix = 0, 
                               p_ch4_nomix = 0, p_ch4_mix = 0, 
                               p_rto_nomix = a_rto, p_rto_mix = a_rto, 
                               p_no2_nomix = p_rto_nomix * p_nox_nomix, 
                               p_no2_mix = p_rto_mix * p_nox_mix, 
                               mlht = curr_mlht)
    
    return(p_last)
}

