
# use receptor UTC time to locate the given TCCON obs and wgt profiles
# grabbing AKs, pressure weighting (integration operator), surface press and hgt, and measured water vapor column

# ppt - hf
# ppb - xco, xn2o, 
# ppm - xco2, xch4, xh2o, 
get.tccon.info = function(tccon.fn, receptor, 
                          tccon.species = c('CH4', 'CO', 'CO2', 'H2O', 'N2O')) {
    
    library(ncdf4); library(reshape2)
    dat = nc_open(tccon.fn)
    
    # load time in UTC
    yr = ncvar_get(dat, 'year')
    doy = ncvar_get(dat, 'day')     # days of the year 1 - 366
    hr_dec = ncvar_get(dat, 'hour')     # decimal_hour
    date = as.Date(doy - 1, origin = paste0(yr, '-01-01'))
    time = as.POSIXct(hr_dec * 3600, 'UTC', date)
    all_timestr = format(time, tz = 'UTC', format = '%Y%m%d%H')

    # select time before grabbing AK to reduce the computational time 
    time_indx = which(time == receptor$run_time)
    if (length(time_indx) == 0) cat('NO TCCON sampling found...\n')

    # geometric altitude in m and atmos pressure in hPa
    zasl_m = ncvar_get(dat, 'zobs')[time_indx] * 1e3          
    pres_hpa = ncvar_get(dat, 'pout')[time_indx]         
    xh2o_ppm = ncvar_get(dat, 'xh2o')[time_indx]
    x_gas = ncvar_get(dat, paste0('x', tolower(tccon.species)))[time_indx]
    xap_gas = ncvar_get(dat, paste0('prior_x', tolower(tccon.species)))[time_indx]

    # integration operator: A vector that, when the dot product is taken with a wet mole fraction profile, applies the TCCON column-average integration. This does NOT include the averaging kernel, those must be applied in addition to this vector. The relates_to attribute indicates which Xgas variable to use this operator for when trying to compare against.
    int_op = ncvar_get(dat, 'integration_operator')[, time_indx]

    # now for aprior [levels, timstamps]
    ap_gas_wet = ncvar_get(dat, paste0('prior_', tolower(tccon.species)))[, time_indx]
    ap_h2o_wet = ncvar_get(dat, 'prior_h2o')[, time_indx]

    # now for AKs [levels, timestamps]
    # Median pressure for the column averaging kernels vertical grid, hPa
    ak_pres = ncvar_get(dat, 'ak_pressure')
    ak_varname = paste0('ak_x', tolower(tccon.species))
    ak_norm = ncvar_get(dat, ak_varname)[, time_indx]
    nc_close(dat)
    
    # store 2D array of int operator and AKs --------------------
    wgt_df = data.frame(ak_pres, ak_norm, ap_gas_wet = ap_gas_wet, 
                        ap_h2o_wet = ap_h2o_wet, int_op)
    info = list(species = tccon.species, zasl_m = zasl_m, psfc_hpa = pres_hpa, 
                xh2o_ppm = xh2o_ppm, x_gas = x_gas, xap_gas = xap_gas, 
                wgt = wgt_df)
    return(info)
}




load_tccon_ak_customize = function(tccon.fn, norm_pres_fn, 
                                   timestr = '2020080120') {

    library(ncdf4); library(reshape2)
    dat = nc_open(tccon.fn)
    alt = ncvar_get(dat, 'alt')
    sza = ncvar_get(dat, 'sza')
    
    # now for AKs [sza, alt]
    ethane_ak = ncvar_get(dat, 'ethane_ak')
    propane_ak = ncvar_get(dat, 'propane_ak')
    nc_close(dat)

    # store 2D array of int operator and AKs --------------------
    dimnames(ethane_ak) = dimnames(propane_ak) = list(sza, alt)
    var_name = list(c('ethane_ak', 'propane_ak'))
    var_list = list(ethane_ak, propane_ak)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    wgt_df = melt(var_array) %>% 
             dcast(Var1 + Var2~Var3, value.var = 'value') %>%
             rename(sza = Var1, alt = Var2) 

    # pressure and pressure weighting --------------------
    pres_df = read.csv('/home/dienwu/X-STILT/data/NormalizedPressure.csv', 
                       header = T)
    pwf_df = data.frame(alt = alt[-length(alt)],   # alt serves as a lower bound
                        pwf = abs(diff(pres_df$Normalized.Pressure)))

    wgt_df = wgt_df %>% left_join(pwf_df, by = 'alt')

    if (F) {

        sz = 0.7; library(viridis)
        w1 = ggplot(data = wgt_df) + labs(x = 'weighting profiles') + 
             geom_point(aes(ethane_ak, alt, color = sza, shape = 'ethane'), size = sz) + 
             geom_point(aes(propane_ak, alt, color = sza, shape = 'propane'), size = sz) + 
             geom_point(aes(pwf * 50, alt), shape = 1, color = 'orange', size = sz) + 
             scale_color_viridis(name = 'SZA', direction = -1) + ylim(c(0, 3))

    }

    return(wgt_df)
}

