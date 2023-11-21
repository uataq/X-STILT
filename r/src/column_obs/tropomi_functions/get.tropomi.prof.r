# function to find out the AK profiles from TROPOMI, 
# will find the most adjacent TROPOMI grid that receptor falls into 
# DW, 08/25/2020 

# add other species like NO2, DW, 09/02/2020
# if tropomi.fn exists, no need to provide tropomi.path or search for tropomi data

# add CH4, DW, 11/09/2020
get.tropomi.prof = function(receptor, 
                            tropomi.species = c('CO', 'CH4', 'NO2')[1], 
                            tropomi.path = NULL, 
                            tropomi.fn = NULL) {

  if (is.null(tropomi.fn)) {
    timestr = strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
    lon_lat = data.frame(minlon = receptor$long - 1, 
                         maxlon = receptor$long + 1, 
                         minlat = receptor$lati - 1, 
                         maxlat = receptor$lati + 1)
    tropomi.info = find.tropomi(tropomi.path, substr(timestr, 1, 8), lon_lat)
    if (length(tropomi.info$fn) == 0) stop('No TROPOMI file found for this timestr\n')
    tropomi.fn = tropomi.info$fn
  } # end if


  # --------------------- load all TROPOMI grid center lat/lon
  if ( grepl('s5p_l2_ch4', tropomi.fn) ) {
    
    dat = nc_open(tropomi.fn)
    lat = ncvar_get(dat, 'instrument/latitude_center')
    lon = ncvar_get(dat, 'instrument/longitude_center')
    xch4 = ncvar_get(dat, 'target_product/xch4'); id = 1:length(xch4)
    xch4_bc = ncvar_get(dat, 'target_product/xch4_corrected')
    xch4_uncert = ncvar_get(dat, 'target_product/xch4_precision')
    qa = ncvar_get(dat, 'diagnostics/qa_value')
    time_mtrx = ncvar_get(dat, 'instrument/time')
    hsfc = ncvar_get(dat, 'meteo/surface_altitude') # m
    psfc = ncvar_get(dat, 'meteo/surface_pressure') # hPa
    h2o_vcd = ncvar_get(dat, 'side_product/h2o_column')     # molec cm-2 
    dp = ncvar_get(dat, 'meteo/dp')  # hPa
    
    time_df = as.data.frame(t(time_mtrx)) %>% 
              mutate(time_utc = paste0(V1, formatC(V2, width = 2, flag = 0),
                                        formatC(V3, width = 2, flag = 0), 
                                        formatC(V4, width = 2, flag = 0),
                                        formatC(V5, width = 2, flag = 0), 
                                        formatC(V6, width = 2, flag = 0)))

    var_df = data.frame(id = id, lon = lon, lat = lat, hsfc = hsfc, psfc = psfc,
                        time_utc = time_df$time_utc, h2o_vcd = h2o_vcd, 
                        qa = qa, xch4 = xch4, xch4_bc = xch4_bc, 
                        xch4_uncert = xch4_uncert, dp = dp) %>% 
             filter(lon <= 1e36, xch4 < 1e36) %>% 
             mutate(h2o_vcd = h2o_vcd / 6.023e23 * 1e4) # to mole m-2

  } else {

    dat = nc_open(tropomi.fn)
    indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
    lat  = ncvar_get(dat, 'PRODUCT/latitude')   # [across, along]
    lon  = ncvar_get(dat, 'PRODUCT/longitude')

    # [corner, across, along]
    lats = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds')
    dimnames(lat) = dimnames(lon) = list(indx_across_track, indx_along_track)

    var_array = array(  
      data = do.call(cbind, list(lon, lat)), dim = c(dim(lon), 2), 
                     dimnames = c(dimnames(lon), list(c('lon', 'lat')))
      )   # assuming all matrices in the list have equal dimensions
    
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')

  } # end if reading file
 

  # ----------- locate the TROPOMI grid that is most close to X-STILT receptor 
  # average observations according to coordinates of limited receptors
  # calculate the spherical distance between a given obs and all receptors
  # and locate the one receptor that has the shortest distance
  tmp.dist = geosphere::distCosine(var_df[, c('lon', 'lat')], 
                                   data.frame(lon = receptor$long, 
                                              lat = receptor$lati)) # in meters
  cat(paste('get.tropomi.prof(): find a TROPOMI grid;', 
            signif(min(tmp.dist) / 1E3, 3), 'km away from receptor\n'))
  tmp.grid = var_df[tmp.dist == min(tmp.dist), ]


  # now locate the one sounding --------------
  if ( grepl('s5p_l2_ch4', tropomi.fn) ) {
    tmp_id = tmp.grid$id

    # [level, nobs] or [corner, nobs]
    air_vcd = sum(ncvar_get(dat, 'meteo/dry_air_subcolumns')[, tmp_id]) / 6.023e23 * 1e4    # total dry air column from molec cm-2 to mol m-2
    ak = ncvar_get(dat, 'target_product/xch4_column_averaging_kernel')[, tmp_id]
    zlvl = ncvar_get(dat, 'meteo/altitude_levels')[, tmp_id]  # in m
    zmid = zoo::rollmean(zlvl, 2)    # middle point for a layer
    find_lats = ncvar_get(dat, 'instrument/latitude_corners')[, tmp_id]
    find_lons = ncvar_get(dat, 'instrument/longitude_corners')[, tmp_id]
    
    # AK, zlvl, zmid, from TOA to surface
    upper_pres = rev(seq(tmp.grid$psfc - tmp.grid$dp, 
                         tmp.grid$psfc - tmp.grid$dp * length(zmid), 
                         by = -tmp.grid$dp))
    lower_pres = rev(seq(tmp.grid$psfc, 
                         tmp.grid$psfc - tmp.grid$dp * (length(zmid) - 1), 
                         by = -tmp.grid$dp))
    nc_close(dat)

    all_info = list(tropomi_lat = tmp.grid$lat, tropomi_lon = tmp.grid$lon, 
                    tropomi_lats = find_lats, tropomi_lons = find_lons, 
                    ak = ak, lower_pres = lower_pres, upper_pres = upper_pres, 
                    zagl = zlvl, zmid = zmid, tropomi_zsfc = tmp.grid$hsfc, 
                    tropomi_psfc = tmp.grid$psfc, xch4 = tmp.grid$xch4, 
                    xch4_uncert = tmp.grid$xch4_uncert, 
                    h2o_vcd = tmp.grid$h2o_vcd, air_vcd = air_vcd)
   return(all_info)
  } # end if for extracing AK from CH4 files downloaded from SNOR Remote repo
  

  along_indx = which(indx_along_track == tmp.grid$indx_along_track)
  across_indx = which(indx_across_track == tmp.grid$indx_across_track)
  time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
  timestr = format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), 
                    format = '%Y%m%d%H%M%S')[along_indx]
  
  # these matrics all have dims of [across, along]
  qa   = ncvar_get(dat, 'PRODUCT/qa_value')[across_indx, along_indx] 
  hsfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude')[across_indx, along_indx]
  psfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure')[across_indx, along_indx] / 100 # Pa to hPa

  ### combine all OCO-2 vertical profiles and other 1D variables
  find_lat  = lat[across_indx, along_indx]
  find_lon  = lon[across_indx, along_indx]
  find_lats = lats[, across_indx, along_indx]
  find_lons = lons[, across_indx, along_indx]
  
  # also need to add total or tropospheric column density of dry air in mol m-2 
  g = 9.8        	         # m s-2
  Mair = 29 / 1E3        	 # kg mol-1
  air_vcd = 1/ g / Mair * psfc * 100

  # ------------- grabbing column concentrations of mol m-2 and AK profiles 
  if (tropomi.species == 'CO') {
    
    co_vcd = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column')[across_indx, along_indx] # mol m-2 
    co_vcd_uncert = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column_precision')[across_indx, along_indx]
    h2o_vcd = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_total_column')[across_indx, along_indx]

    # Pressure at bottom of layer in hPa, from TOA to the surface
    lower_pres = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/pressure_levels')[, across_indx, along_indx] / 100 
    
    # compute pressure for upper level in hPa, press from TOA to sfc
    upper_pres = c(lower_pres[1] * lower_pres[1] / lower_pres[2], 
                   lower_pres[1 : length(lower_pres) - 1] )	

    # zagl in meters
    zagl = as.numeric(ncvar_get(dat, 'PRODUCT/layer'))  

    # column AK in meter, 
    ak.co = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel')[, across_indx, along_indx]

    # convert column density to mixing ratio in ppb
    xco = co_vcd / air_vcd * 1E9
    xco_uncert = co_vcd_uncert / air_vcd * 1E9

    # store as a list, mol m-2 for VCD or ppb for mixing ratio
    all_info = list(tropomi_lat = find_lat, tropomi_lon = find_lon, 
                    tropomi_lats = find_lats, tropomi_lons = find_lons, 
                    ak = ak.co, lower_pres = lower_pres, 
                    upper_pres = upper_pres, agl = zagl, tropomi_zsfc = hsfc, 
                    tropomi_psfc = psfc, co_vcd = co_vcd, 
                    co_vcd_uncert = co_vcd_uncert, h2o_vcd = h2o_vcd, 
                    air_vcd = air_vcd, xco = xco, xco_uncert = xco_uncert) 

  } else if (tropomi.species == 'CH4') {  # grab CH4
    
    # XCH4 with unit of 1E-9 which is ppb
    xch4 = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')[across_indx, along_indx] # 1E-9 = ppb
    xch4_bc = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')[across_indx, along_indx]
    xch4.se = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_precision')[across_indx, along_indx] # standard error
    h2o_vcd = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_total_column')[across_indx, along_indx]  # mol m-2

    # for pressure coordinate with even p spacing (from Pa to hPa), 12 layers, 13 levels
    dp = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/pressure_interval')[across_indx, along_indx] / 100   # in mb

    psfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure')[across_indx, along_indx] / 100   # in mb
   
    # grab altitudes and their associated AK, from TOA to sfc
    zsfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude')[across_indx, along_indx] # in meter

    zlvl = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/altitude_levels')[, across_indx, along_indx] # in meter

    zmid = zoo::rollmean(zlvl, 2)    # middle point for a layer

    ak = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel')[, across_indx, along_indx]
    
    if ( is.na(dp) ) {
      plvl = lower_pres = upper_pres = NA
    } else {
      plvl = seq(0, psfc, dp)
      lower_pres = plvl[-1]
      upper_pres = plvl[-length(plvl)]
    }

    all_info = list(tropomi_lat = find_lat, tropomi_lon = find_lon, 
                    tropomi_lats = find_lats, tropomi_lons = find_lons, 
                    ak = ak, lower_pres = lower_pres, upper_pres = upper_pres, 
                    zagl = zlvl, zmid = zmid, tropomi_zsfc = zsfc, 
                    tropomi_psfc = psfc, xch4 = xch4, xch4_bc = xch4_bc, 
                    xch4_uncert = xch4.se, h2o_vcd = h2o_vcd, air_vcd = air_vcd)

  } else if (tropomi.species == 'NO2') {  # grab NO2

    # column density of tropospheric NO2 in mol m-2
    no2_vcd_tropo = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')[across_indx, along_indx]
    no2_vcd_tropo_uncert = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision')[across_indx, along_indx]

    # read air mass factor 
    amf_tot  = ncvar_get(dat, 'PRODUCT/air_mass_factor_total')[across_indx, along_indx]
    amf_tropo = ncvar_get(dat, 'PRODUCT/air_mass_factor_troposphere')[across_indx, along_indx]
    
    # Water vapor and water liquid slant column density, in mol m-2
    h2ov_vcd_slant = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_slant_column_density')[across_indx, along_indx]
    h2ol_vcd_slant = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_liquid_slant_column_density')[across_indx, along_indx]

    ### compute 34 TROPOMI pressure levels using hybrid_sigma_pressure_coordinate
    # p(t, k, j, i, l) = ap(k, l) + b(k, l) *  ps(t, j, i); i, j for locations
    # ap: tm5_constant_a; b: tm5_constant_b; ps: surface air pressure
    # l=0 for base of layer, l=1 for top of layer.

    # TM5 hybrid A coefficient at upper and lower interface levels in Pa, dim of [34, 2]
    a = ncvar_get(dat, 'PRODUCT/tm5_constant_a') / 100   # convert to hPa now
    
    # TM5 hybrid B coefficient at upper and lower interface levels, unitless, dim of [34, 2]
    b = ncvar_get(dat, 'PRODUCT/tm5_constant_b')   
    
    # k from surface to top of atmosphere, 0 to 33 
    k = ncvar_get(dat, 'PRODUCT/layer')      
    
    # "TM5 layer index of the highest layer in the tropopause"
    tropo_indx = ncvar_get(dat, 'PRODUCT/tm5_tropopause_layer_index')[across_indx, along_indx]

    # vertices: lower or upper pressure level boundaries, v = 0 for lower level
    v = ncvar_get(dat, 'PRODUCT/vertices')  

    # finally compute pressures for lower and upper level boundaries, in hPa now
    pres = (a + b * psfc); dimnames(pres) = list(c('lower', 'upper'), k)  
    pres_df = as.data.frame(t(pres)) %>% mutate(k = as.numeric(k))

    # calculate the pressures for tropopause, matrix of [34, 2] 
    # for pressure of upper [smaller p] and lower [higher p] boundary
    tropo_lower_pres = (pres_df %>% filter(k == tropo_indx))$lower
    #print(tropo_lower_pres)

    # so for tropospheric atmosphere, use both lower bound of tropopause layer
    # and the lower bound for the surface layer
    # *** DW, 03/09/2021, no VCD of H2O vapor, only slant column density, 
    # ignore H2Ov, which is not a perfect idea, need to fix it in the future ***
    air_vcd_tropo = 1 / g / Mair * (psfc - tropo_lower_pres) * 100
    
    # NO2 AK from TOA to surface, unitless, no need for any conversion
    ak_no2 = ncvar_get(dat, 'PRODUCT/averaging_kernel')[, across_indx, along_indx]
    
    # The tropospheric averaging kernel can be obtained by scaling the kernel by
    # M/Mtrop (see [RD20]) and setting all elements of the kernel to zero above 
    # the tropopause layer, i.e. for l > l_TM5_tp, based on user manual, 2019
    ak_no2_tropo = c( ak_no2[k <= tropo_indx] * amf_tot / amf_tropo, 
                      rep(0, length(k[k > tropo_indx])) )

    # convert column density to mixing ratio in ppb
    no2_tropo = no2_vcd_tropo / air_vcd_tropo * 1E9
    no2_tropo_uncert = no2_vcd_tropo_uncert / air_vcd_tropo * 1E9

    # store as a list
    # reverse TROPOMI NO2 pressure levels, now from TOA to sfc
    all_info = list(tropomi_lat = find_lat, tropomi_lon = find_lon, 
                    tropomi_lats = find_lats, tropomi_lons = find_lons, 
                    ak = rev(ak_no2), ak_tropo = rev(ak_no2_tropo), 
                    lower_pres = rev(pres_df$lower), 
                    upper_pres = rev(pres_df$upper), 
                    tropomi_zsfc = hsfc, tropomi_psfc = psfc, 
                    no2_vcd_tropo = no2_vcd_tropo, 
                    no2_vcd_tropo_uncert = no2_vcd_tropo_uncert, 
                    h2ov_vcd_slant = h2ov_vcd_slant, 
                    h2ol_vcd_slant = h2ol_vcd_slant, 
                    air_vcd = air_vcd, air_vcd_tropo = air_vcd_tropo, 
                    tno2 = no2_tropo, tno2_uncert = no2_tropo_uncert)

  } else stop('get.tropomi.prof(): please add codes for grabbing other TROPOMI species other than CO and NO2\n')

  nc_close(dat)
  all_info      # return both profiles and other retrivals
  
} # end of subroutine
