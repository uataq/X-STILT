

### Source R scripts for Dien's X-STILT, 05/23/2018, DW
# rsc = dir('/home/dienwu/postdoc_proj/NOx/chem_lifetime', 
#            pattern = '.*\\.r$', full.names = T, recursive = T)
# invisible(lapply(rsc, source))
library(raster); library(dplyr); library(ncdf4); library(reshape2)
library(GeoLight)


# -----------------------------------------------------------------------------
load_wrfout_var = function(wrfout_fns = NULL, name = 'no2', level = 1, 
                           ext = NULL, wrfout_dir = NULL, tdTF = F, td = 70) {

    if (is.null(wrfout_fns))
        wrfout_fns = list.files(wrfout_dir, 'wrfout_d0', full.names = TRUE)[-1]
        #wrfo = wrfout_fns[1]; wrf_summary(wrfo, c('no2', 'T', 'U', 'V'))
    if (is.null(wrfout_fns) & is.null(wrfout_dir)) 
        stop('NO info about WRF output...\n')

    # loop over each file except for the first one hour to read gases
    cat(paste('load_wrfout_var(): reading variable', name, 
              'from wrfout...it takes a while\n'))
    
    var_stk = wrf_get_mod(wrfout_fns[1], name = name, as_raster = T, 
                          raster_lev = level)
    for (wrfo in wrfout_fns[-1])  # will remove the first WRF output
        var_stk = stack(var_stk, wrf_get_mod(wrfo, name, as_raster = T, 
                                             raster_lev = level))

    var_hrs_ut = as.POSIXct(gsub(paste0(name, '_'), '', names(var_stk)), 'UTC', 
                            format = '%Y.%m.%d_%H.%M.%S')
    #var_hrs_lt = as.character(var_hrs_ut - 8 * 3600)
    names(var_stk) = var_hrs_ut 

    # convert values from ppm to ppb for no, no2, no3, chem_no, chem_no2, 
    # chem_no3, hono, hno3, hno4, n2o5, PANs and organic nitrate ANs...
    if ( grepl('no', name) | name %in% c('n2o5', 'onit', 'pan', 'o3') ) 
        var_stk = var_stk * 1000  

    # neglect high values,  max limit in ppb
    #td_day = 20; atd = seq(0, td_day, td_day / 10)
    if (tdTF) var_stk[var_stk > td] = td
    if (!is.null(ext)) var_stk = crop(var_stk, ext)
    
    return(var_stk)
}


# -----------------------------------------------------------------------------
# bug fix, DW, 01/20/2020: when calculating lifetime, 
# divide hourly mean NOx at the previous hour (e.g., 00 UTC) 
# over the changes in NOx reported at the current hour (e.g., 01 UTC)
# add air temperature
if (F) {
    min_level = 1
    max_level = 15
    domain = 'd01'
    ext = NULL
}

extract_nox_tend = function(wrfout_dir, min_level = 1, max_level = 8, 
                            domain = 'd01', storeTF = T) {
    
    source('/home/dienwu/postdoc_proj/NOx/wrfout_functions.r')
    library(raster); library(dplyr); library(ncdf4); library(reshape2)
    wrfout_fns = list.files(wrfout_dir, paste0('wrfout_', domain), 
                            full.names = TRUE)
    if (length(wrfout_fns) == 0) stop('mising WRF-chem files...\n')
    
    timestr = substr(basename(wrfout_dir), 1, 10); print(timestr)
    rds_fn = file.path(wrfout_dir, paste0('wrfout_nox_no2_L', min_level, '-L', 
                                          max_level, '_', timestr, '.rds'))
    
    ext = NULL 
    comb_df = NULL
    for (level in min_level : max_level) {

        print(level)
        
        # (1 FOR LAND, 0 FOR WATER)
        lmsk = load_wrfout_var(wrfout_fns[1], 'LANDMASK', level, ext)
        no2 = load_wrfout_var(wrfout_fns, 'no2', level, ext)
        no  = load_wrfout_var(wrfout_fns, 'no', level, ext)
        o3  = load_wrfout_var(wrfout_fns, 'o3', level, ext) # already in ppb
        nox = no2 + no

        # d* are accumulated changes in an hour, with unit in ppb / hr
        # accumulative chemical changes in NO2, NO, NOx, and O3
        # initial hour = zero (no change in NOx)
        # the first rasterlayer will be removed later
        dno2_acc = load_wrfout_var(wrfout_fns, 'chem_no2', level, ext)
        dno_acc  = load_wrfout_var(wrfout_fns, 'chem_no', level, ext)
        do3_acc  = load_wrfout_var(wrfout_fns, 'chem_o3', level, ext)
        dnox_acc = dno2_acc + dno_acc 
        
        # grab total pressure = perturbation pres + state pres in Pa -----------
        pp = load_wrfout_var(wrfout_fns, 'P', level, ext)
        pb = load_wrfout_var(wrfout_fns, 'PB', level, ext)
        pres = (pp + pb) / 100        # convert to mb

        # add sfc and altitude-varying temp in degC
        tsfc = load_wrfout_var(wrfout_fns, 'T2', ext) - 273.15

        # grab tempature, T for perturbation theta, thus theta = T + 300 in K
        #t00 = load_wrfout_var(wrfout_fns, 'T00', level, ext)
        theta = load_wrfout_var(wrfout_fns, 'T', level, ext) + 300

        # since theta = T * (P0 / P)^(R/cp), where P0 is 1000 mb in WRF
        talt = theta * (pres / 1000)^0.286 - 273.15        # now in degC

        # since these rate or delta changes are cumulative, --------------------
        # we calculate difference between each two raster layer 
        # also drop the last hour for the same dim
        no = no[[-nlayers(no)]]    
        o3 = o3[[-nlayers(o3)]]
        nox = nox[[-nlayers(nox)]]
        no2 = no2[[-nlayers(no2)]]
        lmsk = lmsk[[-nlayers(lmsk)]]
        tsfc = tsfc[[-nlayers(tsfc)]]
        talt = talt[[-nlayers(talt)]]
        pres = pres[[-nlayers(pres)]]

        # calculate the diff in NOx tendency ----------------------------------
        # as changes are accumulative, calc diff between each 2 rasterlayers
        rnox = dnox_acc     # just for initialization
        rno2 = dno2_acc
        rno = dno_acc
        ro3 = do3_acc       

        # for diff in accumulative changes between two hours (i.e., within an hour, [ppb hr-1]), start with the first hour interval (which is the second hour boundary)
        for (r in 2 : nlayers(dnox_acc)) {  
            rnox[[r]] = dnox_acc[[r]] - dnox_acc[[r - 1]]
            rno2[[r]] = dno2_acc[[r]] - dno2_acc[[r - 1]]
            rno[[r]] = dno_acc[[r]] - dno_acc[[r - 1]]
            ro3[[r]] = do3_acc[[r]] - do3_acc[[r - 1]]
        }   # end for r

        # since the first hour is skipped, remove the first layer
        rnox = rnox[[-1]]; rno2 = rno2[[-1]]; rno = rno[[-1]]; ro3 = ro3[[-1]]
        names(rnox) = names(rno2) = names(rno) = names(ro3) = names(nox)  
        # --------------------------------------------------------------------


        # --------------------------------------------------------------------
        # merge all data frames, mr for mixing ratio
        mr_no2_df = melt(as.data.frame(no2, xy = T), id.vars = c('x', 'y')) 
        mr_nox_df = melt(as.data.frame(nox, xy = T), id.vars = c('x', 'y')) 
        mr_no_df = melt(as.data.frame(no, xy = T), id.vars = c('x', 'y')) 
        mr_o3_df = melt(as.data.frame(o3, xy = T), id.vars = c('x', 'y')) 

        # rate of changes in ppb hr-1
        rnox_df = melt(as.data.frame(rnox, xy = T), id.vars = c('x', 'y'))
        rno2_df = melt(as.data.frame(rno2, xy = T), id.vars = c('x', 'y'))
        rno_df  = melt(as.data.frame(rno,  xy = T), id.vars = c('x', 'y'))
        ro3_df  = melt(as.data.frame(ro3,  xy = T), id.vars = c('x', 'y'))

        # grab temp sfc, temp at each altitude
        tsfc_df = melt(as.data.frame(tsfc, xy = T), id.vars = c('x', 'y'))
        talt_df = melt(as.data.frame(talt, xy = T), id.vars = c('x', 'y'))
        pres_df = melt(as.data.frame(pres, xy = T), id.vars = c('x', 'y'))
        lmsk_df = melt(as.data.frame(lmsk, xy = T), id.vars = c('x', 'y'))

        ### TO DO 
        colnms = c('x', 'y', 'variable')
        tmp_df = full_join(mr_no2_df %>% rename(no2 = value), 
                           mr_no_df  %>% rename(no = value),  by = colnms) %>% 
                 left_join(mr_nox_df %>% rename(nox = value), by = colnms) %>%
                 left_join(mr_o3_df %>% rename(o3 = value), by = colnms) %>%

                 left_join(rnox_df %>% rename(rate_nox = value), by = colnms)%>%
                 left_join(rno2_df %>% rename(rate_no2 = value), by = colnms)%>%
                 left_join(rno_df %>% rename(rate_no = value), by = colnms) %>%
                 left_join(ro3_df %>% rename(rate_o3 = value), by = colnms) %>%

                 left_join(tsfc_df %>% rename(tsfc = value), by = colnms) %>% 
                 left_join(talt_df %>% rename(temp = value), by = colnms) %>%
                 left_join(pres_df %>% rename(pres = value), by = colnms) %>%
                 left_join(lmsk_df %>% rename(lmsk = value), by = colnms) %>%

                 # also calculate the NO2-NOx ratio
                 mutate(utc = as.numeric(substr(as.character(variable), 13,14)),
                        level = level, ratio = no2 / nox) %>% 
                 rename(date = variable)  

        comb_df = rbind(comb_df, tmp_df)
    }   # end for level


    # rate_C = NET changes in the concentration, ppb hr-1
    # freq_C = frequency in NET changes of [C], ppb / hr / ppb = hr-1
    # tau_C = NET lifetime in hr, always positive for NET LOSS

    # by solving the derivative for r = [NO2] / [NOx], DW, YH 09/23/2021
    # dr / dt / r = d[NO2] / dt / [NO2] - d[NOx] / dt / [NOx]
    # freq[r] = freq[NO2] - freq[NOx]
    comb_df = comb_df %>% mutate(freq_nox = rate_nox / nox,  # hr-1
                                 freq_no2 = rate_no2 / no2,    
                                 freq_o3 = rate_o3 / o3, 
                                 ts_nox = 1 / -freq_nox,   # hr
                                 datestr = gsub('X', '', as.character(date)), 
                                 date = as.POSIXct(datestr, 'UTC', 
                                                  format ='%Y.%m.%d.%H.%M.%S')) 
    
    if (storeTF) saveRDS(comb_df, rds_fn)
    return(comb_df)
}



# -----------------------------------------------------------------------------
extract_voc = function(wrfout_dir, min_level = 1, max_level = 8, 
                       domain = 'd01', storeTF = T) {
    
    source('/home/dienwu/postdoc_proj/NOx/wrfout_functions.r')
    library(raster); library(dplyr); library(ncdf4); library(reshape2)
    wrfout_fns = list.files(wrfout_dir, paste0('wrfout_', domain), 
                            full.names = TRUE)
    if (length(wrfout_fns) == 0) stop('mising WRF-chem files...\n')
    
    # remove the first hour, since concentrations at initial time is IC 
    wrfout_fns = wrfout_fns[-1]
    timestr = substr(basename(wrfout_dir), 1, 10); print(timestr)
    voc_vars = c('eth', 'hc3', 'hc5', 'hc8',        # Alkanes
                 'ol2', 'olt', 'oli', 'iso',        # Alkenes
                 'tol', 'csl', 'xyl',               # Aromatic
                 'hcho', 'ald', 'ket', 'gly', 'mgly', 'dcb',    # Carbonyls
                 'pan', 'tpan', 'onit',             # Org nitrogen
                 'op1', 'op2', 'paa')               # Org peroxides
                 #'ora1', 'ora2')                    # Org acids

    # voc concentrations all in ppmv, temps in K, pressure in mb
    colnms = c('x', 'y', 'date')
    voc_mini = function(fn, varname, level, ext, append_df, colnms) {
        
        rt = load_wrfout_var(wrfout_fns, name = varname, level, ext)
        df = melt(as.data.frame(rt, xy = T), id.vars = c('x', 'y')) 
        colnames(df)[colnames(df) == 'value'] = varname
        colnames(df)[colnames(df) == 'variable'] = 'date'
        append_df = append_df %>% left_join(df, by = colnms)
        return(append_df)
    }   # end function

    for (level in min_level : max_level) {

        print(level)
        rds_fn = file.path(wrfout_dir, paste0('wrfout_voc_L', level, 
                                              '_', timestr, '.rds'))
    
        # grab total pressure = perturbation pres + state pres in Pa -----------
        pp = load_wrfout_var(wrfout_fns, 'P', level, ext = NULL)
        pb = load_wrfout_var(wrfout_fns, 'PB', level, ext = NULL)
        pres = (pp + pb) / 100        # convert to mb
        pres_df = melt(as.data.frame(pres, xy = T), id.vars = c('x', 'y')) %>% 
                  rename(pres = value, date = variable)

        # grab altitude-varying temp in K using perturbation theta
        # since theta = T * (P0 / P)^(R/cp), where P0 is 1000 mb in WRF
        theta = load_wrfout_var(wrfout_fns, 'T', level, ext = NULL) + 300
        talt  = theta * (pres / 1000)^0.286      # in K
        talt_df = melt(as.data.frame(talt, xy = T), id.vars = c('x', 'y')) %>% 
                  rename(temp = value, date = variable)
        tmp_df = talt_df %>% left_join(pres_df, by = colnms)

        # grab VOCs
        voc_df = tmp_df
        for (varname in voc_vars) 
            voc_df = voc_mini(fn, varname, level, ext = NULL, voc_df, colnms)

        voc_df = voc_df %>% 
                 mutate(level = level, 
                        utc = as.numeric(substr(as.character(date), 13, 14)),
                        datestr = gsub('X', '', as.character(date)), 
                        date = as.POSIXct(datestr, 'UTC', 
                                          format ='%Y.%m.%d.%H.%M.%S')) 

        if (storeTF) {
            saveRDS(voc_df, rds_fn)
            cat('extract_voc(): RDS file stored...\n\n')
        }
    }   # end for level

}



extract_eno_noxmr = function(wrfout_dir, min_level = 1, max_level = 20, 
                             domain = 'd01', storeTF = T) {
    
    source('/home/dienwu/postdoc_proj/NOx/wrfout_functions.r')
    library(raster); library(dplyr); library(ncdf4); library(reshape2)
    wrfout_fns = list.files(wrfout_dir, paste0('wrfout_', domain), 
                            full.names = TRUE)
    if (length(wrfout_fns) == 0) stop('mising WRF-chem files...\n')
    
    timestr = substr(basename(wrfout_dir), 1, 10); print(timestr)
    rds_fn = file.path(wrfout_dir, paste0('wrfout_eno_noxmr_L', min_level, 
                                         '-L', max_level, '_', timestr, '.rds'))
    
    ext = NULL 
    eno = load_wrfout_var(wrfout_fns[2], 'E_NO', level, ext)
    lmsk = load_wrfout_var(wrfout_fns[2], 'LANDMASK', level, ext) # 1LAND,0WATER
    lmsk_df = melt(as.data.frame(lmsk, xy = T), id.vars = c('x', 'y'))
    eno_df = melt(as.data.frame(eno, xy = T), id.vars = c('x', 'y')) 

    # --------------------------------------------------------------------
    comb_df = NULL
    for (level in min_level : max_level) {

        print(level)
        no2 = load_wrfout_var(wrfout_fns, 'no2', level, ext)
        no  = load_wrfout_var(wrfout_fns, 'no', level, ext)
        o3  = load_wrfout_var(wrfout_fns, 'o3', level, ext) # already in ppb
        nox = no2 + no
        
        # grab total pressure = perturbation pres + state pres in Pa -----------
        pp = load_wrfout_var(wrfout_fns, 'P', level, ext)
        pb = load_wrfout_var(wrfout_fns, 'PB', level, ext)
        pres = (pp + pb) / 100        # convert to mb
        names(nox) = names(pres) = names(no2)

        # merge all data frames, mr for mixing ratio
        mr_no2_df = melt(as.data.frame(no2, xy = T), id.vars = c('x', 'y')) 
        mr_nox_df = melt(as.data.frame(nox, xy = T), id.vars = c('x', 'y')) 
        mr_no_df = melt(as.data.frame(no, xy = T), id.vars = c('x', 'y')) 
        mr_o3_df = melt(as.data.frame(o3, xy = T), id.vars = c('x', 'y')) 
        pres_df = melt(as.data.frame(pres, xy = T), id.vars = c('x', 'y'))
        
        colnms = c('x', 'y', 'variable')
        tmp_df = full_join(mr_no2_df %>% rename(no2 = value), 
                           mr_no_df  %>% rename(no = value),  by = colnms) %>% 
                 left_join(mr_nox_df %>% rename(nox = value), by = colnms) %>%
                 left_join(mr_o3_df %>% rename(o3 = value), by = colnms) %>%
                 left_join(pres_df %>% rename(pres = value), by = colnms) %>%

                 # also calculate the NO2-NOx ratio
                 mutate(utc = as.numeric(substr(as.character(variable), 13,14)),
                        level = level, ratio = no2 / nox) %>% 
                 rename(date = variable)  

        comb_df = rbind(comb_df, tmp_df)
    }   # end for level

    comb_df = comb_df %>% 
              left_join(lmsk_df %>% rename(lmsk = value) %>% 
                        dplyr::select(-variable), by = c('x', 'y')) %>%
              left_join(eno_df %>% rename(eno = value) %>% 
                        dplyr::select(-variable), by = c('x', 'y')) %>%
              mutate(datestr = gsub('X', '', as.character(date)), 
                     date = as.POSIXct(datestr, 'UTC', 
                                       format = '%Y.%m.%d.%H.%M.%S')) 

    if (storeTF) saveRDS(comb_df, rds_fn)
    return(comb_df)
}


# -----------------------------------------------------------------------------
# other choice for colorTheme: RdBuTheme, PuOrTheme
plot_wrfout = function(var, max_var, min_var, d_var = 5, var_string, 
                       layout = c(4, 6), parTheme = infernoTheme, 
                       names.attr = NULL) {

    at = seq(min_var, max_var, d_var)
    at_lab = seq(min_var, max_var, d_var * 2)

    if (length(which(getValues(var) > max_var)) > 0) {
        var[var > max_var] = max_var
        labels = c(at_lab[-length(at_lab)], paste('>=', max_var))
    } else labels = at_lab

    if (is.null(names.attr)) names.attr = names(var)

    v1 = rasterVis::levelplot(var, maxpixels = 3E5 * nlayers(var), at = at, 
                              par.settings = parTheme, layout = layout, 
                              main = list(var_string, fontface = 1), 
                              xlab = 'LONGITUDE', ylab = 'LATITUDE', 
                              names.attr = names.attr, 
                              colorkey = list(at = at, 
                                              labels = list(at = at_lab, labels = labels))) 
                                              #+ 
         #latticeExtra::layer(sp.lines(map_border, col = 'white', lwd = 0.5))
   
    return(v1)

}

# -----------------------------------------------------------------------------
wrf_get_mod = function (file = file.choose(), name = NA, as_raster = FALSE,
                        raster_crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                        raster_lev = 1, verbose = FALSE)
{

    if (!is.na(name)) {
        if (name == 'time') {
            wrfchem <- ncdf4::nc_open(file)
            if (verbose)
                cat(paste0('reading Times from ', file, '\n'))
            TIME <- ncvar_get(wrfchem, 'Times')
            TIME <- as.POSIXlt(TIME, tz = 'UTC', format = '%Y-%m-%d_%H:%M:%OS',
                optional = FALSE)
            if (verbose)
                cat('returning Times in POSIXct\n')
            return(TIME)
        }
    }

    if (verbose) cat(paste0('reading ', name, ' from ', file, '\n'))
    
    wrfchem <- ncdf4::nc_open(file)
    if (is.na(name)) {
        name <- menu(names(wrfchem$var), title = 'Choose the variable:')
        POL <- ncdf4::ncvar_get(wrfchem, names(wrfchem$var)[name])
        name <- names(wrfchem$var)[name]
    } else POL <- ncvar_get(wrfchem, name)

    if (as_raster) {

        if (length(dim(POL)) >= 5) stop('images with 5D or more not suported')
        if (length(dim(POL)) == 4) {
            cat(paste0('4D images not supported, making a 3D RasterBrick using level ', raster_lev, ' of the file\n'))
            POL <- POL[, , raster_lev, , drop = TRUE]
        }
        if (length(dim(POL)) == 3) POL <- POL[, , raster_lev, drop = TRUE]
        
        lat <- ncdf4::ncvar_get(wrfchem, varid = 'XLAT')
        lon <- ncdf4::ncvar_get(wrfchem, varid = 'XLONG')
        time <- ncdf4::ncvar_get(wrfchem, varid = 'Times')
        r.lat <- range(lat)
        r.lon <- range(lon)
        n.lat <- ncdf4::ncatt_get(wrfchem, varid = 0, attname = 'SOUTH-NORTH_PATCH_END_UNSTAG')$value
        n.lon <- ncdf4::ncatt_get(wrfchem, varid = 0, attname = 'WEST-EAST_PATCH_END_UNSTAG')$value
        n <- length(time)

        if (n == 1) {
            r <- raster::raster(x = t(POL), xmn = r.lon[1], xmx = r.lon[2],
                ymn = r.lat[1], ymx = r.lat[2])
            r <- raster::flip(r, 2)
        }

        if (n > 1) {
            r <- raster::brick(x = aperm(POL, c(2, 1, 3)), xmn = r.lon[1],
                xmx = r.lon[2], ymn = r.lat[1], ymx = r.lat[2])
            r <- raster::flip(r, 2)
        }
        
        raster::crs(r) <- sp::CRS(raster_crs)
        names(r) <- paste(name, time, sep = '_')
        ncdf4::nc_close(wrfchem)
        return(r)

    }   else {
        ncdf4::nc_close(wrfchem)
        return(POL)
    }

}


# -----------------------------------------------------------------------------
# modified from eixport::to_wrf by DW, allow 3D array [lon, lat, time] as input for POL
# *** no emissions with multiple vertical level allowed!! 
to_wrf_mod = function (POL, file = file.choose(), name) {
    
    wrf <- ncdf4::nc_open(file)
    g_atributos <- ncdf4::ncatt_get(wrf, 0)
    VAR <- array(0, c(g_atributos$`WEST-EAST_PATCH_END_UNSTAG`,
                      g_atributos$`SOUTH-NORTH_PATCH_END_UNSTAG`, 
                      1,    # force vertical level as one
                      dim(POL)[3]))
    VAR[,,1,] = POL 

    ncdf4::nc_close(wrf)
    wrf_put(file, name = name, POL = VAR)

}


wrf_put_mod = function (file = file.choose(), name = NA, POL, mult = NA, verbose = FALSE)
{

    if (class(POL[1]) == 'POSIXlt' || class(POL[1]) == 'POSIXt') {
        cat('converting POSIXlt to string\n')
        POL <- format(POL, '%Y-%m-%d_%H:%M:%OS')
        if (name == 'time')
            name <- 'Times'
    }

    if (verbose) {
        if (missing(mult)) {
            cat(paste0('writing ', name, ' to   ', file, '\n'))
        }
        else {
            cat(paste0('writing ', name, ' to   ', file, ' multiplier ',
                mult, '\n'))
        }
    }
    
    wrfchem <- ncdf4::nc_open(file, write = TRUE)
    if (missing(mult)) {
        ncdf4::ncvar_put(wrfchem, varid = name, POL)
    }
    else {
        ncdf4::ncvar_put(wrfchem, varid = name, unlist(lapply(seq_along(mult),
            function(i) {
                POL * mult[i]
            })))
    }
    ncdf4::nc_close(wrfchem)
}
