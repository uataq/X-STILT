# default settings for transport error configuration, DW, 07/03/2021 

config_trans_err = function(namelist, site, lon_lat, timestr, xstilt_wd, 
                            simstep_namelist = NULL) {
    
    # Transport and dispersion settings, use default setting in STILT-R v2
    if (!is.null(simstep_namelist)) {
        capemin     = -1
        cmass       = 0
        conage      = 48
        cpack       = 1
        delt        = 2
        dxf         = 1
        dyf         = 1
        dzf         = 0.01
        efile       = ''
        emisshrs    = 0.01
        frhmax      = 3
        frhs        = 1
        frme        = 0.1
        frmr        = 0
        frts        = 0.1
        frvs        = 0.01
        hscale      = 10800
        ichem       = 8
        idsp        = 2
        initd       = 0
        k10m        = 1
        kagl        = 1
        kbls        = 1; if (toupper(namelist$met) == 'NARR') kbls = 2
        kblt        = 5
        kdef        = 0
        khinp       = 0
        khmax       = 9999
        kmix0       = 250
        kmixd       = 3
        kmsl        = 0
        kpuff       = 0
        krand       = 4
        krnd        = 6
        kspl        = 1
        kwet        = 1
        kzmix       = 0
        maxdim      = 1
        maxpar      = 20000
        mgmin       = 10
        mhrs        = 9999
        nbptyp      = 1
        ncycl       = 0
        ndump       = 0
        ninit       = 1
        nstr        = 0
        nturb       = 0
        nver        = 0
        outdt       = 0
        p10f        = 1
        pinbc       = ''
        pinpf       = ''
        poutf       = ''
        qcycle      = 0
        rhb         = 80
        rht         = 60
        splitf      = 1
        tkerd       = 0.18
        tkern       = 0.18
        tlfrac      = 0.1
        tout        = 0
        tratio      = 0.75
        tvmix       = 1
        veght       = 0.5
        vscale      = 200
        vscaleu     = 200
        vscales     = -1
        wbbh        = 0
        wbwf        = 0
        wbwr        = 0
        wvert       = FALSE
        w_option    = 0
        zicontroltf = 0
        z_top       = 25000
        if (toupper(namelist$met) == 'NARR') z_top = 15000 

        # Aggregate STILT/HYSPLIT namelist
        simstep_namelist = list(capemin = capemin, cmass = cmass, conage = conage,
                                cpack = cpack, delt = delt, dxf = dxf, dyf = dyf,
                                dzf = dzf, efile = efile, frhmax = frhmax, 
                                frhs = frhs, frme = frme, frmr = frmr, frts = frts, 
                                frvs = frvs, hnf_plume = hnf_plume, hscale = hscale, 
                                ichem = ichem, idsp = idsp, initd = initd, 
                                k10m = k10m, kagl = kagl, kbls = kbls, kblt = kblt, 
                                kdef = kdef, khinp = khinp, khmax = khmax, 
                                kmix0 = kmix0, kmixd = kmixd, kmsl = kmsl,
                                kpuff = kpuff, krand = krand, krnd = krnd, 
                                kspl = kspl, kwet = kwet, kzmix = kzmix, 
                                maxdim = maxdim, maxpar = maxpar, mgmin = mgmin, 
                                ncycl = ncycl, ndump = ndump, ninit = ninit,
                                nstr = nstr, nturb = nturb, numpar = numpar, 
                                nver = nver, outdt = outdt, p10f = p10f, 
                                pinbc = pinbc, pinpf = pinpf, poutf = poutf, 
                                qcycle = qcycle, rhb = rhb, rht = rht,
                                splitf = splitf, tkerd = tkerd, tkern = tkern, 
                                tlfrac = tlfrac, tout = tout, tratio = tratio, 
                                tvmix = tvmix, varsiwant = varsiwant, veght = veght,
                                vscale = vscale, vscaleu = vscaleu, vscales = vscales, 
                                wbbh = wbbh, wbwf = wbwf, wbwr = wbwr, winderrtf = 0, 
                                wvert = wvert, zicontroltf = zicontroltf)
    }   # end if


    # Transport and dispersion settings, use default setting in STILT-R v2
    # Calculating error stats for horizontal and vertical transport errors
    hor_err = get.uverr(namelist$run_hor_err, site, timestr, xstilt_wd, 
                        run_wind_err = namelist$run_wind_err, 
                        simstep_namelist, 
                        raob.path = namelist$raob_path, 
                        raob.format = 'fsl', 
                        nhrs = namelist$nhrs, 
                        met = namelist$met, 
                        met_path = namelist$met_path, 
                        met_file_format = namelist$met_file_format, 
                        lon.lat = lon_lat, 
                        agl = namelist$agl, 
                        err.path = file.path(namelist$store_path, 'wind_err'))
    print(hor_err)
    pbl_err = get.zierr(namelist$run_ver_err, 
                        nhrs.zisf = abs(namelist$nhrs), 
                        const.zisf = namelist$zisf)

    # and prepare ODIAC based on footprint domain 
    if (namelist$run_hor_err & !is.null(namelist$odiac_path)) {
        foot.ext = extent(lon_lat$minlon, lon_lat$maxlon, lon_lat$minlat, lon_lat$max.at)
        emiss_fn = tif2nc.odiacv3(site, timestr, vname = namelist$odiac_ver, 
                                  workdir = xstilt_wd, foot.ext = foot.ext, 
                                  tiff.path = namelist$odiac_path, gzTF = F)
    } else emiss_fn = NA

    return(list(hor_err = hor_err, pbl_err = pbl_err, emiss_fn = emiss_fn))
}

