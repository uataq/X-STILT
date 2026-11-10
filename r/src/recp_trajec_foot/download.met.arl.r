
# fix timestring, always convert it to character format
# add a line to check if any met files exist; if so, skip it, DW, 10/20/2024
# add n_met_min back, in case users prefer their own fields, DW, 10/22/2024
# add @param download_method, 'ftp' from ARL FTP or 'aws' from NOAA S3 bucket 
#   s3://noaa-oar-arl-hysplit-pds via AWS CLI (no credentials needed), DW, 10/20/2026
download.met.arl = function(timestr, met_file_format, nhrs, met_path, 
                            met = c('wrf27km', 'hrrr', 'gfs0p25')[2], 
                            met_file_tres = c('1 hour', '6 hours', '1 day')[2], 
                            run_trajec = F, n_met_min = NA, selfTF = FALSE, 
                            download_method = c('ftp', 'aws')[1]) {

    # quickly check if met fields exist 
    for (tmp_timestr in timestr) {
        
        if ( nchar(tmp_timestr) == 8) ff = '%Y%m%d'
        if ( nchar(tmp_timestr) == 10) ff = '%Y%m%d%H'
        t_start = as.POSIXct(as.character(tmp_timestr), 'GMT', format = ff)
        met_fns = find_met_files(t_start, n_hours = nhrs, met_path, 
                                 met_file_format, met_file_tres)

        # according to ARL, 1hr for WRF27km, 6hr for HRRR, 1d for GFS
        if (!selfTF) {
            t_end = t_start + nhrs * 3600
            t_min = min(c(t_start, t_end))
            t_max = max(c(t_start, t_end))
            t_hrs = seq(t_min, t_max, by = 'hour')
            t_dys = unique(trunc(t_hrs, units = 'days')) + c(-1, 1) * 3600 * 24
            if (met == 'wrf27km') arl_format = '%Y%m%d%H'
            if (met == 'hrrr') 
                arl_format = ifelse(tmp_timestr < '2019070100', 
                                    '%Y%m%d.%Hz', '%Y%m%d_%H') 
        
            if (met == 'gfs0p25') arl_format = '%Y%m%d'
            t_bound = seq(min(t_dys), max(t_dys), by = met_file_tres)
            t_find_time = unique(t_bound[findInterval(t_hrs, t_bound)])
            t_met = strftime(t_find_time, format = arl_format, tz = 'GMT')
            if (is.na(n_met_min)) n_met_min = length(t_met)
            met_fns = do.call(c, lapply(t_met, 
                                    function(x) met_fns[grepl(x, met_fns)] ))
        }   # end if for downloading files from ARL
        
        
        # downloading met fields if not found or incomplete ---------------
        if ( run_trajec & !selfTF ) {

            arl_path = file.path('ftp://arlftp.arlhq.noaa.gov/archives', met)
            if (met == 'wrf27km') 
                urls = file.path(arl_path, 'avg', substr(t_met, 1, 4), 
                                paste0('wrfout_d01_', t_met, '.ARL'))
            
            if (met == 'gfs0p25') {
                if (tmp_timestr <= '2019050100' & tmp_timestr >= '2016051300') 
                arl_path = gsub('gfs0p25', 'gfs0p25.v1', arl_path)
                urls = file.path(arl_path, paste0(t_met, '_', met))
            }

            if (met == 'hrrr') {
                # ftp://arlftp.arlhq.noaa.gov/archives/hrrr.v1 for data before July 2019
                if ( tmp_timestr < '2019070100' ) {
                    arl_path = gsub('hrrr', 'hrrr.v1', arl_path)
                    urls = file.path(arl_path, paste0('hysplit.', t_met, '.hrrra'))
                } else {
                    hh_end = formatC(as.numeric(substr(t_met, 10, 11)) + 5, 
                                    width = 2, flag = 0) 
                    urls = file.path(arl_path, paste0(t_met, '-', hh_end, '_', met))
                }
            }   # end if

            # NOAA S3 bucket has the same file names as ARL FTP but is organized 
            # as <met>/YYYY/MM/<file>, (wrf27km under avg/YYYY/MM/)
            ftp_urls = urls
            if (download_method == 'aws') {
                if (Sys.which('aws') == '') 
                    stop('download.met.arl(): AWS CLI not found, install it or use download_method = \'ftp\'\n')
                s3_dir = basename(arl_path)

                if (s3_dir == 'gfs0p25.v1') 
                    stop('download.met.arl(): gfs0p25.v1 (before 2019-05) is not on S3, use download_method = \'ftp\'\n')
                if (met == 'wrf27km') s3_dir = file.path(s3_dir, 'avg')
                urls = file.path('s3://noaa-oar-arl-hysplit-pds', s3_dir, 
                                 substr(t_met, 1, 4), substr(t_met, 5, 6), 
                                 basename(urls))
            }

            # only count files whose size matches ARL, partial downloads
            # from a previous interrupted run will be re-downloaded, DW, 10/20/2026
            dest_fns = file.path(met_path, basename(urls))
            okTF = sapply(seq_along(urls), function(i)
                                is.met.complete(dest_fns[i], c(urls[i], ftp_urls[i])))
            met_fns = dest_fns[okTF]

            if ( length(met_fns) < n_met_min ) {
                cat('download.met.arl(): missing or incomplete metfiles...downloading now...\n')

                for (i in which(!okTF)) {
                    
                    # download to a temp file first so an interrupted download
                    # never leaves a partial file under the final name
                    tmp_fn = paste0(dest_fns[i], '.part')
                    if (download_method == 'aws') {
                        status = system2('aws', c('s3', 'cp', shQuote(urls[i]), 
                                                  shQuote(tmp_fn), '--no-sign-request'))
                    } else 
                        status = tryCatch(download.file(urls[i], 
                                          destfile = tmp_fn,
                                          mode = 'wb'),
                                          error = function(e) { print(e); -1 })
                    if (status == 0 & is.met.complete(tmp_fn, c(urls[i], ftp_urls[i]))) {
                        file.rename(tmp_fn, dest_fns[i])
                    } else {
                        unlink(tmp_fn)
                        stop(paste('download.met.arl(): failed to download', urls[i]))
                    }
                }
                met_fns = dest_fns
            }
            
        } else if ( length(met_fns) < n_met_min & run_trajec ) {
            stop('Since you chose your own met fields, insufficient met fields found...\n')
        } # end if

        print(met_fns)
    }

    return(n_met_min)
}


# check if a local met file is complete by comparing its size with the remote
# file size (via a header-only curl request for ftp, or aws s3api for s3), DW, 10/20/2026
# @param urls can have multiple sources of the same file (e.g., s3 then ftp), 
#   the first one that returns a size is used
# if no remote size can be obtained (e.g., offline), fall back to file.exists
is.met.complete = function(fn, urls) {
    if (!file.exists(fn)) return(FALSE)
    for (url in urls) {
        remote_size = get.remote.size(url)
        if (!is.na(remote_size)) return(file.size(fn) == remote_size)
    }
    warning(paste('is.met.complete(): cannot get remote size from', 
                  paste(urls, collapse = ' or '), '; assuming', basename(fn), 
                  'is complete'))
    TRUE
}

# return remote file size in bytes, or NA (with the error message) if it fails
get.remote.size = function(url) {
    if (grepl('^s3://', url)) {
        bucket = sub('^s3://([^/]+)/.*', '\\1', url)
        key = sub('^s3://[^/]+/', '', url)
        args = c('s3api', 'head-object', '--bucket', bucket, 
                 '--key', shQuote(key), '--no-sign-request', '--query', 'ContentLength', '--output', 'text')
        out = suppressWarnings(
                tryCatch(system2('aws', args, stdout = TRUE, stderr = TRUE), 
                         error = function(e) conditionMessage(e)))
        len = grep('^[0-9]+$', trimws(out), value = TRUE)

    } else {
        out = suppressWarnings(
                tryCatch(system2('curl', c('-sSI', '--max-time', '60',
                                           shQuote(url)),
                                 stdout = TRUE, stderr = TRUE),
                        error = function(e) conditionMessage(e)))
        len = grep('^Content-Length:', out, ignore.case = TRUE, value = TRUE)
    }

    if (length(len) == 0) {
        cat('get.remote.size(): failed for', url, '\n', 
            paste(out, collapse = '\n'), '\n')
        return(NA)
    }

    as.numeric(gsub('[^0-9]', '', len[1]))
}
