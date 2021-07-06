#### subroutine to readin ODIAC emissions and then couple with STILT footprint
#' as footprints have already been weighted by AK and PW,
#' just multiple emission with 2D footprint map, originated from 'foot.odiacv3'
#' @author Dien Wu, 09/13/2016

#' @Updates:
#' ADD ODIACv2016, flag 1 for using v2015a, flag 2 for v2016, DW 02/14/2017
#' note that v2015a does not have emission for year 2015
#' ADD PRD ODIAC emissions, DW 03/08/2017
#' ADD TIMES hourly scaling factors for ODIACv2016, DW 03/08/2017
#' Get rid of variable 'odiac.vname',
#'       always preprocess and read ODIAC emission before call this function...

#' version 2 modify based on Ben's code, DW
#'       can work with multiple receptors at a time, now, DW, 06/05/2018
#' fix footprint lat/lon to lower lefts, as Ben uses centered lat/lon
# '      use raster rather than nc_open, DW, 06/19/2018
#' store output contribution map into the same by-d directory,
# '      remove store.path, DW, 07/26/2018
#' 
#' remove foot.path, use full path as foot.file, DW, 07/26/2018
#' fix a minor bug in interpreting footprint filename, DW, 10/11/2018 
#' accommodate the diff in foot.nc filename, using two versions, DW, 01/25/2019
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# simplify txtfile name, DW, 06/28/2020 

# max.yr for max year available 
# only ODIAC is allowed
run.xco2ff.sim = function(site, timestr = '2014100920', vname = 2019, 
                          tiff.path, output_wd, foot.res, xstilt_wd, store.path,
                          nhrs, oco.sensor = NA, oco.ver = NA, met, run_emiss_err, 
                          edgar.file = NA, ffdas.file = NA, max.yr = NULL, 
                          overwriteTF = F){

  library(rgdal)
  
  # grab footprint files and get footprint domain
  foot.path = file.path(output_wd, 'by-id')
  foot.patt = 'X_foot.nc'
  foot.files = list.files(foot.path, foot.patt, recursive = T, full.names = T)
  foot.indx  = grep(signif(foot.res, 3), basename(foot.files))

  # in a previous version of X-STILT, foot.res is shown on the foot nc file, 
  # but, after the refactoring, no more foot.res, 
  # make this if statement, if foot was generated using a previous version
  # DW, 01/25/2019
  if (length(foot.indx) > 0) {
    foot.file = foot.files[foot.indx]
  } else foot.file = foot.files
  
  if (length(foot.file) == 0) {
    stop('run.xco2ff.sim(): NO footprint found...please check by-id...\n')
    return() }

  tmp.foot = stack(foot.file[1])
  foot.extent = extent(tmp.foot)

  # call tif2nc.odiacv2() to subset and get emiss file name
  # get cropped ODIAC emission for the overpass month
  # moved from main script to this subroutine, DW, 10/21/2018
  if (!is.null(max.yr)) max.yr = as.numeric(vname - 1)
  if (as.numeric(substr(timestr, 1, 4)) > max.yr) {
    tmp.timestr = paste0(max.yr, substr(timestr, 5, nchar(timestr)))
    cat(paste('run.xco2ff.sim(): NO data available for', timestr, 
              'use emission in', max.yr, 'instead\n'))
  } else tmp.timestr = timestr 

  emiss.file = tif2nc.odiacv3(site, timestr = tmp.timestr, vname, xstilt_wd, 
                              foot.extent, tiff.path, gzTF = F)


  # txt file name for outputting model results
  txtfile = file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_', 
                                         abs(nhrs), 'hrs_', oco.sensor, oco.ver, 
                                         '_', met, '_odiac', vname, '.txt'))
  if (is.na(oco.sensor)) 
    txtfile = file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_', 
                                           abs(nhrs), 'hrs_',  met, '_odiac', 
                                           vname, '.txt'))

  # add emission error file with absolute emission uncertainty and txtfile
  if (run_emiss_err) {

    # get emission files for ODIAC, FFDAS, EDGAR
    # use year 2008 emissions to calculate absolute emission errors
    odiac.file.2008 = tif2nc.odiacv3(site, timestr = '20081229', vname, xstilt_wd,
                                     foot.extent, tiff.path, gzTF = F)

    # get absolute emission eror that can further be convolved with footprints
    # **** in 0.1 deg resolution, need to generate 0.1 deg res of footprint
    emiss.file = cal.emiss.err(site, timestr, odiac.file.2008, edgar.file,
                               ffdas.file, emiss.file, overwrite = F)

    # txt file name for outputting model results
    txtfile = file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_emiss_err_', 
                                           abs(nhrs), 'hrs_', oco.sensor, oco.ver, 
                                           '_', met, '_odiac', vname, '.txt'))
    
    if (is.na(oco.sensor)) 
      txtfile = file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_emiss_err_', 
                                             abs(nhrs), 'hrs_',  met, '_odiac', 
                                             vname, '.txt'))
  }  # end if run_emiss_err
                                        

  # if cannot find the correct format of nc file for emissions given selected area
  # and return ODIAC file name with path in front
  if (length(emiss.file) == 0) {
    cat('run.xco2.sim(): NO nc file found, check ODIAC tiff file...\n')
    return()

  } else {  ## read in emissions
    emiss.dat = raster(emiss.file)
    emiss.res = res(emiss.dat)[1]
  }  # end if emiss.file

  # from foot.file, get receptor info
  # fix a minor bug in interpreting footprint filename, DW, 10/11/2018 
  nbin = str_count(basename(foot.file[1]), '_') + 1
  receptor = unlist(strsplit(basename(foot.file), '_'))
  receptor = as.data.frame(matrix(receptor, byrow = T, ncol = nbin),
                            stringsAsFactors = F) %>%
              dplyr::select('V1', 'V2', 'V3') %>% 
              
              # mutate_all() convert character to numberic
              mutate_all(funs(as.numeric), colnames(receptor)) %>% 
              rename(timestr = V1, lon = V2, lat = V3)

  order.index = order(receptor$lat)
  receptor  = receptor[order.index, ]
  foot.file = foot.file[order.index]
  if (run_emiss_err) { receptor$xco2.ff.err = NA } else receptor$xco2.ff = NA


  # then loop over each receptor
  for (r in 1 : nrow(receptor)) {

    outfile = file.path(dirname(foot.file[r]),
                        gsub('foot', 'foot_emiss', basename(foot.file[r])))
    if (file.exists(outfile) & overwriteTF == F) {
      xco2.ff.sp = stack(outfile)

    } else {

      # read in footprint
      foot.dat = stack(foot.file[r])
      if (nlayers(foot.dat) > 1) foot.dat = sum(foot.dat)
      #crs(foot.dat) = '+proj=longlat'

      # NOW, foot and emiss should have the same dimension,
      # multiple them to get contribution map of CO2 enhancements
      xco2.ff.sp = raster::overlay(x = emiss.dat, y = foot.dat,
                                    fun = function(x, y){return(x * y)}) # spatial xco2.ff
      
      ### store emission * column footprint = XCO2 contribution grid into .nc file
      # store into the same workding dir
      crs(xco2.ff.sp) = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
      longname = 'XCO2 enhancemnets due to ODIAC emission'; varname = 'XCO2'
      if (run_emiss_err) {
        longname = 'XCO2 error due to ODIAC emission error'; varname = 'XCO2 error'
      }

      # write raster in nc file
      writeRaster(xco2.ff.sp, outfile, overwrite = TRUE, format = 'CDF',
                  varname = varname, varunit = 'PPM', xname = 'lon', yname = 'lat',
                  longname = longname)
      print(outfile)
    } # end if

    # sum the map to get the XCO2 enhancements,
    # note that AK and PW have been incorporated in footprint
    tmp.xco2.ff = sum(getValues(xco2.ff.sp), na.rm = T); print(tmp.xco2.ff)

    if (run_emiss_err) {
      receptor$xco2.ff.err[r] = tmp.xco2.ff
    } else receptor$xco2.ff[r] = tmp.xco2.ff
  }  # end for r

  # finally, write in a txt file
  #if (writeTF) 
  write.table(x = receptor, file = txtfile, sep = ',', row.names = F, quote = F)

  return(txtfile)
} # end of subroutine
