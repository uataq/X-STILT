# script to sim XCO2.ff

cal.sim.xco2 <- function(foot.file, output, foot.info, particle, emiss.file,
  dmassTF = F, stilt.ver, txtfile, lon.lat = NULL) {

  if (stilt.ver == 1) {

    # reform particle to match Trajecfoot
    ensemble <- particle %>%
      dplyr::select(time = time, lat = lati, lon = long, agl = zagl,
                    zi = zsfc, foot = foot, index = indx) %>% as.matrix()
    foottimes <- c(0, abs(foot.info$nhrs))
    foot <- Trajecfoot(ident = NULL, part = ensemble, foottimes = foottimes,
                       dmassTF = dmassTF, lon.ll = foot.info$xmn,
                       lat.ll = foot.info$ymn, lon.res = foot.info$xres,
                       lat.res = foot.info$yres, coarse = 1, 
                       numpix.x = (foot.info$xmx - foot.info$xmn) / foot.info$xres,
                       numpix.y = (foot.info$ymx - foot.info$ymn) / foot.info$yres)

    # change dims to [lon, lat] for write_footprint()
    foot <- aperm(foot, c(2, 1, 3))

    # Set footprint grid
    glong <- head(seq(foot.info$xmn, foot.info$xmx, by = foot.info$xres), -1)
    glati <- head(seq(foot.info$ymn, foot.info$ymx, by = foot.info$yres), -1)

    # Set footprint metadata and write to file
    write_footprint(foot, output = foot.file, glong = glong, glati = glati,
      projection = '+proj=longlat', xres = foot.info$xres, yres = foot.info$yres,
      time_out = as.numeric(as.POSIXct(output$receptor$run_time, tz = 'UTC')))

  } else {

    # Calculate near-field dilution height based on gaussian plume width
    # approximation and recalculate footprint sensitivity for cases when the
    # plume height is less than the PBL height scaled by veght
    if (hnf_plume)
      particle <- calc_plume_dilution(particle, numpar = max(particle$indx),
        r_zagl = output$receptor$zagl, veght = 0.5)

    # Produce footprint --------------------------------------------------------
    # Aggregate the particle trajectory into surface influence footprints. This
    # outputs a .rds file, which can be read with readRDS() containing the
    # resultant footprint and various attributes
    foot <- calc_footprint(particle, output = foot.file,
                           r_run_time = output$receptor$run_time,
                           smooth_factor = 1, time_integrate = T,
                           xmn = foot.info$xmn, xmx = foot.info$xmx,
                           xres = foot.info$xres, yres = foot.info$yres,
                           ymn = foot.info$ymn, ymx = foot.info$ymx)

  }

  xco2.portion <- foot.odiacv3(foot.file, emiss.file, workdir, txtfile,
    lon.lat, plotTF = F)

  return(xco2.portion)
}
