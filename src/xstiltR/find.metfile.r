# subroutine to find met files


find.metfile<-function(timestr, nhrs, metpath, met, trajecTF){

  pattern <- substr(timestr, 1, 6)

  # for WRF--met files for generating trajec
  if(trajecTF & met=="1km"){
    # requires all nested met fields and nhrs, for WRF
    metfile <- rev(list.files(path=metpath, pattern=pattern))[-6]
  }else{
    # met files for trajwind() to interpolate ground heights,
    # requires receptor lat, lon
    filename <- paste("wrfout_d04_", tolower(site), "_", pattern, sep="")
    if(met == "1km") metfile <- list.files(path=metpath, pattern=filename)
  }

  # for 1deg or 0.5deg GDAS--
  filename <- paste(met, "_", pattern, sep="")
  if(met != "1km") metfile <- list.files(path=metpath, pattern=filename)

  return(metfile)
}
