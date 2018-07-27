# subroutine to find met files
# by DW

# updates:
# no need to "cat" met files together,
# simply grab all met files given time strings

find.metfile<-function(timestr, nhrs, metpath, met.format){

  # calculate timestring, given receptor time and nhrs
  # then use timestring to search for met files
  # if simply for interpolating ground height, nhrs = +/- 1
  yr <- as.numeric(substr(timestr, 1, 4))
  mon <- as.numeric(substr(timestr, 5, 6))
  day <- as.numeric(substr(timestr, 7, 8))
  hr <- as.numeric(substr(timestr, 9, 10))

  range.time <- data.frame(weekdayhr(yr=yr, mon=mon, day=day, hr=hr,
                                     runtt=seq(sign(nhrs), nhrs)*60))

  range.timestr <- paste(formatC(range.time$yr, width=4, flag=0), "-",
                         formatC(range.time$mon, width=2,flag=0), "-",
                         formatC(range.time$day, width=2,flag=0), " ",
                         formatC(range.time$hr, width=2,flag=0), ":00:00",
                         sep="")
  range.date <- as.POSIXlt(range.timestr, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  met.pattern <- unique(sort(strftime(range.date, format=met.format)))

  #metfile <- file.path(path=metpath, pattern=met.pattern)

  return(met.pattern)
}
