# dependencies for Trajecmulti()

julian<-function(m, d, y, origin.){
#returns day since 1/1/1960
#
#  $Id: julian.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   only.origin <- all(missing(m), missing(d), missing(y))
   if(only.origin)
           m <- d <- y <- NULL
   # return days since origin
   nonnumeric.p <- !c(is.numeric(m), is.numeric(d), is.numeric(y))
   if(any(nonnumeric.p) && !only.origin) {
           badarg <- paste(c("m", "d", "y")[nonnumeric.p], collapse = ", "
                   )
           stop(paste("Arguments:", badarg, "are not numeric"))
   }
   if(missing(origin.))
           if(is.null(origin. <- .Options$chron.origin))
                   origin. <- c(month = 1, day = 1, year = 1960)
   nms <- names(d)
   max.len <- max(length(m), length(d), length(y))
   #
   # prepend new origin value and rep out to common max. length:
   m <- c(origin.[1], rep(m, length = max.len))
   d <- c(origin.[2], rep(d, length = max.len))
   y <- c(origin.[3], rep(y, length = max.len))
   #
   # code from julian date in the S book (p.269)
   #
   y <- y + ifelse(m > 2, 0, -1)
   m <- m + ifelse(m > 2, -3, 9)
   c <- y %/% 100
   ya <- y - 100 * c
   out <- (146097 * c) %/% 4 + (1461 * ya) %/% 4 + (153 * m + 2) %/% 5 +
           d + 1721119
   #
   # now subtract the new origin from all dates
   #
   if(!only.origin) {
           if(all(origin. == 0))
                   out <- out[-1]
           else out <- out[-1] - out[1]
   }
   names(out) <- nms
   out
}



month.day.year<-function(jul, origin.){
#returns month, day, year
#
#  $Id: month.day.year.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   if(missing(origin.) || is.null(origin.))
           if(is.null(origin. <- .Options$chron.origin))
                   origin. <- c(month = 1, day = 1, year = 1960)
   if(all(origin. == 0))
           shift <- 0
   else shift <- julian(origin = origin.)
   # relative origin
   # "absolute" origin
   j <- jul + shift
   j <- j - 1721119
   y <- (4 * j - 1) %/% 146097
   j <- 4 * j - 1 - 146097 * y
   d <- j %/% 4
   j <- (4 * d + 3) %/% 1461
   d <- 4 * d + 3 - 1461 * j
   d <- (d + 4) %/% 4
   m <- (5 * d - 3) %/% 153
   d <- 5 * d - 3 - 153 * m
   d <- (d + 5) %/% 5
   y <- 100 * y + j
   y <- y + ifelse(m < 10, 0, 1)
   m <- m + ifelse(m < 10, 3, -9)
   list(month = m, day = d, year = y)
}



#---------------------------------------------------------------------------------------------------
#  $Id: unix.r,v 1.4 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

unix<-function(command, intern = TRUE, ignore.stderr = FALSE){
#redefines unix as system() with intern=T as default
if(nchar(Sys.getenv("COMSPEC"))==0) return(system(command = command, intern = intern, ignore.stderr = ignore.stderr)) #unix/linux
if(nchar(Sys.getenv("COMSPEC"))>0) return(system(command = command, intern = intern)) #PC
}


#---------------------------------------------------------------------------------------------------
#  $Id: unix.shell.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

unix.shell<-function(command, shell = "/bin/sh", ...){
        tfile <- tempfile("Sshell")
        on.exit(unlink(tfile))
        cat(file = tfile, command, "\n", sep = "")
        command <- paste(shell, tfile)
        unix(command, ...)
}



day.of.week<-function(month, day, year){
#returns day of week as number (0: Sun, 6: Sat)
#
#  $Id: day.of.week.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

        ix <- year + trunc((month - 14)/12)
        jx <- trunc((13 * (month + 10 - (month + 10) %/% 13 * 12) - 1)/5) +
                day + 77 + (5 * (ix - (ix %/% 100) * 100)) %/% 4 + ix %/% 400 -
                (ix %/% 100) * 2
        jx %% 7
}
