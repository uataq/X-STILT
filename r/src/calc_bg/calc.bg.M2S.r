#' subroutine to compute M2S background based on Silva et al., 
#' @author: Dien Wu, 01/25/2019

#' require: 
#' continent info, oco2 path

calc.bg.M2S <- function(site, all.timestr, oco2.path, oco2.ver, txtfile) {
    
    library(MASS)

    # 4x4 deg box around hotsplot, from Silva and Arellano, 2017
    lon.lat <- get.lon.lat(site, dlat = 2, dlon = 2)

    silva.bg <- NULL
    for (t in 1:length(all.timestr)) {

        # grab observations and calculate the mean and SD
        obs <- grab.oco2(oco2.path, all.timestr[t], lon.lat, oco2.ver) %>% 
               filter(qf == 0)
        tmp.bg <- as.numeric(fitdistr(obs$xco2, 'normal')$estimate[1]) -
                  as.numeric(fitdistr(obs$xco2, 'normal')$estimate[2])
        silva.bg <- c(silva.bg, tmp.bg)
    }  # end for t
    
    silva.bg <- data.frame(timestr = all.timestr, silva.bg = silva.bg)
    write.table(silva.bg, quote = F, row.names = F, sep =',', file = txtfile)

    return(silva.bg)
}
