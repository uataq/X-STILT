# Subroutine to calculate the CO2 contribution from a priori CO2
# what needed are the combined a priori CO2 profiles and the adjusted combine
# AK *PW profiles from "combine.profile" from main script

# Written by Dien WU, 09/16/2016
# clean codes, DW, 06.12/2018

oco2.apriori <- function(combine.profile = combine.profile){

  ak.norm <- combine.profile$ak.norm
  pw <- combine.profile$pw

  # the CO2 contribution from a priori portion is calculated as the product of
  # (I-AK)* CO2.apriori; first grab normalized AK profiles and apriori CO2 profile
  apriori <- combine.profile$oco2.prior
  co2.apriori <- (rep(1, length(ak.norm)) - ak.norm) * apriori
  xco2.apriori <- sum(co2.apriori * pw)

  return(xco2.apriori)
} # end of subroutine
