# script to define spatial domain, DW, 07/05/2018

get.lon.lat <- function(site, dlon = NULL, dlat = NULL) {

  # spatial domains placing receptors and city center, help select OCO-2 data
  # in form of 'lon.lat <- c(minlon, maxlon, minlat, maxlat, city.lon, city.lat)'

  # add region codes,  'ME', 'AS', 'CN', 'AF', 'EU', 'US', 'AM'
  reg.code <- c('ME', 'AS', 'CN', 'AF', 'EU', 'US', 'AM')

  # 1. Middle East
  if (site == 'Riyadh')   lon.lat <- c(46, 48,   23.5, 26.0, 46.72, 24.63, 1)
  if (site == 'Medina')   lon.lat <- c(38, 41,   23.5, 25.5, 39.60, 24.46, 1)
  if (site == 'Mecca')    lon.lat <- c(38, 40.5, 20.5, 22.5, 39.82, 21.42, 1)
  if (site == 'Cairo')    lon.lat <- c(30, 32,   29.0, 32.0, 31.23, 30.05, 1)
  if (site == 'Jerusalem')lon.lat <- c(34, 36,   31.0, 33.0, 35.22, 31.78, 1)
  if (site == 'Jeddah')   lon.lat <- c(38.5, 40,   21,   22,   39.19, 21.49, 1)

  if (site == 'Karachi')  lon.lat <- c(66,   68,   24,   26,   67.01, 24.86, 1)
  if (site == 'Tehran')   lon.lat <- c(50.5, 52.5, 35,   37,   51.39, 35.69, 1)
  if (site == 'Istanbul') lon.lat <- c(28.5, 29.5, 40.5, 41.5, 28.98, 41.01, 1) # turkey
  if (site == 'Baghdad')  lon.lat <- c(43,   45,   32,   34,   44.36, 33.31, 1) # iraq
  if (site == 'Ankara')   lon.lat <- c(32,   34,   39,   41,   32.86, 39.93, 1) # turkey

  # India and southeast Asia
  if (site == 'Mumbai')   lon.lat <- c(72,   74.5, 17.5, 20.0, 72.83, 18.98, 2)
  if (site == 'Delhi')    lon.lat <- c(76.5, 77.5, 28,   30,   77.10, 28.70, 2)
  if (site == 'Bangalore')lon.lat <- c(77,   78,   12.5, 14,   77.59, 12.97, 2)
  if (site == 'Hyderabad')lon.lat <- c(77.5, 79.5, 16.5, 18.5, 78.49, 17.39, 2)
  if (site == 'Ahmedabad')lon.lat <- c(71.5, 73.5, 23.5, 25.5, 72.57, 23.02, 2)

  if (site == 'Manila')     lon.lat <- c(120, 121.5, 14.0, 15.0, 120.98, 14.60, 2)
  if (site == 'Jakarta')    lon.lat <- c(106, 107.5, -7.0, -5.5, 106.86, -6.18, 2)
  if (site == 'Bangkok')    lon.lat <- c(100, 101,   13.5, 15.0, 100.50, 13.76, 2)
  if (site == 'HoChiMinh')  lon.lat <- c(106, 107,   10.0, 11.5, 106.63, 10.82, 2)
  if (site == 'Singapore')  lon.lat <- c(103, 104.5, 1.0,   2.0, 103.87,  1.35, 2)
  if (site == 'KualaLumpur')lon.lat <- c(101, 102,   2.5,   4.0, 101.69,  3.14, 2)

  # China
  if (site == 'YRD')     lon.lat <- c(117,   122,   30,   33, 120,    31.5, 3)
  if (site == 'Shanghai')lon.lat <- c(120,   122,   30.5, 33, 121.52, 31.27, 3)
  if (site == 'Nanjing') lon.lat <- c(118,   120,   31,   33, 118.78, 32.06, 3)
  if (site == 'Suzhou')  lon.lat <- c(119.5, 121.5, 30,   32, 120.65, 31.14, 3)

  if (site == 'PRD')    lon.lat <- c(110,   118,   21,   27, 114.11, 22.40, 3)
  if (site == 'Beijing')lon.lat <- c(115.5, 117.5, 39,   41, 116.41, 39.90, 3)
  if (site == 'Tianjin')lon.lat <- c(116.5, 118.5, 38.5, 41, 117.71, 39.00, 3)
  if (site == 'JJJ')    lon.lat <- c(115,   119,   38,   41, 117,    39.50, 3)

  if (site == 'Xian')   lon.lat <- c(107.5, 110.5, 33,   35.5, 108.90, 34.27, 3)
  if (site == 'Lanzhou')lon.lat <- c(102.5, 105,   35,   37.5, 103.80, 36.04, 3)
  if (site == 'Chengdu')lon.lat <- c(102.5, 105,   29,   32,   104.07, 30.57, 3)
  if (site == 'Zhengzhou')lon.lat <- c(112, 116,   33,   36,   113.63, 34.75, 3)

  # Japan and Korean
  if (site == 'Nagoya') lon.lat <- c(136, 138, 34.5, 36.0, 136.91, 35.18, 2)
  if (site == 'Tokyo-Yokohama')lon.lat <- c(139, 141, 35, 36.5, 139.69, 35.69, 2)
  if (site == 'Osaka-Kobe-Kyoto')lon.lat <- c(134.5, 136, 34, 35.5, 135.50, 34.69, 2)
  if (site == 'Seoul') lon.lat <- c(126.5, 127.5, 37, 38.0, 126.98, 37.57, 2)
  if (site == 'Busan') lon.lat <- c(128.5, 129.5, 35, 35.5, 129.08, 35.18, 2)

  # Africa
  if (site == 'Lagos')       lon.lat <- c(2.5, 4.5,     6,    7,   3.38,   6.52, 4)
  if (site == 'Luanda')      lon.lat <- c(12.5, 14,   -10,   -8,  13.29,  -8.84, 4)
  if (site == 'Kinshasa')    lon.lat <- c(14.5, 16,    -5,   -4,  15.27,  -4.44, 4)
  if (site == 'Johannesburg')lon.lat <- c(27,   29,   -27,  -25,  28.05, -26.20, 4)
  if (site == 'CapeTown')    lon.lat <- c(18, 19.5, -33.5, -34.5, 18.42, -33.92, 4)

  # Europe
  if (site == 'Moscow')   lon.lat <- c(36.5, 38.5, 55,   56.5, 37.62, 55.76, 5)
  if (site == 'Paris')    lon.lat <- c( 1.0,  3.5, 48,   50,    2.35, 48.86, 5)
  if (site == 'London')   lon.lat <- c(-1.0,  1.0, 51,   52.5, -0.13, 51.51, 5)
  if (site == 'Madrid')   lon.lat <- c(-4.5, -3.0, 39.5, 41.5, -3.70, 40.42, 5)
  if (site == 'Barcelona')lon.lat <- c( 1.0 , 3.0, 41,   42,    2.17, 41.39, 5)

  if (site == 'Rome')  lon.lat <- c(12.0, 13.0, 41.5, 42.5, 12.50, 41.90, 5)
  if (site == 'Berlin')lon.lat <- c(12.5, 13.5, 52.0, 53.0, 13.40, 52.52, 5)
  if (site == 'Milan') lon.lat <- c( 8.5,  9.5, 45.0, 46.0,  9.19, 45.46, 5)
  if (site == 'Athens')lon.lat <- c(23.5, 24.0, 37.0, 39.0, 23.73, 37.98, 5)
  if (site == 'StPetersburg')lon.lat <- c(29, 31, 59, 61,   30.34, 59.93, 5)

  # US
  if (site == 'Phoenix')lon.lat <- c(-113, -110, 32,   35,  -112.07, 33.45, 6)
  if (site == 'SLC')    lon.lat <- c(-114, -111, 37,   43,  -111.88, 40.75, 6)
  if (site == 'Denver') lon.lat <- c(-109, -101, 37,   43,  -104.88, 39.76, 6)
  if (site == 'LA')     lon.lat <- c(-122, -115, 32,   38,  -118.41, 34.05, 6)
  if (site == 'LV')     lon.lat <- c(-116, -114, 35,   37,  -115.14, 36.17, 6)
  if (site == 'Seattle')lon.lat <- c(-125, -119, 45,   50,  -122.35, 47.61, 6)

  if (site == 'Dallas')     lon.lat <- c(-98,  -96,  31,   34,   -96.80,  32.78, 6)
  if (site == 'Houston')    lon.lat <- c(-98,  -96,  28.5, 30.5, -95.37,  29.76, 6)
  if (site == 'Albuquerque')lon.lat <- c(-108, -105, 32,   36,   -106.61, 35.11, 6)

  # us with more bio
  if (site == 'NY')     lon.lat <- c(-75, -73,   40, 41.5, -74.01, 40.71, 6)
  if (site == 'DC')     lon.lat <- c(-78, -76,   38, 41,   -77.04, 38.91, 6)
  if (site == 'Indy')   lon.lat <- c(-90, -82,   38, 43,   -86.15, 39.77, 6)
  if (site == 'Chicago')lon.lat <- c(-90, -82,   38, 43,   -86.15, 39.77, 6)
  if (site == 'Miami')  lon.lat <- c(-81, -80,   25, 27,   -80.19, 25.76, 6)
  if (site == 'Atlanta')lon.lat <- c(-85, -83.5, 33, 35,   -84.39, 33.75, 6)

  # Canada
  if (site == 'Toronto')  lon.lat <- c(-80,   -78.5, 43,   44,  -79.38, 43.65, 7)
  if (site == 'Montreal') lon.lat <- c(-74.5, -72.5, 45,   46,  -73.57, 45.50, 7)
  if (site == 'Vancouver')lon.lat <- c(-123.5, -122, 48.5, 50, -123.12, 49.28, 7)

  # central & south america
  if (site == 'SaoPaulo')lon.lat <- c(-47.5, -46,   -24, -23.5, -46.63, -23.55, 7)
  if (site == 'Lima')    lon.lat <- c(-77.5, -75.5, -13, -11,   -77.04, -12.05, 7)

  # Australia
  if (site == 'Perth')    lon.lat <- c(115, 117,   -32.5, -31,   115.86, -31.95, 2)
  if (site == 'Sydney')   lon.lat <- c(150, 151.5, -34.5, -33,   151.21, -33.87, 2)
  if (site == 'Brisbane') lon.lat <- c(152, 154,   -28.5, -26.5, 153.03, -27.47, 2)
  if (site == 'Melbourne')lon.lat <- c(144, 146,   -38.5, -37,   144.96, -37.81, 2)


  if (!is.null(dlat) & !is.null(dlon)) {
    lon.lat[1] <- lon.lat[5] - dlon
    lon.lat[2] <- lon.lat[5] + dlon

    lon.lat[3] <- lon.lat[6] - dlat
    lon.lat[4] <- lon.lat[6] + dlat
  }

  return(lon.lat)
}
