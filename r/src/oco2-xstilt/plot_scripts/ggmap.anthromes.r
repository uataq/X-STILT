# subroutine to read anthromes data, DW, 07/05/2018

ggmap.anthromes <- function(lon.lat, anthro.path, anthro.file = 'anthro2_a2000.nc', 
                            mm = NULL, site = NULL, picpath = NULL, 
                            font.size = rel(1.0), width = 10, height = 11){

  ### nc open ----
  #anthro.dat <- nc_open(anthro.file)
  #atm <- ncvar_get(anthro.dat, 'anthro2_a2000') # 1/12 * 1/12 deg anthromes grid
  # convert to lower left
  #res <- 1/12  # in degree
  #lat <- ncvar_get(anthro.dat, 'lat') - res/2
  #lon <- ncvar_get(anthro.dat, 'lon') - res/2
  #dimnames(atm) <- list(lon, lat)
  #melt.atm <- melt(atm)

  ### raster open ---
  atm <- raster(file.path(anthro.path, anthro.file))
  melt.atm <- raster::as.data.frame(atm, xy = T)
  colnames(melt.atm) <- c('lon', 'lat', 'class')

  sel.atm <- melt.atm %>% filter(lon >= lon.lat$minlon, lon <= lon.lat$maxlon,
                                 lat >= lon.lat$minlat, lat <= lon.lat$maxlat, 
                                 !is.na(class)) 
  ### anthro biomes legend
  #GRID Values = Anthrome classes**
  #-------------------------------
  #value: Anthrome class
  #11: Urban
  #12: Mixed settlements

  #21: Rice villages
  #22: Irrigated villages
  #23: Rainfed villages
  #24: Pastoral villages

  #31: Residential irrigated croplands
  #32: Residential rainfed croplands
  #33: Populated croplands
  #34: Remote croplands

  #41: Residential rangelands
  #42: Populated rangelands
  #43: Remote rangelands

  #51: Residential woodlands
  #52: Populated woodlands
  #53: Remote woodlands
  #54: Inhabited treeless and barren lands

  #61: Wild woodlands
  #62: Wild treeless and barren lands
  #_______________________________

  if (!is.null(mm)) {
    #** Note that a "LABEL" field with these definitions is attached as a table.
    col <- rgb(
      red = c(168, 255, 0, 0, 169, 255, 0, 230, 255, 255, 230, 255, 255, 56,
              165, 211, 178, 218, 225),
      green = c(0, 0, 112, 169, 0, 115, 255, 230, 255, 255, 152, 211, 235, 168,
                245, 255, 178, 242, 225),
      blue = c(0, 0, 255, 230, 230, 223, 197, 0, 115, 190, 0, 127, 175, 0, 122,
              178, 178, 234, 225),
      maxColorValue = 255)

    col.val <- c(
      "11" = col[1], "12" = col[2], "21" = col[3], "22" = col[4], "23" = col[5],
      "24" = col[6], "31" = col[7], "32" = col[8], "33" = col[9],
      "34" = col[10], "41" = col[11], "42" = col[12], "43" = col[13],
      "51" = col[14], "52" = col[15], "53" = col[16], "54" = col[17],
      "61" = col[18], "62" = col[19])

    lab <- c(
      "11" = "Urban", "12" = "Mixed settlements", "21" = "Rice villages",
      "22" = "Irrigated villages", "23" = "Rainfed villages",
      "24" = "Pastoral villages", "31" = "Residential irrigated croplands",
      "32" = "Residential rainfed croplands", "33" = "Populated croplands",
      "34" = "Remote croplands", "41" = "Residential rangelands",
      "42" = "Populated rangelands", "43" = "Remote rangelands",
      "51" = "Residential woodlands", "52" = "Populated woodlands",
      "53" = "Remote woodlands", "54" = "Inhabited treeless and barren lands",
      "61" = "Wild woodlands", "62" = "Wild treeless and barren lands")

    a1 <- mm[[1]] + theme_bw() + coord_equal(1.1) +
      geom_raster(data = sel.atm, aes(lon + mm[[3]], lat + mm[[2]],
                  fill = as.factor(class)), alpha = 0.7)

    a2 <- a1 +
      scale_fill_manual(name = "Land cover", values = col.val, labels = lab) +
      labs(x = "LONGITUDE", y = "LATITUDE") +
      guides(fill = guide_legend(ncol = 3))
      #,title="Land over map from Anthropogenic Biomesv2 [Ellis et al., 2010]")

    a3 <- a2 + theme(legend.position = 'bottom',
                     legend.text = element_text(size = font.size),
                     legend.key = element_blank(), 
                     legend.key.width = unit(2, 'cm'),
                     legend.key.height = unit(0.6, 'cm'),
                     axis.title.y = element_text(size = font.size, angle = 90),
                     axis.title.x = element_text(size = font.size, angle = 0),
                     axis.text = element_text(size = font.size),
                     axis.ticks = element_line(size = font.size),
                     title = element_text(size = font.size),
                     strip.text.x = element_text(size = font.size))

    picname <- paste0('anthromes_', site, '.png')
    picname <- file.path(picpath, picname)
    print(picname)
    ggsave(a3, file = picname, width = width, height = height)
  } # end if is.null

  return(sel.atm)

}  # end of grab.anthromes()
