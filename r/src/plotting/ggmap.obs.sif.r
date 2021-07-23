# subroutine to read OCO-2/3 SIF, DW, 06/06/2018
# 'thred.count' for the sounding # thredshold

ggmap.obs.sif <- function(site, timestr, oco.sensor, oco.ver, sif.path, lon.lat,
                          plotdir, zoom = 8){

  if (is.null(sif.path)) {cat('NO SIF path defined, returning NA...\n'); return()}                          
  sel.sif <- grab.sif(sif.path, timestr, lon.lat, oco.ver)
  if (length(sel.sif) == 0) {cat('NO SIF file matched, returning NA...\n'); return()}
  if (nrow(sel.sif) == 0) {cat('NO SIF retrieved within lon.lat, returning NA...\n'); return()}

  # plot center
  font.size <- rel(1.0); col <- def.col()

  # plot google map
  m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
                   center.lon = lon.lat$citylon)[[1]] + theme_bw() + 
        labs(x = 'LONGITUDE', y = 'LATITUDE') + 
        theme(legend.position = 'right', legend.key = element_blank(), 
                   legend.text = element_text(size = font.size),
                   legend.key.height = unit(1.2, 'cm'),
                   legend.key.width = unit(0.5, 'cm'),
                   axis.title.y = element_text(size = font.size, angle = 90),
                   axis.title.x = element_text(size = font.size, angle = 0),
                   axis.text = element_text(size = font.size),
                   axis.ticks = element_line(size = font.size),
                   title = element_text(size = font.size),
                   strip.text = element_text(size = font.size))

  melt.sif <- sel.sif %>% 
              dplyr::select('timestr', 'lat', 'lon', 'sif757', 'sif771', 'avg.sif') %>% 
              melt(id.var = c('timestr', 'lat', 'lon'))

  title <- paste(oco.sensor, 'SIF [W/m2/sr/µm] and IGBP for', site, 'on', timestr)
  c1 <- m1 + geom_point(data = melt.sif, aes(lon, lat, colour = value), size = 0.9) +
             facet_wrap(~variable, nrow = 3) + 
             scale_colour_gradientn(name = paste(oco.sensor, 'SIF'), colours = col,
                                    limits = c(-1, max(2.5, max(melt.sif$value))),
                                    breaks = seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5))

  # land cover 
  l1 <- ggmap.igbp(sel.sif, m1) + labs(title = 'IGBP', x = NULL, y = NULL) + 
        theme(legend.position = 'bottom', legend.key = element_blank(), 
              legend.text = element_text(size = font.size),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(1.2, 'cm'),
              axis.title.y = element_text(size = font.size, angle = 90),
              axis.title.x = element_text(size = font.size, angle = 0),
              axis.text = element_text(size = font.size),
              axis.ticks = element_line(size = font.size),
              title = element_text(size = font.size),
              strip.text = element_text(size = font.size))

  cl <- ggarrange(c1, l1, ncol = 2, widths = c(1, 1.4))
  cl <- annotate_figure(cl, top = title)

  picname <- paste0('ggmap_sif_', site, '_', timestr, '.png')
  picfile <- file.path(plotdir, picname); print(picfile)
  ggsave(cl, filename = picfile, width = 11, height = 9)

  return(cl)
}



# ---------------------------------
ggmap.igbp <- function(df, map, font.size = rel(0.8), scale.coord = 1.2, 
                       legend.ncol = 4) {

    get.igbp.col <- function() {

        val  <- seq(1, 17)
        name <- c('ENF', 'EBF',   # Evergreen Needleleaf/Broadleaf Forest
                  'DNF', 'DBF',   # Deciduous Needleleaf/Broadleaf Forest
                  'MF',           # Mixed forest
                  'CSHR', 'OSHR', # Open/Closed shrubland
                  'WSAV', 'SAV', 'GRA',  # (woody) Savannas, Grasslands
                  'WET',  'CRO', 'URB',  # wetland, Croplands
                  'CRO/NVM', 'SI',       # Natural Vegetation Mosaic; snow and ice
                  'BAR',  'WAT')  

        col <- c('#008000', '#00FF00', '#99CC00', '#99FF99', '#339966', 
                 '#993366', '#FFCC99', '#CCFFCC', '#FFCC00', '#FF9900', 
                 '#006699', '#FFFF00', '#FF0000', '#999966', '#FFFFFF', 
                 '#808080', '#000080')
        out <- data.frame(val = val, name = name, col = col, stringsAsFactors = F)
        return(out)
    }

    igbp.col.df <- get.igbp.col() # get igbp colors
    name <- igbp.col.df$name; attributes(name)$names <- igbp.col.df$val
    col <- igbp.col.df$col; attributes(col)$names <- igbp.col.df$val
    igbp.col.list <- list(name = name, val = igbp.col.df$val, col = col)

    i1 <- map + coord_equal(scale.coord) + 
          geom_point(data = df, aes(lon, lat, color = as.factor(igbp)), alpha = 0.6) + 
          scale_color_manual(name = NULL, values = igbp.col.list$col, 
                             breaks = igbp.col.list$val, 
                             labels = igbp.col.list$name) + 
          labs(x = 'LONGITUDE', y = 'LATITUDE') + 
          guides(color = guide_legend(ncol = legend.ncol))

    return(i1)
}
