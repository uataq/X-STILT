# script to plot latitude series, DW, 03/05/2021 

plot.bg.3d = function(site, timestr, obs_df, bg_df) {
  
  # plot screened soundings if there is
  obs_uni = obs_df %>% dplyr::select(lat, lon, val, plmTF, group) %>% unique()
  obs_plm = obs_uni %>% filter(plmTF == TRUE)

  bgs = bg_df$bg.mean 
  bg_df = bg_df %>% filter(bg.mean %in% sort(bgs)[c(1, 2)])
  uni_side = sort(unique(bg_df$bg.side))

  library('scatterplot3d')
  title = paste('X-section of Xgas and defined background [ppm]\nfor', site, 'on', timestr)
  l1 = scatterplot3d(obs_uni$lon, obs_uni$lat, obs_uni$val, color = 'black', box = F, 
                     angle = 20, grid = T, cex.symbols = 0.4, main = title, 
                     xlab = 'LONGITUDE', ylab = 'LATITUDE', zlab = 'Xgas')
  l1$points3d(obs_plm$lon, obs_plm$lat, obs_plm$val, col = 'brown')

  # select obs for available bg side 
  pts.col = ggdef.col(4)[c(2, 3)]
  for (s in 1 : length(uni_side)) {
    bg_tmp = bg_df %>% filter(bg.side == uni_side[s])
    obs_bg = obs_uni %>% filter(lon > bg_tmp$bg.xmn, lon < bg_tmp$bg.xmx, 
                                lat > bg_tmp$bg.ymn, lat < bg_tmp$bg.ymx) 
    
    l1$points3d(obs_bg$lon, obs_bg$lat, obs_bg$val, cex = 0.4, col = pts.col[s])
    l1$plane3d(bg_tmp$bg.mean, 0, 0, draw_polygon = T, draw_lines = F, lty.box = NULL, 
               polygon_args = list(border = 'gray70', col = alpha(pts.col, 0.3)[s]))
  } # end for
  legend('top', legend = toupper(uni_side), pch = 1, box.col = 'gray80', 
                horiz = TRUE, col = pts.col)

  return(l1)
}