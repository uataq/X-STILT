# reading trajectories from rds file and select trajec during overpass duration 
# DW, 03/04/2021


# ------------------------------------------------------------------------------ 
# td for threshold of normalized kernel density
# dx, dy - larger the number, the urban plume appears to be less curvy
fit.kde.plume = function(site, timestr, traj_path, obs_df, sensor, 
                         dx = 0.1, dy = 0.15, td = 0.15) {

  # ---------------------------------------------------------------------------- 
  # load trajec directories
  tlist = load_forward_traj(traj_path, sensor, timestr)
  trajdat = tlist$trajdat
  recp_info = tlist$recp_info
  recp_box = tlist$recp_box

  # subset trajec based on satellite overpass duration
  min.xtime = min(obs_df$datestr) - 1 * 60  # buffer for 1 min, convert to sec
  max.xtime = max(obs_df$datestr) + 1 * 60
  sel_traj = trajdat %>% filter(datestr >= min.xtime & datestr <= max.xtime) 

  if (nrow(sel_traj) == 0) {
    cat(paste('*** NO trajectories exist during the overpass time\n',
               'likely caused by incomplete traj files, please check...\n'))
    return()
  }


  # ---------------------------------------------------------------------------- 
  library(ggplot2)

  # calculate 2D kernel density and Normalized it by the max density
  cat('fit.kde.plume(): calculating kernel density...\n')
  dens = kde2d(sel_traj$lon, sel_traj$lat, h = c(dx, dy), n = 100)
  densf = data.frame(expand.grid(lon = dens$x, lat = dens$y),
                     prob = as.vector(dens$z)) %>% 
          mutate(norm.prob = prob / max(prob)) %>% dplyr::select(-prob) 
  
  # check to see if points with kd of td fall within the receptor box
  # if not, overwrite with the minimal density value, DW, 03/25/2021
  lab.norm = sort(unique(c(td, seq(0.05, 1, 0.05))))
  densf_td = densf %>% filter(norm.prob <= td)
  rpTF = which((sp::point.in.polygon(densf_td$lon, densf_td$lat, 
                                     recp_box$lon, recp_box$lat)) > 0)
  
  if ( length(rpTF) == 0 ) {
    densf_recp = densf[(sp::point.in.polygon(densf$lon, densf$lat, 
                                             recp_box$lon, recp_box$lat)) > 0, ]
    dr = ggplot() + geom_contour(data = densf_recp, 
                                 aes(lon, lat, z = norm.prob, color = ..level..), 
                                 breaks = lab.norm)
    
    # count # of pieces for each level grounp, if multiple pieces -> broken contour (remove)
    dr_df = ggplot_build(dr)$data[[1]] %>% dplyr::select(level, piece) %>% unique() %>% 
            group_by(level) %>% tally() %>% filter(n == 1) %>% arrange(level)

    td = max(td, min(dr_df$level))
    if (td > 0.2) td = 0.2
    cat(paste('fit.kde.plume(): no initial @prarm td value OR density values within the receptor box >', 
              td, '\n-> obtaining @param td as', td, '...\n'))
  } # end if update td


  # ---------------------------------------------------------------------------
  # utilize geom_contour to figure out the contours per kernel density bin
  d1 = ggplot(data = densf) + geom_contour(aes(lon, lat, z = norm.prob), breaks = lab.norm)
  kd_df = ggplot_build(d1)$data[[1]] %>% dplyr::select(x, y, level)
  uni_kd = unique(kd_df$level)

  # assign density levels to each selected particle
  sel_traj = sel_traj %>% mutate(dens.level = 0)
  for (kd in uni_kd) {
    tmp.kd = kd_df %>% filter(level == kd)
    tmp.index = sp::point.in.polygon(sel_traj$lon, sel_traj$lat, tmp.kd$x, tmp.kd$y)
    sel_traj[tmp.index > 0, 'dens.level'] = kd
  } # end for u


  # ---------------------- compute outmost boundary -------------------------
  cat('fit.kde.plume(): obtaining the urban plume and intersection between plume and soundings...\n')

  # get the contour line (i.e., urban plume) the user needs based on td 
  # and only select columns for X, Y, set polygon ID as 1
  # add unique() to the end, to eliminate repeated plume lon/lat, DW, 03/09/2021
  plm_df = kd_df %>% filter(level == td) %>% dplyr::select(X = x, Y = y) %>% unique()
  if (nrow(plm_df) == 0) { cat('fit.kde.plume(): NO urban plume defined...\n'); return() }

  # if diff is too large, meaning there are more than one polygon
  abs.dy = abs(diff(plm_df$Y))
  abs.dx = abs(diff(plm_df$X))

  # ---------------------------------------------------------------------------
  # there could be multiple polygons, so select the one we need, DW, 08/21/2018
  if ( length(which(abs.dy > td)) > 0 | length(which(abs.dx > td)) > 0 ) {
    
    cutoff.yindx = which(abs.dy > td)
    cutoff.xindx = which(abs.dx > td)

    y.vec = c(1, cutoff.yindx, nrow(plm_df)); dy.indx = diff(y.vec)
    x.vec = c(1, cutoff.xindx, nrow(plm_df)); dx.indx = diff(x.vec)

    # now select polygons
    plg.yindx = which(dy.indx == max(dy.indx))  # polygon index
    plg.xindx = which(dx.indx == max(dx.indx))  # polygon index

    # it is possible that there are more than 1 polygon that matches, 
    # use the first one that found for now, DW, 07/01/2020
    if (length(plg.xindx) * length(plg.yindx) > 1) {
        for (xi in plg.xindx) {
        for (yi in plg.yindx) {
            plg.xy = intersect(seq(y.vec[yi], y.vec[yi + 1], 1),
                               seq(x.vec[xi], x.vec[xi + 1], 1))
            if (length(plg.xy) > 0) break   # if meet the criteria then break the loop
        }  # end for yi
        }  # end for xi
    } else plg.xy = intersect(seq(y.vec[plg.yindx], y.vec[plg.yindx + 1], 1),
                              seq(x.vec[plg.xindx], x.vec[plg.xindx + 1], 1))

    plm_df = plm_df[plg.xy, ]
  }  # end if multiple polygons


  # ---------------------------------------------------------------------------- 
  # treat both OCO and TROPOMi grid as polygons, not points
  # select soundings whose polygon centers fall within the urban plume
  # use point.in.polygon, DW, 03/08/2021
  plmTF = (sp::point.in.polygon(obs_df$lon, obs_df$lat, plm_df$X, plm_df$Y)) > 0
  sel_obs = obs_df[plmTF == TRUE, ]
  uni_polygon = unique(sel_obs$polygon)

  # initialize plmTF as FALSE -> meaning no soundings within plume as default
  obs_df$plmTF = FALSE 
  obs_df = obs_df %>% mutate(plmTF = polygon %in% uni_polygon)


  # ---------------------------------------------------------------------------- 
  # store intermediate figure even if there is no intersection
  # check to if there is any valid obs falling within the urban plume
  intersectTF = FALSE; if ('TRUE' %in% obs_df$plmTF) intersectTF = TRUE
  if (intersectTF) { 
    cat('fit.kde.plume(): found satellite soundings within urban plume :)\n')
  } else cat(paste('fit.kde.plume(): *** no intersection with screened obs\n',
                   'Likely no valid soundings on', timestr, '\n'))
  
  # store the data frame for urban plume 
  #fn = file.path(traj_path, paste0('urban_plume_', site, '_', 
  #                                  substr(max(sel_obs$time_utc), 1, 10), '_', 
  #                                  sensor, '.txt'))
  #write.table(plm_df, file = fn, sep = ',', row.names = F)

  # store all outputs in a list for plotting 
  outlist = list(recp_box = recp_box, recp_info = recp_info, sel_traj = sel_traj, 
                 densf = densf, obs_df = obs_df, plm_df = plm_df, td = td, 
                 intersectTF = intersectTF)
  
  return(outlist)
} # end of function 





# ------------------------------------------------------------------------------ 
convert.df2sp = function(df, group_var = c('PID', 'polygon'), 
                         crs = CRS("+proj=longlat +datum=WGS84")) {
  
  # make a list, remove group in the list
  uni_group = unique(df[, colnames(df) == group_var])
  tmp_list = split(df, df[, colnames(df) == group_var])
  tmp_list = lapply(tmp_list, function(x) { x[colnames(x) == group_var] = NULL; x })

  # make data.frame into spatial polygon,
  # cf. http://jwhollister.com/iale_open_science/2015/07/05/03-Spatial-Data-In-R/
  # create SpatialPolygons Object, convert coords to polygon
  # ID for how many individual polygons
  #tmp_ps = sp::Polygons(lapply(tmp_list, sp::Polygon), ID = uni_group) 
  #tmp_ply = SpatialPolygons(list(tmp_ps), proj4string = crs)

  tmp_ps <- lapply(tmp_list, Polygon)
  tmp_ps <- lapply(seq_along(tmp_ps), 
                   function(i) Polygons(list(tmp_ps[[i]]), ID = uni_group[i] ))

  # create SpatialPolygons object
  tmp_ply = SpatialPolygons(tmp_ps, proj4string = crs)
  return(tmp_ply)
  
}




# ---------------------------------------------------------------------------- 
load_forward_traj = function(traj_path, sensor, timestr) {

  # load trajec directories
  cat('reading forward-time trajectories...\n')
  traj_dirs = list.dirs(traj_path, recursive = F)
  traj_dirs = traj_dirs[grepl(sensor, traj_dirs)]

  # timestr in YYYYMMDDHH
  traj_dir = traj_dirs[grepl(timestr, traj_dirs)]
  if (length(traj_dir) == 0)  # relax to the same day which is on YYYYMMDD
  traj_dir = traj_dirs[grepl(substr(timestr, 1, 8), traj_dirs)]

  traj_fns = list.files(traj_dir, '.rds', full.names = T, recursive = T)
  if (length(traj_fns) == 0) { 
  cat('*** NO forward trajec found... please check\n'); return() }
  
  # read in forward trajec
  trajdat = do.call(rbind, lapply(traj_fns, function(x) {
      if (!file.exists(x)) { cat('*** NO forward trajec found... please check\n'); return() }
      readRDS(x) 
  }))

  # ---------------------------------------------------------------------------- 
  # before selecting trajec over overpass duration, grab box receptor
  trajdat_recp = trajdat %>% filter(time == min(time))
  llon = c(min(trajdat_recp$lon), max(trajdat_recp$lon))
  llat = c(min(trajdat_recp$lat), max(trajdat_recp$lat))
  recp_box = data.frame(lon = c(llon[1], llon[2], llon[2], llon[1]),
                        lat = c(llat[1], llat[1], llat[2], llat[2]))

  # get receptor info 
  recp_info = ident.to.info(ident = basename(traj_fns), stilt.ver = 1, aglTF = F)[[1]]
  recp_info$rel.date = as.POSIXlt(recp_info$timestr, 'UTC', format = '%Y%m%d%H%M')

  # get first release time 
  min.datestr = min(recp_info$rel.date); print(min.datestr)
  trajdat = trajdat %>% mutate(datestr = min.datestr + time * 60) 

  return(list(trajdat = trajdat, recp_box = recp_box, recp_info = recp_info))
}











# acutally, no need to calculate the convex hull...
if (F) {

  # ---------------------------------------------------------------------------- 
  # calculate background since intersection between swath and trajec exists ----
  # compute convex hull
  obs_df$plmTF = FALSE  # initialize plmTF as FALSE -> no soundings within plume
  if (intersectTF) {

    hull_indx = chull(x = sel_obs$lon, y = sel_obs$lat)
    hull_indx = c(hull_indx, hull_indx[1])
    hull_obs  = sel_obs[hull_indx, ]
    hull_df   = data.frame(X = hull_obs$lon, Y = hull_obs$lat,
                           PID = rep(2, nrow(hull_obs)), 
                           POS = 1 : nrow(hull_obs))

    # find overlapping region which is also called the convex hull
    plm_bound = PBSmapping::joinPolys(plm_df, hull_df) # in PBS Mapping package
    
    if (!is.null(plm_bound)) { 

      # select obs that fall into the polluted lat range
      plm_indx = sp::point.in.polygon(obs_df$lon, obs_df$lat, 
                                      plm_bound$X, plm_bound$Y)
      obs_df[plm_indx > 0, 'plmTF'] = TRUE
    } else intersectTF = FALSE    # if NULL, meaning no soundings in plume

  } # end if intersectTF



  for (p in unique(obs_df$polygon)) {
    test_df = obs_df %>% filter(polygon == p)
    withinTF = sp::point.in.polygon(test_df$lons, test_df$lats, plm_df$X, plm_df$Y)
    #if (unique(withinTF) == 1) print(p)
  }

    # convert data frame of the urban plume to spatial polygons 
  projcrs = CRS("+proj=longlat +datum=WGS84")
  plm_sf = st_as_sf(x = plm_df %>% dplyr::select('X', 'Y', 'PID'),                         
                    coords = c('X', 'Y'), crs = projcrs)
  obs_sf = st_as_sf(x = obs_df %>% dplyr::select('lons', 'lats', 'polygon'),                         
                    coords = c('lons', 'lats'), crs = projcrs)
  pint = st_intersection(obs_sf, plm_sf)
  
  # see polygon overlaps 
  obs_sp = convert.df2sp(obs_df %>% dplyr::select('lons', 'lats', 'polygon'), 'polygon')
  plm_sp = convert.df2sp(plm_df %>% dplyr::select('X', 'Y', 'PID'), 'PID')
  pint = intersect(obs_sp, plm_sp)
  plot(plm_sp, axes = T); plot(crop(obs_sp, extent(plm_sp) + 0.5), add = T)
  plot(pint, add = T, col='red')


 # sometimes, there could be weird contour line, 
  # add receptor locations with the threshold density before fitting kde
  recp_lons = seq(llon[1], llon[2], 0.1)
  recp_lats = seq(llat[1], llat[2], 0.1)
  recp_part = data.frame(lon = c(recp_lons, rep(llon[2], length(recp_lats)), 
                                 rev(recp_lons), rep(llon[1], length(recp_lats))), 
                         lat = c(rep(llat[1], length(recp_lons)), recp_lats, 
                                 rep(llat[2], length(recp_lons)), rev(recp_lats)))
  densf_add = rbind(densf, recp_part %>% mutate(norm.prob = td))
  densf = densf_add 

}



