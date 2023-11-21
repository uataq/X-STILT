
#' @param theta roation angle in degree, positive for clockwise

rotate_plume = function(theta, pr_df, loc_df) {
    
    library(geosphere)

    # construct rotation matrix, theta + for clockwise rotation around origin
    rot = function(theta) matrix(c( cos(theta), sin(theta),
                                    -sin(theta), cos(theta)), 2)

    # need to normalize lat/lon by origin (i.e., emission source)
    # convert lat/lon to dlon/dlat from emission source
    # always rotate obs plume to match sim
    pr_df = pr_df %>% mutate(dlons = lons - loc_df$lon, 
                             dlats = lats - loc_df$lat, 
                             dlon = lon - loc_df$lon, 
                             dlat = lat - loc_df$lat) 

    pc_mtrx = pr_df %>% dplyr::select(dlon, dlat) %>% as.matrix()
    pr_mtrx = pr_df %>% dplyr::select(dlons, dlats) %>% as.matrix()
    pc_rot = pc_mtrx %*% rot(theta / 180 * pi) %>% as.data.frame() %>% 
             mutate(lon = V1 + loc_df$lon, lat = V2 + loc_df$lat) %>% 
             dplyr::select(lon, lat)

    pr_rot = pr_mtrx %*% rot(theta / 180 * pi) %>% as.data.frame() %>% 
             mutate(lons = V1 + loc_df$lon, lats = V2 + loc_df$lat) %>% 
             dplyr::select(lons, lats)

    # this is the df with rotated lat/lon 
    prr_df = cbind(pc_rot, pr_rot) %>% unique() %>% 
             mutate(corner = rep(seq(1, 4), nrow(pr_rot) / 4), 
                    indx = rep(1 : (nrow(pr_rot) / 4), each = 4)) %>% 
             left_join(pr_df[, c('indx', 'obs_tno2')] %>% unique(),by = 'indx') 

    # now match rotated OBS plume with initial MODELED plume 
    # look over each model polygon and find the closest rotated obs polygon
    pr_new = NULL
    for (p in unique(pr_df$indx)) {
        print(p)

        # initial model polygon 
        pr_tmp = pr_df %>% filter(indx == p) %>% 
                 rename(obs_tno2_init = obs_tno2, 
                        obs_tno2_uncert_init = obs_tno2_uncert)

        # calculate distance between rotated obs center and original mod center
        dist_df = do.call(rbind, lapply(unique(prr_df$indx), function(x) {
            tmp_df = prr_df %>% filter(indx == x) %>% 
                        dplyr::select(lon, lat, indx) %>% unique()
            tmp_dist = distCosine(p1 = c(pr_tmp$lon[1], pr_tmp$lat[1]), 
                                  p2 = c(tmp_df$lon, tmp_df$lat))
            data.frame(indx = x, dist = tmp_dist)
        }))
        find_indx = dist_df[dist_df$dist == min(dist_df$dist), 'indx']
        
        # based on initial model df, add a column for rematched obs
        pr_tmp = pr_tmp %>% 
                 mutate(indx_rot = find_indx, 
                        obs_tno2 = pr_df[pr_df$indx == find_indx, 'obs_tno2'])
        pr_new = rbind(pr_new, pr_tmp)
    }


    return(pr_new)
}


# add rotated modeled x rotated modeled
# now five sets: gfs * obs; hrrr * obs; gfs * gfs; hrrr * hrrr; obs * obs
#' @param ini_df dataframe to be rotated
#' @param pts_df dataframe containing remapped pixels, if not, use @param ini_df
rotate_plumev2 = function(ini_df, dangle = 5, rot_center = NULL, 
                          thetas = NULL, pts_df = NULL) {
    
    library(artKIT); library(dplyr); library(sp)
    ini_df = ini_df %>% rename(x = lons, y = lats, group = indx) %>% 
             mutate(rotation = 0, num.edges = 4)
    
    if (is.null(pts_df)) {
        pts_df = ini_df
    } else pts_df = pts_df %>% rename(x = lons, y = lats, group = indx) %>% 
                    mutate(rotation = 0, num.edges = 4)

    # default rotation center is the center of the image
    if (is.null(rot_center))    
        rot_center = data.frame(x = mean(df$lon), y = mean(df$lat))

    # it will only rotate x, y
    if (is.null(thetas)) thetas = seq(-180, 180, dangle) 
    rads = thetas * pi / 180

    all_df = NULL
    for (tt in 1: length(rads)) {
        
        if (tt %% 5 == 0) 
            print(paste('# --', signif(tt / length(rads) * 100, 3), '% -- #'))

        #' @param rotation unit in radian, not degrees
        rot_df2 = rotate_polygon(ini_df, rotation = rads[tt], 
                                 center.of.rotation = rot_center) %>% 
                  group_by(group) %>% 
                  mutate(cx = mean(x), cy = mean(y)) %>% ungroup()

        # then write a script to match rotated grid with simulation 
        # try using point.in.polygon() 
        # here we rotated simulated plumes to match the fixed observed plumes 
        # loop over each initial pixels
        new_df = NULL
        for (p in unique(pts_df$group)) {
            
            # keep observed the same
            tmp_pts_df = pts_df %>% filter(group == p)

            # pip will create holes, need to use distantce between the centroids
            pts = tmp_pts_df %>% dplyr::select(lon, lat) %>% unique() %>% 
                  rbind(rot_df2 %>% dplyr::select(lon = cx, lat = cy) %>% 
                        unique())
            
            # need the distance between row 1 and row 2+ 
            dist = dist(pts, method = "euclidean")[1:(nrow(pts) - 1)]
            min_dist = min(dist)

            row_indx = which(dist == min_dist) + 1
            sel_rot = rot_df2 %>% filter(group == row_indx - 1)

            # if the coordinates after the rotation are outside the domain of interest, drop rotated grids
            bff = 0.02
            if (max(sel_rot$x) < min(tmp_pts_df$x) - bff | 
                min(sel_rot$x) > max(tmp_pts_df$x) + bff |
                max(sel_rot$y) < min(tmp_pts_df$y) - bff | 
                min(sel_rot$y) > max(tmp_pts_df$y) + bff) {
                next
            }

            tmp_pts_df$sim_tno2_rot = rep(mean(sel_rot$sim_tno2_mix), 
                                          nrow(tmp_pts_df))
            tmp_pts_df$sim_tno2_nochem_rot = rep(mean(sel_rot$sim_tno2_nochem), 
                                                 nrow(tmp_pts_df))
            tmp_pts_df$obs_tno2_rot = rep(mean(sel_rot$obs_tno2), 
                                          nrow(tmp_pts_df))
            new_df = rbind(new_df, tmp_pts_df)
        }

        all_df = rbind(all_df, new_df %>% mutate(radian = rads[tt], 
                                                 angle = radian / pi * 180))
    }

    return(all_df)
}



if (F) {

p1 = p0 + ggtitle('Initial sim') + 
     geom_polygon(data = ini_df, aes(x, y, group = group, fill = sim_tno2_mix), 
                    color = 'white', linewidth = 0.3) + 
     geom_point(data = ini_df, aes(lon, lat, color = group), size = 0.5)

p3 = p0 + ggtitle('Initial obs') + 
     geom_polygon(data = ini_df, aes(x, y, group = group, 
                    fill = obs_tno2), color = 'white', linewidth = 0.3)+
     geom_point(data = ini_df, aes(lon, lat, color = group), size = 0.5)

pp = ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
ggsave(pp, filename = 'Fig4_test_rotate_plume.png', 
          width = 10, height = 11, bg = 'white')



if (F) {    # sanity check
    library(viridis)
    r1 = ggplot() + theme_bw() + scale_fill_viridis() + 
            geom_point(data = rot_center, aes(lon, lat), shape = 4) + 
            geom_polygon(data = ini_df, aes(x, y, group = group, fill = sim_tno2_mix))  
    r2 = ggplot() + theme_bw() + scale_fill_viridis() + 
            geom_point(data = rot_center, aes(lon, lat), shape = 4) + 
            geom_polygon(data = rot_df2, aes(x, y, group = group, fill = sim_tno2_mix))  
    ggsave(ggarrange(r1, r2, ncol = 2), filename = 'test_rot_urban.png')
}

}

