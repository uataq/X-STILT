
#' @param df a data frame of data.frame(lon = , lat = , wgt = )
#' modified from https://towardsdatascience.com/using-weighted-k-means-clustering-to-determine-distribution-centres-locations-2567646fc31d


#haversine distance function
haversine_dist = function(point1, point2) { 
    
    #each argument is a numeric vector with two elements (lon, lat)
    lon1 = point1[1] 
    lat1 = point1[2]
    lon2 = point2[1]
    lat2 = point2[2]
    
    R = 6371000 #earth radius in meters
    phi1 = lat1 * pi / 180 #convert to radian
    phi2 = lat2 * pi / 180 #convert to radian
    delta_phi = (lat2 - lat1) * pi / 180
    delta_lambda = (lon2 - lon1) * pi / 180
    
    a = (sin(delta_phi/2))^2 + cos(phi1) * cos(phi2) * ((sin(delta_lambda/2))^2)
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    
    #haversine distance between point1 and point 2 in meters
    distance = R * c 
    return(round(distance, 2))
}


cluster_tno2_plume = function(df, K) {
    
    #initial centroids by random
    init_centroids_index = sample(nrow(df), K)

    #initiate containers
    distance_matrix = matrix(data = NA, nrow = nrow(df), ncol = K)
    cluster = vector()
    centroid_long = vector()
    centroid_lat = vector()

    #compute distance between cities and initial centroids
    for (k in 1 : K) {
        for (i in 1 : nrow(df)) {
            loc_i = as.numeric(df[i, c('lon', 'lat')])
            centroid_k = as.numeric(df[init_centroids_index[k], 
                                       c('lon', 'lat')])
            distance_matrix[i, k] = haversine_dist(loc_i, centroid_k)
        }
    }

    #initial cluster assignment for each loc
    for (i in c(1:nrow(df))) cluster[i] = which.min(distance_matrix[i,])

    #iteration baseline
    old_cluster = vector(length = length(cluster))
    new_cluster = cluster
    
    #iterations
    while (!all(old_cluster == new_cluster)) {

        #update old cluster assignment
        old_cluster = new_cluster

        #calculate centroids using weighted average
        for (k in 1 : K) {
            cluster_k = which(old_cluster == k) #loc index of cluster k
            centroid_long[k] = weighted.mean(df$lon[cluster_k], 
                                             df$wgt[cluster_k])
            centroid_lat[k] = weighted.mean(df$lat[cluster_k], 
                                            df$wgt[cluster_k])
        }   # end for

        df_centroid = as.data.frame(cbind(centroid_long, centroid_lat))

        #compute distance between cities and centroids
        for (k in 1 : K) {
            for (i in 1 : nrow(df)) {

                loc_i = as.numeric(df[i, c('lon', 'lat')])
                centroid_k = as.numeric(df_centroid[k, ])
                distance_matrix[i, k] = haversine_dist(loc_i, centroid_k)
            }
        }

        #update cluster assignment for each loc
        for (i in 1:nrow(df)) cluster[i] = which.min(distance_matrix[i,])

        #update new_cluster
        new_cluster = cluster
    }

    df$cl = new_cluster 

    return(df)
}