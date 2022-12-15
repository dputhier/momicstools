#' @name compute_dist_to_border
#' @title Given a boarder between two spot groups, infer for each spot a: (i) the distance the
#' closest spot (b) in the other class, (ii) the x/y coordinate of the point lying at the crossing
#' between the border and (a,b) segment, (iii) the distance to the border.
#' @param coord_spot a data.frame with spot coordinates (columns "x" and "y") and a column k (0 or 1) giving the class of the spot.
#' @param border_segments a data.frame as produce by compute_visium_ortho_hull.
#' @param diagnostic_plot Whether to produce a diagnostic diagram. Highly recommanded to visualy inspect the results.
#' @keywords hull, spatial transcriptomics, visium, border, spot.
#' @return a dataframe with the following columns: x (input x coordinate), y (input y coordinate),
#' tgt_name (the name/index of the closest spot), tgt_x (the x coordinate of the closest spot),
#' tgt_y (the y coordinate of the closest spot), dist2tgt (the euclidean distance to the closest spot),
#' x_inter (the x coordinate of the intersection with the border), y_inter (the y coordinate of the
#' intersection with the border), dist2_inter (the distance to the intersection with the border),
#' k (input k).
#' @examples
#' data(coord_spot_brain)
#' cluster_to_show <- 0 # Could be also c(a, b)
#' coord_spot_brain$k <- ifelse(coord_spot_brain$group %in% cluster_to_show, 1, 0)
#' border_segments <- compute_visium_ortho_hull(coord_spot_brain, size_x=3.6, size_y=3.4, delta=0.5)
#' compute_dist_to_border(coord_spot, border_segments)
#' @export compute_dist_to_border
library(spatstat.geom) # to compute segment intersection
library(pdist)
library(tibble)
library(pracma)

# coord_spot <- coord_spot_brain

compute_dist_to_border <- function(coord_spot, border_segments, diagnostic_plot=TRUE){

  rownames(coord_spot) <- 1:nrow(coord_spot)

  # Compute first the distance of spot 0 to the border
  # then of spot 1
  k1 <- 0
  k2 <- 1

  results <- list()

  for(iter in 1:2){

    if(iter == 2){
      k1 <- 1
      k2 <- 0
    }

    source_point_and_info <- coord_spot[coord_spot$k==k1, 1:2]
    source_point_and_info$tgt_name <- NA
    source_point_and_info$tgt_x <- NA
    source_point_and_info$tgt_y <- NA
    source_point_and_info$dist2tgt <- NA
    source_point_and_info$x_inter <- NA
    source_point_and_info$y_inter <- NA
    source_point_and_info$dist2_inter <- NA
    source_point_and_info$k <- NA

    target_point <- coord_spot[coord_spot$k==k2, 1:2]
    dist_to_point_from_other_class <- distmat(as.matrix(source_point_and_info[,c("x", "y")]),
                                              as.matrix(target_point))


    window_delta_x <- 10/100 * (max(coord_spot$x) - min(coord_spot$x))
    window_delta_y <- 10/100 * (max(coord_spot$y) - min(coord_spot$y))

    source_point_and_info$tgt_name <- NA
    source_point_and_info$dist2tgt <- NA

    for(i in 1:nrow(source_point_and_info)){
      pos <- which(dist_to_point_from_other_class[i,] == min(dist_to_point_from_other_class[i,]))[1]
      source_point_and_info[i,c("tgt_name", "tgt_x", "tgt_y", "dist2tgt", "k")] <- c(as.integer(names(pos)),
                                                                                coord_spot$x[as.integer(names(pos))],
                                                                                coord_spot$y[as.integer(names(pos))],
                                                                                dist_to_point_from_other_class[i,pos],
                                                                                iter-1)
    }


    src_to_tgt_segment_psp <- psp(x0=coord_spot[coord_spot$k==k1,]$x,
                                  y0=coord_spot[coord_spot$k==k1,]$y,
                                  x1=coord_spot[source_point_and_info$tgt_name,]$x,
                                  y1=coord_spot[source_point_and_info$tgt_name,]$y,
                                  window=owin(xrange=c(min(coord_spot$x - window_delta_x),
                                                       max(coord_spot$x + window_delta_x)),
                                              yrange=c(min(coord_spot$y - window_delta_y),
                                                       max(coord_spot$y + window_delta_y))))

    border_segments_psp <- psp(x0=border_segments$x1,
                               y0=border_segments$y1,
                               x1=border_segments$x2,
                               y1=border_segments$y2,
                               window=owin(xrange=c(min(coord_spot$x - window_delta_x),
                                                    max(coord_spot$x + window_delta_x)),
                                           yrange=c(min(coord_spot$y - window_delta_y),
                                                    max(coord_spot$y + window_delta_y))))
    for(i in 1:nrow(source_point_and_info)){
      tmp <- crossing.psp(border_segments_psp, src_to_tgt_segment_psp[i,])
      source_point_and_info[i,c("x_inter", "y_inter")] <- c(tmp$x, tmp$y)
    }

    for(i in 1:nrow(source_point_and_info)){
      source_point_and_info$dist2_inter[i] <- dist(matrix(c(source_point_and_info$x[i],
                                                         source_point_and_info$y[i],
                                                         source_point_and_info$x_inter[i],
                                                         source_point_and_info$y_inter[i]),
                                                       ncol=2,
                                                       byrow = TRUE))
    }


    src_to_intersect_segment_psp <- psp(x0=source_point_and_info$x,
                                        y0=source_point_and_info$y,
                                        x1=source_point_and_info$x_inter,
                                        y1=source_point_and_info$y_inter,
                                        window=owin(xrange=c(min(coord_spot$x - window_delta_x),
                                                             max(coord_spot$x + window_delta_x)),
                                                    yrange=c(min(coord_spot$y - window_delta_y),
                                                             max(coord_spot$y + window_delta_y))))
    results[[iter]] <- source_point_and_info
  }

  results <- rbind(results[[1]], results[[2]])
  results$k <- as.factor(results$k)

  if(diagnostic_plot){

      p <- ggplot(data=results, mapping=aes(x=x, y=y)) + geom_point() +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            geom_segment(data=results, aes(x=x, y=y, xend=tgt_x, yend=tgt_y, color=k), size=0.3) +
            geom_segment(data=border_segments,
                       mapping=aes(x=x1, y=y1, xend=x2, yend=y2), size=0.8, col="black") +
            geom_segment(data=results,
                   mapping=aes(x=x, y=y, xend=x_inter, yend=y_inter), col="green", size=0.3)
     print(p)

  }

  return(results)
}

