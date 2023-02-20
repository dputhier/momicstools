#' @name visium_hull
#' @title Internal function. Compute location of segments used to create an orthogonal 
#' hull around points related to a particular class.
#' @param data a data.frame.
#' @param x The column name storing the x coord
#' @param y The column name storing the y coord
#' @param k The column name storing the classes (0 not part of the class of interest, 1 part of the class of interest)
#' @param size_y The size of the square (y axis)
#' @param size_x The size of the square (x axis)
#' @param step_y The distance between two points on the y axis.
#' @param step_x The distance between two points on the x axis.
#' @param delta Add more or less flexibility to search for neighbor points
#' @keywords hull, spatial transcriptomics, visium
#' @return A dataframe with coordinates x, y, xend, yend (see geom_segments).
#' @examples
#' if (! "stxBrain" %in% InstalledData()$Dataset) InstallData('stxBrain')
#' LoadData("stxBrain", type = "anterior1")
#' # Preprocess the data
#' anterior1_preprocessed <- preprocess_seurat_data(anterior1, normalization_method = "NormalizeData", approx=F)
#' ## Retrieve x/y coordinates and group from ggplot object
#' coord_st_data <- getFlippedTissueCoordinates(anterior1_preprocessed, as_data_frame=TRUE)
#' # Create a hull around cluster 0
#' coord_st_data$k <- ifelse(Idents(anterior1_preprocessed)==0, 1, 0)
#' ## Compute the segments of the hull.
#' path <- visium_hull(coord_st_data, size_x=3.2, size_y=3.6, delta=0.5)
#' ## Add the segments to the ggplot diagram
#' spatial_plot <- SpatialDimPlot(anterior1_preprocessed, label = TRUE, label.size = 3, pt.size.factor = 1.5)
#' spatial_plot +
#'   theme_bw() +
#'   geom_segment(data=path,
#'                mapping=aes(x=x1,
#'                            y=y1,
#'                            xend=x2,
#'                            yend=y2),
#'                inherit.aes = F,
#'                color="white",
#'                size=0.7)
#' @export visium_hull
visium_hull <- function(data,
                        x="x",
                        y="y",
                        k="k",
                        size_x=4.6,
                        size_y=5,
                        step_x=2.6,
                        step_y=2.4,
                        delta=0.3){

    center_and_rotate <- function(data, center_x, center_y, angle, x="x", y="y"){

      data[, x] <- data[, x] - center_x
      data[, y] <- data[, y] - center_y
      angle_rad <- pi * angle / 180

      # center
      data[p, x] <- data[p, x] - center_x
      data[p, y] <- data[p, y] - center_y

      # size_y and size_x will be half the size:
      size_y <- size_y / 2
      size_x <- size_x / 2

      # rotate
      rotation_mat <- matrix(c(cos(angle_rad),
                               -sin(angle_rad),
                               sin(angle_rad),
                               cos(angle_rad)),
                             ncol=2, byrow = T)

      rotated_mat <- as.matrix(data[, c(x,y)]) %*% rotation_mat
      colnames(rotated_mat) <- c(x,y)
      data[,c(x,y)] <- rotated_mat[,c(x,y)]

      return(data)
  }


  get_neighbor_class <- function(data_cr){
    test_x <- data_cr$x > 0 & data_cr$x < step_x * 2 + (step_x * 2) * delta
    test_y <- data_cr$y > 0 - (step_y/2)  & data_cr$y < 0 + (step_y/2)
    hit <- data[test_x & test_y, ]

    if(nrow(hit) > 0){
      neighbor_class <- hit[, k]
    }else{
      neighbor_class <- NULL
    }

    return(neighbor_class)
  }

  # Create a dataframe to store
  # square segment

  print_msg("Creating a dataframe to store output", msg_type="INFO")

  df_coord <- data.frame(matrix(NA, ncol=4, nrow=6))
  rownames(df_coord) <- sapply("region_name",
                               paste0, "_",
                               c("north_west", "north_east",
                                 "south_east", "south_west",
                                 "west", "east"))
  colnames(df_coord) <- c("x1", "x2", "y1", "y2")

  print_msg("Looping over the points.", msg_type="INFO")

  for(p in 1:nrow(data)){

    x_p <- data[p, x]
    y_p <- data[p, y]
    pt_class <- data[p, k]

    p_name <- paste0(p, "_", x_p, "_", y_p, "_")

    print_msg(paste0("Processing: ", p_name), msg_type="DEBUG")


    if(pt_class == 1){
      print_msg(paste0("x_p, y_p", x_p, " ", y_p), msg_type="DEBUG")

      ## West
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=180, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)

      if(neighbor_class == 0 || is.null(neighbor_class)){
          print_msg(paste0(p_name, " West"), msg_type="DEBUG")
          print_msg(as.character(neighbor_class), msg_type="DEBUG")

        df_coord[paste0(p_name, "west"),] <- c(x1 = x_p - size_x,
                                               x2 = x_p - size_x,
                                               y1 = y_p -  size_y,
                                               y2 = y_p +  size_y
        )
      }


      ## Check_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=0, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)


      if(neighbor_class == 0 || is.null(neighbor_class)){
        
        print_msg(paste0(p_name, " East"), msg_type="DEBUG")
        print_msg(as.character(neighbor_class), msg_type="DEBUG")
        
        df_coord[paste0(p_name, "east"),] <- c(x1 = x_p + size_x,
                                               x2 = x_p + size_x,
                                               y1 = y_p -  size_y,
                                               y2 = y_p +  size_y
        )
      }

      ## North_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=60, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)

      if(neighbor_class == 0 || is.null(neighbor_class)){

        print_msg(paste0(p_name, " North_east"), msg_type="DEBUG")
        print_msg(as.character(neighbor_class), msg_type="DEBUG")

        df_coord[paste0(p_name, "north_east"),] <- c(x1 = x_p,
                                                     x2 = x_p + size_x,
                                                     y1 = y_p +  size_y,
                                                     y2 = y_p +  size_y
        )
      }

      ## North_west
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=120, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)

      if(neighbor_class == 0 || is.null(neighbor_class)){
        
        print_msg(paste0(p_name, " North_west"), msg_type="DEBUG")
        print_msg(as.character(neighbor_class), msg_type="DEBUG")
        
        df_coord[paste0(p_name, "north_west"),] <- c(x1 = x_p -  size_x,
                                                     x2 = x_p ,
                                                     y1 = y_p +  size_y,
                                                     y2 = y_p +  size_y
        )
      }

      ## South_west
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=240, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)

      if(neighbor_class == 0 || is.null(neighbor_class)){
        
        print_msg(paste0(p_name, " South_west"), msg_type="DEBUG")
        print_msg(as.character(neighbor_class), msg_type="DEBUG")

        df_coord[paste0(p_name, "south_west"),] <- c(x1 = x_p - size_x,
                                                     x2 = x_p,
                                                     y1 = y_p - size_y,
                                                     y2 = y_p - size_y
        )
      }

      ## South_east
      data_cr <- center_and_rotate(data, center_x = x_p, center_y = y_p, angle=300, x=x, y=y)
      neighbor_class <- get_neighbor_class(data_cr)

      if(neighbor_class == 0 || is.null(neighbor_class)){
        
        print_msg(paste0(p_name, " South_east"), msg_type="DEBUG")
        print_msg(as.character(neighbor_class), msg_type="DEBUG")

        df_coord[paste0(p_name, "south_east"),] <- c(x1 = x_p,
                                                     x2 = x_p +  size_x,
                                                     y1 = y_p - size_y,
                                                     y2 = y_p - size_y
        )
      }
    }
  }


  return(stats::na.omit(df_coord))
}
