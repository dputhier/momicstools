#################################################################
##    set_verbosity
#################################################################
#' Set the verbosity level for the scomicstools package
#'
#' This function sets the verbosity level for the scomicstools package,
#' which controls the amount of information that is printed to the console by
#' the \code{\link{print_msg}} function. The verbosity level can be set to
#'  any non-negative integer, with higher values indicating more detailed output.
#'  By default, the verbosity level is set to 1.
#'
#' @param verbosity_value A non-negative integer indicating the verbosity level to be set.
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' # Set verbosity level to 2
#' set_verbosity(2)
#'
#' # Set verbosity level to 0
#' set_verbosity(0)

# 0 : No message
# 1 : Display only INFO type message
# 2 : Display both INFO and DEBUG type message

set_verbosity <- function(verbosity_value) {
  if (!is.null(verbosity_value) &
      verbosity_value >= 0 & is.numeric(verbosity_value)) {
    options(scomicstools_verbosity = verbosity_value)
  }
}

#################################################################
##    get_verbosity()
#################################################################
#' Get the current verbosity level.
#'
#' This function get the verbosity level of the scomicstools package which
#' controls the amount of information that is printed to the console by
#' the \code{\link{print_msg}} function.
#'
#'
#' @return A vector
#'
#' @export
#'
#' @examples
#' get_verbosity()
#'

get_verbosity <- function() {
  if (is.null(unlist(options()["scomicstools_verbosity"]))) {
    options(scomicstools_verbosity = 1)
  }
  return(options()$scomicstools_verbosity)
}

#################################################################
##    print_msg
#################################################################
#' Print a message based on the level of verbosity
#'
#' @param msg The message to be printed
#' @param msg_type The type of message, one of "INFO", "DEBUG", or "WARNING"
#'
#' @return None
#'
#' @export
#' @examples
#' set_verbosity(1)
#' print_msg("Hello world!", "INFO")
#' set_verbosity(3)
#' print_msg("Debugging message", "DEBUG")
#' set_verbosity(0)
#' print_msg("Hello world!", "INFO")
#' print_msg("Debugging message", "DEBUG")
#' options(warn=0)
#' print_msg("Warning message", "WARNING")
#' options(warn=-1)
#' print_msg("A warning message not displayed", "WARNING")
#' options(warn=opt_warn)
print_msg <-
  function(msg,
           msg_type = c("INFO", "DEBUG", "WARNING", "STOP")) {
    if (is.null(unlist(options()["scomicstools_verbosity"]))) {
      options(scomicstools_verbosity = 1)
    }
    if (msg_type == "INFO")
      if (unlist(options()["scomicstools_verbosity"]) > 0)
        cat(paste("|-- INFO : ", msg, "\n"))
    if (msg_type == "DEBUG")
      if (unlist(options()["scomicstools_verbosity"]) > 1)
        cat(paste("|-- DEBUG : ", msg, "\n"))
    if (msg_type == "WARNING")
      warning("|-- WARNING : ", msg, call. = FALSE)
    if (msg_type == "STOP")
      stop(paste0("|-- STOP : ", msg), call. = FALSE)
  }

#################################################################
##    print_stat
#################################################################
#' Mostly a debugging function that will print some summary
#' statistics about a numeric vector, matrix or dataframe.
#'
#' @param msg The message to users.
#' @param data  The vector (numeric) for which the stats are to be produced
#' (a vector, )
#' @param msg_type The type of message, one of "INFO", "DEBUG", or "WARNING"
#' @param round_val Round the values in its first argument to the specified number
#'  of decimal. Set argument to -1 for no rounding
#' @return None
#'
#' @export
#' @examples
#' print_stat("My data", 1:10, msg_type="INFO")
#' set_verbosity(3)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' set_verbosity(0)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' (opt_warn <- options()$warn)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="WARNING")
#' options(warn=-1)
#' print_stat("My data", matrix(rnorm(10), nc=2), msg_type="WARNING")
#' options(warn=opt_warn)
print_stat <-
  function(msg,
           data,
           round_val = 2,
           msg_type = c("INFO", "DEBUG", "WARNING")) {
    
    msg_type <- match.arg(arg = msg_type, c("DEBUG", "WARNING", "INFO"))
    
    if (inherits(data, "data.frame")) {
      data <- as.matrix(data)
    }
    
    data <- as.vector(data)
    
    if (!is.numeric(data)) {
      
      print_msg("Can't print stats from numeric object", msg_type = "WARNING")
      stats="No Statistics"
    }else{
      stats <- summary(data)
      names(stats) <- c("Min", "Q1", "Med", "Mean", "Q3", "Max")
      
      if (round_val > 0 & is.numeric(round_val)) {
        stats <- round(stats, round_val)
      }
      
      stats <- paste(names(stats), stats, sep = ":", collapse = " ")
    }
    
    print_msg(paste0(msg, ": ", stats), msg_type = msg_type)
    
  }

#################################################################
##    print_stat
#################################################################
#' Generate a vector of colors for a gradient
#'
#' This function generates a vector of colors for a gradient, given 
#' a specified palette name.
#'
#' @param palette A character vector specifying the palette to use. One of: "Je1", 
#' "Seurat_like", "Ju1", "De1",  "De2", "De3", "De4", "De5", "De6", "De7", "De8".
#' @return A character vector of color codes.
#' @export colors_for_gradient
#' @examples
#' colors_for_gradient()
#' colors_for_gradient(palette = "Seurat_like")
#' 
colors_for_gradient <- function(palette=c("Je1", "Seurat_like", "Ju1", "De1", 
                                          "De2", "De3", "De4", "De5",
                                          "De6", "De7", "De8")){
  if(palette == "Seurat_like"){
    return(c("#5D50A3", "#9FD7A4", "#FBFDBA", "#FEB163", "#A80B44"))
  }else if(palette == "Ju1"){
    return(c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"))
  }else if(palette == "De1"){
    return(c("#d73027","#fc8d59","#fee090","#e0f3f8","#91bfdb","#253494"))
  }else if(palette == "De2"){
    return(c("#FFF7FB","#ECE2F0","#D0D1E6","#A6BDDB","#67A9CF","#3690C0","#02818A","#016450"))
  }else if(palette == "De3"){
    return(c("#1A1835","#15464E","#2B6F39","#757B33","#C17A70","#D490C6","#C3C1F2","#CFEBEF'"))
  }else if(palette == "De4"){
    return(c("#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000"))
  }else if(palette == "De5"){
    return(c("#0000AA","#0000FF","#00FFFF","#80FF80","#FFFF00","#FF0000","#AA0000"))
  }else if(palette == "De6"){
    return(c("#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027"))
  }else if(palette == "De7"){
    return(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))
  }else if(palette == "De7"){
    return(c("#2b83ba","#abdda4","#fdae61","#d7191c"))
  }else if(palette == "De8"){
    return(c("#0000BF","#0000FF","#0080FF","#00FFFF","#40FFBF","#80FF80","#BFFF40","#FFFF00","#FF8000","#FF0000","#BF0000"))
  }else if(palette == "Je1"){
    c("#27408B", "#3A5FCD", "#3288BD", "#66C2A5","#ABDDA4", "#E6F598","#FEE08B", "#FDAE61","#F46D43","#D53E4F","#8B2323")
  }
}

