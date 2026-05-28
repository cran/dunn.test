################################################################################
# Returns s truncated to length width, and, if necessary padded with spaces to 
# the left so that the string is exactly width wide.
# Author: Alexis Dinno
# Date: May 20, 2026
# Takes: A string s, and a positive integer width.
pad.left <- function(s,width) {
  width <- max(0, width)
  len.s <- nchar(s)
  if (len.s < width) {
    s <- paste0(pad.spaces(width - len.s),s)
    }
  return(substr(s,1,width))
  }
