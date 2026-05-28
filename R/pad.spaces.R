################################################################################
# pad.spaces returns a string containing the number of spaces given as an 
# argument.
# Author: Alexis Dinno
# Date: May 20, 2026
# Takes: A single positive integer
pad.spaces <- function(n) {
  n <- max(0, n)
  return(paste0(rep(" ",n),collapse=""))
  }