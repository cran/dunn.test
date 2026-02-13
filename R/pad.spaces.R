################################################################################
# pad.spaces returns a string containing the number of spaces given as an 
# argument.
# Author: Alexis Dinno
# Date: September 6, 2024
# Takes: A single positive integer
pad.spaces <- function(n) {
  return(paste0(rep(" ",n),collapse=""))
  }