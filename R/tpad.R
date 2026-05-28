# tpad: returns multiples of a string, defaulting to spaces. Takes:
#   n: a number of replications,
#   pad: string to replicate
# Returns: a formatted string
# Date: March 31, 2026
tpad <- function(n=1, pad=" ") {
  if (n <= 0) {
    return("")
    }
  return( paste0(rep(x=pad,times=n),collapse="") )
  }