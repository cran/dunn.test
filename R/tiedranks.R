# tiedranks: enumerates tied values from a vector of ranks provided:
#   ranks: a vector of rank values
# Returns: ties, a vector of values that are tied
# Date: March 31, 2026
tiedranks <- function(ranks) {
  ranks <- sort(ranks)
  #enumerate tied values
  ties <- c()
  for (i in 2:length(ranks)) {
    if (ranks[i-1] == ranks[i]) {
      if (length(ties) > 0) {
        if (ranks[i-1] != tail(ties,n=1)) {
          ties <- c(ties,ranks[i-1])
          }
        }
       else {
        ties <- c(ranks[i-1])
        }
      }
    }
  return(ties)
  }