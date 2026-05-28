# get_order: takes a vector of values or labels v1, and
# a vector of observations of those values or labels v2
# and returns a vector of the ordering of v1 in v2.
# E.g., v1==c("b","a"), and v2==c("b","a","a",b") will
# return c(1, 3, 4, 2)
# Date: March 31, 2026
get_order <- function(v1,v2) {
  n1 <- length(v1)
  n2 <- length(v2)
  out <- rep(NA, times=n2)
  count <- 1
  for (i in 1:n1) {
    for (j in 1:n2) {
      if (v2[j] == v1[i]) {
        out[j] <- count
        count <- count + 1
        }
      }
    }
  return(out)
  }