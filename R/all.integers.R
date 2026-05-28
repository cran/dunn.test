# all.integers: robustly tests whether all elements of a vector are integers
# Date: March 31, 2026
all.integers <- function(x) {
  for (i in 1:length(x)) {
    if (is.na(x[i]) | is.list(x[i]) | length(x[i]) > 1 | !is.numeric(x[i])) {
      return(FALSE)
      }
    if (x[i]%%1!=0) {
      return(FALSE)
      }
    }
  return(TRUE)
  }

