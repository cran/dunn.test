# alphabetize.factor: alphabetizes a factor
# this is very quick and dirty with no checking.
# Date: March 31, 2026
alphabetize.factor <- function(x) {
  return(factor(sort(c(as.character(x)))))
  }