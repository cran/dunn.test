# centertext: centers a string within a specific width, taking:
#   text: a string,
#   width: an integer indicating the desired width of the centered string,
#   lmargin: an integer indicating the left margin in spaces,
#   rmargin: an integer indicating the right margin in spaces,
# NOTE: text may be truncated depending on the values of width, lmargin, and rmargin
# Returns: a formatted string
# Date: April 1, 2026
centertext <- function(text, width=80, lmargin=2, rmargin=2) {
  textwidth <- nchar(text)
  if (textwidth > width - lmargin - rmargin) {
    text <- substr(text, 1, width - lmargin - rmargin)
    textwidth <- width - lmargin - rmargin
    }
  # buff is the number of spaces needed to add to the possibly truncated text
  # in order to center it. buff is divided into leftspace and rightspace.
  buff <- (width - textwidth)
  leftspace <- ""
  rightspace <- ""
  # Use floor() and ceiling() to add the extra space on the left when buff is
  # odd. (leftspace equals rightspace if buff is even.)
  if (buff > 0) {
    leftspace <- paste0(rep(" ", ceiling(buff/2)), collapse="")
    rightspace <- paste(rep(" ", floor(buff/2)), collapse="")
    }
  out <- paste0(leftspace, text, rightspace, collapse="", sep="")
  return(out)
  }