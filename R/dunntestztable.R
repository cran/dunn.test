# dunntestztable displays Dunn's test z values, taking:
#   groupvar: group variable
#   index: index value of the group variable
#   Z: real valued Z test statistic
#   colstart: integer indicating starting column
#   colstop: integer indicating stopping column
# Outputs: formatted Z statistics table row message for dunn.test()
# Date: March 31, 2026
dunntestztable <- function(groupvar, index, Z, colstart, colstop) {
  groupvalues  <- levels(factor(groupvar))

  # Row headers
  vallab  <- substr(groupvalues[index], 1, 8)
  pad     <- 8-nchar(vallab)
  rowhead <- paste0(tpad(n=pad), vallab)
  z_row_head <- paste0(rowhead, " \U2502")

  # Table z entries
  z_row_tail <- ""
  for (i in colstart:colstop) {
    z <- Z[(((index-2)*(index-1))/2) + i]
    z_row_tail <- paste0(z_row_tail, "  ", zformat(z), sep="")
    }
  rlang::inform(message=paste0(z_row_head, z_row_tail, sep=""))
  }