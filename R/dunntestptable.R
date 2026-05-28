# dunntestptable: displays Dunn's test p values
#   groupvar: group variable
#   index: index value of the group variable
#   P: real valued p-values (adjusted or not as appropriate)
#   colstart: integer indicating starting column
#   colstop: integer indicating stopping column
#   Reject: 0/1 indicators of whether test is rejected given alpha and method
#   last: 0/1 indicator of whether this is the last row
# Outputs: formatted p-values table row message for dunn.test()
# Date: April 1, 2026
dunntestptable <- function(index, P, colstart, colstop, Reject, last) {
  # p row header
  p_row_head <- "         \U2502"

  # Table p entries
  p_row_tail = c(" ")
  for (i in colstart:colstop) {
    p <- P[(((index-2)*(index-1))/2) + i]
    if ( Reject[(((index-2)*(index-1))/2) + i] == 0) {
      p_row_tail <- paste0(p_row_tail, "    ", sprintf("%1.4f",p), " ", sep="")
      }
     else {
      p_row_tail <- paste0(p_row_tail, "    ", sprintf("%1.4f",p), "*", sep="")
      }
    }
  rlang::inform(message=paste0(p_row_head, p_row_tail, sep=""))

  # Close out with another blank row header
  if (last == 0) {
    rlang::inform(message="         \U2502") 
    }
  }

