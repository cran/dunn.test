# dunntestheader: displays Dunn's test table headers, taking:
#   groupvar: name the group variable
#   colstart: integer indicating column value to start at
#   colstop: integer indicating column value to stop at
#   rmc: boolean indicating what p-value convention to use
# Outputs: table header message for dunn.test()
# Date: March 31, 2026
dunntestheader <- function(groupvar, colstart, colstop, rmc) {
  if (rmc==FALSE) {
    rlang::inform(message="Col Mean-\U2502")
    RMCM_head <- "Row Mean \U2502"
    }
   else {
    rlang::inform(message="Row Mean-\U2502")
    RMCM_head <- "Col Mean \U2502"
    }
  groupvalues  <- levels(factor(groupvar))
  RMCM_tail <- ""
  for (col in colstart:colstop) {
    vallab  <- substr(groupvalues[col], 1, 8)
    pad     <- 8-nchar(vallab)
    colhead <- paste0(tpad(n=pad),vallab)
    RMCM_tail <- paste0(RMCM_tail, "   ", substr(colhead, 1, 8), sep="")
    }
  rlang::inform(message=paste0(RMCM_head, RMCM_tail, sep=""))  
  separatorlength = 10 + 11*(colstop-colstart)+1
  rlang::inform(message=paste0(tpad(n=9, pad="\U2500"), "\U253C", tpad(n=separatorlength, pad="\U2500"), sep=""))
  }