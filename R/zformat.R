# zformat: formats z values for display in table, taking
#   z: a real z-value
# Returns: a formatted string
# Date: May 20, 2026
# Author: Alexis Dinno
zformat <- function(z) {
  if (z < 0) {
    sign_z <- "-"
    }
   else {
     sign_z <- " "
     }
  if (abs(z) >= 1) {
     leftspaces <- max(2,floor(log10(abs(z)))+2)
     leftdigits <- floor(abs(z))
     rightspaces <- 8 - leftspaces
     rightdigits <- substr(paste0(abs(z) - floor(abs(z)),"00000000"),3,rightspaces+2)
     return(paste0(sign_z,leftdigits,".",rightdigits))
     }
  if (abs(z) <= 1) {
    if (z < 0) {
      return(sprintf("%0.6f",z))
      }
    if (z >= 0) {
      return(paste0(" ", sprintf("%0.6f",z)))
      }
    }
  }

