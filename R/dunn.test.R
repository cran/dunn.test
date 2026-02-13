# version 1.3.8 February 11, 2026 by alexis.dinno@pdx.edu
# perform Dunn's test of multiple comparisons using rank sums

p.adjustment.methods <- c("none","bonferroni","sidak","holm","hs","hochberg","bh","by")

dunn.test <- function(
  x=NA, 
  g=NA, 
  method=p.adjustment.methods, 
  kw=TRUE, 
  label=TRUE, 
  wrap=FALSE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.05, 
  altp=FALSE, 
  interpret=TRUE) {

  # FUNCTIONS

  # kwallis.test: a custom Kruskal Wallis test function to support dunn.test
  # Note: does not currently play nicely with missing values.
  kwallis.test <- function(x=NA, g=NA, pass=0) {
    # FUNCTIONS
    # tiedranks: enumerates tied values from a vector of ranks
    # ranks: a vector of rank values
    # returns: ties, a vector of values that are tied
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

    #set up data by lists
    if (is.list(x)) {
      N <- 0
      for (i in 1:length(x)) {
        N <- N + length(x[[i]])
        }
      Data <- matrix(NA,N,4)
      for (i in 1:N) {
        Data[i,1] <- i
        }
      obs <- c()
      group <- c()
      for (i in 1:length(x)) {
        obs <- c(obs,x[[i]])
        group <- c(group,rep(i,length(x[[i]])))
        }
      Data[,2] <- obs
      if (length(g) > 1) {
        Data[,3] <- g
        }
       else {
        Data[,3] <- group        
        }
      Data[,4] <- rank(Data[,2], ties.method="average", na.last=NA)
      }
  
    #set up data by groups
    if (!is.list(x)) {
      N <- length(x)
      Data <- matrix(NA,length(x),4)
      Data[,1] <- 1:length(x)
      Data[,2] <- x
      Data[,3] <- g
      Data[,4] <- rank(Data[,2], ties.method="average", na.last=NA)
      }
    k <- length(unique(Data[,3]))
  
    #calculate ties adjustment
    ranks <- Data[,4]
    ranks <- ranks[!is.na(ranks)]
    ties <- tiedranks(ranks)
    r <- length(ties)
    tiesadj <- 0
    if (r > 0) {
      for (s in 1:r) {
        tau <- sum(ranks==ties[s])
        tiesadj <- tiesadj + (tau^{3} - tau)
        }
      }
    tiesadj <- 1-(tiesadj/((N^3) - N))
  
    #calculate H
    ranksum <- 0
    for (i in 1:k) {
      ranksum <- ranksum + ((sum(Data[,4][Data[,3]==i]))^2)/(sum(Data[,3]==i))
      }
    H  <- ((12/(N*(N+1)))*ranksum - 3*(N+1))/tiesadj
    df <- k-1
    p  <- pchisq(H,k-1,lower.tail=FALSE)
  
    #present output
    output <- paste("Kruskal-Wallis chi-squared = ",round(H,digits=4),", df = ",df,", p-value = ",round(p,digits=2),"\n" ,sep="")
  
    invisible(list(output=output,H=H,df=df,p=p,N=N,Data=Data,data.name=data))
    }

  # all.integers robustly tests whether all elements of a vector are integers
  all.integers <- function(x) {
    for (i in length(x)) {
      if (is.na(x[i]) | is.list(x[i]) | length(x[i]) > 1 | !is.numeric(x[i])) {
        return(FALSE)
        }
      if (x[i]%%1!=0) {
        return(FALSE)
        }
      }
    return(TRUE)
    }

  # tiedranks: enumerates tied values from a vector of ranks
  # ranks: a vector of rank values
  # returns: ties, a vector of values that are tied
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
  
  # tpad: returns the pad string, n times, defaulting to spaces
  # n: a number of replications; pad: string to replicate
  tpad <- function(n=1, pad=" ") {
    if (n == 0) {
      return("")
      }
    return( paste0(rep(x=pad,times=n),collapse="") )
    }

  # zformat: formats z values for display in table
  # z: a real z-value
  # returns: a formatted string
  zformat <- function(z) {
    if (z < 0) {
      sign <- "-"
      }
     else {
       sign <- " "
       }
     leftspaces <- max(2,floor(log10(abs(z)))+2)
     leftdigits <- floor(abs(z))
     rightspaces <- 8 - leftspaces
     rightdigits <- substr(paste0(abs(z) - floor(abs(z)),"00000000"),3,rightspaces+2)
     return(paste0(sign,leftdigits,".",rightdigits))
    }

  # centertext: centers a string within a specific width
  centertext <- function(text, width=80, lmargin=2, rmargin=2) {
    textwidth <- nchar(text)
    if (textwidth <= width-lmargin-rmargin) {
      text <- substr(text, 1, width-lmargin-rmargin)
      }
    buff <- (width-lmargin-rmargin-textwidth)
    if (buff%%2 == 0) {
      return(paste(paste(rep(" ", buff/2), collapse=""), text, paste(rep(" ", buff/2), collapse=""), collapse=""))
      }
     else {
      return(paste(paste(rep(" ", 1+(buff-1)/2), collapse=""), text, paste(rep(" ", (buff-1)/2), collapse=""), collapse=""))
      }
    }

  # dunntestheader: displays Dunn's test table headers.
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

  # dunntestztable displays Dunn's test z values
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

  # dunntestptable: displays Dunn's test p values
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
    
  # alphabetize.factor: alphabetizes a factor
  # this is very quick and dirty with no checking.
  alphabetize.factor <- function(x) {
    return(factor(sort(c(as.character(x)))))
    }


  # VALIDATIONS & PREPARATIONS
  # names for output
  xname <- paste(if(is.name(substitute(x))) {deparse(substitute(x))} else {"x"})
  gname <- paste(if(is.name(substitute(g))) {deparse(substitute(g))} else {"group"})
  xaslist <- is.list(x)
  tempg <- 1:length(x)
  
  # casewise deletion of missing data
  # Start dealing with the fact that g may be numeric, string, or factor data
  if (length(g) > 1) {  
    if (label==TRUE) {
      if (is.factor(g)) {
        glevels <- levels(g)
        }
       else {
       	g <- factor(g)
        glevels <- levels(g)
        }
      }
     else {
      glevels <- names(table(addNA(g, ifany = TRUE)))
      }
    Data <- data.frame(x,g)[order(levels(g)[c(g)]),]
    Data <- Data[!is.na(unlist(Data$x)),]
    Data <- Data[!is.na(Data$g),]
    x <- Data[,1]
	g <- alphabetize.factor(Data$g)
    }
   else {
     g <- c()
     for (i in 1:length(x)) {
       for (j in 1:length(x[[i]])) {
         g <- c(g,i)
         }
       }
    x <- as.numeric(unlist(x)[!is.na(unlist(x))])
    }

  # validate method
  if (length(method) > 1) {
    method <- "none"
    }
    
  if (tolower(method) != "none" & tolower(method) != "bonferroni" & tolower(method) != "sidak" & tolower(method) != "holm" & tolower(method) != "hs" & tolower(method) != "hochberg" & tolower(method) != "bh" & tolower(method) != "by") {
    rlang::abort(message="method must be one of: none, bonferroni, sidak, hs, bh, or by")
    }

  if (tolower(method)=="none") {
    Name <- "(No adjustment)"
    }
  if (tolower(method)=="bonferroni") {
    Name <- "(Bonferroni)"
    }
  if (tolower(method)=="sidak") {
    Name <- "(\u0160id\u00E1k)"
    }
  if (tolower(method)=="holm") {
    Name <- "(Holm)"
    }
  if (tolower(method)=="hs") {
    Name <- "(Holm-\u0160id\u00E1k)"
    }
  if (tolower(method)=="hochberg") {
    Name <- "(Hochberg)"
    }
  if (tolower(method)=="bh") {
    Name <- "(Benjamini-Hochberg)"
    }
  if (tolower(method)=="by") {
    Name <- "(Benjamini-Yekuteili)"
    }

  # validate that x is longer than 1
  if (length(unlist(x))==1) {
    rlang::abort(message="too few observations in x.")
    }
  # validate that x is a numeric vector of data values, or a list of numeric 
  # data vectors, and is not NA
  if (TRUE %in% is.na(unlist(x)) | !is.vector(unlist(x), mode="numeric") ) {
    rlang::abort(message="x must contain a numeric vector of data values, or a list of numeric data vectors.")
    }
  # validate that g is not missing if x is a vector
  if (!xaslist & TRUE %in% is.na(g) ) {
    rlang::abort(message="when specifying x as a vector, you must include g.")
    }
  # validate that g is not NA
  if (length(g) > 1 & TRUE %in% is.na(g)) {
    rlang::abort(message="g must have no missing values.")
    }
  # validate that g is factor or vector.
  if (length(g) > 1 & ( !is.factor(g) & !is.vector(g) ) ) {
    
    rlang::abort(message="g must be a factor, character vector, or integer vector.")
    }
  # validate that g is a vector of mode = character or mode = integer.
  if (length(g) > 1 & is.vector(g) ) {
    if ( !is.vector(g, mode="character") & !all.integers(g) ) {
      rlang::abort(message="g must be a factor, character vector, or integer vector.")
      }
    }

  # CALCULATIONS
  out <- NULL
  if (xaslist & length(g)==1) {
    kwallis.test(x, 1:length(x), pass=1) -> out
    }
  if (length(g)>1 & xaslist) {
    kwallis.test(x, g, pass=1) -> out
    }
  if (length(g)>1 & !xaslist) {
    kwallis.test(x, g, pass=0) -> out
    }

  if (kw==TRUE) {
    rlang::inform(message=paste0("  Kruskal-Wallis rank sum test\n\ndata: ", xname, " and ", gname, sep=""))
    rlang::inform(message=out$output)
    }
    chi2 <- out$H
    df   <- out$df
    p    <- out$p
    N    <- out$N
    Data <- out$Data
    k    <- df+1
    m    <- k*(k-1)/2
    Z    <- rep(NA,m)
    P    <- rep(NA,m)

  #calculate ties adjustment to be used in pooled variance estimate later
  ranks <- Data[,4]
  ties <- tiedranks(ranks)
  r <- length(ties)
  tiesadj <- 0
  if (r > 0) {
    for (s in 1:r) {
      tau <- sum(ranks==ties[s])
      tiesadj <- tiesadj + (tau^(3) - tau)
      }
    }
  tiesadj <- tiesadj/(12*(N-1))

  # Generate differences in mean ranks, standard deviation of same, and z statistic
  Y      <- rep(NA,m)
  Sigma  <- rep(NA,m)
  row    <- c()
  col    <- c()
  for (i in 1:(k-1)) {
    row <- c(row,rep(i,(k-i)))
    col <- c(col,(i+1):k)
    }

  # Calculate approximate Z test statistics
  Z <- rep(0,m)
  for (i in 2:k) {
    for (j in 1:(i-1)) {
      ni <- sum(Data[,3]==i)
      nj <- sum(Data[,3]==j)
      meanranki <- mean(Data[,4][Data[,3]==i])
      meanrankj <- mean(Data[,4][Data[,3]==j])
      if (rmc==TRUE) {
        z <- (meanranki - meanrankj) / sqrt( ((N*(N+1)/12) - tiesadj) * ((1/nj) + (1/ni)) )
        }
       else {
        z <- (meanrankj - meanranki) / sqrt( ((N*(N+1)/12) - tiesadj) * ((1/nj) + (1/ni)) )
        }
      index <- ((i-2)*(i-1)/2) + j
      Z[index] <- z
      }
    }
    
  # Calculate p-values for Z statistics, and adjust as needed
  # If p = P(Z >= |z|)
  if (altp==FALSE) {
    P <- pnorm(abs(Z), lower.tail=FALSE)
    }
   # Otherwise, if p = P(|Z| >= |z|)
   else {
   	P <- 2*pnorm(abs(Z), lower.tail=FALSE)
   	}
  #calculate adjusted p-values based on method argument
  Reject <- rep(0,m)
  # No adjustment for multiple comparisons
  if (tolower(c(method))=="none") {
    P.adjust <- P
    # If p = P(Z >= |z|)
    if (altp==FALSE) {
      Reject <- P.adjust <= alpha/2
      }
     # Otherwise, if p= P(|Z| >= |z|)
     else {
      Reject <- P.adjust <= alpha
      }
    }
  # Control FWER using (Dunn's) Bonferroni
  if (tolower(c(method))=="bonferroni") {
    P.adjust <- pmin(1,P*m)
    # If p = P(Z >= |z|)
    if (altp==FALSE) {
      Reject <- P.adjust <= alpha/2
      }
     # Otherwise, if p= P(|Z| >= |z|)
     else {
      Reject <- P.adjust <= alpha
      }
    }
  # Control FWER using Šidák
  if (tolower(c(method))=="sidak") {
    P.adjust <- pmin(1,1 - (1-P)^m)
    # If p = P(Z >= |z|)
    if (altp==FALSE) {
      Reject <- P.adjust <= alpha/2
      }
     # Otherwise, if p= P(|Z| >= |z|)
     else {
      Reject <- P.adjust <= alpha
      }
    }
  # Control FWER using Holm(-Bonferroni)
  if (tolower(c(method))=="holm") {
    Psort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    Psort <- Psort[,order(Psort[1,])]
    for (i in 1:m) {
      adjust <- m+1-i
      Psort[1,i] <- pmin(1, Psort[1,i]*adjust)
      # If p = P(Z >= |z|)
      if (altp==FALSE) {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha/2
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha/2) & Psort[3,i-1] != 0)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|)
       else {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha) & Psort[3,i-1] != 0)
          }
        }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }
  # Control FWER using Holm-Šidák
  if (tolower(c(method))=="hs") {
    Psort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    Psort <- Psort[,order(Psort[1,])]
    for (i in 1:m) {
      adjust <- m+1-i
      Psort[1,i] <- pmin(1, (1 - ((1 - Psort[1,i])^adjust)))
      # If p = P(Z >= |z|)
      if (altp==FALSE) {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha/2
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha/2) & Psort[3,i-1] != 0)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|)
       else {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha) & Psort[3,i-1] != 0)
          }
        }      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }
  # Control FWER using Hochberg
  if (tolower(c(method))=="hochberg") {
    Psort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- i
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      # If p = P(Z >= |z|)
      if (altp==FALSE) {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha/2
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|)
       else {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha) | Psort[3,i-1] == 1)
          }
        }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }
  # Control FDR using Benjamini-Hochberg
  if (tolower(c(method))=="bh") {
    Psort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- (m/(m+1-i))
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      # If p = P(Z >= |z|)
      if (altp==FALSE) {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha/2
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|)
       else {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha) | Psort[3,i-1] == 1)
          }
        }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }
  # Control FDR using Benjamini-Yekuteili
  if (tolower(c(method))=="by") {
    Psort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- (m/(m+1-i))*sum(1/(1:m))
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      # If p = P(Z >= |z|)
      if (altp==FALSE) {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha/2
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|)
       else {
        if (i==1) {
          Psort[3,i] <- Psort[1,i] <= alpha
          }
         else {
          Psort[3,i] <- ((Psort[1,i] <= alpha) | Psort[3,i-1] == 1)
          }
        }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # OUTPUT
  rlang::inform(message="")
  if (table==TRUE | list==TRUE) {
    if ((TRUE %in% is.na(g))) {
      title <- paste("Dunn's Pairwise Comparison of ",xname," across ",k," groups",sep="")
      }
     else {
      title <- paste("Dunn's Pairwise Comparison of ",xname," by ",gname,sep="")
      }
    rlang::inform(message=centertext(title))
    rlang::inform(message=paste0(centertext(Name), sep=""))
    }

  if (table==TRUE) {
    # Need to determine how many tables (reps) to output
    reps      <- floor((k-1)/6)
    laststart <- (reps*6) + 1
    kminusone <- k - 1
    if (label==FALSE) {
      g <- as.numeric(g)
      }
    if (length(g)==1) {
      g <- 1:k
      }
  
    # Replication loop for >7 groups, no wrap
    if (wrap==FALSE) {
      if (k > 7) {
        for (rep in 1:reps) {
          colstart <- (6*rep)-5
          colstop  <- 6*rep
          dunntestheader(g,colstart,colstop,rmc)
          # Table body
          for (i in (colstart+1):k) {
            colstop <- min(i-1,6*rep)
            dunntestztable(g,i,Z,colstart,colstop)
            if (i < k) {
              dunntestptable(i,P.adjust,colstart,colstop,Reject,0)
              }
             else {
              dunntestptable(i,P.adjust,colstart,colstop,Reject,1)
              }
            }
          }
        # End of table
        if (laststart < k) {
          dunntestheader(g,laststart,kminusone,rmc)
          # Table body
          for (i in (laststart+1):k) {
          	colstop <- i-1
            dunntestztable(g,i,Z,laststart,colstop) 
            if (i < k) {
              dunntestptable(i,P.adjust,laststart,colstop,Reject,0)
              }
             else {
              dunntestptable(i,P.adjust,laststart,colstop,Reject,1)
              }
            }
          }
        }
  
      # Replication loop for <=7 groups
      if (k <= 7) {
        dunntestheader(g,1,kminusone,rmc)
        # Table body
        for (i in 2:k) {
          colstop <- i-1
          dunntestztable(g,i,Z,1,colstop) 
          if (i < k) {
            dunntestptable(i,P.adjust,1,colstop,Reject,0)
            }
           else {
            dunntestptable(i,P.adjust,1,colstop,Reject,1)
            }
          }
        }
      }
  
    # Replication loop for >7 groups, with wrap
    if (wrap==TRUE) {
      dunntestheader(g,1,kminusone,rmc)
      # Table body
      for (i in 2:k) {
        colstop <- i-1
        dunntestztable(g,i,Z,1,colstop) 
        if (i < k) {
          dunntestptable(i,P.adjust,1,colstop,Reject,0)
          }
         else {
          dunntestptable(i,P.adjust,1,colstop,Reject,1)
          }
        }
      }
      
    rlang::inform(message="")
    }

  # Output pairwise comparisons as a list if requested.
  if (list==TRUE) {
    groupvalues  <- levels(factor(g))
    # get the lengths of each group name (whether explicitly labeled or not)
    sort(unlist(lapply(groupvalues,nchar))) -> lengths
    # stringlength will be the sum of the two largest values in lengths plus 6
    stringlength <- lengths[length(lengths)] + lengths[length(lengths)-1] + 6

    # Output list header
    if (tolower(method)=="none") {
      rlang::inform(message="List of pairwise comparisons: Z statistic (p-value)")
      }
     else {
      rlang::inform(message="List of pairwise comparisons: Z statistic (adjusted p-value)")
      }
    list_head_length <- max(nchar(groupvalues)) + tail(sort(nchar(unique(groupvalues))), n=2)[1] + 4
    rlang::inform(message=paste0(paste0(rep("\U2500", list_head_length), collapse=""), "\U252C", paste0(rep("\U2500", 19+max(Reject)), collapse=""), collapse="", sep=""))
    index <- 0
    for (i in 2:k) {
      for (j in 1:(i-1)) {
        index <- index + 1
        buffer <- list_head_length - (nchar(groupvalues[i]) + nchar(groupvalues[j]) + 4)
        if ( Reject[index] == 0) {
          pformatted <- paste0("(", sprintf("%1.4f", P.adjust[index]), ")", sep="")
          }
         else {
         	pformatted <- paste0("(", sprintf("%1.4f", P.adjust[index]), ")", "*", sep="")
         	}
        if (rmc==FALSE) {
          rlang::inform(message=paste0(groupvalues[j]," - ",groupvalues[i],paste0(rep(" ",buffer), collapse="")," \U2502 ",zformat(Z[index]), " ",pformatted, sep=""))
          }
        if (rmc==TRUE) {
          rlang::inform(message=paste0(groupvalues[i]," - ",groupvalues[j],paste0(rep(" ",buffer), collapse="")," \U2502 ",zformat(Z[index]), " ",pformatted, sep=""))
          }
        }
      }
  if (interpret==TRUE) {
      rlang::inform(message="")
    }
  }
  symbol <- "\U03B1"
  procedure <- ""
  if (method == "holm" | method == "sidak" | method == "hs" | method == "hochberg") symbol <- "FWER" 
  if (method == "bh" | method == "by") symbol <- "FDR"
  if (method == "holm" | method == "hs" | method == "hochberg" | method == "bh" | method == "by") procedure <- " with stopping rule" 
  if (interpret==TRUE) {
    rlang::inform(message=paste0(symbol, " = ",alpha, sep=""))
    # If p = P(Z >= |z|)
    if (altp==FALSE) {
      rlang::inform(message=paste0("Reject Ho if p \U2264 ", symbol, "/2", procedure, ", where p = Pr(Z \U2265 |z|)"))
      }
     # Otherwise, if p = P(|Z| <= |z|)
     else {
      rlang::inform(message=paste0("Reject Ho if p \U2264 ", symbol, procedure, ", where p = Pr(|Z| \U2265 |z|)"))
      }
    }

  # Create comparisons variable for returned values (whether the list option
  # is TRUE or FALSE
  comparisons <- rep(NA,(k*(k-1)/2))
  groupvalues  <- levels(factor(g))
  index <- 0
  for (i in 2:k) {
    for (j in 1:(i-1)) {
      index <- index + 1
      if (rmc==FALSE) {
        comparisons[index] <- paste0(groupvalues[j]," - ",groupvalues[i])
        }
      if (rmc==TRUE) {
        comparisons[index] <- paste0(groupvalues[i]," - ",groupvalues[j])
        }
      }
    }
   
  # If p = P(Z <= |z|)
  if (altp==FALSE) {
   invisible(list(chi2=chi2, Z=Z, P=P, P.adjusted=P.adjust, comparisons=comparisons))
   }
   # Otherwise, if p = P(|Z| <= |z|)
   else {
   invisible(list(chi2=chi2, Z=Z, altP=P, altP.adjusted=P.adjust, comparisons=comparisons))
   }
  }
