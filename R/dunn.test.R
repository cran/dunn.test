# Perform Dunn's test of multiple comparisons using rank sums
# Version: 1.4.0
# Author: Alexis Dinno
# Date: May 23, 2026

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
    interpret=TRUE, 
    single.contrast=NA
    ) {

  ##########
  # VALIDATIONS & PREPARATIONS
  # names for output
  xname <- paste(if(is.name(substitute(x))) {deparse(substitute(x))} else {"x"})
  gname <- paste(if(is.name(substitute(g))) {deparse(substitute(g))} else {"group"})
  xaslist <- is.list(x)
  tempg <- 1:length(x)

  # casewise deletion of missing data
  # Start dealing with the fact that g may be numeric, string, or factor data
  if (length(g) > 1) {  

    # Case-wise deletion of incomplete data
    x_miss <- is.na(x)
    g_miss <- is.na(g)
    all_miss <- (x_miss | g_miss)
    x <- x[!all_miss]
    g <- g[!all_miss]
    if (label==TRUE) {
      if (is.factor(g)) {
        glevels <- levels(factor(g, levels=unique(factor(g))))
        }
       else {
       	g <- factor(g, levels=unique(factor(g)))
        glevels <- levels(g)
        }
      }
     else {
      glevels <- names(table(addNA(g, ifany = TRUE)))
      }
    Data <- data.frame(x,g)[get_order(sort(glevels), g),]
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
  # Check if single.contrast is not the empty string, and if it is a number, or 
  # should be treated as a label string, prepare a variable to hold either the 
  # value or the label, and change the value of single.contrast to indicate 
  # which
  single.index <- NA
  single.contrast.type <- ""
  single.contrast.error.message <- "the single.contrast argument must contain either a numeric value within g, a \nlevel of g if g is a factor, a string value within g, or numerically indicate \nwhich single group within x will be the contrast if g is omitted."
  if (!is.na(single.contrast)) {
    # Validate that single.contrast is a single value
    if (length(single.contrast) > 1) {
      rlang::abort(single.contrast.error.message)
      }
    # Validate that single.contrast is a single value within g, or that it 
    # indexes a group within x
    if (!(FALSE %in% !is.na(g)) & !(single.contrast %in% g)) {
      rlang::abort(single.contrast.error.message)
      }
    if ((FALSE %in% !is.na(g)) & (single.contrast > length(x) | single.contrast < 1)) {
      rlang::abort(single.contrast.error.message)
      }
    # Assign single.value if single.contrast is numeric
    if (is_numeric_like(single.contrast)) {
      single.value <- as.numeric(single.contrast)
      single.contrast.type <- "single.value"
      if (!(FALSE %in% !is.na(g))) {
        single.label <- paste0(gname, "=", single.contrast, sep="") 
        }
      }
     else {
      single.label = single.contrast 
      single.contrast.type = "single.label"
      }
    }        # Close if (!is.na(single.contrast)) 
  # Index single.value in g or x if single.contrast is numeric
  if (single.contrast.type == "single.value") {
    if (!(FALSE %in% !is.na(g))) {
      single.index <- match(single.contrast, unique(g))
      k = length(g)
      }
    if ((FALSE %in% !is.na(g))) {
      single.index <- single.value
      k = length(x)
      }
    }        # Close if (single.contrast.type == "single.value")

  # Index single.label in g if single.contrast is string
  if (single.contrast.type == "single.label") {
    single.index <- match(single.contrast, unique(g))
    k = length(g)
    }        # Close if (single.contrast.type == "single.label")

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
  if (!xaslist & TRUE %in% (FALSE %in% !is.na(g)) ) {
    rlang::abort(message="when specifying x as a vector, you must include g.")
    }
  # validate that g is not NA
  if (length(g) > 1 & TRUE %in% (FALSE %in% !is.na(g))) {
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

  ##########
  # CALCULATIONS
  
  out <- NULL
  if (xaslist & length(g)==1) {
    kwallis.test(x, 1:length(x)) -> out
    }
  if (length(g)>1 & xaslist) {
    kwallis.test(x, g) -> out
    }
  if (length(g)>1 & !xaslist) {
    kwallis.test(x, g) -> out
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
    if (is.na(single.contrast)) {
      m    <- k*(k-1)/2
      }
    if (!is.na(single.contrast)) {
      m    <- k-1
      }
    Z    <- rep(NA,m)
    P    <- rep(NA,m)

  # Prepare to modify behavior if single.contrast has been specified
  if (!is.na(single.contrast)) {
    table <- FALSE
    list <- TRUE
    rmc <- FALSE
    }

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
  if (is.na(single.contrast)) {
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
      }        # Close for (i in 2:k)
    }        # Close if (is.na(single.contrast))
   else {
    # Calculate comparisons only between the group in single.comparison and
    # each remaining group
    for (j in 1:k) {
      if (single.index != j) {
        jindex <- j
        if (single.index < j) {
          jindex <- j - 1
          }
        ni <- sum(Data[,3]==single.index)
        nj <- sum(Data[,3]==j)
        meanranki <- mean(Data[,4][Data[,3]==single.index])
        meanrankj <- mean(Data[,4][Data[,3]==j])
        z <- (meanranki - meanrankj) / sqrt( ((N*(N+1)/12) - tiesadj) * ((1/nj) + (1/ni)) )
        if (jindex < k) {
          Z[jindex] <- z
          }
        }        # Close if (single.index != j)
      }        # Close for (j in 1:k)
    }        # Close else
    
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
  MCA <- multiple_comparisons_adjustment(P=P, m=m, method=method, alpha=alpha, altp=altp)
  P.adjust <- MCA$P.adjust
  Reject <- MCA$Reject
  
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
  # If single.contrast is specified then:
  # suppress table
  # output list
  # set rmc is FALSE
  # notify user
  if (!is.na(single.contrast)) {
    if ((TRUE %in% is.na(g))) {
      single.contrast.title <- paste("\nComparisons limited to group number ",single.contrast," vs every other group",sep="")
      }
     else {
      single.contrast.title <- paste("\nComparisons limited to ",gname,"=",single.contrast," vs every other group",sep="")
      }
    rlang::inform(message=single.contrast.title)
    }
  if (table==TRUE) {
    rlang::inform(message="")
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
    }

  # Prepare to modify behavior if single.contrast has been specified
  if (is.na(single.index)) {
    kindex    <- k
    liststart <- 2
    offset    <- 1
    }
   else {
    kindex    <- k - 1
    liststart <- 1
    offset    <- 0
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
      rlang::inform(message="\nList of pairwise comparisons: Z statistic (p-value)")
      }
     else {
      rlang::inform(message="\nList of pairwise comparisons: Z statistic (adjusted p-value)")
      }
    list_head_length <- max(nchar(groupvalues)) + tail(sort(nchar(unique(groupvalues))), n=2)[1] + 4
    rlang::inform(message=paste0(paste0(rep("\U2500", list_head_length), collapse=""), "\U252C", paste0(rep("\U2500", 19+max(Reject)), collapse=""), collapse="", sep=""))
    index <- 0
    for (i in liststart:k) {
      for (j in 1:(i-offset)) {
        # The below conditional allows this nested loop to work when 
        # single.contrast has been specified
        if ((single.index != i & i == j) | is.na(single.index)) {
          index <- index + 1
          labeli <- groupvalues[i]
          if (is.na(single.index)) {
            labelj <- groupvalues[j]
            }
          else {
            labelj <- groupvalues[single.index]
            }
          buffer <- list_head_length - (nchar(labeli) + nchar(labelj) + 4)
          if ( Reject[index] == 0) {
            pformatted <- paste0("(", sprintf("%1.4f", P.adjust[index]), ")", sep="")
            }
           else {
           	pformatted <- paste0("(", sprintf("%1.4f", P.adjust[index]), ")", "*", sep="")
           	}
          if (rmc==FALSE) {
            rlang::inform(message=paste0(labelj," - ",labeli,paste0(rep(" ",buffer), collapse="")," \U2502 ",zformat(Z[index]), " ",pformatted, sep=""))
            }
          if (rmc==TRUE) {
            rlang::inform(message=paste0(labeli," - ",labelj,paste0(rep(" ",buffer), collapse="")," \U2502 ",zformat(Z[index]), " ",pformatted, sep=""))
            }
          }
        }
      }
    }

  # Manage interpretation output
  symbol <- "\U03B1"
  procedure <- ""
  unadj <- ""
  adjust <- ""
  if (method == "bonferroni" | method == "sidak" | method == "holm" | method == "sidak" | method == "hs" | method == "hochberg") {
    symbol <- "FWER" 
    adjust <- "adjusted "
    unadj <- "(unadjusted) "
    }
  if (method == "bh" | method == "by") {
    symbol <- "FDR"
    adjust <- "adjusted "
    unadj <- "(unadjusted) "
    }
  if (method == "holm" | method == "hs" | method == "hochberg" | method == "bh" | method == "by") procedure <- " with stopping rule" 
  if (interpret==TRUE) {
    rlang::inform(message=paste0("\n",symbol, " = ",alpha, sep=""))
    # If p = P(Z >= |z|)
    if (altp==FALSE) {
      rlang::inform(message=paste0("Reject Ho if ", adjust, "p \U2264 ", symbol, "/2", procedure, ", where ", unadj, "p = Pr(Z \U2265 |z|)"))
      }
     # Otherwise, if p = P(|Z| <= |z|)
     else {
      rlang::inform(message=paste0("Reject Ho if ", adjust, "p \U2264 ", symbol, procedure, ", where ", unadj, "p = Pr(|Z| \U2265 |z|)"))
      }
    }

  # Create comparisons variable for returned values (whether the list option
  # is TRUE or FALSE
  comparisons <- rep(NA,m)
  groupvalues  <- levels(factor(g))
    index <- 0
    for (i in liststart:k) {
      for (j in 1:(i-offset)) {
        # The below conditional allows this nested loop to work when 
        # single.contrast has been specified
        if ((single.index != i & i == j) | is.na(single.index)) {
          index <- index + 1
          labeli <- groupvalues[i]
          if (is.na(single.index)) {
            labelj <- groupvalues[j]
            }
          else {
            labelj <- groupvalues[single.index]
            }
          if (rmc==FALSE) {
            comparisons[index] <- paste0(labelj," - ",labeli)
            }
          if (rmc==TRUE) {
            comparisons[index] <- paste0(labeli," - ",labelj)
            }
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
