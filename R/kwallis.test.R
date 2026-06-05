# kwallis.test: a custom Kruskal Wallis test function to support dunn.test
# Date: June 5, 2026
kwallis.test <- function(x=NA, g=NA) {

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
  for (i in unique(Data[,3])) {
    ranksum <- ranksum + ((sum(as.numeric(Data[,4][Data[,3]==i])))^2)/(sum(as.numeric(Data[,3]==i)))
    }
  H  <- ((12/(N*(N+1)))*ranksum - 3*(N+1))/tiesadj
  df <- k-1
  p  <- pchisq(H,k-1,lower.tail=FALSE)

  #present output
  output <- paste("Kruskal-Wallis chi-squared = ",round(H,digits=4),", df = ",df,", p-value = ",round(p,digits=2),"\n" ,sep="")
  
  invisible(list(output=output,H=H,df=df,p=p,N=N,Data=Data))
  }