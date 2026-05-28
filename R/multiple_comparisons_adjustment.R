# multiple_comparisons_adjustment: adjusts p-values for multiple comparisons 
# appropriate to the provided:
#   vector of p-values P,
#   number of pairwise comparisons m,
#   multiple adjustment procedure method,
#   Type I error rate alpha, and
#   p-value convention altp
# Returns: vector p-values and vector rejection decisions
# Date: March 31, 2026

multiple_comparisons_adjustment <- function(P, m, method, alpha, altp) {
  # No adjustment for multiple comparisons
  if (tolower(c(method))=="none") P.adjust <- P

  # Control FWER using (Dunn's) Bonferroni
  if (tolower(c(method))=="bonferroni") P.adjust <- pmin(1,P*m)

  # Control FWER using Šidák
  if (tolower(c(method))=="sidak") P.adjust <- pmin(1,1 - (1-P)^m)

  # Rejection decision for non-stepped methods
  if (tolower(c(method)) == "none" | tolower(c(method)) == "bonferroni" | tolower(c(method)) == "sidak") {
    # If p = P(Z >= |z|)
    if (altp == FALSE) {
      Reject <- P.adjust <= alpha/2
      }
     # Otherwise, if p= P(|Z| >= |z|)
     else {
      Reject <- P.adjust <= alpha
      }
    }        # Close if method none, bonferroni, or sidak

  # Rejection decision preparation for step-down methods
  if (tolower(c(method)) == "holm" | tolower(c(method)) == "hs") {
    PSort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    PSort <- PSort[,order(PSort[1,])]
    for (i in 1:m) {
      adjust <- m+1-i
      if (tolower(c(method)) == "holm") {
        # Holm adjustment
        PSort[1,i] <- pmin(1, PSort[1,i]*adjust)
        }
       else {
        # Holm-Sidak adjustment
        PSort[1,i] <- pmin(1, (1 - ((1 - PSort[1,i])^adjust)))
        }
      # If p = P(Z >= |z|) via step-up
      if (altp==FALSE) {
        if (i==1) {
          PSort[3,i] <- PSort[1,i] <= alpha/2
          }
         else {
          PSort[3,i] <- ((PSort[1,i] <= alpha/2) & PSort[3,i-1] != 0)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|) via step-up
       else {
        if (i==1) {
          PSort[3,i] <- PSort[1,i] <= alpha
          }
         else {
          PSort[3,i] <- ((PSort[1,i] <= alpha) & PSort[3,i-1] != 0)
          }
        }
      }        # Close for(i in 1:m)
    }        # Close if step-down methods
  # Rejection decision preparation for step-up methods
  if (tolower(c(method)) == "hochberg" | tolower(c(method)) == "bh" | tolower(c(method)) == "by") {
    PSort <- matrix(c(P,1:m,rep(0,m)), 3, m, byrow=TRUE)
    PSort <- PSort[,order(PSort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      if (tolower(c(method))=="hochberg") {
        adjust <- i
        PSort[1,i] <- min(1,PSort[1,i]*adjust)
        }
      if (tolower(c(method))=="bh") {
        adjust <- (m/(m+1-i))
        PSort[1,i] <- min(1,PSort[1,i]*adjust)
        }
      if (tolower(c(method))=="by") {
        adjust <- (m/(m+1-i))*sum(1/(1:m))
        PSort[1,i] <- min(1,PSort[1,i]*adjust)
        }
      # If p = P(Z >= |z|) via step-down
      if (altp==FALSE) {
        if (i==1) {
          PSort[3,i] <- PSort[1,i] <= alpha/2
          }
         else {
          PSort[3,i] <- ((PSort[1,i] <= alpha/2) | PSort[3,i-1] == 1)
          }
        }
       # Otherwise, if p = P(|Z| >= |z|) via step-down
       else {
        if (i==1) {
          PSort[3,i] <- PSort[1,i] <= alpha
          }
         else {
          PSort[3,i] <- ((PSort[1,i] <= alpha) | PSort[3,i-1] == 1)
          }
        }
      }        # Close if(i in 1:m)
    }        # Close if step-up methods
  # Rejection decisions for stepped methods
  if (tolower(c(method)) == "holm" | tolower(c(method)) == "hs" | tolower(c(method)) == "hochberg" | tolower(c(method)) == "bh" | tolower(c(method)) == "by") {
    PSort <- PSort[,order(PSort[2,])]
    P.adjust <- PSort[1,]
    Reject <- PSort[3,]
    }        # Close if stepped methods
  invisible(list(Reject=Reject, P.adjust=P.adjust))
  }        # Close multiple_comparisons_adjustment <- function()