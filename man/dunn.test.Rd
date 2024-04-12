\encoding{utf8}
\name{dunn.test}
\alias{Dunn's test}
\alias{dunn.test}
\alias{p.adjustment.methods}
\concept{Dunn's test}
\concept{nonparametric pairwise multiple comparisons}
\title{Dunn's Test}
\description{Performs Dunn's test of multiple comparisons using rank sums}
\usage{
dunn.test  (x, g=NA, method=p.adjustment.methods, kw=TRUE, label=TRUE, 
      wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

p.adjustment.methods
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")
}
\arguments{
 \item{x}{a numeric vector, or a list of numeric vectors.  Missing values are ignored.  If the former, then groups \emph{must} be specified using \code{g}.}
 \item{g}{a factor variable, numeric vector, or character vector indicating group.  Missing values are ignored.}
 \item{method}{adjusts the \emph{p}-value for multiple comparisons using the Bonferroni, \enc{Šidák}{Sidak}, Holm, Holm-\enc{Šidák}{Sidak}, Hochberg, Benjamini-Hochberg, or Benjamini-Yekutieli adjustment (see Details).  The default is no adjustment for multiple comparisons.}
 \item{kw}{if \option{TRUE} then the results of the Kruskal-Wallis test are reported.}
 \item{label}{if \option{TRUE} then the factor labels are used in the output table.}
 \item{wrap}{does not break up tables to maintain nicely formatted output.  If \option{FALSE} then output of large tables is broken up across multiple pages.}
 \item{table}{outputs results of Dunn's test in a table format, as qualified by the \code{label} and \code{wrap} options.}
 \item{list}{outputs results of Dunn's test in a list format.}
 \item{rmc}{if \option{TRUE} then the reported test statistics and table are based on row minus column, rather than the default column minus row (i.e. the signs of the test statistic are flipped).}
 \item{alpha}{the nominal level of significance used in the step-up/step-down multiple comparisons procedures (Holm, Holm-\enc{Šidák}{Sidak}, Hochberg, Benjamini-Hochberg, and Benjamini-Yekutieli).}
 \item{altp}{if \option{TRUE} then express \emph{p}-values in alternative format.  The default is to express \emph{p}-value = P(\emph{Z} \eqn{\ge} |\emph{z}|), and reject Ho if \emph{p} \eqn{\le} \eqn{\alpha}/2.  When the \code{altp} option is used, \emph{p}-values are instead expressed as \emph{p}-value = P(|\emph{Z}| \eqn{\ge} |\emph{z}|), and reject Ho if \emph{p} \eqn{\le} \eqn{\alpha}.  \strong{These two expressions give identical test results.} Use of \code{altp} is therefore merely a semantic choice.}
}
\details{\code{dunn.test} computes Dunn's test (1964) for stochastic dominance and reports the results among multiple pairwise comparisons after a Kruskal-Wallis test for stochastic dominance among \emph{k} groups (Kruskal and Wallis, 1952).  The interpretation of stochastic dominance requires an assumption that the \acronym{CDF} of one group does not cross the \acronym{CDF} of the other.  \code{dunn.test} makes \emph{m} = \emph{k}(\emph{k}-1)/2 multiple pairwise comparisons based on Dunn's \emph{z}-test-statistic approximations to the actual rank statistics.  The null hypothesis for each pairwise comparison is that the probability of observing a randomly selected value from the first group that is larger than a randomly selected value from the second group equals one half; this null hypothesis corresponds to that of the Wilcoxon-Mann-Whitney rank-sum test.  Like the rank-sum test, if the data can be assumed to be continuous, and the distributions are assumed identical except for a difference in location, Dunn's test may be understood as a test for median difference.  \code{dunn.test} accounts for tied ranks.  

\code{dunn.test} outputs both \emph{z}-test-statistics for each pairwise comparison and the \emph{p}-value = P(\emph{Z} \eqn{\ge} |\emph{z}|) for each. Reject Ho based on \emph{p} \eqn{\le} \eqn{\alpha}/2 (and in combination with \emph{p}-value ordering for stepwise \code{method} options).  If you prefer to work with \emph{p}-values expressed as \emph{p}-value = P(|\emph{Z}| \eqn{\ge} |\emph{z}|) use the \code{altp=TRUE} option, and reject Ho based on \emph{p} \eqn{\le} \eqn{\alpha} (and in combination with \emph{p}-value ordering for stepwise \code{method} options).  These are exactly equivalent rejection decisions.

Several options are available to adjust \emph{p}-values for multiple comparisons, including methods to control the family-wise error rate (\acronym{FWER}) and methods to control the false discovery rate (\acronym{FDR}):

\tabular{rl}{
\option{none} \tab no adjustment is made.  Those comparisons rejected without adjustment at the \eqn{\alpha} \cr
\tab level (two-sided test) are starred in the output table, and starred in the list when \cr
\tab using the \code{list=TRUE} option.\cr\cr
\option{bonferroni} \tab the \acronym{FWER} is controlled using Dunn's (1961) Bonferroni adjustment, and adjusted \cr
\tab \emph{p}-values = max(1, \emph{pm}).  Those comparisons rejected with the Bonferroni \cr
\tab adjustment at the \eqn{\alpha} level (two-sided test) are starred in the output table, and \cr
\tab starred in the list when using the \code{list=TRUE} option.\cr\cr
\option{sidak}\tab the \acronym{FWER} is controlled using Šidák's (1967) adjustment, and adjusted \cr
\tab \emph{p}-values = max(1, 1 - (1 - \emph{p})^\emph{m}).  Those comparisons rejected with the Šidák \cr
\tab adjustment at the \eqn{\alpha} level (two-sided test) are starred in the output table, and \cr
\tab starred in the list when using the \code{list=TRUE} option.\cr\cr
\option{holm}\tab the \acronym{FWER} controlled using Holm's (1979) progressive step-up procedure to relax \cr
\tab control on subsequent tests.  \emph{p} values are ordered from smallest to largest, and \cr
\tab adjusted \emph{p}-values = max[1, \emph{p}(\emph{m}+1-\emph{i})], where \emph{i} indexes the ordering.  All tests after \cr
\tab and including the first test to not be rejected are also not rejected.\cr\cr
\option{hs}\tab the \acronym{FWER} is controlled using the Holm-Šidák adjustment (Holm, 1979): another \cr
\tab progressive step-up procedure but assuming dependence between tests.  \emph{p} values \cr
\tab are ordered from smallest to largest, and adjusted \cr
\tab \emph{p}-values = max[1, 1 - (1 - \emph{p})^(\emph{m}+1-\emph{i})], where \emph{i} indexes the ordering.  All tests after \cr
\tab and including the first test to not be rejected are also not rejected.\cr\cr
\option{hochberg}\tab the \acronym{FWER} is controlled using Hochberg's (1988) progressive step-down \cr
\tab procedure to increase control on successive tests.  \emph{p} values are ordered from largest–\cr
\tab smallest, and adjusted \emph{p}-values = max[1, \emph{p}*\emph{i}], where \emph{i} indexes the ordering.  All \cr
\tab tests after and including the first to be rejected are also rejected.\cr\cr
\option{bh}\tab the \acronym{FDR} is controlled using the Benjamini-Hochberg adjustment (1995), a step-\cr
\tab down procedure appropriate to independent tests or tests that are positively \cr
\tab dependent.  \emph{p}-values are ordered from largest to smallest, and adjusted \cr
\tab \emph{p}-values = max[1, \emph{pm}/(\emph{m}+1-\emph{i})], where \emph{i} indexes the ordering.  All tests after and \cr
\tab including the first to be rejected are also rejected.\cr\cr
\option{by}\tab the \acronym{FDR} is controlled using the Benjamini-Yekutieli adjustment (2011), a step-\cr
\tab down procedure appropriate to dependent tests.  \emph{p}-values are ordered from largest to \cr
\tab smallest, and adjusted \emph{p}-values = max[1, \emph{pmC}/(\emph{m}+1-\emph{i})], where \emph{i} indexes the \cr
\tab ordering, and the constant \emph{C} = 1 + 1/2 + .  .  .  + 1/\emph{m}.  All tests after and including the \cr
\tab first to be rejected are also rejected.\cr
}
Because the sequential step-up/step-down tests rejection decisions depend on both the \emph{p}-values and their ordering, those tests rejected using \code{"holm"}, \code{"hs"}, \code{"hochberg"}, \code{"bh"}, or \code{"by"} at the indicated \eqn{\alpha} level are starred in the output table, and starred in the list when using the \code{list=TRUE} option.  
}
\value{
\code{dunn.test} returns:
 \item{chi2}{a scalar of the Kruskal-Wallis test statistic adjusted for ties.}
 \item{Z}{a vector of all \emph{m} of Dunn \emph{z} test statistics.}
 \item{P}{a vector of \emph{p}-values corresponding to \code{Z}. --OR--}
 \item{altP}{a vector of \emph{p}-values corresponding to \code{Z} when using the \code{altp=TRUE} option.}
 \item{P.adjust}{a vector of \emph{p}-values corresponding to \code{Z}, but adjusted for multiple comparisons as per \code{method} (\code{P} = \code{P.adjust} if \code{method}=\code{"none"}). --OR--}
 \item{altP.adjust}{a vector of \emph{p}-values corresponding to \code{Z}, but adjusted for multiple comparisons as per \code{method} (\code{P} = \code{P.adjust} if \code{method}=\code{"none"}) when using the \code{altp=TRUE} option.}
 \item{comparisons}{a vector of strings labeling each pairwise comparison, as qualified by the \code{rmc} option, using either the variable values, or the factor labels or (or factor values if unlabeled). These labels match the corresponding position in the \code{Z}, \code{P}, and \code{P.adjust} vectors.}
 }
\author{
Alexis Dinno (\email{alexis.dinno@pdx.edu})

Please contact me with any questions, bug reports or suggestions for improvement.  Fixing bugs will be facilitated by sending along:
\tabular{rl}{
[1] \tab a copy of the data (de-labeled or anonymized is fine),\cr
[2] \tab a copy of the command syntax used, and\cr
[3] \tab a copy of the exact output of the command.\cr
}
}
\references{
Benjamini, Y. and Hochberg, Y.  (1995)  Controlling the false discovery rate: A practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society.  Series B (Methodological)}.  \bold{57}, 289--300.  <\doi{10.1111/j.2517-6161.1995.tb02031.x}>.

Benjamini, Y.  and Yekutieli, D.  (2001)  The control of the false discovery rate in multiple testing under dependency. \emph{Annals of Statistics}.  \bold{29}, 1165--1188.  <\doi{10.1214/aos/1013699998}>.

Dunn, O.  J.  (1961)  Multiple comparisons among means.  \emph{Journal of the American Statistical Association}.  \bold{56}, 52--64.  <\doi{10.1080/01621459.1961.10482090}>.

Dunn, O.  J.  (1964)  Multiple comparisons using rank sums.  \emph{Technometrics}.  \bold{6}, 241--252.  <\doi{10.1080/00401706.1964.10490181}>.

Hochberg, Y.  (1988)  A sharper Bonferroni procedure for multiple tests of significance.  \emph{Biometrika}.  \bold{75}, 800--802.  <\doi{10.1093/biomet/75.4.800}>.

Holm, S.  (1979)  \href{https://www.jstor.org/stable/4615733}{A simple sequentially rejective multiple test procedure}.  \emph{Scandinavian Journal of Statistics}.  \bold{6}, 65--70.

Kruskal, W.  H. and Wallis, A.  (1952)  Use of ranks in one-criterion variance analysis.  \emph{Journal of the American Statistical Association}.  \bold{47}, 583--621.  <\doi{10.1080/01621459.1952.10483441}>.

\enc{Šidák}{Sidak}, Z.  (1967)  Rectangular confidence regions for the means of multivariate normal distributions.  \emph{Journal of the American Statistical Association}.  \bold{62}, 626--633.  <\doi{10.1080/01621459.1967.10482935}>.
}
\examples{
## Example cribbed and modified from the kruskal.test documentation
## Hollander & Wolfe (1973), 116.
## Mucociliary efficiency from the rate of removal of dust in normal
##  subjects, subjects with obstructive airway disease, and subjects
##  with asbestosis.  
x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
dunn.test(x=list(x,y,z))

x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)),
            labels = c("Normal",
                       "COPD",
                       "Asbestosis"))
dunn.test(x, g)

## Example based on home care data from Dunn (1964)
data(homecare)
attach(homecare)
dunn.test(occupation, eligibility, method="hs", list=TRUE)

## Air quality data set illustrates differences in different
## multiple comparisons adjustments
attach(airquality)
dunn.test(Ozone, Month, kw=FALSE, method="bonferroni")
dunn.test(Ozone, Month, kw=FALSE, method="hs")
dunn.test(Ozone, Month, kw=FALSE, method="bh")
detach(airquality)
}
\keyword{htest}
\keyword{nonparametric}
\keyword{stats}
