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
dunn.test(x, g=NA, method=p.adjustment.methods, kw=TRUE, label=TRUE, 
      wrap=FALSE, alpha=0.05)

p.adjustment.methods
# c("none", "bonferroni", "sidak", "holm", "hs", "bh", "by")
}
\arguments{
 \item{x}{a numeric vector, or a list of numeric vectors. If the former, then groups \emph{must} be specified using \code{g}.}
 \item{g}{a factor variable, numeric vector, or character vector indicating group.}
 \item{method}{adjusts the \emph{p}-value for multiple comparisons using the Bonferroni, \enc{Šidák}{Sidak}, Holm, Holm-\enc{Šidák}{Sidak}, Benjamini-Hochberg, or Benjamini-YekuTeili adjustment (see Details). The default is no adjustment for multiple comparisons.}
 \item{kw}{if \code{TRUE} then the results of the Kruskal-Wallis test are reported.}
 \item{label}{if \code{TRUE} then the factor labels are used in the output table.}
 \item{wrap}{does not break up tables to maintain nicely formatted output. If \code{FALSE} then output of large tables is broken up across multiple pages.}
 \item{alpha}{the nominal level of significance used in the step-up/step-down multiple comparisons procedures (Holm, Holm-\enc{Šidák}{Sidak}, and Benjamini-Hochberg).}
}
\details{\code{dunn.test} computes Dunn's test (1964) for stochastic dominance and reports the results among multiple pairwise comparisons after a Kruskal-Wallis test for stochastic dominance among \emph{k} groups (Kruskal and Wallis, 1952). The interpretation of stochastic dominance requires an assumption that the \acronym{CDF} of one group does not cross the \acronym{CDF} of the other. \code{dunn.test }makes \emph{k}(\emph{k}-1)/2 multiple pairwise comparisons based on Dunn's \emph{z}-test-statistic approximations to the actual rank statistics. The null hypothesis for each pairwise comparison is that the probability of observing a randomly selected value from the first group that is larger than a randomly selected value from the second group equals one half; this null hypothesis corresponds to that of the Wilcoxon-Mann-Whitney rank-sum test. Like the rank-sum test, if the data can be assumed to be continuous, and the distributions are assumed identical except for a difference in location, Dunn's test may be understood as a test for median difference. \code{dunn.test} accounts for tied ranks.

The choices of multiple comparison adjustment available for \code{p.adjustment.methods} are \code{"none"} for no adjustment, several familywise error rate methods include \code{"bonferroni"} for Dunn's (1961) so-called Bonferroni adjustment, \code{"sidak"} for the Šidák (1967) adjustment for dependent tests, and \code{"holm"} and \code{"hs"} for the step-up Holm (1979) and Holm-\enc{Šidák}{Sidak} (Holm, 1979) adjustments, respectively, which progressively relax the adjustments. Finally, the \code{"bh"} option for the step-down Benjamini-Hochberg (1995) adjustment controls the false discovery rate for independent or strictly positively dependent tests, while the \code{"by"} option for the step-down Benjamini-Yekuteili (2001) adjustment controls the false discovery rate for all dependent tests.

Because the sequential step-up/step-down tests rejection decisions depend on both the \emph{p}-values and their ordering, those tests rejected using \code{"holm"}, \code{"hs"}, \code{"bh"}, or \code{"by"} at the \code{alpha} level (two-sided test) have an asterisk placed immediately after the adjusted \emph{p}-value in the output.
}
\value{
\code{dunn.test} returns:
 \item{Z}{a vector of all \emph{m} of Dunn's \emph{z} test statistics}
 \item{P}{a vector of \emph{p}-values corresponding to \code{Z}.}
 \item{P.adjust}{a vector of \emph{p}-values corresponding to \code{Z} but adjusted for multiple comparisons as per \code{method}.}
}
\author{
Alexis Dinno (\email{alexis.dinno@pdx.edu})
}
\references{
Benjamini, Y. and Hochberg, Y. 1995. \dQuote{Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.} \emph{Journal of the Royal Statistical Society. Series B (Methodological)}. 57: 289--300.

Benjamini, Y. and Yekutieli, D. 2001. \dQuote{The control of the false discovery rate in multiple testing under dependency.} \emph{Annals of Statistics}, 29(4):1165--1188.

Dunn, O. J. 1961. \dQuote{Multiple comparisons among means.} \emph{Journal of the American Statistical Association}. 56: 52--64.

Dunn, O. J. 1964. \dQuote{Multiple comparisons using rank sums}. \emph{Technometrics}. 6: 241--252.

Holm, S. 1979. \dQuote{A simple sequentially rejective multiple test procedure}. \emph{Scandinavian Journal of Statistics}. 6: 65--70.

Kruskal, W. H. and Wallis, A. 1952. \dQuote{Use of ranks in one-criterion variance analysis}. \emph{Journal of the American Statistical Association}. 47: 583--621.

\enc{Šidák}{Sidak}, Z. 1967. \dQuote{Rectangular confidence regions for the means of multivariate normal distributions}. \emph{Journal of the American Statistical Association}. 62: 626--633.
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
dunn.test(occupation,eligibility, method="hs")

}
\keyword{htest}
\keyword{nonparametric}
\keyword{stats}
