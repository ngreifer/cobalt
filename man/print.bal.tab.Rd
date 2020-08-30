\name{print.bal.tab}
\alias{print.bal.tab}
\alias{print.bal.tab.cluster}
\alias{print.bal.tab.subclass}
\alias{print.bal.tab.imp}
\alias{print.bal.tab.multi}
\alias{print.bal.tab.msm}

\title{
Print Results of a Call to \code{bal.tab()}
}
\description{
Prints \code{bal.tab()} output in a clean way. Provides options for printing.}
\usage{
\method{print}{bal.tab}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)

\method{print}{bal.tab.subclass}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      disp.subclass = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)

\method{print}{bal.tab.cluster}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      which.cluster, 
      cluster.summary = "as.is", 
      cluster.fun = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)

\method{print}{bal.tab.imp}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      which.imp, 
      imp.summary = "as.is", 
      imp.fun = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)
    
\method{print}{bal.tab.multi}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      which.treat, 
      multi.summary = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)
    
\method{print}{bal.tab.msm}(x, 
      imbalanced.only = "as.is", 
      un = "as.is", 
      disp.bal.tab = "as.is", 
      stats = "as.is", 
      disp.thresholds = "as.is", 
      disp = "as.is", 
      which.time, 
      msm.summary = "as.is", 
      digits = max(3, getOption("digits") - 3), 
      ...)
}
\arguments{
  \item{x}{
a \code{bal.tab} object; the output of a call to \code{\link{bal.tab}}.
}
  \item{imbalanced.only}{
whether to display only the covariates that failed to meet at least one of balance thresholds. Depends only on whether threshold were initial set in the call to \code{bal.tab()} and not on any arguments to \code{print()} (except \code{disp.bal.tab}).
}
    \item{un}{
whether to display balance values for the unadjusted sample. Ignored (and set to \code{TRUE}) if no conditioning was performed.
}
  \item{disp.bal.tab}{
whether to display the table of balance statistics. If \code{FALSE}, only other values (e.g., the call, sample sizes, balance tallies, and maximum imbalances) will be presented.
}
  \item{stats}{
\code{character}; which statistic(s) should be reported. For binary or multi-category treatments, the options are "mean.diffs" for mean differences (standardized or not according the selected \code{bal.tab} options), "variance.ratios" for variance ratios, and "ks.statistics" for Kolmogorov-Smirnov statistics. "mean.diffs" is the default. For continuous treatments, the only option is "correlations" for treatment-covariate correlations. Multiple options are allowed. Abbreviations allowed. Statistics that weren't requested in the original call to \code{bal.tab()} cannot be requested with \code{print} unless \code{quick = FALSE} in the original call.
}
  \item{disp.thresholds}{
\code{logical}; whether to display thresholds for each statistic for which thresholds were originally requested in the call to \code{bal.tab}. Should be a named logical vector with names corresponding to the thresholds. For example, if thresholds for mean differences were requested in \code{bal.tab}, set \code{disp.thresholds = c(m = FALSE)} to prevent them from being printed. If a statistic was prevented from being displayed by another argument to \code{print}, the thresholds will not be displayed.
}
  \item{disp}{
\code{character}; which distribution summary statistics to display. Allowable options include "means" and "sds". Statistics that weren't requested in the original call to \code{bal.tab()} cannot be requested with \code{print} unless \code{quick = FALSE} in the original call.
}
  \item{disp.subclass}{
whether to display balance information for individual subclasses if subclassification is used in conditioning.
}
  \item{which.cluster}{
which cluster(s) to display. If \code{NULL}, all clusters will be displayed. If \code{NA}, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. To display the clusters requested in the original call to \code{bal.tab()}, omit this argument, as specifying \code{"as.is"} will request a cluster called "as.is.".
}
  \item{cluster.summary}{
whether to display the cluster summary table. If \code{which.cluster} is \code{NA}, \code{cluster.summary} will be set to \code{TRUE}.
}
  \item{cluster.fun}{
a character vector of functions of balance statistics to display when displaying balance across clusters. Can be "mean", "min", or "max". More than one are allowed.
}
  \item{which.imp}{
which imputation(s) to display. If \code{NULL}, all imputations will be displayed. If \code{NA}, no imputations will be displayed. Otherwise, can be a vector of imputations numbers for which to display balance. To display the imputations requested in the original call to \code{bal.tab()}, omit this argument, or enter \code{"as.is"}.
}
  \item{imp.summary}{
whether to display the imputation summary table. If \code{which.imp} is \code{NA}, \code{cluster.summary} will be set to \code{TRUE}.
}
  \item{imp.fun}{
a character vector of functions of balance statistics to display when displaying balance across imputations. Can be "mean", "min", or "max". More than one are allowed.
}
  \item{which.treat}{
which treatments to display when multi-category treatments are used. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.
}
 \item{multi.summary}{
\code{logical}; whether to display the balance summary table across pairwise comparisons when multi-category treatments are used. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.
}
  \item{which.time}{
which time periods to display if longitudinal treatments are used. See \code{\link[=class-bal.tab.msm]{bal.tab.msm}} for details.
}
  \item{msm.summary}{
\code{logical}; whether to display the balance summary table across time periods when longitudinal treatments are used. See \code{\link[=class-bal.tab.msm]{bal.tab.msm}} for details.
}
  \item{digits}{
the number of digits to display.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
Simply calling \code{bal.tab()} will print its results, but it can be useful to store the results into an object and print them again later, possibly with different print options specified. The \code{print()} function automatically dispatches the correct method for the \code{bal.tab} object given. For balance tables generated from using weighting, matching, or no adjustment with binary or continuous treatments, \code{print.bal.tab()} will be used. For balance tables generated from using weighting, matching, or no adjustment with clusters, \code{print.bal.tab.cluster()} will be used. For balance tables generated from using subclassification, \code{print.bal.tab.subclass()} will be used. For balance tables generated with multi-category treatments, \code{print.bal.tab.multi()} will be used. For balance tables generated with sequential treatments, \code{print.bal.tab.msm()} will be used. For balance tables generated with multiply imputed data, \code{print.bal.tab.imp()} will be used.

For all parameters except \code{which.cluster}, \code{which.imp}, \code{which.treat}, and \code{which.time}, either omitting the argument or setting it to \code{"as.is"} will use the corresponding print option stored in the \code{bal.tab} object, which results from the original call to \code{bal.tab()}. For these other arguments, to retain the display option of the original call to \code{bal.tab()}, the argument must omitted, as using \code{"as.is"} would cause \code{print()} to attempt to display balance for, e.g., a cluster called \code{"as.is"}. If such a cluster existed and it was desired, it would otherwise be impossible to display it.

Any parameter used in \code{bal.tab()} for calculations, such as \code{int}, \code{addl}, or \code{distance}, cannot be used with \code{print()}; only those parameters listed above, those that solely determine printing options, can be used. To change computation options, a new call to \code{bal.tab()} must be performed.

Prior versions of \code{print()} used different arguments, such as \code{m.threshold} to display thresholds and \code{disp.means} to display group means, and these arguments still work but are not documented here. See \code{\link{display_options}}.
}
\note{
Unless \code{quick = FALSE} in the original call to \code{bal.tab()}, some values may not be calculated, in which case using \code{print()} will not display these values even when requested. For example, if \code{stats = "m"} and \code{quick = TRUE} in the original call to \code{bal.tab()} (the default for both), setting \code{stats = "ks"} in \code{print()} will not print the KS statistics because they were not calculated.
}
\author{
Noah Greifer 
}

\seealso{
\code{\link[base]{print}}, \code{\link{bal.tab}}
\link[=options-display]{display options} for further information on some of these options.
}
\examples{
data(lalonde)
library(WeightIt)
w.out <- weightit(treat ~ age + educ + married + race + re74 + re75, 
                  data = lalonde)
b <- bal.tab(w.out, stats = c("m", "v", "ks"), 
             un = TRUE, v.threshold = 2)
print(b, un = FALSE, stats = c("m", "v"),
      disp.thresholds = c(v = FALSE))
}
