\name{bal.tab.cem.match}
\alias{bal.tab.cem.match}
\title{
Balance Statistics for \code{cem} Objects
}
\description{
Generates balance statistics for \code{cem.match} objects from \pkg{cem}.
}
\usage{
\method{bal.tab}{cem.match}(x, 
        data, 
        stats,
        int = FALSE, 
        poly = 1, 
        distance = NULL, 
        addl = NULL, 
        continuous,  
        binary, 
        s.d.denom, 
        thresholds = NULL,
        weights = NULL,
        cluster = NULL, 
        imp = NULL,
        pairwise = TRUE,
        abs = FALSE,
        subset = NULL,
        quick = TRUE, 
        ...)
}

\arguments{
  \item{x}{
a \code{cem.match} or \code{cem.match.list} object; the output of a call to \code{cem::cem()}.
}
  \item{data}{
a data frame containing the treatment, covariates, and variables that might be named in arguments to \code{distance}, \code{addl}, \code{cluster}, and \code{imp}. Can also be \code{mids} object, the output of a call to \pkgfun{mice}{mice} from the \pkg{mice} package, containing multiply imputed data sets. In this case, \code{imp} is automatically supplied using the imputation variable created from processing the \code{mids} object. An argument to \code{data} is required. It must be the same data used in the call to \code{cem()} or a \code{mids} object from which the data supplied to \code{datalist} in the {cem()} call originated.
}
  \item{stats}{
\code{character}; which statistic(s) should be reported. See \code{\link[=balance.stats]{stats}} for allowable options. For binary and multi-category treatments, "mean.diffs" (i.e., mean differences) is the default. Multiple options are allowed.
}
  \item{int}{
\code{logical} or \code{numeric}; whether or not to include 2-way interactions of covariates included in \code{covs} and in \code{addl}. If \code{numeric}, will be passed to \code{poly} as well.
}
  \item{poly}{
\code{numeric}; the highest polynomial of each continuous covariate to display. For example, if 2, squares of each continuous covariate will be displayed (in addition to the covariate itself); if 3, squares and cubes of each continuous covariate will be displayed, etc. If 1, the default, only the base covariate will be displayed. If \code{int} is numeric, \code{poly} will take on the value of \code{int}.
}
  \item{distance}{
an optional formula or data frame containing distance values (e.g., propensity scores) or a character vector containing their names. If a formula or variable names are specified, \code{bal.tab()} will look in the argument to \code{data}, if specified.
}
  \item{addl}{
an optional formula or data frame containing additional covariates for which to present balance or a character vector containing their names. If a formula or variable names are specified, \code{bal.tab()} will look in the arguments to \code{covs} and \code{data}, if specified.
}
  \item{continuous}{
whether mean differences for continuous variables should be standardized (\code{"std"}) or raw (\code{"raw"}). Default \code{"std"}. Abbreviations allowed. This option can be set globally using \fun{set.cobalt.options}. This argument is used to set \code{std} in \fun{col_w_smd}.
}
  \item{binary}{
whether mean differences for binary variables (i.e., difference in proportion) should be standardized (\code{"std"}) or raw (\code{"raw"}). Default \code{"raw"}. Abbreviations allowed. This option can be set globally using \fun{set.cobalt.options}. This argument is used to set \code{std} in \fun{col_w_smd}.
}
  \item{s.d.denom}{
\code{character}; how the denominator for standardized mean differences should be calculated, if requested. See \fun{col_w_smd} for allowable options. The default is "treated", where the treated group corresponds to the \code{baseline.group} in the call to \code{cem()}. Abbreviations allowed.
}
  \item{thresholds}{
a named vector of balance thresholds, where the name corresponds to the statistic (i.e., in \code{stats}) that the threshold applies to. For example, to request thresholds on mean differences and variance ratios, one can set \code{thresholds = c(m = .05, v = 2)}. Requesting a threshold automatically requests the display of that statistic. See Details.
}
  \item{weights}{
a named list containing additional weights on which to assess balance. Each entry can be a vector of weights, the name of a variable in \code{data} that contains weights, or an object with a \fun{get.w} method.
}
  \item{cluster}{
either a vector containing cluster membership for each unit or a string containing the name of the cluster membership variable in \code{data}. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
}
  \item{imp}{
either a vector containing imputation indices for each unit or a string containing the name of the imputation index variable in \code{data}. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details. Not necessary if \code{data} is a \code{mids} object.
}
  \item{pairwise}{
whether balance should be computed for pairs of treatments or for each treatment against all groups combined. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details. This can also be used with a binary treatment to assess balance with respect to the full sample.
}
  \item{abs}{
\code{logical}; whether displayed balance statistics should be in absolute value or not. 
}
  \item{subset}{
a \code{logical} or \code{numeric} vector denoting whether each observation should be included or which observations should be included. If \code{logical}, it should be the same length as the variables in the original call to \code{cem()}. \code{NA}s will be treated as \code{FALSE}. This can be used as an alternative to \code{cluster} to examine balance on subsets of the data.
}
  \item{quick}{
\code{logical}; if \code{TRUE}, will not compute any values that will not be displayed. Set to \code{FALSE} if computed values not displayed will be used later.
}
  \item{...}{
further arguments to control display of output. See \link[=options-display]{display options} for details.
}
}
\details{
\code{bal.tab.cem.match()} generates a list of balance summaries for the \code{cem.match} object given, and functions similarly to \code{cem::imbalance()}.

The \code{threshold} argument controls whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure.

}
\value{
If clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{cem.match} object. See \fun{bal.tab} for details.

If imputations are specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details.

If \code{cem()} is used with multi-category treatments, an object of class \code{"bal.tab.multi"} containing balance summaries for each pairwise treatment comparison. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
}

\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab} for details of calculations.
}
\examples{\dontshow{if (requireNamespace("cem", quietly = TRUE)) \{}
library(cem); data("lalonde", package = "cobalt")

## Coarsened exact matching
cem.out <- cem("treat", data = lalonde, drop = "re78")

bal.tab(cem.out, data = lalonde, un = TRUE, 
        stats = c("m", "k"))
\dontshow{\}}}
\keyword{tables}
