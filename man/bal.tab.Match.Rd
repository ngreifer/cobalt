\name{bal.tab.Match}
\alias{bal.tab.Match}
\alias{bal.tab.optmatch}
\alias{bal.tab.ebalance}
\alias{bal.tab.designmatch}

\title{
Balance Statistics for \code{Matching}, \code{optmatch}, \code{ebal}, and \code{designmatch} Objects
}
\description{
Generates balance statistics for output objects from \pkg{Matching}, \pkg{optmatch}, \pkg{ebal}, and \pkg{designmatch}. Note that several arguments that used to be documented here are now documented in \link[=options-display]{display options}. They are still available.
}
\usage{
\method{bal.tab}{Match}(M, 
    formula = NULL, 
    data = NULL, 
    treat = NULL, 
    covs = NULL, 
    int = FALSE, 
    poly = 1, 
    distance = NULL, 
    addl = NULL, 
    continuous,  
    binary,
    s.d.denom, 
    m.threshold = NULL, 
    v.threshold = NULL, 
    ks.threshold = NULL, 
    cluster = NULL, 
    abs = FALSE,
    subset = NULL,
    quick = TRUE, 
    ...)
    
\method{bal.tab}{optmatch}(optmatch, ...)

\method{bal.tab}{ebalance}(ebal, ...)

\method{bal.tab}{designmatch}(dm, ...)
}

\arguments{
  \item{M}{
a \code{Match} object; the output of a call to \code{Match()} from the \pkg{Matching} package.
}
  \item{optmatch}{
an \code{optmatch} object; the output of a call to \code{pairmatch()} or \code{fullmatch()} from the \pkg{optmatch} package. This should be a factor vector containing the matching stratum membership for each unit.
}
  \item{ebal}{
an \code{ebalance} object; the output of a call to \code{ebalance()} or \code{ebalance.trim()} from the \pkg{ebal} package.
}
  \item{dm}{
the output of a call to \code{bmatch()}, \code{nmatch()}, or related wrapper functions from the \pkg{designmatch} package. This should be a list containing the IDs of the matched treated and control units.
}
  \item{formula}{
a \code{formula} with the treatment variable as the response and the covariates for which balance is to be assessed as the predictors. All named variables must be in \code{data}. See Details.
}
  \item{data}{
a data frame containing all the variables named in \code{formula}. See Details.
}
  \item{treat}{
a vector of treatment statuses. See Details.
}
  \item{covs}{
a data frame of covariate values for which to check balance. See Details.
}
  \item{int}{
\code{logical} or \code{numeric}; whether or not to include 2-way interactions of covariates included in \code{covs} and in \code{addl}. If \code{numeric}, will be passed to \code{poly} as well. In older versions of \pkg{cobalt}, setting \code{int = TRUE} displayed squares of covariates; to replicate this behavior, set \code{int = 2}.
}
  \item{poly}{
\code{numeric}; the highest polynomial of each continuous covariate to display. For example, if 2, squares of each continuous covariate will be displayed (in addition to the covariate itself); if 3, squares and cubes of each continuous covariate will be displayed, etc. If 1, the default, only the base covariate will be displayed. If \code{int} is numeric, \code{poly} will take on the value of \code{int}.
}
  \item{distance}{
optional; either a vector or data.frame containing distance values (e.g., propensity scores) for each unit or a string containing the name of the distance variable in \code{data}.
}
  \item{addl}{
a data frame of additional covariates for which to present balance. These may be covariates included in the original dataset but not included in \code{formula} or \code{covs}. In general, it makes more sense to include all desired variables in \code{formula} or \code{covs} than in \code{addl}. See note in Details for using \code{addl}.
}
  \item{continuous}{
whether mean differences for continuous variables should be standardized ("std") or raw ("raw"). Default "std". Abbreviations allowed. This option can be set globally using \code{\link{set.cobalt.options}}.
}
  \item{binary}{
whether mean differences for binary variables (i.e., difference in proportion) should be standardized ("std") or raw ("raw"). Default "raw". Abbreviations allowed. This option can be set globally using \code{\link{set.cobalt.options}}.
}
  \item{s.d.denom}{
whether the denominator for standardized differences (if any are calculated) should be the standard deviation of the treated group ("treated"), the standard deviation of the control group ("control"), or the pooled standard deviation ("pooled"), computed as the square root of the mean of the group variances. Abbreviations allowed. If not specified, for \code{Match} objects, \code{bal.tab()} will use "treated" if the estimand of the call to \code{Match()} is the ATT, "pooled" if the estimand is the ATE, and "control" if the estimand is the ATC; for \code{optmatch}, \code{ebal}, and \code{designmatch} objects, \code{bal.tab()} will use "treated".
}
  \item{m.threshold}{
a numeric value for the threshold for mean differences. .1 is recommended. 
}
  \item{v.threshold}{
a numeric value for the threshold for variance ratios. Will automatically convert to the inverse if less than 1.
}
  \item{ks.threshold}{
a numeric value for the threshold for Kolmogorov-Smirnov statistics. Must be between 0 and 1. 
}
  \item{cluster}{
either a vector containing cluster membership for each unit or a string containing the name of the cluster membership variable in \code{data} or the CBPS object. See \code{\link{bal.tab.cluster}} for details.
}
  \item{abs}{
\code{logical}; whether displayed balance statistics should be in absolute value or not. 
}
  \item{subset}{
a \code{logical} vector denoting whether each observation should be included. It should be the same length as the variables in the original call to the conditioning function. \code{NA}s will be treated as \code{FALSE}. This can be used as an alternative to \code{cluster} to examine balance on subsets of the data.
}
  \item{quick}{
\code{logical}; if \code{TRUE}, will not compute any values that will not be displayed. Set to \code{FALSE} if computed values not displayed will be used later.
}
  \item{...}{
for \code{bal.tab.optmatch()}, \code{bal.tab.ebalance()}, and \code{bal.tab.designmatch()}, the same arguments as those passed to \code{bal.tab.Match()}. Otherwise, further arguments to control display of output. See \link[=options-display]{display options} for details.
}
}
\details{
\code{bal.tab()} generates a list of balance summaries for the object given, and function similarly to \code{MatchBalance()} in \pkg{Matching} and \code{meantab()} in \pkg{designmatch}. Note that output objects from \pkg{designmatch} do not have their own class; \code{bal.tab()} first checks whether the object meets the criteria to be treated as a \code{designmatch} object before dispatching the correct method. In particular, renaming or removing items from the output object can create unintended consequences.

The input to \code{bal.tab.Match()}, \code{bal.tab.optmatch()}, \code{bal.tab.ebalance()}, and \code{bal.tab.designmatch()} must include either both \code{formula} and \code{data} or both \code{covs} and \code{treat}. Using the \code{formula} + \code{data} inputs mirrors how \code{MatchBalance()} is used in \pkg{Matching}, and using the \code{covs} + \code{treat} input mirrors how \code{meantab()} is used in \pkg{designmatch}. (Note that to see identical results to \code{meantb()}, \code{s.d.denom} must be set to \code{"pooled"}, though this is not recommended.)

All balance statistics are calculated whether they are displayed by print or not, unless \code{quick = TRUE}. The threshold values (\code{m.threshold}, \code{v.threshold}, and \code{ks.threshold}) control whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure.

The inputs (if any) to \code{covs} must be a data frame; if more than one variable is included, this is straightforward (i.e., because \code{data[,c("v1", "v2")]} is already a data frame), but if only one variable is used (e.g., \code{data[,"v1"]}), R will coerce it to a vector, thus making it unfit for input. To avoid this, simply wrap the input to \code{covs} in \code{data.frame()} or use \code{subset()} if only one variable is to be added. Again, when more than one variable is included, the input is general already a data frame and nothing needs to be done.
}
\value{
For point treatments, if clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the given object. See \code{\link{bal.tab}} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link{bal.tab.cluster}} for details.

}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}} for details of calculations.
}
\examples{
########## Matching ##########
library(Matching); data("lalonde", package = "cobalt")

p.score <- glm(treat ~ age + educ + race + 
            married + nodegree + re74 + re75, 
            data = lalonde, family = "binomial")$fitted.values
Match.out <- Match(Tr = lalonde$treat, X = p.score)

## Using formula and data
bal.tab(Match.out, treat ~ age + educ + race + 
        married + nodegree + re74 + re75, data = lalonde)

########## optmatch ##########
library("optmatch"); data("lalonde", package = "cobalt")

lalonde$prop.score <- glm(treat ~ age + educ + race + 
            married + nodegree + re74 + re75, 
            data = lalonde, family = binomial)$fitted.values
pm <- pairmatch(treat ~ prop.score, data = lalonde)

## Using formula and data
bal.tab(pm, treat ~ age + educ + race + 
        married + nodegree + re74 + re75, data = lalonde,
        distance = "prop.score")

########## ebal ##########
library("ebal"); data("lalonde", package = "cobalt")

covariates <- subset(lalonde, select = -c(re78, treat, race))
e.out <- ebalance(lalonde$treat, covariates)

## Using treat and covs
bal.tab(e.out, treat = lalonde$treat, covs = covariates)

\dontrun{
########## designmatch ##########
library("designmatch"); data("lalonde", package = "cobalt")

covariates <- as.matrix(lalonde[c("age", "educ", "re74", "re75")])
dmout <- bmatch(lalonde$treat,
                total_groups = sum(lalonde$treat == 1),
                mom = list(covs = covariates,
                           tols = absstddif(covs, treat, .05))
                )

## Using treat and covs
bal.tab(dmout, treat = lalonde$treat, covs = covariates)
}
}