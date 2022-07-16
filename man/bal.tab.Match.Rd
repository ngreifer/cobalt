\name{bal.tab.Match}
\alias{bal.tab.Match}
\alias{bal.tab.optmatch}
\alias{bal.tab.ebalance}
\alias{bal.tab.designmatch}

\title{
Balance Statistics for \code{Matching}, \code{optmatch}, \code{ebal}, and \code{designmatch} Objects
}
\description{
Generates balance statistics for output objects from \pkg{Matching}, \pkg{optmatch}, \pkg{ebal}, and \pkg{designmatch}.
}
\usage{
\method{bal.tab}{Match}(x, 
        formula = NULL,
        data = NULL,
        treat = NULL,
        covs = NULL,
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
        s.weights = NULL,
        abs = FALSE,
        subset = NULL,
        quick = TRUE,
        ...)
}

\arguments{
  \item{x}{
either a \code{Match} object (the output of a call to \pkgfun{Matching}{Match} or \pkgfun{Matching}{Matchby}), an \code{optmatch} object (the output of a call to \pkgfun{optmatch}{pairmatch} or \pkgfun{optmatch}{fullmatch}), an \code{ebalance} object (the output of a call to \pkgfun{ebal}{ebalance} or \pkgfun{ebal}{ebalance.trim}), or the output of a call to \pkgfun{designmatch}{bmatch} or related wrapper functions from the \pkg{designmatch} package.
}
  \item{formula}{
a \code{formula} with the treatment variable as the response and the covariates for which balance is to be assessed as the predictors. All named variables must be in \code{data}. See Details.
}
  \item{data}{
a data frame containing variables named in \code{formula}, if supplied, and other arguments.
}
  \item{treat}{
a vector of treatment statuses. See Details.
}
  \item{covs}{
a data frame of covariate values for which to check balance. See Details.
}
  \item{stats, int, poly, distance, addl, continuous, binary, thresholds, weights, cluster, imp, pairwise, s.weights, abs, subset, quick, \dots}{
see \fun{bal.tab} for details.

See below for a special note on the \code{s.d.denom} argument.
}

The following argument has a special notes when used with these input objects:
  \item{s.d.denom}{
if not specified, for \code{Match} objects, \code{bal.tab()} will use "treated" if the estimand of the call to \code{Match()} is the ATT, "pooled" if the estimand is the ATE, and "control" if the estimand is the ATC; for \code{optmatch}, \code{ebal}, and \code{designmatch} objects, \code{bal.tab()} will determine which value makes the most sense based on the input. Abbreviations allowed.
}

}
\details{
\code{bal.tab()} generates a list of balance summaries for the object given, and function similarly to \pkgfun{Matching}{MatchBalance} and \pkgfun{designmatch}{meantab}. Note that output objects from \pkg{designmatch} do not have their own class; \code{bal.tab()} first checks whether the object meets the criteria to be treated as a \code{designmatch} object before dispatching the correct method. Renaming or removing items from the output object can create unintended consequences.

The input to \code{bal.tab.Match()}, \code{bal.tab.optmatch()}, \code{bal.tab.ebalance()}, and \code{bal.tab.designmatch()} must include either both \code{formula} and \code{data} or both \code{covs} and \code{treat}. Using the \code{formula} + \code{data} inputs mirrors how \pkgfun{Matching}{MatchBalance} is used, and using the \code{covs} + \code{treat} input mirrors how \pkgfun{designmatch}{meantab} is used. (Note that to see identical results to \code{meantab()}, \code{s.d.denom} must be set to \code{"pooled"}, though this is not recommended.) For \code{optmatch} output objects, specifying a treatment is not required.
}
\value{
For point treatments, if clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the given object. See \fun{bal.tab} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.

}

\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab} for details of calculations.
}
\examples{\dontshow{if (requireNamespace("Matching", quietly = TRUE)) force(\{ # examplesIf}
########## Matching ##########

library(Matching); data("lalonde", package = "cobalt")

p.score <- glm(treat ~ age + educ + race + 
            married + nodegree + re74 + re75, 
            data = lalonde, family = "binomial")$fitted.values
Match.out <- Match(Tr = lalonde$treat, X = p.score)

## Using formula and data
bal.tab(Match.out, formula = treat ~ age + educ + race + 
        married + nodegree + re74 + re75, data = lalonde)
\dontshow{\}) # examplesIf} \dontshow{if (requireNamespace("optmatch", quietly = TRUE)) force(\{ # examplesIf}
########## optmatch ##########

library("optmatch"); data("lalonde", package = "cobalt")

lalonde$prop.score <- glm(treat ~ age + educ + race + 
            married + nodegree + re74 + re75, 
            data = lalonde, family = binomial)$fitted.values
pm <- pairmatch(treat ~ prop.score, data = lalonde)

## Using formula and data
bal.tab(pm, formula = treat ~ age + educ + race + 
        married + nodegree + re74 + re75, data = lalonde,
        distance = "prop.score")
\dontshow{\}) # examplesIf} \dontshow{if (requireNamespace("ebal", quietly = TRUE)) force(\{ # examplesIf}
########## ebal ##########

library("ebal"); data("lalonde", package = "cobalt")

covariates <- subset(lalonde, select = -c(re78, treat, race))
e.out <- ebalance(lalonde$treat, covariates)

## Using treat and covs
bal.tab(e.out, treat = lalonde$treat, covs = covariates)
\dontshow{\}) # examplesIf} \dontshow{if (requireNamespace("designmatch", quietly = TRUE)) force(\{ # examplesIf}
########## designmatch ##########
\donttest{
library("designmatch"); data("lalonde", package = "cobalt")

covariates <- as.matrix(lalonde[c("age", "educ", "re74", "re75")])
treat <- lalonde$treat
dmout <- bmatch(treat,
                total_groups = sum(treat == 1),
                mom = list(covs = covariates,
                           tols = absstddif(covariates, 
                                            treat, .05))
                )

## Using treat and covs
bal.tab(dmout, treat = treat, covs = covariates)
}
\dontshow{\}) # examplesIf}}
\keyword{tables}
