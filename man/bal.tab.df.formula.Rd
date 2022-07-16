\name{bal.tab.df.formula}
\alias{bal.tab.data.frame}
\alias{bal.tab.formula}
\title{
Balance Statistics for Data Sets
}
\description{
Generates balance statistics for unadjusted, matched, weighted, or stratified data using either a \code{data.frame} or formula interface.
}
\usage{
\method{bal.tab}{data.frame}(x, 
        treat,
        stats,
        int = FALSE,
        poly = 1,
        distance = NULL,
        addl = NULL,
        data = NULL,
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
        subclass = NULL,
        match.strata = NULL,
        method,
        estimand = NULL,
        focal = NULL,
        ...)

\method{bal.tab}{formula}(x, 
        data = NULL,
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
        subclass = NULL,
        match.strata = NULL,
        method,
        estimand = NULL,
        focal = NULL,
        ...)
}
\arguments{
  \item{x}{
either a \code{data.frame} containing covariate values for each unit or a \code{formula} with the treatment variable as the response and the covariates for which balance is to be assessed as the terms. If a formula is supplied, all terms must be present as variable names in \code{data} or the global environment.
}
  \item{treat}{
either a vector containing treatment status values for each unit or a string containing the name of the treatment variable in \code{data}. Required for the \code{data.frame} method.
}
  \item{stats, int, poly, distance, addl, data, continuous, binary, thresholds, weights, cluster, imp, pairwise, s.weights, abs, subset, quick, \dots}{
see \fun{bal.tab} for details. See below for a special note on the \code{s.d.denom} argument.
}
  \item{subclass}{
optional; either a vector containing subclass membership for each unit or a string containing the name of the subclass variable in \code{data}.
}
\item{match.strata}{
optional; either a vector containing matching stratum membership for each unit or a string containing the name of the matching stratum variable in \code{data}. See Details.
}
  \item{method}{
\code{character}; the method of adjustment, if any. If \code{weights} are specified, the user can specify either "matching" or "weighting"; "weighting" is the default. If multiple sets of weights are used, each must have a corresponding value for \code{method}, but if they are all of the same type, only one value is required. If \code{subclass} is specified, "subclassification" is the default. Abbreviations allowed. The only distinction between "matching" and "weighting" is how sample sizes are displayed.
}
  \item{estimand}{
\code{character}; whether the desired estimand is the "ATT", "ATC", or "ATE" for each set of weights. This argument can be used in place of \code{s.d.denom} to specify how standardized differences are calculated.
}
  \item{focal}{
the name of the focal treatment when multi-category treatments are used. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.
}

The following argument has a special note when used with \code{data.frame} or \code{formula} input objects:
  \item{s.d.denom}{
if weights are supplied, each set of weights should have a corresponding entry to \code{s.d.denom}; a single entry will be recycled to all sets of weights. If left blank and one of \code{weights}, \code{subclass}, or \code{match.strata} are supplied, \code{bal.tab()} will figure out which one is best based on \code{estimand}, if given (for ATT, "treated"; for ATC, "control"; otherwise "pooled") and other clues if not.
}
}
\details{
\code{bal.tab.data.frame()} generates a list of balance summaries for the covariates and treatment status values given. \code{bal.tab.formula()} does the same but uses a formula interface instead.  When the formula interface is used, the formula and data are reshaped into a treatment vector and \code{data.frame} of covariates and then simply passed through the \code{data.frame} method.  

If \code{weights}, \code{subclass} and \code{match.strata} are all \code{NULL}, balance information will be presented only for the unadjusted sample.

The argument to \code{match.strata} corresponds to a factor vector containing the name or index of each pair/stratum for units conditioned through matching, for example, using the \pkg{optmatch} package. If more than one of \code{weights}, \code{subclass}, or \code{match.strata} are specified, \code{bal.tab()} will attempt to figure out which one to apply. Currently only one of these can be applied ta a time. \code{bal.tab()} behaves differently depending on whether subclasses are used in conditioning or not. If they are used, \code{bal.tab()} creates balance statistics for each subclass and for the sample in aggregate. See \code{\link[=class-bal.tab.subclass]{bal.tab.subclass}} for more information.

Multiple sets of weights can be supplied simultaneously by entering a \code{data.frame} or a character vector containing the names of weight variables found in \code{data} or a list of weights vectors or names. The arguments to \code{method}, \code{s.d.denom}, and \code{estimand}, if any, must be either the same length as the number of sets of weights or of length one, where the sole entry is applied to all sets. When standardized differences are computed for the unadjusted group, they are done using the first entry to \code{s.d.denom} or \code{estimand}. When only one set of weights is supplied, the output for the adjusted group will simply be called \code{"Adj"}, but otherwise will be named after each corresponding set of weights. Specifying multiple sets of weights will also add components to other outputs of \code{bal.tab()}.

}
\value{
For point treatments, if clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the specified treatment and covariates. See \fun{bal.tab} for details.

If imputations are specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details.

If multi-category treatments are used, an object of class \code{"bal.tab.multi"} containing balance summaries for each pairwise treatment comparison. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
}

\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab} for output and details of calculations.
\code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for more information on clustered data.
\code{\link[=class-bal.tab.imp]{bal.tab.imp}} for more information on multiply imputed data.
\code{\link[=class-bal.tab.multi]{bal.tab.multi}} for more information on multi-category treatments.
}
\examples{
data("lalonde", package = "cobalt")
lalonde$p.score <- glm(treat ~ age + educ + race, data = lalonde, 
                       family = "binomial")$fitted.values
covariates <- subset(lalonde, select = c(age, educ, race))
                     
## Propensity score weighting using IPTW
lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
                               1/lalonde$p.score, 
                               1/(1-lalonde$p.score))

# data frame interface:
bal.tab(covariates, treat = "treat", data = lalonde, 
        weights = "iptw.weights", s.d.denom = "pooled")

# Formula interface:
bal.tab(treat ~ age + educ + race, data = lalonde, 
        weights = "iptw.weights", s.d.denom = "pooled")
      
## Propensity score subclassification
lalonde$subclass <- findInterval(lalonde$p.score, 
                        quantile(lalonde$p.score, 
                        (0:6)/6), all.inside = TRUE)

# data frame interface:
bal.tab(covariates, treat = "treat", data = lalonde, 
        subclass = "subclass", disp.subclass = TRUE, 
        s.d.denom = "pooled")

# Formula interface:
bal.tab(treat ~ age + educ + race, data = lalonde, 
      subclass = "subclass", disp.subclass = TRUE, 
      s.d.denom = "pooled")
}
\keyword{tables}
