\name{bal.tab.default}
\alias{bal.tab.default}
\title{
Balance Statistics for Other Objects
}
\description{
Generates balance statistics using an object for which there is not a defined method. 
}
\usage{
\method{bal.tab}{default}(x,
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
        ...)

}
\arguments{
  \item{x}{
An object containing information about conditioning. See Details.
}
  \item{stats, int, poly, distance, addl, data, continuous, binary, s.d.denom, thresholds, weights, cluster, imp, pairwise, s.weights, abs, subset, quick}{
see \fun{bal.tab} for details.
}
  \item{\dots}{
other arguments that would be passed to \fun{bal.tab.formula}, \fun{bal.tab.data.frame}, or \fun{bal.tab.time.list}. See Details.
}
}
\details{
\code{bal.tab.default()} processes its input and attempt to extract enough information from it to display covariate balance for \code{x}. The purpose of this method is to allow users who have created their own objects containing conditioning information (i.e., weights, subclasses, treatments, covariates, etc.) to access the capabilities of \code{bal.tab()} without having a special method written for them. By including the correct items in \code{x}, \code{bal.tab.default()} can present balance tables as if the input was the output of one of the specifically supported packages (e.g., \pkg{MatchIt}, \pkg{twang}, etc.).

The function will search \code{x} for the following named items and attempt to process them:
\describe{
\item{\code{treat}}{A vector (\code{numeric}, \code{character}, \code{factor}) containing the values of the treatment for each unit or the name of the column in \code{data} containing them. Essentially the same input to \code{treat} in \fun{bal.tab.data.frame}.
}
\item{\code{treat.list}}{A list of vectors (\code{numeric}, \code{character}, \code{factor}) containing, for each time point, the values of the treatment for each unit or the name of the column in \code{data} containing them. Essentially the same input to \code{treat.list} in \fun{bal.tab.time.list}.
}
\item{\code{covs}}{A \code{data.frame} containing the values of the covariates for each unit. Essentially the same input to \code{covs} in \fun{bal.tab.data.frame}.
}
\item{\code{covs.list}}{A list of \code{data.frame}s containing, for each time point, the values of the covariates for each unit. Essentially the same input to \code{covs.list} in \fun{bal.tab.time.list}.
}
\item{\code{formula}}{A \code{formula} with the treatment variable as the response and the covariates for which balance is to be assessed as the terms. Essentially the same input to \code{formula} in \fun{bal.tab.formula}.
}
\item{\code{formula.list}}{A list of \code{formula}s with, for each time point, the treatment variable as the response and the covariates for which balance is to be assessed as the terms. Essentially the same input to \code{formula.list} in \fun{bal.tab.time.list}.
}
\item{\code{data}}{A \code{data.frame} containing variables with the names used in other arguments and components (e.g., \code{formula}, \code{weights}, etc.). Essentially the same input to \code{data} in \fun{bal.tab.formula}, \fun{bal.tab.data.frame}, or \fun{bal.tab.time.list}.
}
\item{\code{weights}}{A vector, list, or \code{data.frame} containing weights for each unit or a string containing the names of the weights variables in \code{data}. Essentially the same input to \code{weights} in \fun{bal.tab.data.frame} or \fun{bal.tab.time.list}.
}
\item{\code{distance}}{
A vector, formula, or data frame containing distance values (e.g., propensity scores) or a character vector containing their names. If a formula or variable names are specified, \code{bal.tab()} will look in the argument to \code{data}, if specified. Essentially the same input to \code{distance} in \fun{bal.tab.data.frame}.
}
\item{\code{formula.list}}{A list of vectors or \code{data.frame}s containing, for each time point, distance values (e.g., propensity scores) for each unit or a string containing the name of the distance variable in \code{data}. Essentially the same input to \code{distance.list} in \fun{bal.tab.time.list}.
}
\item{\code{subclass}}{A vector containing subclass membership for each unit or a string containing the name of the subclass variable in \code{data}. Essentially the same input to \code{subclass} in \fun{bal.tab.data.frame}.
}
\item{\code{match.strata}}{A vector containing matching stratum membership for each unit or a string containing the name of the matching stratum variable in \code{data}. Essentially the same input to \code{match.strata} in \fun{bal.tab.data.frame}.
}
\item{\code{estimand}}{A \code{character} vector; whether the desired estimand is the "ATT", "ATC", or "ATE" for each set of weights. Essentially the same input to \code{estimand} in \fun{bal.tab.data.frame}.
}
\item{\code{s.weights}}{A vector containing sampling weights for each unit or a string containing the name of the sampling weight variable in \code{data}. Essentially the same input to \code{s.weights} in \fun{bal.tab.data.frame} or \fun{bal.tab.time.list}.
}
\item{\code{focal}}{The name of the focal treatment when multi-category treatments are used. Essentially the same input to \code{focal} in \fun{bal.tab.data.frame}.
}
\item{\code{call}}{A \code{call} object containing the function call, usually generated by using \fun{match.call} inside the function that created \code{x}.
}
}
Any of these items can also be supplied directly to \code{bal.tab.default}, e.g., \code{bal.tab.default(x, formula = treat ~ x1 + x2)}. If supplied, it will override the object with the same role in \code{x}. In addition, any arguments to \fun{bal.tab.formula}, \fun{bal.tab.data.frame}, and \fun{bal.tab.time.list} are allowed and perform the same function.

At least some inputs containing information to create the treatment and covariates are required (e.g., \code{formula} and \code{data} or \code{covs} and \code{treat}). All other arguments are optional and have the same defaults as those in \fun{bal.tab.data.frame} or \fun{bal.tab.time.list}. If \code{treat.list}, \code{covs.list}, or \code{formula.list} are supplied in \code{x} or as an argument to \code{bal.tab.default()}, the function will proceed considering a longitudinal treatment. Otherwise, it will proceed considering a point treatment.

\code{bal.tab.default()}, like other \code{bal.tab} methods, is just a shortcut to supply arguments to \code{bal.tab.data.frame()} or \code{bal.tab.time.list()}. Therefore, any matters regarding argument priority or function are described in the documentation for these methods.

}
\value{
For point treatments, if clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the specified treatment and covariates. See \fun{bal.tab} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.

If imputations are specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations, just as with clusters. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details.

If multi-category treatments are used, an object of class \code{"bal.tab.multi"} containing balance summaries for each pairwise treatment comparison and a summary of balance across pairwise comparisons. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.

If longitudinal treatments are used, an object of class \code{"bal.tab.msm"} containing balance summaries at each time point. Each balance summary is its own \code{bal.tab} object. See \code{\link[=class-bal.tab.msm]{bal.tab.msm}} for more details.
}


\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab.data.frame} and \fun{bal.tab.time.list} for additional arguments to be supplied.
\fun{bal.tab} for output and details of calculations.
\code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for more information on clustered data.
\code{\link[=class-bal.tab.imp]{bal.tab.imp}} for more information on multiply imputed data.
\code{\link[=class-bal.tab.multi]{bal.tab.multi}} for more information on multi-category treatments.
}
\examples{
data("lalonde", package = "cobalt")
covs <- subset(lalonde,  select = -c(treat, re78))

##Writing a function the produces output for direct
##use in bal.tab.default

ate.weights <- function(treat, covs) {
    data <- data.frame(treat, covs)
    formula <- formula(data)
    ps <- glm(formula, data = data, 
              family = "binomial")$fitted.values
    weights <- treat/ps + (1-treat)/(1-ps)
    call <- match.call()
    out <- list(treat = treat,
                covs = covs,
                distance = ps,
                weights = weights,
                estimand = "ATE",
                call = call)
    return(out)
}

out <- ate.weights(lalonde$treat, covs)

bal.tab(out, un = TRUE)
}
\keyword{tables}
