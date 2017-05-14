\name{bal.plot}
\alias{bal.plot}

\title{
Generate Univariate Balance Plots
}
\description{
Generates density plots, bar graphs, or scatterplots displaying distributional balance between treatment and covariates using \pkg{ggplot2}.
}
\usage{
bal.plot(obj, var.name, ..., un = FALSE, which.sub = NULL, 
         cluster = NULL, which.cluster = NULL, imp = NULL, 
         which.imp = NULL)
}
\arguments{
  \item{obj}{
the object for which balance is to be assessed; can be a \code{matchit} object from \pkg{MatchIt}, a \code{ps} object from \pkg{twang}, a \code{Match} object from \pkg{Matching}, an \code{optmatch} object from \pkg{optmatch}, a \code{CBPS} or \code{CBPSContinuous} object from \pkg{CBPS}, an \code{ebalance} object from \pkg{ebal}, or a \code{formula} or data frame containing covariate information.
}
  \item{var.name}{
\code{character}; the name of the variable whose values are to be plotted. To view distrbutions of the distance measure (e.g., propensity score), if any, use \code{"distance"} as the argument unless the distance variable has been named. If there are duplicate variable names across inputs, \code{bal.plot()} will first look in the covariate data frame from \code{obj}, followed by \code{addl}, and then \code{distance}, if any.
}
  \item{\dots}{
other arguments to define the variable, treatment, and weights. Some inputs are required depending on the method. See Additional Arguments.
}
  \item{un}{
  \code{logical}; whether to display distributional balance for the unadjusted sample or not. The default is to display balance for the adjusted sample.
 }
  \item{which.sub}{
\code{numeric}; if subclassification is used, a vector corresponding to the subclass(es) for which the distributions are to be displayed. If \code{NULL} (the default), distributions from all subclasses are displayed in a grid.
}
  \item{cluster}{
optional; a vector of cluster membership, or the name of a variable in an available data set passed to \code{bal.plot()} that contains cluster membership.
}
  \item{which.cluster}{
if clusters are used, which cluster(s) to display. Can be cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. If \code{NULL} (the default), all clusters are displayed. If \code{NA}, clusters are ignored.
}
  \item{imp}{
optional; a vector of imputation indices, or the name of a variable in an available data set passed to \code{bal.plot()} that contains imputation indices.
}
  \item{which.imp}{
if imputations are used, which imputations(s) to display. Must be numerical indices for which to display balance. If \code{NULL} (the default), all imputations are displayed. If \code{NA}, imputations are ignored.
}
}
\section{Additional Arguments}{
\code{bal.plot()} works like \code{bal.tab()} in that it can take a variety of types of inputs and yield the same output for each. Depending on what kind of input is given, different additional paramaters are required in \code{\dots}. For details on what is required and allowed for each additional input and their defaults, see the help file for the \code{\link{bal.tab}} method associated with the input. The following are the required additional arguments based on each input type (optional arguments in parentheses):

For \code{matchit} objects: None

For \code{ps} objects: (\code{full.stop.method}; see \link[=bal.tab.ps]{defaults}).

For \code{Match} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

For \code{optmatch} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

For \code{CBPS} objects: None

For \code{ebalance} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

For \code{formula}s: \code{data}, \code{weights}, (\code{distance}), (\code{method}; see \link[=bal.tab.formula]{defaults})

For data frames: \code{treat}, (\code{data}), \code{weights}, (\code{distance}), (\code{method}; see \link[=bal.tab.data.frame]{defaults})
}
\details{
\code{bal.plot()} uses \code{ggplot()} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. For categorical treatments with continuous covariates or continuous treatments with categorical covariates, density plots are created using \pkg{ggplot2}'s \code{geom_density()} method; for categorical treatments with categorical covariates, bar graphs are created using \pkg{ggplot2}'s \code{geom_bar()} method; for continuous treatments with continuous covariates, scatterplots are created using \pkg{ggplot2}'s \code{geom_point()} method.

For continuous treatments with continuous covariates, three additional lines are presented for aid in balance assessment. The dashed blue line is the linear fit line. The solid blue line is a Loess fit curve generated with \pkg{ggplot2}'s \code{geom_smooth(method = "loess")}. The solid black line is a horizontal reference line intercepting the (weighted) treatment mean. Proximity of the fit lines to the reference line is indicative of independence between the covariate and treatment variable.

When multiple plots are to be displayed (i.e., when requesting subclass balance, cluster balance, or imputation balance), the plots will be displayed in a grid using \pkg{ggplot2}'s \code{facet_grid()}. Subclassification cannot be used with clusters or multiply imputed data. If both clusters and multiple imputation are used (and the arguments to \code{which.cluster} and \code{which.imp} are not both \code{NA}), at least one of \code{which.cluster} and \code{which.imp} must correspond to a single (non-\code{NULL}) value. That is, plots can be displayed for one cluster and multiple imputations or multiple clusters and only one imputation.
}
\value{
A \code{"ggplot"} object, returned invisbly.
}

\author{
Noah Greifer \email{noah@unc.edu}
}


\seealso{
\code{\link{bal.tab}}
}
\examples{
library(MatchIt); data("lalonde", package = "cobalt")

#Nearest Neighbor Matching
m.out <- matchit(treat ~ age + educ + race + 
                 married + nodegree + re74 + re75, 
                 data = lalonde)

bal.plot(m.out, "age")
bal.plot(m.out, "race")
}
\keyword{design}