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
         cluster = NULL, which.cluster = NULL)
}
\arguments{
  \item{obj}{
the object for which balance is to be assessed; can be a \code{matchit} object from MatchIt, a \code{ps} object from twang, a \code{Match} object from Matching, a \code{CBPS} or \code{CBPSContinuous} object from CBPS, or a \code{formula} or data frame containing covariate information.
}
  \item{var.name}{
\code{character}; the name of the variable whose values are to be plotted. To view distrbutions of the distance measure (e.g., propensity score), if any, use \code{".distance"} as the argument.
}
  \item{\dots}{
other arguments to define the variable, treatment, and weights. Some inputs are required depending on the method. See Additional Arguments.
}
  \item{un}{
  \code{logical}; whether to display distributional balance for the unadjusted sample or not. The default is to display balance for the adjusted sample.
 }
  \item{which.sub}{
if subclassification is used, a single number corresponding to the subclass for which the distributions are to be displayed.
}
  \item{cluster}{
optional; a vector of cluster membership.
}
  \item{which.cluster}{
if clusters are used, which cluster to display. \code{bal.plot()} can only display balance for one cluster at a time. Can be a cluster name or numerical index for which to display balance. Indices correspond to the alphabetical order of cluster names.
}
}
\section{Additional Arguments}{
\code{bal.plot()} works like \code{bal.tab()} in that it can take a variety of types of inputs and yield the same output for each. Depending on what kind of input is given, different additional paramaters are required in \code{\dots}. For details on what is required and allowed for each additional input and their defaults, see the help file for the \code{\link{bal.tab}} method associated with the input. The following are the required and optional additional arguments based on each input type:

For \code{matchit} objects: None

For \code{ps} objects: \code{full.stop.method} (optional; see defaults).

For \code{Match} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

For \code{CBPS} objects: None

For \code{formula}s: \code{data}, \code{weights}, \code{distance} (optional), \code{method} (optional; see defaults)

For data frames: \code{treat}, \code{data} (optional), \code{weights}, \code{distance} (optional), \code{method} (optional; see defaults)
}
\details{
\code{bal.plot()} uses \code{ggplot()} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. This mean that users can edit aspects of the plot using \pkg{ggplot2} syntax. For binary treatments with continuous covariates or continuous treatments with categorical or binary covariates, density plots are created using \pkg{ggplot2}'s \code{geom_density()} method; for binary treatments with binary covariates, bar graphs are created using \pkg{ggplot2}'s \code{geom_bar()} method; for continuous treatments with continuous covariates, scatterplots are created using \pkg{ggplot2}'s \code{geom_point()} method.

For continuous treatments with continuous covariates, three additional lines are presented for aid in balance assessment. The dashed blue line is the linear fit line. The solid blue line is a Loess fit curve generated with \pkg{ggplot2}'s default \code{geom_smooth()} method. The solid black line is a horizontal reference line intercepting the (weighted) treatment mean. Proximity of the fit lines to the reference line is indicative of independence between the covariate and treatment variable.
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
m.out <- matchit(treat ~ age + educ + black + hispan + 
                 married + nodegree + re74 + re75, 
                 data = lalonde)

bal.plot(m.out, "age")
bal.plot(m.out, "married")
}
\keyword{design}