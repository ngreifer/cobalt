\name{bal.plot}
\alias{bal.plot}

\title{
Generate Univariate Balance Plots
}
\description{
Generates density plots (for continuous variables) and histograms (for binary variables) displaying distributional balance between the treated and control groups using \pkg{ggplot2}.
}
\usage{
bal.plot(obj, var.name, ...)
}
\arguments{
  \item{obj}{
the object for which balance is to be assessed; can be a \code{matchit} object from MatchIt, a \code{ps} object from twang, a \code{Match} object from Matching, or a \code{formula} or data frame containing covariate information.
}
  \item{var.name}{
\code{character}; the name of the variable whose values are to be plotted.
}
  \item{\dots}{
other arguments to define the variable, treatment, and weights. Some inputs are required depending on the method. See Additional Arguments.
}
}
\section{Additional Arguments}{
\code{bal.plot()} works like \code{bal.tab()} in that it can take a variety of types of inputs and yield the same output for each. Depending on what kind of input is given, different additional paramaters are required in \code{\dots}. For details on what is required and allowed for each additional input and their defaults, see the help file for the \code{\link{bal.tab}} method associated with the input. The following are the required an optiona ladditiona larguments based on each input type:

For \code{matchit} objects: None

For \code{ps} objects: \code{full.stop.method} (optional; see defaults).

For \code{Match} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

For \code{CBPS} objects: None

For \code{formula}s: \code{data}, \code{weights}, \code{distance} (optional), \code{method} (optional; see defaults)

For data frames: \code{treat}, \code{data} (optional), \code{weights}, \code{distance} (optional), \code{method} (optional; see defaults)

}
\details{
\code{bal.plot()} uses \code{ggplot()} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. This mean that users can edit aspects of the plot using \pkg{ggplot2} syntax. For continuous variables, density plots are created using \pkg{ggplot2}'s \code{geom_density()} method; for binary variables, bar graphs are created using \pkg{ggplot2}'s \code{geom_bar()} method.
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