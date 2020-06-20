\name{bal.plot}
\alias{bal.plot}

\title{
Generate Univariate Balance Plots
}
\description{
Generates density plots, bar graphs, or scatterplots displaying distributional balance between treatment and covariates using \pkg{ggplot2}.
}
\usage{
bal.plot(x, 
    var.name, 
    ..., 
    which, 
    which.sub = NULL, 
    cluster = NULL,
    which.cluster = NULL, 
    imp = NULL,
    which.imp = NULL, 
    which.treat = NULL,
    which.time = NULL,
    mirror = FALSE, 
    type = "density", 
    colors = NULL,
    grid = FALSE,
    sample.names,
    position = "right", 
    facet.formula = NULL, 
    disp.means = getOption("cobalt_disp.means", FALSE),
    alpha.weight = TRUE)
}
\arguments{
  \item{x}{
the object for which balance is to be assessed; can be any object for which there is support in \code{bal.tab}.
}
  \item{var.name}{
\code{character}; the name of the variable whose values are to be plotted. To view distributions of the distance measure (e.g., propensity score), if any, use \code{"distance"} as the argument unless the distance variable has been named. If there are duplicate variable names across inputs, \code{bal.plot()} will first look in the covariate \code{data.frame} from \code{obj}, followed by \code{addl}, and then \code{distance}, if any. If not specified, will use the first covariate available with a warning.
}
  \item{\dots}{
other arguments to define the variable, treatment, and weights. Some inputs are required depending on the method. See Additional Arguments. Can also be used to supply the \code{bw}, \code{adjust}, \code{kernel}, and \code{n} arguments for \code{\link[ggplot2]{geom_density}} and the \code{bins} argument for \code{\link[ggplot2]{geom_histogram}}.
}
  \item{which}{
whether to display distributional balance for the adjusted (\code{"adjusted"}) or unadjusted sample (\code{"unadjusted"}) or both at the same time (\code{"both"}). When multiple weights are present, the names of the weights can be supplied, too. The default is to display balance for the adjusted sample only unless no weights, subclasses, or matching strata are specified. Multiple values and abbreviations allowed.
 }
  \item{which.sub}{
\code{numeric}; if subclassification is used, a vector corresponding to the subclass(es) for which the distributions are to be displayed. If \code{NULL} (the default), distributions from all subclasses are displayed in a grid.
}
  \item{cluster}{
optional; a vector of cluster membership, or the name of a variable in an available data set passed to \code{bal.plot()} that contains cluster membership.
}
  \item{which.cluster}{
if clusters are used, which cluster(s) to display. Can be cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. If \code{NULL} (the default) or \code{NA}, all clusters are displayed.
}
  \item{imp}{
optional; a vector of imputation indices, or the name of a variable in an available data set passed to \code{bal.plot()} that contains imputation indices.
}
  \item{which.imp}{
if imputations are used, which imputations(s) to display. Must be numerical indices for which to display balance. If \code{NULL} (the default) or \code{NA}, all imputations are displayed.
}
  \item{which.treat}{
which treatment groups to display. If \code{NULL} (the default) or \code{NA}, all treatment groups are displayed.
}
  \item{which.time}{
for longitudinal treatments, which time points to display. Can be treatment names or time period indices. If \code{NULL} (the default) or \code{NA}, all time points are displayed.
}
  \item{mirror}{
\code{logical}; if the treatment is binary, the covariate is continuous, and densities or histograms are requested, whether to display mirrored densities/histograms or overlapping densities/histograms. Ignored otherwise. 
}
  \item{type}{
\code{character}; for binary and multi-category treatments with a continuous covariate, whether to display densities (\code{"density"}), histograms  (\code{"histogram"}), or empirical cumulative density function plots (\code{"ecdf"}). The default is to display densities. Abbreviations are allowed.
}
  \item{colors}{
a vector of colors for the plotted densities/histograms. See 'Color Specification' at \code{\link[graphics]{par}}. Defaults to the default \pkg{ggplot2} colors.
}
  \item{grid}{
\code{logical}; whether gridlines should be shown on the plot. Default is \code{TRUE}.
}
  \item{sample.names}{
\code{character}; new names to be given to the samples (i.e., in place of "Unadjusted Sample" and "Adjusted Sample"). For example, when matching it used, it may be useful to enter \code{c("Unmatched", "Matched")}.
}
  \item{position}{
the position of the legend. This can be any value that would be appropriate as an argument to \code{legend.position} in \code{\link[ggplot2]{theme}}.
}
  \item{facet.formula}{
a \code{formula} designating which facets should be on the rows and columns. This should be of the "historical" formula interface to \code{\link[ggplot2]{facet_grid}}. If of the form \code{a ~ b}, \code{a} will be faceted on the rows and \code{b} on the columns. To only facet on the rows, provide a one-sided formula with an empty left-hand side. To only facet on the columns, the formula should be of the form \code{a ~ .} (i.e., with only \code{.} on the right-hand side). The allowable facets depend on which arguments have been supplied to \code{bal.plot}; possible values include \code{which}, \code{cluster}, \code{imp}, and (for longitudinal treatments) \code{time}. If \code{NULL}, \code{bal.plot} will decide what looks best; this argument exists in case you disagree with its choice.
}
  \item{disp.means}{
\code{logical}; for a categorical treatment with a continuous covariate, whether a line should be drawn for each treatment level denoting the (weighted) mean of the covariate. Ignored if \code{type} is not "density" or "histogram". Default is \code{FALSE}.
}
  \item{alpha.weight}{
\code{logical}; if both the treatment and the covariate are continuous, whether points should be shaded according to their weight. Fainter points are those that have smaller weights. Default is \code{TRUE}.
}

}
\section{Additional Arguments}{
\code{bal.plot()} works like \code{bal.tab()} in that it can take a variety of types of inputs and yield the same output for each. Depending on what kind of input is given, different additional parameters are required in \code{\dots}. For details on what is required and allowed for each additional input and their defaults, see the help file for the \code{\link{bal.tab}} method associated with the input. The following are the required additional arguments based on each input type:
\itemize{
\item For \code{matchit} objects: None

\item For \code{weightit} objects: None

\item For \code{ps}, \code{ps.cont}, \code{mnps}, and \code{iptw} objects: (\code{stop.method}; see \link[=bal.tab.ps]{defaults}).

\item For \code{Match} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

\item For \code{optmatch} objects: \code{formula} and \code{data} or \code{covs} (\code{treat} is not required).

\item For \code{CBPS} objects: None

\item For \code{ebalance} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

\item For \code{formula}s: \code{data}

\item For \code{data.frame}s: \code{treat}

\item For \code{designmatch} objects: \code{formula} and \code{data} or \code{covs} and \code{treat}.

\item For \code{sbw} objects: None

\item For \code{mimids} and \code{wimids} objects: None, but an argument to \code{which.imp} should be specified.

\item For other objects processed through \code{bal.tab}'s default method, whichever arguments are required to identify treatment, variables, and a conditioning method (if any).
}
}
\details{
\code{bal.plot()} uses \code{ggplot()} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. For categorical treatments with continuous covariates or continuous treatments with categorical covariates, density plots are created using \code{\link[ggplot2]{geom_density}}, histograms are created using \code{\link[ggplot2]{geom_histogram}}, and empirical CDF plots are created using \code{\link[ggplot2:geom_path]{geom_step}}; for categorical treatments with categorical covariates, bar graphs are created using \code{\link[ggplot2]{geom_bar}}; for continuous treatments with continuous covariates, scatterplots are created using \code{\link[ggplot2]{geom_point}}.

For continuous treatments with continuous covariates, four additional lines are presented for aid in balance assessment. The red line is the linear fit line. The blue line is a smoothing curve generated with \pkg{ggplot2}'s \code{\link[ggplot2]{geom_smooth}} with \code{method = "auto"}. The horizontal black line is a horizontal reference line intercepting the (unweighted) treatment mean. The vertical black line is a reference line intercepting the (unweighted) treatment mean. Balance is indicated by the flatness of both fit lines and whether they pass through the intersection of the two black reference lines.

When multiple plots are to be displayed (i.e., when requesting subclass balance, cluster balance, or imputation balance, or when multiple sets of weights are provided or \code{which = "both"}, or when treatment is longitudinal), the plots will be displayed in a grid using \pkg{ggplot2}'s \code{\link[ggplot2]{facet_grid}}. Subclassification cannot be used with clusters or multiply imputed data.
}
\value{
A \code{"ggplot"} object, returned invisibly.
}

\author{
Noah Greifer 
}


\seealso{
\code{\link{bal.tab}}
}
\examples{
data("lalonde", package = "cobalt")

#Nearest Neighbor Matching
library(MatchIt)
m.out <- matchit(treat ~ age + educ + race + 
                 married + nodegree + re74 + re75, 
                 data = lalonde)

bal.plot(m.out, "age", which = "both")
bal.plot(m.out, "re74", which = "both", type = "ecdf")
bal.plot(m.out, "race", which = "both")
bal.plot(m.out, "distance", which = "both", mirror = TRUE,
         type = "histogram", colors = c("white", "black"))
         
#PS weighting with a continuous treatment
library(WeightIt)
w.out <- weightit(re75 ~ age + I(age^2) + educ + 
                  race + married + nodegree,
                  data = lalonde)
                  
bal.plot(w.out, "age", which = "both")
bal.plot(w.out, "married", which = "both")

}
\keyword{plots}