\name{love.plot}
\alias{love.plot}
\alias{plot.bal.tab}
\title{
Generate Balance Plots for Publication
}
\description{
Generates a "Love" plot graphically displaying covariate balance before and after adjusting.
}
\usage{
love.plot(x, 
    stat = c("mean.diffs", "variance.ratios", 
             "ks.statistics"), 
    threshold = NULL, 
    abs = TRUE, 
    var.order = NULL, 
    no.missing = TRUE, 
    var.names = NULL, 
    drop.distance = FALSE, 
    agg.fun = c("mean", "median", "max", "range"), 
    colors = NULL, 
    shapes = NULL, 
    line = FALSE, 
    ...)
}

\arguments{
  \item{x}{
a \code{bal.tab} object; the output of a call to \code{bal.tab()}. \code{m.threshold}, \code{v.threshold}, and \code{r.threshold} can be used in \code{bal.tab()} instead of \code{love.plot()}'s \code{threshold} parameter. 
}
  \item{stat}{
\code{character}; which statistic should be reported if treatment is binary. The options are "mean.diffs" for mean differences (standardized or not according the options selected in \code{bal.tab} object), "variance.ratios" for variance ratios, and "ks.statistics" for Kolmogorov-Smirnov statistics. "mean.diffs" is the default. Abbreviations allowed.
}
  \item{threshold}{
\code{numeric}; an optional value to be used as a threshold marker in the plot. Overrides the threshold set in the \code{bal.tab} object.
}
  \item{abs}{
\code{logical}; whether to present the statistic in absolute value or not if \code{stat = "mean.diffs"} or the treatment variable is continuous. See Details.
}
  \item{var.order}{
\code{character}; how to order the variables in the plot. If \code{NULL}, they will be displayed in the same order as in the call to \code{bal.tab()}, which is the order of the underlying data set. If "alphabetical", they will be displayed in alphabetical order. If "adjusted", they will be ordered by the balance statistic of the adjusted sample. if "unadjusted", they will be ordered by the balance statistic of the unadjusted sample. "unadjusted" looks the nicest, but \code{NULL} or "alphabetical" should be used if comparing variables across data sets to maintain variable order. If multiple plots are produced simultaneously (i.e., for individual clusters or imputations), \code{var.order} cannot be "unadjusted" or "adjusted".
}
  \item{no.missing}{
\code{logical}; whether to drop rows for variables for which the statistic has a value of \code{NA}, for example, variance ratios for binary variables. If \code{FALSE}, there will be rows for these variables but no points representing their value. Default is \code{TRUE}, so that variables with missing values are absent.
}
  \item{var.names}{
an optional object providing alternate names for the variables in the plot, which will otherwise be the variable names as they are stored. This may be useful when variables have ugly names. See Details on how to specify \code{var.names}.
}
  \item{drop.distance}{
\code{logical}; whether to ignore the distance measure (if there are any) in plotting.
}
  \item{agg.fun}{
if balance is to be displayed across clusters or imputations rather than within a single cluster or imputation, which summarizing function (mean, median, max, or range) of the balance statistics sould be used. If "range" is entered, \code{love.plot()} will display a line from the min to the max with a point at the mean for each covariate; it can only be used if \code{quick = FALSE} in the \code{bal.tab()} call. Abbreviations allowed; "mean" is default.
}
  \item{colors}{
  The colors of the points on the plot. See 'Color Specification' at \code{\link[graphics]{par}}. The first value corresponds to the color for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to "red" and "blue".
  }
  \item{shapes}{
  \code{numeric}; the shapes of the points on the plot. See Example at \code{\link[ggplot2]{shape}}. Must be one or two numbers between 21 and 25. The first value corresponds to the shape for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to 21 (circle).
  }
  \item{line}{
  \code{logical}; whether to display a line connecting the points for each sample.
  }
  \item{...}{
Options for display of the plot. The following arguments are currently accepted:
\describe{
    \item{\code{size}}{\code{numeric}; the size of the points on the plot. Defaults to 1.
    }
    \item{\code{title}}{\code{character}; the title of the plot.}
    \item{\code{subtitle}}{\code{character}; the subtitle of the plot.}
    \item{\code{sample.names}}{\code{character}; new names to be given to the samples (i.e., in place of "Unadjusted" and "Adjusted").}
    \item{\code{limits}}{\code{numeric}; the bounds for the x-axis of the plot. Must a vector of length 2 in ascending order.}
    \item{\code{which.cluster}}{which clusters to display. Overrides the \code{which.cluster} option in the original \code{bal.tab} object.}
    \item{\code{which.imp}}{which imputations to display. Overrides the \code{which.imp} option in the original \code{bal.tab} object.}
    \item{\code{which.treat}}{which treatment groups to display. Overrides the \code{which.treat} option in the original \code{bal.tab} object.}
    \item{\code{disp.subclass}}{whether to display individual subclasses. Overrides the \code{disp.subclass} option in the original \code{bal.tab} object.}
  }
}
}
\details{
\code{love.plot()} uses \code{ggplot} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. This means that users can edit aspects of the plot using ggplot2 syntax.

The default in \code{love.plot()} is to present variables as they are named in the output of the call to \code{bal.tab()}, so it is important to know this output before specifying alternate variable names when using \code{var.names}, as the displayed variable names may differ from those in the original data.

There are several ways to specify alternate names for presentation in the displayed plot using the \code{var.names} argument by specifying a list of old and new variable names, pairing the old name with the new name. You can do this in three ways: 1) use a vector or list of new variable names, with the \code{names} of the values the old variable names; 2) use a data frame with exactly one column containing the new variable names and the row names containing the old variable names; or 3) use a data frame with two columns, the first (or the one named "old") containing the old variable names and the second (or the one named "new") containing the new variable names. If a variable in the output from \code{bal.tab()} is not provided in the list of old variable names, \code{love.plot()} will use the original old variable name.

\code{love.plot()} can be used with clusters, imuputations, or both as well. The cluster or imputation arguments must be specified in the call to \code{bal.tab()}. Several types of plots can be requested: a plot aggregating over all imputations across all clusters, a plot displaying individual clusters aggregating over imputations (if any), a plot displaying individual imputations across clusters, or a plot displaying individual clusters within one imputation (if any) or individual imputations for one cluster. The choice of these displays is controlled by the arguments to \code{which.cluster} and \code{which.imp}. If either of them are \code{NA}, the plot will aggregate over that collection. If either are individual values, the plot will display the values for those specific clusters or imputations. If either of them are \code{NULL}, the plot will display the values for all clusters or imputations. If both clusters and imputations are specified, at least one of \code{which.cluster} or \code{which.imp} must be a single value, or \code{NA}. When aggregating, an argument should be specified to \code{agg.fun} referring to whether the mean, minimum ("min"), or maximum ("max") balance statistic or range ("range") of balance statistics for each covariate should be presented in the plot.

With subclasses, balance will be displayed for the unadjusted sample and the aggregated subclassified sample. If \code{disp.subclass} is \code{TRUE}, each subclass will be displayed additionally as a number on the plot. 
If no aggregation (e.g., over clusters or imputations) is to take place, \code{abs} will be \code{TRUE} if it is \code{TRUE} in the call to \code{bal.tab()} or in the call to \code{love.plot()}. Otherwise, if balance is requested with multinomial or longitudinal treatments, \code{abs} is set to \code{TRUE}, and if balance is requested with clustered or multiply imputed data, \code{abs} will take on the value of \code{abs} set in the original call to \code{bal.tab()}.

}
\value{
A \code{"ggplot"} object, returned invisbly.
}
\note{
\code{love.plot()} can also be called by using \code{plot()} on a \code{bal.tab} object. That is, the "love." prefix is optional.
}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}}
}
\examples{
library(MatchIt); data("lalonde", package = "cobalt")

## Nearest Neighbor matching
m.out1 <- matchit(treat ~ age + educ + race + 
                  married + nodegree + re74 + re75, 
                  data = lalonde)

love.plot(bal.tab(m.out1), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted")

## Using alternate variable names
v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
                        "race_white", "married", "nodegree", "re74", 
                        "re75"),
                new = c("Age", "Years of Education", "Black", 
                        "Hispanic", "White", "Married", "No Degree", 
                        "Earnings 1974", "Earnings 1975"))
                
love.plot(bal.tab(m.out1), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", var.names = v)
          
#Changing colors and shapes and adding lines
love.plot(bal.tab(m.out1), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", var.names = v, abs = TRUE,
          shapes = c(22, 25), colors = c("darkblue", "lightblue"),
          line = TRUE)

}
