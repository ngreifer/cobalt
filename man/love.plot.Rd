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
    agg.fun = c("range", "max", "mean"), 
    colors = NULL, 
    shapes = NULL, 
    line = FALSE, 
    stars = "none",
    ...)
}

\arguments{
  \item{x}{
  the valid input to a call to \code{\link{bal.tab}()} (e.g., the output of a preprocessing function). Other arguments that would be supplied to \code{bal.tab()} can be entered with \code{...}. Can also be a \code{bal.tab} object, the output of a call to \code{bal.tab()}. See Examples. If \code{x} is not a \code{bal.tab} object, \code{love.plot()} calls \code{bal.tab()} with the arguments supplied.
}
  \item{stat}{
\code{character}; which statistic(s) should be reported. For binary or multinomial treatments, the options are "mean.diffs" for mean differences (standardized or not according the selected \code{bal.tab} options), "variance.ratios" for variance ratios, and "ks.statistics" for Kolmogorov-Smirnov statistics. "mean.diffs" is the default. For continuous treatments, the only option is "correlations" for treatment-covariate correlations. Multiple options are allowed. Abbreviations allowed.
}
  \item{threshold}{
\code{numeric}; an optional value to be used as a threshold marker in the plot. Can be either an unnamed vector with the same length as \code{stat}, where \code{NA} indicates no threshold, or a named vector containing the thresholds desired. For example, if \code{stat = c("mean.diffs", "ks.statistics")} one could set \code{threshold = c(.1, NA)} or \code{threshold = c(mean.diffs = .1)} to set a threshold only for mean differences. The names can be abbreviated (e.g., \code{threshold = c(m = .1)} for mean differences). If \code{x} is a \code{bal.tab} object and a theshold was set in it, its threshold will be used unless overridden using the \code{threshold} argument in \code{love.plot}.
}
  \item{abs}{
\code{logical}; whether to present the statistic in absolute value or not. For variance ratios, this will force all ratios to be greater than or equal to 1. If \code{x} is a \code{bal.tab} object, \code{love.plot()} might ignore \code{abs} depending on the original \code{bal.tab()} call.
}
  \item{var.order}{
\code{character}; how to order the variables in the plot. If \code{NULL}, they will be displayed in the same order as in the call to \code{bal.tab}, which is the order of the underlying data set. If "alphabetical", they will be displayed in alphabetical order. If "adjusted", they will be ordered by the balance statistic of the adjusted sample. if "unadjusted", they will be ordered by the balance statistic of the unadjusted sample. "unadjusted" looks the nicest, but \code{NULL} or "alphabetical" should be used if comparing variables across data sets to maintain variable order. If multiple \code{stat}s are requested, the order will be determined by the first \code{stat} (e.g., if "mean.diffs" and "ks.statistics" are requested, and \code{var.order = "unadjusted"}, the variables iwll be displayed in order of the unadjusted mean differences for both plots). If multiple plots are produced simultaneously (i.e., for individual clusters or imputations), \code{var.order} cannot be "unadjusted" or "adjusted".
}
  \item{no.missing}{
\code{logical}; whether to drop rows for variables for which the statistic has a value of \code{NA}, for example, variance ratios for binary variables. If \code{FALSE}, there will be rows for these variables but no points representing their value. Default is \code{TRUE}, so that variables with missing balance statistics are absent. When multiple \code{stat}s are requested, only variables with \code{NA}s for all \code{stat}s will be dropped if \code{no.missing = TRUE}.
}
  \item{var.names}{
an optional object providing alternate names for the variables in the plot, which will otherwise be the variable names as they are stored. This may be useful when variables have ugly names. See Details on how to specify \code{var.names}. \code{link{var.names}} can be a useful tool for extracting and editing the names from the \code{bal.tab} object.
}
  \item{drop.distance}{
\code{logical}; whether to ignore the distance measure (if there are any) in plotting.
}
  \item{agg.fun}{
if balance is to be displayed across clusters or imputations rather than within a single cluster or imputation, which summarizing function ("mean", "max", or "range") of the balance statistics should be used. If "range" is entered, \code{love.plot} will display a line from the min to the max with a point at the mean for each covariate. Abbreviations allowed; "range" is default. Remember to set \code{which.<ARG> = .none} (where \code{<ARG>} is the grouping argument, such as \code{cluster} or \code{imp}) to use \code{agg.fun}.
}
  \item{colors}{
the colors of the points on the plot. See 'Color Specification' at \code{\link[graphics]{par}} or the \code{ggplot2} \href{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#colour-and-fill}{aesthetic specifications} page. The first value corresponds to the color for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to the default \pkg{ggplot2} colors.
  }
  \item{shapes}{
the shapes of the points on the plot. Must be one or two numbers between 1 and 25 or the name of a valid shape. See the \code{ggplot2} \href{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}{aesthetic specifications} page for valid options. Values 15 to 25 are recommended. The first value corresponds to the shape for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to 19 (\code{"circle filled"}).
  }
  \item{line}{
  \code{logical}; whether to display a line connecting the points for each sample.
  }
  \item{stars}{
when mean differences are to be displayed, which variable names should have a star (i.e., an asterisk) next to them. Allowable values are "none", "std" (for variables with mean differences that have been standardized), or "raw" (for variables with mean differences that have not been standardized). If "raw", the x-axis title will be "Standardized Mean Differences". Otherwise, it will be "Mean Differences". Ignored when mean difference are not displayed. See Details for an explanation of the purpose of this option.
  }
  \item{...}{
additional arguments. These can be arguments that are passed to \code{bal.tab()} or options for display of the plot. The following display-related arguments are currently accepted:
\describe{
    \item{\code{size}}{\code{numeric}; the size of the points on the plot. Defaults to 1.}
    \item{\code{alpha}}{\code{numeric}; the transperancy of the points. See \code{\link[ggplot2]{scale_alpha}}.}
    \item{\code{grid}}{\code{logical}; whether gridlines should be shown on the plot. Default to \code{TRUE}.}
    \item{\code{title}}{\code{character}; the title of the plot.}
    \item{\code{subtitle}}{\code{character}; the subtitle of the plot.}
    \item{\code{sample.names}}{\code{character}; new names to be given to the samples (i.e., in place of "Unadjusted" and "Adjusted").}
    \item{\code{limits}}{\code{numeric}; the bounds for the x-axis of the plot. Must a vector of length 2 in ascending order.}
    \item{\code{wrap}}{\code{numeric}; the number of characters at which to wrap labels to the next line. See \code{\link[ggplot2]{label_wrap_gen}}; this value is supplied to the \code{width} argument. Defaults to 20.}
    \item{\code{disp.subclass}}{whether to display individual subclasses if subclassification is used. Overrides the \code{disp.subclass} option in the original \code{bal.tab()} call if \code{x} is a \code{bal.tab} object.}
  }
Additionally, any of the \code{which.} arguments used with clustered or multiply imputed data or longitudinal or multinomial treatments can be specified to display balance on selected groupings. Set to \code{.none} to aggregate across groups (in which \code{agg.fun} comes into effect) and set to \code{.all} to view all groups. See \code{display_options} for options, and see the vignette "Appendix 2: Using cobalt with Clustered and Multiply Imputed Data" for details and examples.
}
}
\details{
\code{love.plot} uses \code{ggplot} from the \pkg{ggplot2} package, and (invisibly) returns a \code{"ggplot"} object. This means that users can edit aspects of the plot using ggplot2 syntax.

The default in \code{love.plot} is to present variables as they are named in the output of the call to \code{bal.tab}, so it is important to know this output before specifying alternate variable names when using \code{var.names}, as the displayed variable names may differ from those in the original data.

There are several ways to specify alternate names for presentation in the displayed plot using the \code{var.names} argument by specifying a list of old and new variable names, pairing the old name with the new name. You can do this in three ways: 1) use a vector or list of new variable names, with the \code{names} of the values the old variable names; 2) use a data frame with exactly one column containing the new variable names and the row names containing the old variable names; or 3) use a data frame with two columns, the first (or the one named "old") containing the old variable names and the second (or the one named "new") containing the new variable names. If a variable in the output from \code{bal.tab} is not provided in the list of old variable names, \code{love.plot} will use the original old variable name.

\code{love.plot} can replace old variables names with new ones based on exact matching for the name strings or matching using the variable name components. For example, if a factor variable \code{"X"} with levels \code{"a"}, \code{"b"}, and \code{"c"} is displayed with \code{love.plot}, the variables \code{"X_a"}, \code{"X_b"}, and \code{"X_c"} will be displayed. You can enter replacement names for all three variables individually with \code{var.names}, or you can simply specify a replacement name for \code{"X"}, and \code{"X"} will be replaced by the given name in all instances it appears, including not just factor expansions, but also polynomials and interactions in \code{int = TRUE} in the original \code{bal.tab} call. In an interaction with another variable, say \code{"Y"}, there are several ways to replace the name of the interaction term \code{"X_a * Y"}. If the entire string (\code{"X_a * Y"}) is included in \code{var.names}, the entire string will be replaced. If \code{"X_a"} is included in \code{var.name}, only it will be replaced (and it will be replaced everywhere else it appears). If \code{"X"} is included in \code{var.name}, only it will be replaced (and it will be replaced everywhere else it appears). See example at \code{\link{var.names}}.

\code{love.plot} can be used with clusters, imputations, or both as well. The cluster or imputation arguments must be specified in the call to \code{bal.tab}. Several types of plots can be requested: a plot aggregating over all imputations across all clusters, a plot displaying individual clusters aggregating over imputations (if any), a plot displaying individual imputations across clusters, or a plot displaying individual clusters within one imputation (if any) or individual imputations for one cluster. The choice of these displays is controlled by the arguments to \code{which.cluster} and \code{which.imp}. If either of them are \code{.none}, the plot will aggregate over that collection. If either are individual values, the plot will display the values for those specific clusters or imputations. If either of them are \code{.all}, the plot will display the values for all clusters or imputations. If both clusters and imputations are specified, at least one of \code{which.cluster} or \code{which.imp} must be a single value, or \code{.none}. When aggregating, an argument should be specified to \code{agg.fun} referring to whether the mean, minimum ("min"), or maximum ("max") balance statistic or range ("range") of balance statistics for each covariate should be presented in the plot.

With subclasses, balance will be displayed for the unadjusted sample and the aggregated subclassified sample. If \code{disp.subclass} is \code{TRUE}, each subclass will be displayed additionally as a number on the plot. 

\bold{Stars and the x-axis label with mean differences}

When mean differences are to be displayed, \code{love.plot()} attempts to figure out the appropriate label for the x-axis. If all mean differences are standardized, the x-axis label will be "Standardized Mean Differences". If all mean differences are raw (i.e., unstandardized), the x-axis label will be "Mean Differences". Otherwise, \code{love.plot()} turns to the \code{stars} argument. If "raw", the x-axis label will be "Standardized Mean Differences" (i.e., because un-starred variables have standardized mean differences displayed). If "std", the x-axis label will be "Mean Differences" (i.e., because un-starred variables have raw mean differences displayed). If "none", the x-axis label will be "Mean Differences" and a warning will be issued recommending the use of \code{stars}. 

The default is to display standardized mean differences for continuous variables, raw mean differences for binary variables, and no stars, so this warning will be issued in most default uses of \code{love.plot()}. The purpose of this is to correct behavior of previous versions of \pkg{cobalt} in which the default x-axis label was "Mean Differences", even when standardized mean differences were displayed, yielding a potentially misleading plot. This warning requires the user to think about what values are being displayed. The idea of using \code{stars} is that the user can, in a caption for the plot, explain that variables with an asterisk have standardized (or raw) mean diffences display, in contrast to un-starred variables.

}
\value{
A \code{"ggplot"} object, returned invisibly.
}
\note{
\code{love.plot} can also be called by using \code{plot} on a \code{bal.tab} object. That is, the \code{"love."} prefix is optional, unless it is being called directly on the preprocessing output objects.
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

love.plot(m.out1, stat = "mean.diffs", 
          threshold = .1, var.order = "unadjusted")

#This gives the same results using older syntax
love.plot(bal.tab(m.out1), stat = "mean.diffs", 
          threshold = .1, var.order = "unadjusted")

## Using alternate variable names
v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
                        "race_white", "married", "nodegree", "re74", 
                        "re75"),
                new = c("Age", "Years of Education", "Black", 
                        "Hispanic", "White", "Married", "No Degree", 
                        "Earnings 1974", "Earnings 1975"))
                
love.plot(m.out1, stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", var.names = v)
          
#Changing colors and shapes and adding lines
love.plot(m.out1, stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", var.names = v, abs = TRUE,
          shapes = c("triangle filled", "circle"), 
          colors = c("red", "blue"), line = TRUE)

}
