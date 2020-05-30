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
          stats, 
          abs, 
          agg.fun = NULL, 
          var.order = NULL, 
          drop.missing = TRUE, 
          drop.distance = FALSE, 
          thresholds = NULL, 
          line = FALSE, 
          stars = "none", 
          grid = FALSE, 
          limits = NULL, 
          colors = NULL, 
          shapes = NULL, 
          alpha = 1, 
          size = 3, 
          wrap = 30, 
          var.names = NULL, 
          title, 
          sample.names, 
          labels = FALSE,
          position = "right", 
          themes = NULL,
          ...)
}

\arguments{
  \item{x}{
  the valid input to a call to \code{\link{bal.tab}()} (e.g., the output of a preprocessing function). Other arguments that would be supplied to \code{bal.tab()} can be entered with \code{...}. Can also be a \code{bal.tab} object, i.e., the output of a call to \code{bal.tab()}. See Examples. If \code{x} is not a \code{bal.tab} object, \code{love.plot()} calls \code{bal.tab()} with the arguments supplied.
}
  \item{stats}{
\code{character}; which statistic(s) should be reported. See \code{\link[=balance.stats]{stats}} for allowable options. For binary and multi-category treatments, "mean.diffs" (i.e., mean differences) is the default. For continuous treatments, "correlations" (i.e., treatment-covariate Pearson correlations) is the default. Multiple options are allowed.
}
  \item{abs}{
\code{logical}; whether to present the statistic in absolute value or not. For variance ratios, this will force all ratios to be greater than or equal to 1. If \code{x} is a \code{bal.tab} object, \code{love.plot()} might ignore \code{abs} depending on the original \code{bal.tab()} call. If unspecified, uses whatever was used in the call to \code{bal.tab()}.
}
  \item{agg.fun}{
if balance is to be displayed across clusters or imputations rather than within a single cluster or imputation, which summarizing function ("mean", "max", or "range") of the balance statistics should be used. If "range" is entered, \code{love.plot} will display a line from the min to the max with a point at the mean for each covariate. Abbreviations allowed; "range" is default. Remember to set \code{which.<ARG> = .none} (where \code{<ARG>} is the grouping argument, such as \code{cluster} or \code{imp}) to use \code{agg.fun}. See Details.
}
  \item{var.order}{
a \code{character} or \code{love.plot} object; how to order the variables in the plot. See Details. 
}
  \item{drop.missing}{
\code{logical}; whether to drop rows for variables for which the statistic has a value of \code{NA}, for example, variance ratios for binary variables. If \code{FALSE}, there will be rows for these variables but no points representing their value. Default is \code{TRUE}, so that variables with missing balance statistics are absent. When multiple \code{stats} are requested, only variables with \code{NA}s for all \code{stats} will be dropped if \code{drop.missing = TRUE}. This argument used to be called \code{no.missing}, and that name still works (but has been deprecated).
}
  \item{drop.distance}{
\code{logical}; whether to ignore the distance measure (if there are any) in plotting.
}
  \item{thresholds}{
\code{numeric}; an optional value to be used as a threshold marker in the plot. Should be a named vector where each name corresponds to the statistic for which the threshold is to be applied. See example at \code{\link[=balance.stats]{stats}}. If \code{x} is a \code{bal.tab} object and a threshold was set in it (e.g., with \code{thresholds}), its threshold will be used unless overridden using the \code{threshold} argument in \code{love.plot}.
}
  \item{line}{
\code{logical}; whether to display a line connecting the points for each sample.
}
  \item{stars}{
when mean differences are to be displayed, which variable names should have a star (i.e., an asterisk) next to them. Allowable values are "none", "std" (for variables with mean differences that have been standardized), or "raw" (for variables with mean differences that have not been standardized). If "raw", the x-axis title will be "Standardized Mean Differences". Otherwise, it will be "Mean Differences". Ignored when mean difference are not displayed. See Details for an explanation of the purpose of this option.
}
  \item{grid}{
\code{logical}; whether gridlines should be shown on the plot. Default is \code{FALSE}.
}
  \item{limits}{
\code{numeric}; the bounds for the x-axis of the plot. Must a (named) list of vectors of length 2 in ascending order, one for each value of \code{stats} that is to have limits; e.g., \code{list(m = c(-.2, .2))}. If values exceed the limits, they will be plotted at the edge.
}
  \item{colors}{
the colors of the points on the plot. See 'Color Specification' at \code{\link[graphics]{par}} or the \code{ggplot2} \href{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#colour-and-fill}{aesthetic specifications} page. The first value corresponds to the color for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to the default \pkg{ggplot2} colors.
}
  \item{shapes}{
the shapes of the points on the plot. Must be one or two numbers between 1 and 25 or the name of a valid shape. See the \code{ggplot2} \href{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}{aesthetic specifications} page for valid options. Values 15 to 25 are recommended. The first value corresponds to the shape for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to 19 (\code{"circle filled"}).
}
  \item{alpha}{
\code{numeric}; the transparency of the points. See \code{\link[ggplot2]{scale_alpha}}.
}
  \item{size}{
\code{numeric}; the size of the points on the plot. Defaults to 3. In previous versions, the size was scaled by a factor of 3. Now \code{size} corresponds directly to the \code{size} aesthetic in \code{\link[ggplot2]{geom_point}}.
}
  \item{wrap}{
\code{numeric}; the number of characters at which to wrap axis labels to the next line. Defaults to 30. Decrease this if the axis labels are excessively long.
}
  \item{var.names}{
an optional object providing alternate names for the variables in the plot, which will otherwise be the variable names as they are stored. This may be useful when variables have ugly names. See Details on how to specify \code{var.names}. \code{link{var.names}} can be a useful tool for extracting and editing the names from the \code{bal.tab} object.
}
  \item{title}{
\code{character}; the title of the plot.
}
  \item{sample.names}{
\code{character}; new names to be given to the samples (i.e., in place of "Unadjusted" and "Adjusted"). For example, when matching it used, it may be useful to enter \code{c("Unmatched", "Matched")}.
}
  \item{labels}{
\code{logical} or \code{character}; labels to give the plots when multiple \code{stats} are requested. If \code{TRUE}, the labels will be capital letters. Otherwise, must be a string with the same length as \code{stats}. This can be useful when the plots are to be used in an article.
}
  \item{position}{
the position of the legend. When \code{stats} has length 1, this can be any value that would be appropriate as an argument to \code{legend.position} in \code{\link[ggplot2]{theme}}. When \code{stat} has length greater than 1, can be one of "none", "left", "right", "bottom", or "top".
}
  \item{themes}{
an optional list of \code{theme} objects to append to each individual plot. Each entry should be the output of a call to \code{\link[ggplot2]{theme}} in \pkg{ggplot2}. This is a way to customize the individual plots when multiple \code{stats} are requested since the final output is not a manipulable \code{ggplot} object. It can be used with length-1 \code{stats}, but it probably makes more sense to just add the \code{theme()} call after \code{love.plot()}.
}
  \item{...}{
additional arguments passed to \code{bal.tab()} or options for display of the plot. The following related arguments are currently accepted:
\describe{
  \item{\code{use.grid}}{whether to use \code{\link[gridExtra]{arrangeGrob}} in \code{gridExtra} to make the plot when \code{stats} has length 1. See section Value.}
  \item{\code{disp.subclass}}{whether to display individual subclasses if subclassification is used. Overrides the \code{disp.subclass} option in the original \code{bal.tab()} call if \code{x} is a \code{bal.tab} object.}
  \item{\code{star_char}}{\code{character}; when \code{stars} are used, the character that should be the "star" next to the starred variables. The default is \code{"*"}. \code{"†"} or \code{"\u2020"} (i.e., dagger) might be appealing as well.}
  }
Additionally, any of the \code{which.} arguments used with clustered or multiply imputed data or longitudinal or multi-category treatments can be specified to display balance on selected groupings. Set to \code{.none} to aggregate across groups (in which \code{agg.fun} comes into effect) and set to \code{.all} to view all groups. See \code{display_options} for options, and see the vignette "Appendix 2: Using cobalt with Clustered, Multiply Imputed, and Other Segmented Data" for details and examples.
}
}
\details{

\code{love.plot} can be used with clusters, imputations, or both in addition to the standard case. The cluster or imputation arguments must be specified in the call to \code{bal.tab}. Several types of plots can be requested: a plot aggregating over all imputations across all clusters, a plot displaying individual clusters aggregating over imputations (if any), a plot displaying individual imputations across clusters, or a plot displaying individual clusters within one imputation (if any) or individual imputations for one cluster. The choice of these displays is controlled by the arguments to \code{which.cluster} and \code{which.imp}. If either of them are \code{.none}, the plot will aggregate over that collection. If either are individual values, the plot will display the values for those specific clusters or imputations. If either of them are \code{.all}, the plot will display the values for all clusters or imputations. If both clusters and imputations are specified, at least one of \code{which.cluster} or \code{which.imp} must be a single value, or \code{.none}. When aggregating, an argument should be specified to \code{agg.fun} referring to whether the mean, minimum ("min"), or maximum ("max") balance statistic or range ("range", the default) of balance statistics for each covariate should be presented in the plot. See the vignette "Appendix 2: Using cobalt with Clustered, Multiply Imputed, and Other Segmented Data" for examples.

With subclasses, balance will be displayed for the unadjusted sample and the aggregated subclassified sample. If \code{disp.subclass} is \code{TRUE}, each subclass will be displayed additionally as a number on the plot. 

\bold{Variable order using \code{var.order}}

The order that the variables are presented in depends on the argument to \code{var.order}. If \code{NULL}, the default, they will be displayed in the same order as in the call to \code{bal.tab}, which is the order of the underlying data set. If "alphabetical", they will be displayed in alphabetical order. If "unadjusted", they will be ordered by the balance statistic of the unadjusted sample. To order by the values of the adjusted sample, "adjusted" can be supplied if only one set of weights (or subclasses) are specified; otherwise, the name of the set of weights should be specified.

If multiple \code{stats} are requested, the order will be determined by the first entry to \code{stats} (e.g., if both "mean.diffs" and "ks.statistics" are requested, and \code{var.order = "unadjusted"}, the variables will be displayed in order of the unadjusted mean differences for both plots). If multiple plots are produced simultaneously (i.e., for individual clusters or imputations), \code{var.order} can only be \code{NULL} or "alphabetical".

If a \code{love.plot} object is supplied, the plot being drawn will use the variable order in the supplied \code{love.plot} object. This can be useful when making more than one plot and the variable order should be the same across plots.

\bold{Variable names using \code{var.names}}

The default in \code{love.plot} is to present variables as they are named in the output of the call to \code{bal.tab}, so it is important to know this output before specifying alternate variable names when using \code{var.names}, as the displayed variable names may differ from those in the original data.

There are several ways to specify alternate names for presentation in the displayed plot using the \code{var.names} argument by specifying a list of old and new variable names, pairing the old name with the new name. You can do this in three ways: 1) use a vector or list of new variable names, with the \code{names} of the values the old variable names; 2) use a data frame with exactly one column containing the new variable names and the row names containing the old variable names; or 3) use a data frame with two columns, the first (or the one named "old") containing the old variable names and the second (or the one named "new") containing the new variable names. If a variable in the output from \code{bal.tab} is not provided in the list of old variable names, \code{love.plot} will use the original old variable name.

\code{love.plot} can replace old variables names with new ones based on exact matching for the name strings or matching using the variable name components. For example, if a factor variable \code{"X"} with levels \code{"a"}, \code{"b"}, and \code{"c"} is displayed with \code{love.plot}, the variables \code{"X_a"}, \code{"X_b"}, and \code{"X_c"} will be displayed. You can enter replacement names for all three variables individually with \code{var.names}, or you can simply specify a replacement name for \code{"X"}, and \code{"X"} will be replaced by the given name in all instances it appears, including not just factor expansions, but also polynomials and interactions in \code{int = TRUE} in the original \code{bal.tab} call. In an interaction with another variable, say \code{"Y"}, there are several ways to replace the name of the interaction term \code{"X_a * Y"}. If the entire string (\code{"X_a * Y"}) is included in \code{var.names}, the entire string will be replaced. If \code{"X_a"} is included in \code{var.name}, only it will be replaced (and it will be replaced everywhere else it appears). If \code{"X"} is included in \code{var.name}, only it will be replaced (and it will be replaced everywhere else it appears). See example at \code{\link{var.names}}.

\bold{Stars and the x-axis label with mean differences}

When mean differences are to be displayed, \code{love.plot()} attempts to figure out the appropriate label for the x-axis. If all mean differences are standardized, the x-axis label will be "Standardized Mean Differences". If all mean differences are raw (i.e., unstandardized), the x-axis label will be "Mean Differences". Otherwise, \code{love.plot()} turns to the \code{stars} argument. If "raw", the x-axis label will be "Standardized Mean Differences" (i.e., because un-starred variables have standardized mean differences displayed). If "std", the x-axis label will be "Mean Differences" (i.e., because un-starred variables have raw mean differences displayed). If "none", the x-axis label will be "Mean Differences" and a warning will be issued recommending the use of \code{stars}. 

The default is to display standardized mean differences for continuous variables, raw mean differences for binary variables, and no stars, so this warning will be issued in most default uses of \code{love.plot()}. The purpose of this is to correct behavior of previous versions of \pkg{cobalt} in which the default x-axis label was "Mean Differences", even when standardized mean differences were displayed, yielding a potentially misleading plot. This warning requires the user to think about what values are being displayed. The idea of using \code{stars} is that the user can, in a caption for the plot, explain that variables with an asterisk have standardized (or raw) mean differences display, in contrast to un-starred variables.

}
\value{
When only one type of balance statistic is requested, the returned object is a standard \code{ggplot} object that can be manipulated using \pkg{ggplot2} syntax. This facilitates changing fonts, background colors, and features of the legend outside of what \code{love.plot} provides automatically. 

When more than one type of balance statistic is requested, the plot is constructed using \code{\link[gridExtra]{arrangeGrob}} in \code{gridExtra}, which arranges multiple plots and their shared legend into one plot. Because the output of \code{arrangeGrob} is a \code{gtable} object, its features cannot be manipulated in the standard way. Use the \code{themes} argument to change theme elements of the component plots. The original plots are stored in the \code{"plots"} attribute of the output object.
}
\note{
\code{love.plot} can also be called by using \code{plot} or \code{autoplot} on a \code{bal.tab} object. If used in this way, some messages may appear twice. It is recommended that you just use \code{love.plot} instead.
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

love.plot(m.out1, thresholds = c(m = .1), var.order = "unadjusted")

## Using alternate variable names
v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
                        "race_white", "married", "nodegree", "re74", 
                        "re75", "distance"),
                new = c("Age", "Years of Education", "Black", 
                        "Hispanic", "White", "Married", "No Degree", 
                        "Earnings 1974", "Earnings 1975", 
                        "Propensity Score"))
                
love.plot(m.out1, stats = "m", threshold = .1, 
          var.order = "unadjusted", var.names = v)
          
#Using multiple stats
love.plot(m.out1, stats = c("m", "ks"), 
          thresholds = c(m = .1, ks = .05), 
          var.order = "unadjusted", var.names = v, stars = "raw",
          position = "bottom", wrap = 20)
          
#Changing visual elements
love.plot(m.out1, thresholds = c(m = .1), 
          var.order = "unadjusted", var.names = v, abs = TRUE,
          shapes = c("triangle filled", "circle"), 
          colors = c("red", "blue"), line = TRUE,
          grid = FALSE, sample.names = c("Unmatched", "Matched"),
          stars = "raw", position = "top")

}
\keyword{plots}