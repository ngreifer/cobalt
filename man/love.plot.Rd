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
love.plot(b, stat = c("mean.diffs", "variance.ratios"), 
          threshold = NULL, abs = FALSE, var.order = NULL, 
          no.missing = FALSE, var.names = NULL,
          drop.distance = TRUE)
}
\arguments{
  \item{b}{
a \code{bal.tab} object. Print options are ignored, but \code{m.threshold} and \code{v.threshold} can be used in \code{bal.tab()} instead of \code{love.plot()}'s \code{threshold} parameter. 
}
  \item{stat}{
\code{character}; which statistic should be reported. The options are "mean.diffs" for mean differences (standardized or not according the options selected in \code{bal.tab} object) and "variance.ratios" for variance ratios. "mean.diffs" is default. Abbreviations allowed.
}
  \item{threshold}{
\code{numeric}; an optional value to be used as a threshold marker in the plot. Overrides the threshold set in the \code{bal.tab} object.
}
  \item{abs}{
\code{logical}; whether to present the statistic in absolute value or not if \code{stat="mean.diffs"}.
}
  \item{var.order}{
\code{character}; how to order the variables in the plot. If \code{NULL}, they will be displayed in alphabetical order. If "adjusted", they will be ordered by the balance statistic of the adjusted sample. if "unadjusted", they will be ordered by the balance statistic of the unadjusted sample. "unadjusted" looks the nicest, but \code{NULL} should be used if comparing variables across data sets to maintain variable order.
}
  \item{no.missing}{
\code{logical}; whether to drop rows for variables for which the statistic has a value of \code{NA}, for example, variance ratios for binary variables. If \code{FALSE}, there will be rows for these variables but no points representing their value, and a warning message from ggplot2 will appear.
}
  \item{var.names}{
an optional object providing alternate names for the variables in the plot, which will otherwise be the variable names as they are stored. This may be useful when variables have ugly names. If \code{var.order} is \code{NULL}, the variables will be placed in alphabetical order of the new variable names. See Details on how to specify \code{var.names}.
}
  \item{drop.distance}{
\code{logical}; whether to ignore the distance measure (if there is one) in plotting. Because balance on the covariates is primary goal of conditioning, including balance on the distance measure can be misleading; therefore, the default is \code{TRUE}, so that the distance measure is not displayed.
}
}
\details{
\code{love.plot} uses \code{ggplot} from the ggplot2 package, and (invisibly) returns a \code{"ggplot"} object. This means that users can edit aspects of the plot using ggplot2 syntax.

The default in \code{love.plot()} is to present variables as they are named in the output of the call to \code{bal.tab()}, so it is important to know this output before specifying alternate variable names when using \code{var.names}, as the displayed variable names may differ from those in the original data. To see a vector of the variable names from \code{bal.tab()}, use \code{row.names(bal.tab(args)$Balance)} (or \code{row.names(bal.tab(args)$Balance.Across.Subclass} if subclassification is used). Note that if \code{drop.distance = TRUE}, which the default, the distance measure, if any, will not count as a variable below; otherwise, it will count.

There are several ways to specify alternate names for presentation in the displayed plot using the \code{var.names} argument. You can use a vector of alternate names the same length as the variable list output from \code{bal.tab()}, and \code{love.plot()} will use these names instead. To leave a variable name as is, enter \code{""} or \code{NA} in the position of that variable. Another way is to specify a list of old and new variable names, pairing the old name with the new name. You can do this in three ways: 1) use a vector of new variable names, with the \code{names} of the values the old variable names; 2) use a data frame with exactly one column containing the new variable names and the row names containing the old variable names; or 3) use a data frame with two columns, the first containing the old variable names and the second containing the new variable names. This third method is the safest because the coersion rules in \R are least likely to affect the input. If a variable in the output from \code{bal.tab()} is not provided in the list of old variable names, \code{love.plot()} will use the original old variable name.
}
\value{
A \code{"ggplot"} object, returned invisbly.
}

\author{
Noah Greifer \email{noah@unc.edu}
}
\references{
Love, T. (2004), "Graphical Display of Covariate Balance," \url{http://chrp.org/love/JSM2004RoundTableHandout.pdf}
}
\seealso{
\code{\link{bal.tab}}
}
\examples{
library(MatchIt); data("lalonde", package = "cobalt")

## Nearest Neighbor matching
m.out1 <- matchit(treat ~ age + educ + black + hispan + 
                  married + nodegree + re74 + re75, 
                  data = lalonde)

love.plot(bal.tab(m.out1), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted")

## Using alternate variable names
v <- data.frame(old = c("age", "educ", "black", "hispan", 
                        "married", "nodegree", "re74", "re75"),
                new = c("Age", "Years of Education", "Black", 
                        "Hispanic", "Married", "No Degree", 
                        "Earnings 1974", "Earnings 1975"))
                
love.plot(bal.tab(m.out1), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", var.names = v)

}
\keyword{design}