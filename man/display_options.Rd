\name{options-display}
\alias{options-display}

\title{Options for Displaying bal.tab Output}

\description{
Several additional arguments can be passed to \code{bal.tab} that control the display of the output but don't affect computations. These arguments are documented here.

Not all arguments are applicable to all uses of \code{bal.tab}; for example, \code{disp.subclass}, which controls whether individual subclasses are displayed when subclassification is used, won't do anything when subclassification is not used.

Note that when \code{quick = TRUE} is set in the call to \code{bal.tab}, setting any of these arguments to \code{FALSE} can prevent some values from being computed, which can have unintended effects. For example, setting \code{quick = TRUE} and \code{un = FALSE} in \code{bal.tab} and then using \code{love.plot} on the output will yield an error because \code{love.plot} depends on the computation of the unadjusted balance statistics, which are not computed when \code{quick = TRUE} and \code{un = FALSE}.
}
\arguments{
\item{imbalanced.only}{\code{logical}; whether to display only the covariates that failed to meet at least one of balance thresholds. Default is \code{FALSE}, so all covariates are displayed.
}
\item{un}{\code{logical}; whether to print statistics for the unadjusted sample as well as for the adjusted sample. Default is \code{FALSE}, so only the statistics for the adjusted sample are displayed.
}
\item{disp.means}{\code{logical}; whether to print the group means in balance output. Default is \code{FALSE}, so means are not displayed.}
\item{disp.sds}{\code{logical}; whether to print the group standard deviations in balance output. Default is \code{FALSE}, so standard deviations are not displayed.}
\item{disp.v.ratio}{\code{logical}; whether to display variance ratios in balance output for binary and multinomial treatments. Default is \code{FALSE}, so variance ratios are not displayed.
}
\item{disp.ks}{\code{logical}; whether to display Kolmogorov-Smirnov (KS) statistics in balance output for binary and multinomial treatments. Default is \code{FALSE}, so KS statistics are not displayed.
}
\item{disp.subclass}{\code{logical}; whether to display balance information for individual subclasses if subclassification is used in conditioning. See \code{\link{bal.tab.subclass}} for details. Default is \code{FALSE}, so individual subclasses are not displayed.
}
\item{disp.bal.tab}{\code{logical}; whether to display the table of balance statistics. Default is \code{TRUE}, so the balance table is displayed.
}

}
\details{
These arguments used to be named arguments to \code{bal.tab}, but were relegated to arguments passed through \code{...} to consolidate documentation and make it easier to add new features and options. For users of previous versions of \pkg{cobalt}, nothing should change as far as how the arguments are used, but new options may be added.

In addition to being able to be specified as arguments, if you find you frequently set a display option to something other than its default, you can set that as a permanent option (for the present R session) using \code{\link{options}} and retreive it using \code{\link{getOption}}. The name of each option is \code{"cobalt_arg.name"}, where \code{arg.name} is replaced by the name of the argument. For example, to set the global option so that \code{un = TRUE} unless otherwise specified (which is opposite from the default), use \code{options(cobalt_un = TRUE)}. To determine whether the global option for \code{un} has been set and what its value is, use \code{getOption("cobalt_un")}.
}
\seealso{
\link{bal.tab}
}