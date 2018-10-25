\name{options-display}
\alias{options-display}

\title{Options for Displaying \code{bal.tab} Output}

\description{
Several additional arguments can be passed to \code{bal.tab()} that control the display of the output but don't affect computations. These arguments are documented here.

Not all arguments are applicable to all uses of \code{bal.tab()}; for example, \code{disp.subclass}, which controls whether individual subclasses are displayed when subclassification is used, won't do anything when subclassification is not used.

Note that when \code{quick = TRUE} is set in the call to \code{bal.tab}, setting any of these arguments to \code{FALSE} can prevent some values from being computed, which can have unintended effects. For example, setting \code{quick = TRUE} and \code{un = FALSE} in \code{bal.tab} and then using \code{love.plot} on the output will yield an error because \code{love.plot} depends on the computation of the unadjusted balance statistics, which are not computed when \code{quick = TRUE} and \code{un = FALSE}.
}
\arguments{
\item{imbalanced.only}{\code{logical}; whether to display only the covariates that failed to meet at least one of balance thresholds. Default is \code{FALSE}, so all covariates are displayed.
}
\item{un}{\code{logical}; whether to print statistics for the unadjusted sample as well as for the adjusted sample. Default is \code{FALSE}, so only the statistics for the adjusted sample are displayed.
}
\item{disp.means}{\code{logical}; whether to print the group means in balance output. Default is \code{FALSE}, so means are not displayed.
}
\item{disp.sds}{\code{logical}; whether to print the group standard deviations in balance output. Default is \code{FALSE}, so standard deviations are not displayed.
}
\item{disp.v.ratio}{\code{logical}; whether to display variance ratios in balance output for binary and multinomial treatments. Default is \code{FALSE}, so variance ratios are not displayed.
}
\item{disp.ks}{\code{logical}; whether to display Kolmogorov-Smirnov (KS) statistics in balance output for binary and multinomial treatments. Default is \code{FALSE}, so KS statistics are not displayed.
}
\item{disp.bal.tab}{\code{logical}; whether to display the table of balance statistics. Default is \code{TRUE}, so the balance table is displayed.
}
\item{factor_sep}{\code{character}; the string used to separate factor variables from their levels when variable names are printed. Default is \code{"_"}.
}
\item{int_sep}{\code{character}; the string used to separate two variables involved in an interaction when variable names are printed. Default is \code{" * "}. Older versions of \pkg{cobalt} used \code{"_"}.
}
\bold{When subclassification is used}
\item{disp.subclass}{\code{logical}; whether to display balance information for individual subclasses if subclassification is used in conditioning. See \code{\link{bal.tab.subclass}} for details. Default is \code{FALSE}, so individual subclasses are not displayed.
}
\bold{When the treatment is multinomial}
\item{which.treat}{For which treatments or treatment combinations balance tables should be displayed. If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatments that both appear in the input will be displayed. For example, setting  \code{which.treat = "A"} will display "A vs. B" and "A vs. C", while setting \code{which.treat = c("A", "B")} will only display "A vs. B". \code{NA} indicates no treatment comparisons will be displayed, and \code{NULL} indicates all treatment comparisons will be displayed. Default is \code{NA}. See \code{\link{bal.tab.multi}}.
}
\item{multi.summary}{\code{logical}; whether to display the balance summary across all treatment pairs. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. Default is \code{TRUE}, and if \code{which.treat} is \code{NA}, it will automatically be set to \code{TRUE}. See \code{\link{bal.tab.multi}}.
}
\bold{When clusters are present}
\item{which.cluster}{For which clusters balance tables should be displayed. If \code{NULL} (the default), all clusters in \code{cluster} will be displayed. If \code{NA}, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor). See \code{\link{bal.tab.cluster}}.
}
\item{cluster.summary}{\code{logical}; whether to display the balance summary across clusters. Default is \code{TRUE}, and if \code{which.cluster} is \code{NA}, it will automatically be set to \code{TRUE}. See \code{\link{bal.tab.cluster}}.
}
\item{cluster.fun}{Which function is used in the across-cluster summary to combine results across clusters. Can be "min", "mean", or "max". For example, if \code{cluster.fun = "mean"} the mean balance statistic across clusters will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic. See \code{\link{bal.tab.cluster}}.
}
\bold{When multiple imputations are present}
\item{which.imp}{For which imputations balance tables should be displayed. If \code{NULL}, all imputations in \code{imp} will be displayed. If \code{NA} (the default), no imputations will be displayed. Otherwise, can be a vector of imputation indices for which to display balance. See \code{\link{bal.tab.imp}}.
}
\item{imp.summary}{\code{logical}; whether to display the balance summary across imputations. Default is \code{TRUE}, and if \code{which.cluster} is \code{NA}, it will automatically be set to \code{TRUE}. See \code{\link{bal.tab.imp}}.
}
\item{imp.fun}{Which function is used in the across-imputation summary to combine results across imputations. Can be "min", "mean", or "max". For example, if \code{imp.fun = "mean"} the mean balance statistic across imputations will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic. See \code{\link{bal.tab.imp}}.
}
\bold{When the treatment is longitudinal}
\item{which.time}{For which time points balance tables should be displayed. If \code{NULL}, all time points will be displayed. If \code{NA}, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance. Default is \code{NULL}. See \code{\link{bal.tab.msm}}.}
\item{msm.summary}{\code{logical}; whether to display the balance summary across time points. Default is \code{TRUE}, and if \code{which.time} is \code{NA}, it will automatically be set to \code{TRUE}. See \code{\link{bal.tab.msm}}.
}
}

\details{
These arguments used to be named arguments to \code{bal.tab}, but were relegated to arguments passed through \code{...} to consolidate documentation and make it easier to add new features and options. For users of previous versions of \pkg{cobalt}, nothing should change as far as how the arguments are used, but new options may be added.

In addition to being able to be specified as arguments, if you find you frequently set a display option to something other than its default, you can set that as a permanent option (for the present R session) using \code{\link{options}} and retrieve it using \code{\link{getOption}}. The name of each option is \code{"cobalt_argname"}, where \code{argname} is replaced by the name of the argument. For example, to set the global option so that \code{un = TRUE} unless otherwise specified (which is opposite from the default), use \code{options(cobalt_un = TRUE)}. To determine whether the global option for \code{un} has been set and what its value is, use \code{getOption("cobalt_un")}. Note that global options cannot be set for \code{which.cluster}, \code{which.imp}, \code{which.treat}, or \code{which.time}.
}

\seealso{
\code{\link{bal.tab}}, \code{\link{print.bal.tab}}
}
