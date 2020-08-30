\name{Display Options}
\alias{options-display}
\alias{display_options}

\title{Options for Displaying \code{bal.tab} Output}

\description{
Several additional arguments can be passed to \code{bal.tab()} that control the display of the output but don't affect computations. These arguments are documented here.

Not all arguments are applicable to all uses of \code{bal.tab()}; for example, \code{disp.subclass}, which controls whether individual subclasses are displayed when subclassification is used, won't do anything when subclassification is not used.

Note that when \code{quick = TRUE} is set in the call to \code{bal.tab}, setting any of these arguments to \code{FALSE} can prevent some values from being computed, which can have unintended effects.
}
\arguments{
\item{disp.bal.tab}{\code{logical}; whether to display the table of balance statistics. Default is \code{TRUE}, so the balance table is displayed.
}
\item{imbalanced.only}{\code{logical}; whether to display only the covariates that failed to meet at least one of balance thresholds. Default is \code{FALSE}, so all covariates are displayed.
}
\item{un}{\code{logical}; whether to print statistics for the unadjusted sample as well as for the adjusted sample. Default is \code{FALSE}, so only the statistics for the adjusted sample are displayed.
}
\item{disp}{\code{character}; which distribution summary statistic(s) should be reported. Allowable options include \code{"means"} and \code{"sds"}. Multiple options are allowed. Abbreviations allowed.
}
\item{stats}{\code{character}; which statistic(s) should be reported. See \code{\link[=balance.stats]{stats}} to see which options are available. Multiple options are allowed. Abbreviations allowed. For binary and multi-category treatments, the default is \code{"mean.diffs"} (i.e., [standardized] mean differences), and for continuous treatments, the default is \code{"correlations"} (i.e., treatment-covariate Pearson correlations).
}
\item{factor_sep}{\code{character}; the string used to separate factor variables from their levels when variable names are printed. Default is \code{"_"}.
}
\item{int_sep}{\code{character}; the string used to separate two variables involved in an interaction when variable names are printed. Default is \code{" * "}. Older versions of \pkg{cobalt} used \code{"_"}.
}
\bold{When subclassification is used}
\item{disp.subclass}{\code{logical}; whether to display balance information for individual subclasses if subclassification is used in conditioning. See \code{\link[=class-bal.tab.subclass]{bal.tab.subclass}} for details. Default is \code{FALSE}, so individual subclasses are not displayed.
}
\bold{When the treatment is multi-category}
\item{which.treat}{For which treatments or treatment combinations balance tables should be displayed. If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatments that both appear in the input will be displayed. For example, setting  \code{which.treat = "A"} will display "A vs. B" and "A vs. C", while setting \code{which.treat = c("A", "B")} will only display "A vs. B". \code{.none} indicates no treatment comparisons will be displayed, and \code{.all} indicates all treatment comparisons will be displayed. Default is \code{.none}. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}}.
}
\item{multi.summary}{\code{logical}; whether to display the balance summary across all treatment pairs. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. Default is \code{TRUE} when \code{which.treat} is \code{.none} and \code{FALSE} otherwise. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}}.
}
\bold{When clusters are present}
\item{which.cluster}{For which clusters balance tables should be displayed. If \code{.all}, all clusters in \code{cluster} will be displayed. If \code{.none}, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor). Default is \code{.all}. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}}.
}
\item{cluster.summary}{\code{logical}; whether to display the balance summary across clusters. Default is \code{TRUE} when \code{which.cluster} is \code{.none} and \code{FALSE} otherwise (note the default for \code{which.cluster} is \code{.all}). See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}}.
}
\item{cluster.fun}{Which function is used in the across-cluster summary to combine results across clusters. Can be "min", "mean", or "max". For example, if \code{cluster.fun = "mean"} the mean balance statistic across clusters will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}}.
}
\bold{When multiple imputations are present}
\item{which.imp}{For which imputations balance tables should be displayed. If \code{.all}, all imputations in \code{imp} will be displayed. If \code{.none}, no imputations will be displayed. Otherwise, can be a vector of imputation indices for which to display balance. Default is \code{.none}. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}}.
}
\item{imp.summary}{\code{logical}; whether to display the balance summary across imputations. Default is \code{TRUE} when \code{which.imp} is \code{.none} and \code{FALSE} otherwise. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}}.
}
\item{imp.fun}{Which function is used in the across-imputation summary to combine results across imputations. Can be "min", "mean", or "max". For example, if \code{imp.fun = "mean"} the mean balance statistic across imputations will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}}.
}
\bold{When the treatment is longitudinal}
\item{which.time}{For which time points balance tables should be displayed. If \code{.all}, all time points will be displayed. If \code{.none}, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance. Default is \code{.none}. See \code{\link[=class-bal.tab.msm]{bal.tab.msm}}.}
\item{msm.summary}{\code{logical}; whether to display the balance summary across time points. Default is \code{TRUE} when \code{which.time} is \code{.none} and \code{FALSE} otherwise. See \code{\link[=class-bal.tab.msm]{bal.tab.msm}}.
}

\emph{Deprecated arguments}

The following arguments are deprecated but still work.

\item{disp.means}{\code{logical}; whether to print the group means in balance output. Default is \code{FALSE}, so means are not displayed. Deprecated; use \code{disp} instead.
}
\item{disp.sds}{\code{logical}; whether to print the group standard deviations in balance output. Default is \code{FALSE}, so standard deviations are not displayed. Deprecated; use \code{disp} instead.
}
\item{disp.diff}{\code{logical}; whether to display (standardized) mean differences in balance output for binary and multi-category treatments. Default is \code{TRUE}, so mean differences are displayed. Deprecated; use \code{stats} instead.
}
\item{disp.v.ratio}{\code{logical}; whether to display variance ratios in balance output for binary and multi-category treatments. Default is \code{FALSE}, so variance ratios are not displayed. Deprecated; use \code{stats} instead.
}
\item{disp.ks}{\code{logical}; whether to display Kolmogorov-Smirnov (KS) statistics in balance output for binary and multi-category treatments. Default is \code{FALSE}, so KS statistics are not displayed. Deprecated; use \code{stats} instead.
}
\item{disp.ovl}{\code{logical}; whether to display overlapping (OVL) coefficients in balance output for binary and multi-category treatments. Default is \code{FALSE}, so OVL coefficients are not displayed. Deprecated; use \code{stats} instead.
}
}

\details{
In addition to being able to be specified as arguments, if you find you frequently set a display option to something other than its default, you can set that as a global option (for the present R session) using \code{\link{set.cobalt.options}} and retrieve it using \code{\link{get.cobalt.options}}. Note that global options cannot be set for \code{which.cluster}, \code{which.imp}, \code{which.treat}, or \code{which.time}.
}
\note{
When calling \code{bal.tab()} using \code{\link{do.call}}, if you are using \code{.all} or \code{.none} as inputs to arguments, you need to use \code{\link{alist}} rather than \code{list} to group the arguments. For example, \code{do.call(bal.tab, list(., which.cluster = .none))} will produce an error, but \code{do.call(bal.tab, alist(., which.cluster = .none))} should work correctly.
}

\seealso{
\code{\link{bal.tab}}, \code{\link{print.bal.tab}}
}
