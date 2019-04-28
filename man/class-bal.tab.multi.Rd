\name{class-bal.tab.multi}
\alias{class-bal.tab.multi}
\alias{bal.tab.multi}

\title{Using \code{bal.tab()} with Multiple Categorical Treatments}

\description{
When using \code{bal.tab()} with multiple categorical treatments, the output will be different from the case with binary or continuous treatments, and there are some options that are common across all \code{bal.tab} methods. This page outlines the outputs and options in this case.

There are two main components of the output of \code{bal.tab()} with multiple categorical treatments: the two-group treatment comparisons and the balance summary. The two-group treatment comparisons are standard binary treatment comparison either for pairs of groups (e.g., for treatments A, B, and C, "A vs. B", "A vs. C", and "B vs. C") or each group against the others.

The balance summary is, for each variable, the greatest imbalance across all two-group comparisons. So, for variable X1, if "A vs. B" had a standardized mean difference of 0.52, "A vs. C" had a standardized mean difference of .17,  and "B vs. C" had a standardized mean difference of .35, the balance summary would have 0.52 for the value of the standardized mean difference for X1. The same goes for other variables and other measures of balance. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. (Note that when \code{s.d.denom = "pooled"}, i.e., when the estimand is the ATE, the pooled standard deviation in the denominator will be the average of the standard deviations across all treatment groups, not just those used in the pairwise comparison.)
}
\arguments{
There are four arguments for each \code{bal.tab()} method that can handle multiple categorical treatments: \code{pairwise}, \code{focal}, \code{which.treat}, and \code{multi.summary}.
\item{pairwise}{Whether to compute the two-group comparisons pairwise or not. If \code{TRUE}, \code{bal.tab()} will compute comparisons for each pair of treatments. This can be valuable if treatments are to be compared with one another (which is often the case). If \code{FALSE}, \code{bal.tab()} will compute balance for each treatment group against all other groups combined. This only makes sense if the ATE is desired. When \code{focal} is specified, \code{pairwise} is automatically set to \code{TRUE}.}
\item{focal}{When one group is to be compared to multiple control groups in an ATT analysis, the group considered "treated" is the focal group. Only comparisons between other groups and the focal group are of interest. By specifying the name or index of the treatment condition considered focal, \code{bal.tab()} will only compute and display pairwise balance for treatment comparisons that include the focal group. For example, if "A" was to be compared with "B" and "C", "A" would be considered focal, and the comparison between groups "B" and "C" would not be computed. In general it is only a good idea to specify \code{focal} if an only if the ATT is sought.}
\item{which.treat}{This is a display option that does not affect computation. When displaying the \code{bal.tab} output, which treatments should be displayed? If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatments that both appear in the input will be displayed. For example, inputting \code{"A"} will display "A vs. B" and "A vs. C", while entering \code{c("A", "B")} will only display "A vs. B". \code{.none} indicates no treatment comparisons will be displayed, and \code{.all} indicates all treatment comparisons will be displayed. \code{.none} is the default.}
\item{multi.summary}{This is a display option that does not affect computation. If \code{TRUE}, the balance summary across all comparisons will be displayed. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. The default is \code{TRUE}, and if \code{which.treat} is \code{.none}, it will automatically be set to \code{TRUE}.}
}
\value{
The output is a \code{bal.tab.multi} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Pair.Balance}{For each pair of treatment groups, a regular \code{bal.tab} object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified. If \code{focal} is specified, only the comparisons involving the focal group are computed. If \code{pairwise} is \code{FALSE}, the comparisons will be between each group and the other groups combined (labeled "Others").}
\item{Balance.Across.Pairs}{The balance summary across two-group comparisons. This will include the greatest (i.e., maximum) absolute mean difference for each covariate across all comparisons computed, and the same for variance ratios and KS statistics if requested. Thresholds can be requested for each balance measure as with binary treatments.}
\item{Observations}{A table of sample sizes or effective sample sizes for each treatment group before and after adjustment.}

As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}}, \code{\link{bal.tab.data.frame}}, \code{\link{print.bal.tab.multi}}
}