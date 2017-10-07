\name{bal.tab.multi}
\alias{bal.tab.multi}

\title{Using \code{bal.tab()} with Multiple Categorical Treatments}

\description{
When using \code{bal.tab()} with multiple categorical treatments, the output will be different from the case with binary or continuous treatments, and there are some options that are common across all \code{bal.tab} methods. This page outlines the outputs and options in this case.

There are two main componenets of the output of \code{bal.tab()} with multiple categorical treatments: the pairwise treatment comparisons and the balance summary. The pairwise treatment comparisons essentially are standard binary treatment balance displays for each possible pair of treatments. So, for example, with treatments "A", "B", and "C", the pairs "A vs. B", "A vs. C", and "B vs. C" comprise the pairs that will be examined. 

The balance summary is, for each variable, the greatest imbalance across all treatment comparisons. So, for variable X1, if "A vs. B" had a standardized mean difference of 0.52, "A vs. C" had a standardized mean difference of .17,  and "B vs. C" had a standardized mean difference of .35, the balance summary would have 0.52 for the value of the standardized mean difference for X1. The same goes for other variables and other measures of balance. If the greatest observed imbalance is tolerable, then all other imabalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall.
}
\arguments{
There are three arguments for each \code{bal.tab()} method that can handle multiple categorical treatments: \code{focal}, \code{which.treat}, and \code{multi.summary}.
\item{focal}{When one group is to be compared to multiple control groups in an ATT analysis, the group considered "treated" is the focal group. Only comparisons between other groups and the focal group are of interest. By specifying the name or index of the treatment condition considered focal, \code{bal.tab()} will only compute and display pairwise balance for treatment comparisons that include the focal group. For example, if "A" was to be compared with "B" and "C", "A"" would be considered focal, and the comparison between groups "B" and "C" would not be computed. In general it is only a good idea to specify \code{focal} if an only if the ATT is saught.}
\item{which.treat}{This is a display option that does not affect computation. When displaying the \code{bal.tab} output, which treatments should be displayed? If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatmenets that both appear in the input will be displayed. For example, inputting \code{"A"} will display "A vs. B" and "A vs. C", while entering \code{c("A", "B")} will only display "A vs. B". \code{NA} indicates no treatment comparisons will be displayed, and \code{NULL} indicates all treatment comparisons will be displayed. \code{NA} is the default.}
\item{multi.summary}{This is a display option that does not affect computation. If \code{TRUE}, the balance summary across all comparisons will be displayed. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. The default is \code{TRUE}, and if \code{which.treat} is \code{NA}, it will automatically be set to \code{TRUE}.}
}
\value{
The output is a \code{bal.tab.multi} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Pair.Balance}{For each pair of treatment groups, a regular \code{bal.tab} object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified. If \code{focal} is specified, only the comparisons involving the focal group are computed.}
\item{Balance.Across.Pairs}{The balance summary across pairs. This will include the greatest (i.e., maximum) absulte mean difference for each covariate across all pairwise comparisons computed, and the same for variance ratios and KS statistics if requested. Thresholds can be requested for each balance measure as with binary treatments.}
\item{Observations}{A table of sample sizes or effective sample sizes for each treatment group before and after adjustment.}

As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\link{bal.tab}, \link{bal.tab.data.frame}, \link{print.bal.tab.multi}
}