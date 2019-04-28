\name{class-bal.tab.cluster}
\alias{class-bal.tab.cluster}
\alias{bal.tab.cluster}

\title{Using \code{bal.tab()} with Clustered Data}

\description{
When using \code{bal.tab()} with clustered data, the output will be different from the case with single-level data, and there are some options that are common across all \code{bal.tab} methods. This page outlines the outputs and options in this case.

There are two main components of the output of \code{bal.tab()} with clustered data: the within-cluster balance summaries and the across-cluster balance summary. The within-cluster balance summaries display balance for units within each cluster separately.

The across-cluster balance summary pools information across the within-cluster balance summaries to simplify balance assessment. It provides a combination (e.g., mean or maximum) of each balance statistic for each covariate across all clusters. This allows you to see how bad the worst imbalance is and what balance looks like on average.
}
\arguments{
There are four arguments for each \code{bal.tab()} method that can handle clustered data: \code{cluster}, \code{which.cluster}, \code{cluster.summary}, and \code{cluster.fun}.
\item{cluster}{A vector of cluster membership. This can be factor, character, or numeric vector. This argument is required to let \code{bal.tab()} know that the data is clustered. If a \code{data} argument is specified, this can also be the name of a variable in \code{data} that contains cluster membership.}
\item{which.cluster}{This is a display option that does not affect computation. If \code{.all} (the default), all clusters in \code{cluster} will be displayed. If \code{.none}, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor).}
\item{cluster.summary}{This is a display option that does not affect computation. If \code{TRUE}, the balance summary across clusters will be displayed. The default is \code{TRUE}, and if \code{which.cluster} is \code{.none}, it will automatically be set to \code{TRUE}.}
\item{cluster.fun}{This is a display option that does not affect computation. Can be "min", "mean", or "max" and corresponds to which function is used in the across-cluster summary to combine results across clusters. For example, if \code{cluster.fun = "mean"} the mean balance statistic across clusters will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic.
}
}
\value{
The output is a \code{bal.tab.cluster} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Cluster.Balance}{For each cluster, a regular \code{bal.tab} object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.}
\item{Cluster.Summary}{The balance summary across clusters. This will include the combination of absolute mean differences for each covariate across all clusters according to the value of \code{cluster.fun}, and the same for variance ratios and KS statistics if requested.}
\item{Observations}{A table of sample sizes or effective sample sizes for each cluster before and after adjustment.}

As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
}
\seealso{
\code{\link{bal.tab}}, \code{\link{bal.tab.data.frame}}, \code{\link{print.bal.tab.cluster}}
}