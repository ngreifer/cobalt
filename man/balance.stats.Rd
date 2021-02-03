\name{Balance Statistics}
\alias{balance.stats}

\title{Balance Statistics in \code{bal.tab} and \code{love.plot}}

\description{
\fun{bal.tab} and \fun{love.plot} display balance statistics for the included covariates. The \code{stats} argument in each of these functions controls which balance statistics are to be displayed. The argument to \code{stats} should be a character vector with the names of the desired balance statistics. This page describes all of the available balance statistics and how to request them. Abbreviations are allowed, so you can use the first few letters of each balance statistics to request it instead of typing out its whole name. That convention is used throughout the documentation. For example, to request mean differences and variance ratios in \code{bal.tab()} or \code{love.plot()}, you could include \code{stats = c("m", "v")}. In addition, the \code{thresholds} argument uses the same naming conventions and can be used to request balance thresholds on each statistic. For example, to request a balance threshold of .1 for mean differences, you could include \code{thresholds = c(m = .1)}.

Below, each allowable entry to \code{stats} and \code{thresholds} are described, along with other details or option that accompany them.

\subsection{Binary/Multi-Category Treatments}{
\describe{
\item{\code{"mean.diffs"}}{Mean differences as computed by \fun{col_w_smd}. Can be abbreviated as \code{"m"}. Setting the arguments \code{continuous} and \code{binary} to either \code{"std"} or \code{"raw"} will determine whether standardized mean differences or raw mean differences are calculated for continuous and categorical variables, respectively. When standardized mean differences are requested, the \code{s.d.denom} argument controls how the standardization occurs. When \code{abs = TRUE}, negative values become positive. Mean differences are requested by default when no entry to \code{stats} is provided.}

\item{\code{"variance.ratios"}}{Variance ratios as computed by \fun{col_w_vr}. Can be abbreviated as \code{"v"}. Will not be computed for binary variables. When \code{abs = TRUE}, values less than 1 will have their inverse taken. When used with \code{love.plot}, the x-axis scaled will be logged so that, e.g., .5 is as far away from 1 as 2 is.}

\item{\code{"ks.statistics"}}{Kolmogorov-Smirnov (KS) statistics as computed by \fun{col_w_ks}.}

\item{\code{"ovl.coefficients"}}{Overlapping (OVL) statistics as computed by \fun{col_w_ovl}. Can be abbreviated as \code{"ovl"}. Additional arguments passed to \code{col_w_ovl()}, such as \code{integrate} or \code{bw}, can be supplied to \code{bal.tab()} or \code{love.plot()}.}
}
}
\subsection{Continuous Treatments}{
\describe{
\item{\code{"correlations"}}{Pearson correlations as computed by \fun{col_w_cov}. Can be abbreviated as \code{"cor"}. Setting the arguments \code{continuous} and \code{binary} to either \code{"std"} or \code{"raw"} will determine whether correlations or covariances are calculated for continuous and categorical variables, respectively (they are both \code{"std"} by default). When correlations are requested, the \code{s.d.denom} argument controls how the standardization occurs. When \code{abs = TRUE}, negative values become positive. Pearson correlations are requested by default when no entry to \code{stats} is provided.}

\item{\code{"spearman.correlations"}}{Spearman correlations as computed by \fun{col_w_cov}. Can be abbreviated as \code{"sp"}. All arguments are the same as those for \code{"correlations"}. When \code{abs = TRUE}, negative values become positive.}

\item{\code{"mean.diffs.target"}}{Mean differences computed between the weighted and unweighted sample to ensure the weighted sample is representative of the original population. Can be abbreviated as \code{"m"}. Setting the arguments \code{continuous} and \code{binary} to either \code{"std"} or \code{"raw"} will determine whether standardized mean differences or raw mean differences are calculated for continuous and categorical variables, respectively. The standardization factor will be computed in the unweighted sample. When \code{abs = TRUE}, negative values become positive. This statistic is only computed for the adjusted samples.}

\item{\code{"ks.statistics.target"}}{KS-statistics computed between the weighted and unweighted sample to ensure the weighted sample is representative of the original population. Can be abbreviated as \code{"ks"}. This statistic is only computed for the adjusted samples.}
}
}

If a statistic is requested in \code{thresholds}, it will automatically be placed in \code{stats}. For example, \code{bal.tab(..., stats = "m", thresholds = c(v = 2))} will display both mean differences and variance ratios, and the variance ratios will have a balance threshold set to 2.
}
\examples{
data(lalonde)

#Binary treatments
bal.tab(treat ~ age + educ + married + re74, data = lalonde,
        stats = c("m", "v", "ks"))
love.plot(treat ~ age + educ + married + re74, data = lalonde,
          stats = c("m", "v", "ks"), binary = "std",
          thresholds = c(m = .1, v = 2))

#Continuous treatments
bal.tab(re75 ~ age + educ + married + re74, data = lalonde,
        stats = c("cor", "sp"))
love.plot(re75 ~ age + educ + married + re74, data = lalonde,
          thresholds = c(cor = .1, sp = .1))
}
