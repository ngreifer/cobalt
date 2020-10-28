\name{Balance Summary}
\alias{col_w_mean}
\alias{col_w_sd}
\alias{col_w_smd}
\alias{col_w_vr}
\alias{col_w_ks}
\alias{col_w_ovl}
\alias{col_w_cov}
\alias{col_w_corr}

\title{
Compute Balance Statistics for Covariates
}
\description{
These functions quickly compute balance statistics for the given covariates. These functions are used in \code{\link{bal.tab}}, but they are available for use in programming without having to call \code{bal.tab()} to get them.
\itemize{
\item{\code{col_w_mean} computes the (weighted) means for a set of covariates and weights and is essentially a weighted version of \code{\link{colMeans}}.}
\item{\code{col_w_sd} computes the (weighted) standard deviations for a set of covariates and weights.}
\item{\code{col_w_smd} computes the (weighted) (absolute) (standardized) difference in means for a set of covariates, a binary treatment, and weights.}
\item{\code{col_w_vr} computes the (weighted) variance ratio for a set of covariates, a binary treatment, and weights.}
\item{\code{col_w_ks} computes the (weighted) Kolmogorov-Smirnov (KS) statistic for a set of covariates, a binary treatment, and weights.}
\item{\code{col_w_ovl} computes the complement of the (weighted) overlapping coefficient for a set of covariates, a binary treatment, and weights (based on Franklin et al, 2014).}
\item{\code{col_w_cov} and \code{col_w_corr} compute the (weighted) (absolute) treatment-covariate covariance or correlation for a set of covariates, a continuous treatment, and weights.}
}
}

\usage{
col_w_mean(mat, weights = NULL, s.weights = NULL, 
           subset = NULL, na.rm = TRUE, ...)
           
col_w_sd(mat, weights = NULL, s.weights = NULL, 
         bin.vars, subset = NULL, na.rm = TRUE, ...)
         
col_w_smd(mat, treat, weights = NULL, std = TRUE, 
          s.d.denom = "pooled", abs = FALSE, 
          s.weights = NULL, bin.vars, 
          subset = NULL, weighted.weights = weights, 
          na.rm = TRUE, ...)
          
col_w_vr(mat, treat, weights = NULL, abs = FALSE, 
         s.weights = NULL, bin.vars, 
         subset = NULL, na.rm = TRUE, ...)
         
col_w_ks(mat, treat, weights = NULL, 
         s.weights = NULL, bin.vars, 
         subset = NULL, na.rm = TRUE, ...)
         
col_w_ovl(mat, treat, weights = NULL, 
         s.weights = NULL, bin.vars, 
         integrate = FALSE, subset = NULL, 
         na.rm = TRUE, ...)
         
col_w_cov(mat, treat, weights = NULL, type = "pearson",
          std = FALSE, s.d.denom = "all", abs = FALSE, 
          s.weights = NULL, bin.vars, 
          subset = NULL, weighted.weights = weights, 
          na.rm = TRUE, ...)
           
col_w_corr(mat, treat, weights = NULL, type = "pearson",
           s.d.denom = "all", abs = FALSE, 
           s.weights = NULL, bin.vars, 
           subset = NULL, weighted.weights = weights, 
           na.rm = TRUE, ...)
}

\arguments{
  \item{mat}{
a numeric matrix or a data.frame containing the covariates for which the statistic is to be computed. If a data.frame, \code{\link{splitfactor}} with \code{drop.first = "if2"} will be called if any character or factor variables are present. This can slow down the function, so it's generally best to supply a numeric matrix. If a numeric vector is supplied, it will be converted to a 1-column matrix first.
}
  \item{weights}{
\code{numeric}; an optional set of weights used to compute the weighted statistics. If sampling weights are supplied through \code{s.weights}, the \code{weights} should not incorporate these weights, as \code{weights} and \code{s.weights} will be multiplied together prior to computing the weighted statistics.
}
  \item{s.weights}{
\code{numeric}; an optional set of sampling weights used to compute the weighted statistics. If weights are supplied through \code{weights}, \code{weights} and \code{s.weights} will be multiplied together prior to computing the weighted statistics. Some functions use \code{s.weights} in a particular way; for others, supplying \code{weights} and \code{s.weights} is equivalent to supplying their product to either \code{weights} or \code{s.weights}. See Details.
}
  \item{subset}{
a \code{logical} vector with length equal to the number of rows of \code{mat} used to subset the data. See Details for notes on its use with \code{col_w_smd}, \code{col_w_cov}, and \code{col_w_corr}.
}
  \item{na.rm}{
\code{logical}; whether \code{NA}s should be ignored or not. If \code{FALSE}, any variable with any \code{NA}s will have its corresponding statistic returned as \code{NA}. If \code{TRUE}, any variable with any \code{NA}s will have its corresponding statistic computed as if the missing value were not there.
}
  \item{treat}{
a vector of treatment status for each individual. For \code{col_w_smd}, \code{col_w_vr}, \code{col_w_ks}, and \code{col_w_ovl}, \code{treat} should have exactly two unique values. For \code{col_w_cov} and \code{col_w_corr}, \code{treat} should be a many-valued numeric vector.
}
  \item{std}{
\code{logical}; for \code{col_w_smd}, whether the computed mean differences for each variable should be standardized; for \code{col_w_cov}, whether treatment-covariate correlations should be computed (\code{TRUE}) rather than covariances (\code{FALSE}). Can be either length 1, whereby all variables will be standardized or not, or length equal to the number of columns of \code{mat}, whereby only variables with a value of \code{TRUE} will be standardized. See Details.
}
  \item{s.d.denom}{
for \code{col_w_smd} and \code{col_w_cov} when \code{std} is \code{TRUE} for some variables, and for \code{col_w_corr}, how the standardization factor should be computed. For \code{col_w_smd} (i.e., when computing standardized mean differences), allowable options include 
\itemize{
\item{\code{"treated"} - uses the standard deviation of the variable in the treated group}
\item{\code{"control"} - uses the standard deviation of the variable in the control group}
\item{\code{"pooled"} - uses the square root of the average of the variances of the variable in the treated and control groups}
\item{\code{"all"} - uses the standard deviation of the variable in the full sample}
\item{\code{"weighted"} - uses the standard deviation of the variable in the full sample weighted by \code{weighted.weights}}
\item{\code{"hedges"} - uses the small-sample corrected version of Hedge's G described in the WWC Procedures Handbook (see References)}
\item{the name of one of the treatment values - uses the standard deviation of the variable in that treatment group.}
}
For \code{col_w_cov} and \code{col_w_corr}, only "all" and "weighted" are allowed. Abbreviations allowed. This can also be supplied as a numeric vector of standard deviations with length equal to the number of columns of \code{mat}; the values will be used as the standardization factors.
}
  \item{abs}{
\code{logical}; for \code{col_w_smd}, \code{col_w_cov}, and \code{col_w_corr}, whether the returned statistics should be in absolute value (\code{TRUE}) or not. For \code{col_w_vr}, whether the ratio should always include the larger variance in the numerator, so that the ratio is always greater than or equal to 1. Default is \code{FALSE}.
}
  \item{bin.vars}{
a vector used to denote whether each variable is binary or not. Can be a \code{logical} vector with length equal to the number of columns of \code{mat} or a vector of numeric indices or character names of the binary variables. If missing (the default), the function will figure out which covariates are binary or not, which can increase computation time. If \code{NULL}, it will be assumed no variables are binary. All functions other than \code{col_w_mean} treat binary variables different from continuous variables. If a factor or character variable is in \code{mat}, all the dummies created will automatically be marked as binary, but it should still receive an entry when \code{bin.vars} is supplied as \code{logical}.
}
  \item{weighted.weights}{
for \code{col_w_smd}, \code{col_w_cov}, and \code{col_w_corr}, when \code{std = TRUE} and \code{s.d.denom = "weighted"}, a vector of weights to be applied to the computation of the denominator standard deviation. If not specified, will use the argument to \code{weights}. When \code{s.d.denom} is not "weighted", this is ignored. The main purpose of this is to allow \code{weights} to be \code{NULL} while weighting the denominator standard deviations for assessing balance in the unweighted sample but using the standard deviations of the weighted sample.
}
  \item{type}{
for \code{col_w_cov} and \code{col_w_corr}, the type of covariance/correlation to be computed. Allowable options include "pearson" and "spearman". When "spearman" is requested, the covariates and treatment are first turned into ranks using \code{\link{rank}} with \code{na.last = "keep"}.
}
  \item{integrate}{
\code{logical}; for \code{col_w_ovl}, whether to use \code{\link{integrate}} to calculate the area of overlap. If \code{FALSE}, a midpoint Riemann sum with 1000 partitions will be used instead. The Riemann sum is a little slower and very slightly imprecise (unnoticibly in most contexts), but the integral can fail sometimes and thus is less stable. The default is to use the Riemann sum.
}
  \item{...}{
for all functions, additional arguments supplied to \code{\link{splitfactor}} when \code{mat} is a data.frame. \code{data}, \code{var.name}, \code{drop.first}, and \code{drop.level} are ignored; \code{drop.first} is automatically set to \code{"if2"}. For \code{col_w_ovl}, other arguments passed to \code{\link{density}} besides \code{x} and \code{weights}. Note that the default value for \code{bw} when unspecified is "nrd" rather than the default in \code{density}, which is "nrd0".
}
}
\details{
\code{col_w_mean} computes column weighted means for a matrix of variables. It is similar to \code{colMeans} but (optionally) incorporates weights. \code{weights} and \code{s.weights} are multiplied together prior to being used, and there is no distinction between them. This could be used to compute the weighted means of each covariate in the general population to examine the degree to which a weighting method has left the weighted samples resembling the original population.

\code{col_w_sd} computes column weighted standard deviations for a matrix of variables. \code{weights} and \code{s.weights} are multiplied together prior to being used, and there is no distinction between them. The variance of binary variables is computed as \eqn{p(1-p)}, where \code{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula; the standard deviation is the square root of this variance.

\code{col_w_smd} computes the mean difference for each covariate between treatment groups defined by \code{treat}. These mean differences can be optionally be weighted, standardized, and/or in absolute value. The standardization factor is computed using the unweighted standard deviation or variance when \code{s.weights} are absent, and is computed using the sampling weighted standard deviation or variance when \code{s.weights} are present, except when \code{s.d.denom = "weighted"}, in which case the product of \code{weighted.weights} and \code{s.weights} (if present) are used to weight the standardization factor. The standardization factor is computed using the whole sample even when \code{subset} is used. Note that unlike \code{bal.tab}, \code{col_w_smd} requires the user to specify whether each individual variable should be standardized using \code{std} rather than relying on \code{continuous} or \code{binary}. The mean difference is computed using the product of \code{weights} and \code{s.weights}, if specified. The variance of binary variables is computed as \eqn{p(1-p)}, where \code{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. 

\code{col_w_vr} computes the variance ratio for each covariate between treatment groups defined by \code{treat}. When \code{abs = TRUE}, \code{max(out, 1/out)} is applied to the output so that the ratio is always greater than or equal to 1. For binary variables, the variance is computed as \eqn{p(1-p)}, where \code{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. Note that in \code{bal.tab}, variance ratios are not computed for binary variables, while here, they are (but likely should not be interpreted). \code{weights} and \code{s.weights} are multiplied together prior to being used, and there is no distinction between them.

\code{col_w_ks} computes the KS statistic for each covariate using the method implemented in \pkg{twang}. The KS statistics can optionally be weighted. For binary variables, the KS statistic is just the difference in proportions. \code{weights} and \code{s.weights} are multiplied together prior to being used, and there is no distinction between them.

\code{col_w_ovl} computes the complement of the overlapping coefficient as described by Franklin et al. (2014). It does so by computing the density of the covariate in the treated and control groups, then finding the area where those density overlap, and subtracting that number from 1, yielding a value between 0 and 1 where 1 indicates complete imbalance, and 0 indicates perfect balance. \code{\link{density}} is used to model the density in each group. The bandwidth of the covariate in the smaller treatment group is used for both groups. The area of overlap can be computed using \code{integrate}, which quickly and accurately computes the integral, or using a midpoint Riemann sum with 1000 partitions, which approximates the area more slowly. A reason to prefer the Riemann sum is that \code{integrate} can fail for unknown reasons, though Riemann sums will fail with some extreme distributions. When either method fails, the resulting value will be \code{NA}. For binary variables, the complement of the overlapping coefficient is just the difference in proportions. \code{weights} and \code{s.weights} are multiplied together prior to being used, and there is no distinction between them. The weights are used to compute the weighted density by supplying them to the \code{weights} argument of \code{density}.

\code{col_w_cov} computes the covariance between a continuous treatment and the covariates to assess balance for continuous treatments as recommended in Austin (2019). These covariance can optionally be weighted or in absolute value or can be requested as correlations (i.e., standardized covariances). The correlations are computed as the covariance between the treatment and covariate divided by a standardization factor, which is equal to the square root of the product of the variance of treatment and the variance of the covariate. The standardization factor is computed using the unweighted variances when \code{s.weights} are absent, and is computed using the sampling weighted variances when \code{s.weights} are present, except when \code{s.d.denom = "weighted"}, in which case the product of \code{weighted.weights} and \code{s.weights} (if present) are used to weight the standardization factor. For this reason, the computed correlation can be greater than 1 or less than -1. The standardization factor is always computed using the whole sample even when \code{subset} is used. The covariance is computed using the product of \code{weights} and \code{s.weights}, if specified. The variance of binary variables is computed as \eqn{p(1-p)}, where \code{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. 

\code{col_w_corr} is a wrapper for \code{col_w_cov} with \code{std} set to \code{TRUE}. 
}
\value{
A vector of balance statistics, one for each variable in \code{mat}. If \code{mat} has column names, the output will be named as well.
}
\references{
Franklin, J. M., Rassen, J. A., Ackermann, D., Bartels, D. B., & Schneeweiss, S. (2014). Metrics for covariate balance in cohort studies of causal effects. \emph{Statistics in Medicine}, 33(10), 1685–1699. \doi{10.1002/sim.6058}

Austin, P. C. (2019). Assessing covariate balance when using the generalized propensity score with quantitative or continuous exposures. \emph{Statistical Methods in Medical Research}, 28(5), 1365–1377. \doi{10.1177/0962280218756159}

What Works Clearinghouse. (2020). WWC Procedures Handbook (Version 4.1). Retrieved from
\url{https://ies.ed.gov/ncee/wwc/Handbooks}

}

\seealso{
\code{\link{bal.tab}}
}
\examples{\dontshow{if (requireNamespace("WeightIt", quietly = TRUE)) \{}
library(WeightIt); data("lalonde", package = "cobalt")

treat <- lalonde$treat
covs <- subset(lalonde, select = -c(treat, re78))
covs <- splitfactor(covs, drop.first = "if2")
bin.vars <- c(FALSE, FALSE, TRUE, TRUE, TRUE,
              TRUE, TRUE, FALSE, FALSE)
W <- weightit(treat ~ covs, method = "ps", 
              estimand = "ATE", tols = .01)
weights <- W$weights

round(data.frame(
    m0 = col_w_mean(covs, weights = weights, subset = treat == 0),
    sd0 = col_w_sd(covs, weights = weights,
                   bin.vars = bin.vars, subset = treat == 0),
    m1 = col_w_mean(covs, weights = weights, subset = treat == 1),
    sd1 = col_w_sd(covs, weights = weights,
                   bin.vars = bin.vars, subset = treat == 1),
    smd = col_w_smd(covs, treat = treat, weights = weights,
                     std = TRUE, bin.vars = bin.vars),
    vr = col_w_vr(covs, treat = treat, weights = weights,
                  bin.vars = bin.vars),
    ks = col_w_ks(covs, treat = treat, weights = weights,
                  bin.vars = bin.vars),
    row.names = colnames(covs)
), 4)

# Compare to bal.tab():
bal.tab(covs, treat, weights = weights, disp.means = TRUE,
        disp.sds = TRUE, disp.v.ratio = TRUE, disp.ks = TRUE, 
        estimand = "ATE", method = "weighting", binary = "std")
\dontshow{\}}}
\keyword{support functions}
