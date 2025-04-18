#' @title Display Balance Statistics in a Table
#' 
#' @description Generates balance statistics on covariates in relation to an observed treatment variable. It is a generic function that dispatches to the method corresponding to the class of the first argument.
#' 
#' @usage
#' bal.tab(x, ...)
#' 
#' ## # Arguments common across all input types:
#' ## bal.tab(x,
#' ##         stats,
#' ##         int = FALSE,
#' ##         poly = 1,
#' ##         distance = NULL,
#' ##         addl = NULL,
#' ##         data = NULL,
#' ##         continuous,
#' ##         binary,
#' ##         s.d.denom,
#' ##         thresholds = NULL,
#' ##         weights = NULL,
#' ##         cluster = NULL,
#' ##         imp = NULL,
#' ##         pairwise = TRUE,
#' ##         s.weights = NULL,
#' ##         abs = FALSE,
#' ##         subset = NULL,
#' ##         quick = TRUE,
#' ##         ...)
#' 
#' @param x an input object on which to assess balance. Can be the output of a call to a balancing function in another package or a formula or data frame. Input to this argument will determine which `bal.tab()` method is used. Each input type has its own documentation page, which is linked in the See Also section below. Some input types require or allow additional arguments to be specified. For inputs with no dedicated method, the [default method][bal.tab.default] will be dispatched. See Details below.
#' @param stats `character`; which statistic(s) should be reported. See [`stats`][balance-statistics] for allowable options. For binary and multi-category treatments, `"mean.diffs"` (i.e., mean differences) is the default. For continuous treatments, `"correlations"` (i.e., treatment-covariate Pearson correlations) is the default. Multiple options are allowed.
#' @param int `logical` or `numeric`; whether or not to include 2-way interactions of covariates included in `covs` and in `addl`. If `numeric`, will be passed to `poly` as well.
#' @param poly `numeric`; the highest polynomial of each continuous covariate to display. For example, if 2, squares of each continuous covariate will be displayed (in addition to the covariate itself); if 3, squares and cubes of each continuous covariate will be displayed, etc. If 1, the default, only the base covariate will be displayed. If `int` is numeric, `poly` will take on the value of `int`.
#' @param distance an optional formula or data frame containing distance values (e.g., propensity scores) or a character vector containing their names. If a formula or variable names are specified, `bal.tab()` will look in the argument to `data`, if specified. For longitudinal treatments, can be a list of allowable arguments, one for each time point.
#' @param addl an optional formula or data frame containing additional covariates for which to present balance or a character vector containing their names. If a formula or variable names are specified, `bal.tab()` will look in the arguments to the input object, `covs`, and `data`, if specified. For longitudinal treatments, can be a list of allowable arguments, one for each time point.
#' @param data an optional data frame containing variables named in other arguments. For some input object types, this is required.
#' @param continuous whether mean differences for continuous variables should be standardized (`"std"`) or raw (`"raw"`). Default `"std"`. Abbreviations allowed. This option can be set globally using [set.cobalt.options()].
#' @param binary whether mean differences for binary variables (i.e., difference in proportion) should be standardized (`"std"`) or raw (`"raw"`). Default `"raw"`. Abbreviations allowed. This option can be set globally using [set.cobalt.options()].
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. If weights are supplied, each set of weights should have a corresponding entry to `s.d.denom`. Abbreviations allowed. If left blank and weights, subclasses, or matching strata are supplied, `bal.tab()` will figure out which one is best based on the `estimand`, if given (for ATT, `"treated"`; for ATC, `"control"`; otherwise `"pooled"`) and other clues if not.
#' @param thresholds a named vector of balance thresholds, where the name corresponds to the statistic (i.e., in `stats`) that the threshold applies to. For example, to request thresholds on mean differences and variance ratios, one can set `thresholds = c(m = .05, v = 2)`. Requesting a threshold automatically requests the display of that statistic. When specified, extra columns are inserted into the Balance table describing whether the requested balance statistics exceeded the threshold or not. Summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure are added to the output.
#' @param weights a vector, list, or `data.frame` containing weights for each unit, or a string containing the names of the weights variables in `data`, or an object with a [get.w()] method or a list thereof. The weights can be, e.g., inverse probability weights or matching weights resulting from a matching algorithm.
#' @param cluster either a vector containing cluster membership for each unit or a string containing the name of the cluster membership variable in `data` or the input object. See [`class-bal.tab.cluster`] for details.
#' @param imp either a vector containing imputation indices for each unit or a string containing the name of the imputation index variable in `data` or the input object. See [`class-bal.tab.imp`] for details. Not necessary if `data` is a `mids` object.
#' @param pairwise whether balance should be computed for pairs of treatments or for each treatment against all groups combined. See [`bal.tab.multi()`][class-bal.tab.multi] for details. This can also be used with a binary treatment to assess balance with respect to the full sample.
#' @param s.weights Optional; either a vector containing sampling weights for each unit or a string containing the name of the sampling weight variable in `data`. These function like regular weights except that both the adjusted and unadjusted samples will be weighted according to these weights if weights are used.
#' @param abs `logical`; whether displayed balance statistics should be in absolute value or not. 
#' @param subset a `logical` or `numeric` vector denoting whether each observation should be included or which observations should be included. If `logical`, it should have length equal to the number of units. `NA`s will be treated as `FALSE`. This can be used as an alternative to `cluster` to examine balance on subsets of the data.
#' @param quick `logical`; if `TRUE`, will not compute any values that will not be displayed. Set to `FALSE` if computed values not displayed will be used later.
#' @param ... for some input types, other arguments that are required or allowed. Otherwise, further arguments to control display of output. See [display options][display-options] for details.
#' 
#' @returns
#' An object of class `"bal.tab"`. The use of continuous treatments, subclasses, clusters, and/or imputations will also cause the object to inherit other classes. The class `"bal.tab"` has its own `print()` method ([print.bal.tab()]), which formats the output nicely and in accordance with print-related options given in the call to `bal.tab()`, and which can be called with its own options.
#' 
#' For scenarios with binary point treatments and no subclasses, imputations, or clusters, the following are the elements of the `bal.tab` object:
#'     
#' \item{Balance}{A data frame containing balance information for each covariate. Balance contains the following columns, with additional columns present when other balance statistics are requested, and some columns omitted when not requested:
#' \itemize{
#' \item{`Type`: Whether the covariate is binary, continuous, or a measure of distance (e.g., the propensity score).}
#' \item{`M.0.Un`: The mean of the control group prior to adjusting.}
#' \item{`SD.0.Un`: The standard deviation of the control group prior to adjusting.}
#' \item{`M.1.Un`: The mean of the treated group prior to adjusting.}
#' \item{`SD.1.Un`: The standard deviation of the treated group prior to adjusting.}
#' \item{`Diff.Un`: The (standardized) difference in means between the two groups prior to adjusting. See the `binary` and `continuous` arguments on the `bal.tab` method pages to determine whether standardized or raw mean differences are being reported. By default, the standardized mean difference is displayed for continuous variables and the raw mean difference (difference in proportion) is displayed for binary variables.}
#' \item{`M.0.Adj`: The mean of the control group after adjusting.}
#' \item{`SD.0.Adj`: The standard deviation of the control group after adjusting.}
#' \item{`M.1.Adj`: The mean of the treated group after adjusting.}
#' \item{`SD.1.Adj`: The standard deviation of the treated group after adjusting.}
#' \item{`Diff.Adj`: The (standardized) difference in means between the two groups after adjusting. See the `binary` and `continuous` arguments on the `bal.tab` method pages to determine whether standardized or raw mean differences are being reported. By default, the standardized mean difference is displayed for continuous variables and the raw mean difference (difference in proportion) is displayed for binary variables.}
#' \item{`M.Threshold`: Whether or not the calculated mean difference after adjusting exceeds or is within the threshold given by `thresholds`.  If a threshold for mean differences is not specified, this column will be `NA`.}
#' }}
#' \item{Balanced.Means}{If a threshold on mean differences is specified, a table tallying the number of variables that exceed or are within the threshold.}
#' \item{Max.Imbalance.Means}{If a threshold on mean differences is specified, a table displaying the variable with the greatest absolute mean difference.}
#' \item{Observations}{A table displaying the sample sizes before and after adjusting. Often the effective sample size (ESS) will be displayed. See Details.}
#' \item{call}{The original function call, if adjustment was performed by a function in another package.}
#' 
#' If the treatment is continuous, instead of producing mean differences, `bal.tab()` will produce correlations between the covariates and the treatment. The default corresponding entries in the output will be `"Corr.Un"`, "`Corr.Adj"`, and `"R.Threshold"` (and accordingly for the balance tally and maximum imbalance tables).
#' 
#' If multiple weights are supplied, `"Adj"` in `Balance` will be replaced by the provided names of the sets of weights, and extra columns will be added for each set of weights. Additional columns and rows for other items in the output will be created as well.
#' 
#' For `bal.tab` output with subclassification, see [`class-bal.tab.subclass`].
#' 
#' @details
#' `bal.tab()` performs various calculations on the the data objects given. This page details the arguments and calculations that are used across `bal.tab()` methods.
#' 
#' ### With Binary Point Treatments
#'     
#' Balance statistics can be requested with the [`stats`][balance-statistics] argument. The default balance statistic for mean differences for continuous variables is the standardized mean difference, which is the difference in the means divided by a measure of spread (i.e., a d-type effect size measure). This is the default because it puts the mean differences on the same scale for comparison with each other and with a given threshold. For binary variables, the default balance statistic is the raw difference in proportion. Although standardized differences in proportion can be computed, raw differences in proportion for binary variables are already on the same scale, and computing the standardized difference in proportion can obscure the true difference in proportion by dividing the difference in proportion by a number that is itself a function of the observed proportions.
#'     
#' Standardized mean differences are calculated using [col_w_smd()] as follows: the numerator is the mean of the treated group minus the mean of the control group, and the denominator is a measure of spread calculated in accordance with the argument to `s.d.denom` or the default of the specific method used. Common approaches in the literature include using the standard deviation of the treated group or using the "pooled" standard deviation (i.e., the square root of the mean of the group variances) in calculating standardized mean differences. The computed spread `bal.tab()` uses is always that of the full, unadjusted sample (i.e., before matching, weighting, or subclassification), as recommended by Stuart (2010).
#'     
#' Prior to computation, all variables are checked for variable type, which allows users to differentiate balance statistic calculations based on type using the arguments to `continuous` and `binary`. First, if a given covariate is numeric and has only 2 levels, it is converted into a binary (0,1) variable. If 0 is a value in the original variable, it retains its value and the other value is converted to 1; otherwise, the lower value is converted to 0 and the other to 1. Next, if the covariate is not numeric or logical (i.e., is a character or factor variable), it will be split into new binary variables, named with the original variable and the value, separated by an underscore. Otherwise, the covariate will be used as is and treated as a continuous variable.
#'     
#' When weighting or matching are used, an "effective sample size" is calculated for each group using the following formula: \eqn{(\sum w)^2 / \sum w^2}. The effective sample size is "approximately the number of observations from a simple random sample that yields an estimate with sampling variation equal to the sampling variation obtained with the weighted comparison observations" (Ridgeway et al., 2016). The calculated number tends to underestimate the true effective sample size of the weighted samples. The number depends on the variability of the weights, so sometimes trimming units with large weights can actually increase the effective sample size, even though units are being down-weighted. When matching is used, an additional "unweighted" sample size will be displayed indicating the total number of units contributing to the weighted sample.
#'     
#' When subclassification is used, the balance tables for each subclass stored in `$Subclass.Balance` use values calculated as described above. For the aggregate balance table stored in `$Balance.Across.Subclass`, the values of each statistic are computed as a weighted average of the statistic across subclasses, weighted by the proportion of units in each subclass. See [`class-bal.tab.subclass`] for more details.
#' 
#' ### With Continuous Point Treatments
#'     
#' When continuous treatment variables are considered, the balance statistic calculated is the Pearson correlation between the covariate and treatment. The correlation after adjustment is computed using [col_w_cov()] as the weighted covariance between the covariate and treatment divided by the product of the standard deviations of the unweighted covariate and treatment, in an analogous way to how how the weighted standardized mean difference uses an unweighted measure of spread in its denominator, with the purpose of avoiding the analogous paradox (i.e., where the covariance decreases but is accompanied by a change in the standard deviations, thereby distorting the actual resulting balance computed using the weighted standard deviations). This can sometimes yield correlations greater than 1 in absolute value; these usually indicate degenerate cases anyway.
#' 
#' ### With Multi-Category Point Treatments
#'     
#' For information on using `bal.tab()` with multi-category treatments, see [`class-bal.tab.multi`]. Essentially, `bal.tab()` compares pairs of treatment groups in a standard way.
#' 
#' ### With Longitudinal Treatments
#'     
#' For information on using `bal.tab()` with longitudinal treatments, see [`class-bal.tab.msm`] and `vignette("longitudinal-treat")`. Essentially, `bal.tab()` summarizes balance at each time point and summarizes across time points.
#' 
#' ### With Clustered or Multiply Imputed Data
#'     
#' For information on using `bal.tab()` with clustered data, see [`class-bal.tab.cluster`]. For information on using `bal.tab()` with multiply imputed data, see [`class-bal.tab.imp`]. 
#' 
#' ### `quick`
#'     
#' Calculations can take some time, especially when there are many variables, interactions, or clusters. When certain values are not printed, by default they are not computed. In particular, summary tables are not computed when their display has not been requested. This can speed up the overall production of the output when these values are not to be used later. However, when they are to be used later, such as when output is to be further examined with `print()` or is to be used in some other way after the original call to `bal.tab()`, it may be useful to compute them even if they are not to be printed initially. To do so, users can set `quick = FALSE`, which will cause `bal.tab()` to calculate all values and components it can. Note that `love.plot()` is fully functional even when `quick = TRUE` and values are requested that are otherwise not computed in `bal.tab()` with `quick = TRUE`.
#' 
#' ### Missing Data
#'     
#' If there is missing data in the covariates (i.e., `NA`s in the covariates provided to `bal.tab()`), a few additional things happen. A warning will appear mentioning that missing values were present in the data set. The computed balance summaries will be for the variables ignoring the missing values. New variables will be created representing missingness indicators for each variable, named `var: <NA>` (with `var` replaced by the actual name of the variable). If `int = TRUE`, balance for the pairwise interactions between the missingness indicators will also be computed. These variables are treated like regular variables once created.
#' 
#' @references 
#' Ridgeway, G., McCaffrey, D., Morral, A., Burgette, L., & Griffin, B. A. (2016). Toolkit for Weighting and Analysis of Nonequivalent Groups: A tutorial for the twang package. R vignette. RAND.
#' 
#' Stuart, E. A. (2010). Matching Methods for Causal Inference: A Review and a Look Forward. Statistical Science, 25(1), 1-21. \doi{10.1214/09-STS313}
#' 
#' @seealso 
#' For information on the use of `bal.tab()` with specific types of objects, use the following links:
#' 
#' * [bal.tab.matchit()] for the method for objects returned by \pkg{MatchIt}.
#' * [bal.tab.weightit()] for the method for `weightit` and `weightitMSM` objects returned by \pkg{WeightIt}.
#' * [bal.tab.ps()] for the method for `ps`, `mnps`, and `iptw` objects returned by \pkg{twang} and for `ps.cont` objects returned by \pkg{twangContinuous}.
#' * [bal.tab.Match()] for the method for objects returned by \pkg{Matching}.
#' * [bal.tab.optmatch()] for the method for objects returned by \pkg{optmatch}.
#' * [bal.tab.cem.match()] for the method for objects returned by \pkg{cem}.
#' * [bal.tab.CBPS()] for the method for objects returned by \pkg{CBPS}.
#' * [bal.tab.ebalance()] for the method for objects returned by \pkg{ebal}.
#' * [bal.tab.designmatch()] for the method for objects returned by \pkg{designmatch}.
#' * [bal.tab.mimids()] for the method for objects returned by \pkg{MatchThem}.
#' * [bal.tab.sbwcau()] for the method for objects returned by \pkg{sbw}.
#' * [bal.tab.formula()] and [bal.tab.data.frame()] for the methods for `formula` and data frame interfaces when the user has covariate values and weights (including matching weights) or subclasses or wants to evaluate balance on an unconditioned data set. For data that corresponds to a longitudinal treatment (i.e., to be analyzed with a marginal structural model), see [bal.tab.time.list()].
#' 
#' See `vignette("faq")` for answers to frequently asked questions about `bal.tab()`.
#'     
#' @examples 
#' ## See individual pages above for examples with
#' ## different inputs, or see `vignette("cobalt")`

#' @rdname bal.tab
#' @export 
bal.tab <- function(x, ...) {
  
  .call <- match.call()
  
  x <- try_chk(force(x), warn = TRUE)
  
  #Replace .all and .none with NULL and NA respectively
  if (!inherits(x, "cobalt.processed.obj")) {
    
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
      .call[.alls] <- expression(NULL)
      .call[.nones] <- expression(NA)
    }
    .call[["x"]] <- process_obj(x)
    
    return(eval.parent(.call))
  }
  
  UseMethod("bal.tab")
}
