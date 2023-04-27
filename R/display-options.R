#' @title Options for Displaying `bal.tab()` Output
#' @name display-options
#' 
#' @description
#' Several additional arguments can be passed to [bal.tab()] that control the display of the output; these arguments are documented here. Not all arguments are applicable to all uses of `bal.tab()`; for example, `which.subclass`, which controls which subclasses are displayed when subclassification is used, won't do anything when subclassification is not used. Note that when `quick = TRUE` is set in the call to `bal.tab()` (which is the default), setting any of these arguments to `FALSE` can prevent some values from being computed, which can have unintended effects.
#' 
#' @section Allowed arguments:
#' 
#' \describe{
#' \item{`disp.bal.tab`}{`logical`; whether to display the table of balance statistics. Default is `TRUE`, so the balance table is displayed.
#' }
#' 
#' \item{`imbalanced.only`}{`logical`; whether to display only the covariates that failed to meet at least one of balance thresholds. Default is `FALSE`, so all covariates are displayed.
#' }
#' 
#' \item{`un`}{`logical`; whether to print statistics for the unadjusted sample as well as for the adjusted sample. Default is `FALSE`, so only the statistics for the adjusted sample are displayed.
#' }
#' 
#' \item{`disp`}{`character`; which distribution summary statistic(s) should be reported. Allowable options include `"means"` and `"sds"`. Multiple options are allowed. Abbreviations allowed.
#' }
#' 
#' \item{`stats`}{`character`; which statistic(s) should be reported. See [`stats`][balance-statistics] to see which options are available. Multiple options are allowed. Abbreviations allowed. For binary and multi-category treatments, the default is `"mean.diffs"` (i.e., \[standardized\] mean differences), and for continuous treatments, the default is `"correlations"` (i.e., treatment-covariate Pearson correlations).
#' }
#' 
#' \item{`factor_sep`}{`character`; the string used to separate factor variables from their levels when variable names are printed. Default is `"_"`.
#' }
#' 
#' \item{`int_sep`}{`character`; the string used to separate two variables involved in an interaction when variable names are printed. Default is `" * "`. Older versions of \pkg{cobalt} used `"_"`.
#' }
#' 
#' \item{`disp.call`}{`logical`; whether to display the function call from the original input object, if present. Default is `FALSE`, so the function call is not displayed.
#' }
#' }
#' 
#' ## When subclassification is used
#' 
#' \describe{
#' \item{`which.subclass`}{Which subclasses (if any) should be displayed. If `.all`, all subclasses will be displayed. If `.none` (the default), no subclasses will be displayed. Otherwise, can be a vector of subclass indices for which to display balance.
#' }
#' 
#' \item{`subclass.summary`}{`logical`; whether to display the balance summary across subclasses. If `TRUE`, the balance summary across subclasses will be displayed. The default is `TRUE`, and if `which.subclass` is `.none`, it will automatically be set to `TRUE`.
#' }
#' }
#' ## When the treatment is multi-category
#' 
#' \describe{
#' \item{`which.treat`}{For which treatments or treatment combinations balance tables should be displayed. If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatments that both appear in the input will be displayed. For example, setting  `which.treat = "A"` will display "A vs. B" and "A vs. C", while setting `which.treat = c("A", "B")` will only display "A vs. B". `.none` indicates no treatment comparisons will be displayed, and `.all` indicates all treatment comparisons will be displayed. Default is `.none`. See [`bal.tab.multi()`][class-bal.tab.multi].
#' }
#' 
#' \item{`multi.summary`}{`logical`; whether to display the balance summary across all treatment pairs. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. Default is `TRUE` when `which.treat` is `.none` and `FALSE` otherwise. See [`bal.tab.multi()`][class-bal.tab.multi].
#' }
#' }
#' 
#' ## When clusters are present
#' 
#' \describe{
#' \item{`which.cluster`}{For which clusters balance tables should be displayed. If `.all`, all clusters in `cluster` will be displayed. If `.none`, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor). Default is `.all`. See [`class-bal.tab.cluster`].
#' }
#' 
#' \item{`cluster.summary`}{`logical`; whether to display the balance summary across clusters. Default is `TRUE` when `which.cluster` is `.none` and `FALSE` otherwise (note the default for `which.cluster` is `.all`). See [`class-bal.tab.cluster`].
#' }
#' 
#' \item{`cluster.fun`}{Which function is used in the across-cluster summary to combine results across clusters. Can be "min", "mean", or "max". For example, if `cluster.fun = "mean"` the mean balance statistic across clusters will be displayed. The default when `abs = FALSE` in the `bal.tab()` call is to display all three. The default when `abs = FALSE` in the `bal.tab()` call is to display just the mean and max balance statistic. See [`class-bal.tab.cluster`].
#' }
#' }
#' 
#' ## When multiple imputations are present
#' 
#' \describe{
#'\item{`which.imp`}{For which imputations balance tables should be displayed. If `.all`, all imputations in `imp` will be displayed. If `.none`, no imputations will be displayed. Otherwise, can be a vector of imputation indices for which to display balance. Default is `.none`. See [`class-bal.tab.imp`].
#' }
#' 
#'\item{`imp.summary`}{`logical`; whether to display the balance summary across imputations. Default is `TRUE` when `which.imp` is `.none` and `FALSE` otherwise. See [`class-bal.tab.imp`].
#' }
#' 
#' \item{`imp.fun`}{Which function is used in the across-imputation summary to combine results across imputations. Can be "min", "mean", or "max". For example, if `imp.fun = "mean"` the mean balance statistic across imputations will be displayed. The default when `abs = FALSE` in the `bal.tab()` call is to display all three. The default when `abs = FALSE` in the `bal.tab()` call is to display just the mean and max balance statistic. See [`class-bal.tab.imp`].
#' }
#' }
#' 
#' ## When the treatment is longitudinal
#' 
#' \describe{
#' \item{`which.time`}{For which time points balance tables should be displayed. If `.all`, all time points will be displayed. If `.none`, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance. Default is `.none`. See [`class-bal.tab.msm`].
#' }
#' 
#' \item{`msm.summary`}{`logical`; whether to display the balance summary across time points. Default is `TRUE` when `which.time` is `.none` and `FALSE` otherwise. See [`class-bal.tab.msm`].
#' }
#' }
#' 
#' @section Setting options globally:
#' 
#' In addition to being able to be specified as arguments, if you find you frequently set a display option to something other than its default, you can set that as a global option (for the present R session) using [set.cobalt.options()] and retrieve it using [get.cobalt.options()]. Note that global options cannot be set for `which.subclass`, `which.cluster`, `which.imp`, `which.treat`, or `which.time`.
#' 
#' @note 
#' When calling `bal.tab()` using [do.call()], if you are using `.all` or `.none` as inputs to arguments, you need to use [alist()] rather than [list()] to group the arguments. For example, `do.call(bal.tab, list(., which.cluster = .none))` will produce an error, but `do.call(bal.tab, alist(., which.cluster = .none))` should work correctly.
#' 
#' 
#' @seealso
#' [bal.tab()], [print.bal.tab()]
#' 
#' 
NULL
