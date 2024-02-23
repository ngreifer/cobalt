#' @keywords internal
#' @description 
#' \if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
#' 
#' A set of tools for assessing covariate balance in observational studies numerically and graphically. The functions provide integration with the major R packages used for balancing covariates, including \pkg{MatchIt}, \pkg{WeightIt}, \pkg{twang}, \pkg{CBPS}, and many others, and support objects not made using these packages. They support binary, multi-category and continuous treatments, point and longitudinal treatments, and clustered and multiply imputed data.
#' 
#' The main functions of \pkg{cobalt} are the following:
#'     
#' * [bal.tab()] - generate tables of balance statistics before and after matching, weighting, or subclassification
#' * [bal.plot()] - generate plots to assess balance visually on one covariate at a time
#' * [love.plot()] - generate plots to summarize and report balance statistics
#' 
#' 
#' Other functions include [get.w()] for extracting weights from objects produced by other packages, [col_w_smd()] (and friends documented on the same page) for computing (weighted) balance statistics outside of `bal.tab()`, [bal.compute()] for computing scalar balance statistics efficiently, and [splitfactor()] for splitting factor variables in a dataset into dummy variables.
#' 
#' \pkg{cobalt} has several vignettes, which can be accessed using `vignette(package = "cobalt")` or visiting the website at <https://ngreifer.github.io/cobalt/>.
"_PACKAGE"

## usethis namespace: start
#' @import stats
#' @importFrom crayon italic
#' @importFrom crayon underline
#' @importFrom ggplot2 .data
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
## usethis namespace: end
NULL
