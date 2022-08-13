\name{cobalt-package}
\alias{cobalt-package}
\alias{cobalt}
\docType{package}
\title{
cobalt: \packageTitle{cobalt}
}
\description{
\if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}

A set of tools for assessing covariate balance in observational studies numerically and graphically. The functions provide integration with the major R packages used for balancing covariates, including \pkg{MatchIt}, \pkg{WeightIt}, \pkg{twang}, \pkg{CBPS}, and many others, and support objects not made using these packages. They support binary, multi-category and continuous treatments, point and longitudinal treatments, and clustered and multiply imputed data.

The main functions of \pkg{cobalt} are the following:
\itemize{
\item{\fun{bal.tab} - generate tables of balance statstics before and after matching, weighting, or subclassification}
\item{\fun{bal.plot} - generate plots to assess balance visually on one covariate at a time}
\item{\fun{love.plot} - generate plots to summarize and report balance statistics}
}

Other functions include \fun{get.w} for extracting weights from objects produced by other packages, \fun{col_w_smd} (and friends documentated on the same page) for computing (weighted) balance statistics outside of \code{bal.tab()}, and \fun{splitfactor} for splitting factor variables in a dataset into dummy variables.

\pkg{cobalt} has several vignettes, which can be accessed using \code{\link[=vignette]{vignette(package = "cobalt")}} or visiting the website at \url{https://ngreifer.github.io/cobalt/}.
}
\author{
Noah Greifer (\href{https://orcid.org/0000-0003-3067-7154}{ORCID})

To contact me about \pkg{cobalt}, please use the GitHub \href{https://github.com/ngreifer/cobalt/issues}{issues page}.
}
\keyword{internal}