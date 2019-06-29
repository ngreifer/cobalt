\docType{package}
\name{cobalt-package}
\alias{cobalt-package}
\alias{cobalt}
\title{
cobalt: Covariate Balance Tables and Plots
}
\description{
A set of tools for assessing, displaying, and reporting covariate balance in observational studies before and after preprocessing through matching, weighting, or subclassification (e.g., using propensity scores). Compatible with many of the major preprocessing packages, including \pkg{MatchIt}, \pkg{twang}, \pkg{Matching}, \pkg{WeightIt}, and others, serving as a replacement or supplement to their balance assessment functions. The focus in \pkg{cobalt} is on flexible output, methodologically recommended practices, and smart defaults.

It's critical that you read some of the documentation to understand what values are being produced. See the links below for the main functions and what output they create:

\itemize{
\item \code{\link{bal.tab}} - Generate balance tables for binary, multi-category, and continuous treatments, longitudinal/sequential treatments, and clustered/subgrouped and multiply imputed data.
\item \code{\link{bal.plot}} - Generate plots to assess distributional balance in a single covariate.
\item \code{\link{love.plot}} - Generate Love plots to display covariate balance graphically for publication.
}

In addition to the main functions, there are also some helper functions that users might find valuable:
\itemize{
\item \code{\link{f.build}} - Create a formula (e.g., \code{Z ~ X1 + X2}) from strings or data frames to reduce programming burden.
\item \code{\link{splitfactor}} - Split factors in data frame into dummy variables with flexibility in how categories are dropped.
\item \code{\link{get.w}} - Extract weights from the output of a preprocessing function from another packahe (e.g., \pkg{MatchIt}, \pkg{twang}, or \pkg{Matching}).
}
}
\section{Citing}{
Please cite \pkg{cobalt} if you use it to produce balance tables or plots in your paper to ensure replicability by others performing the same analysis. Use \code{citation("cobalt")} to generate a current citation.
}
\author{
Noah Greifer \email{noah.greifer@gmail.com}
}

