\name{Balance Summary}
\alias{col_w_mean}
\alias{col_w_smd}
\alias{col_w_ks}
\alias{col_w_ovl}
\alias{col_w_corr}

\title{
Compute balance statistics for covariates
}
\description{
These functions quickly compute balance statistics for the given covariates. These functions are used in \code{\link{bal.tab}}, but they are available for use in programming without having to call \code{bal.tab()} to get them.
}
\usage{
col_w_mean(mat, weights = NULL, s.weights = NULL, na.rm = TRUE, ...)
col_w_smd(mat, treat, weights = NULL, std = TRUE, s.d.denom = "pooled", abs = FALSE, s.weights = NULL, bin.vars = NULL, na.rm = TRUE, ...)
col_w_ks(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, ...)
col_w_ovl(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, integrate = FALSE, ...)
col_w_corr(mat, treat, weights = NULL, abs = FALSE, s.weights = NULL, bin.vars = NULL, na.rm = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{x} here~~
}
  \item{weights}{
%%     ~~Describe \code{x} here~~
}
  \item{s.weights}{
%%     ~~Describe \code{x} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{x} here~~
}
  \item{treat}{
%%     ~~Describe \code{x} here~~
}
  \item{std}{
%%     ~~Describe \code{x} here~~
}
  \item{s.d.denom}{
%%     ~~Describe \code{x} here~~
}
  \item{abs}{
\code{logical}; whether the return statistics should be in absolute value (\code{TRUE}) or not. Default is \code{FALSE}.
}
  \item{bin.vars}{
%%     ~~Describe \code{x} here~~
}
  \item{integrate}{
\code{logical};for \code{col_w_ovl}, whether to use \code{\link{integrate}} to calculate the area over overlap. If \code{FALSE}, a midpoint Riemann sum with 1000 partitions will be used instead. The Riemann sum is a little slower very slightly imprecise (unnoticibly in most contexts), but the integral can fail sometimes and thus is less stable. The default is to use the Riemann sum.
}
  \item{...}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bal.tab}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
