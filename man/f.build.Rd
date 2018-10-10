\name{f.build}
\alias{f.build}
\title{Convenient Formula Generation}
\description{
\code{f.build()} returns a formula of the form \code{y ~ x1 + x2 + ...} from a data frame input. It can be much quicker to use \code{f.build()} than to hand-write the precise formula, which may contain errors. It can be used in place of a formula in, for example, \code{\link{glm}()}, \code{matchit()}, or \code{\link{bal.tab}()}.
}
\usage{
f.build(y, rhs)
}
\arguments{
\item{y}{the quoted name of the response (left hand side) variable in the formula. Only one variable is supported. If missing, \code{NULL}, or the empty string (\code{""}), the formula will have no response variable.}

\item{rhs}{a data frame whose variable names will be the terms on the right hand side of the formula, or a character vector whose values will be the terms on the right hand side of the formula.}
}
\value{
an object of class \code{"formula"}.
}


\examples{
data(lalonde)
covs <- subset(lalonde, select = -c(treat, re78))
lm(f.build("treat", covs), data = lalonde)
}

