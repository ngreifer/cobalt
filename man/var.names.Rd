\name{var.names}
\alias{var.names}
\title{Extract Variable Names from \code{bal.tab} Objects}
\description{
This function extracts variable names from a \code{bal.tab} object for use in specifying alternate variable names in \fun{love.plot}. Optionally, a file can be written for easy editing of names.
}
\usage{
var.names(b, 
          type, 
          file = NULL, 
          minimal = FALSE)
}
\arguments{
\item{b}{a \code{bal.tab} object; the output of a call to \fun{bal.tab}.
}

\item{type}{the type of output desired. Can either be \code{"df"} for a \code{data.frame} or \code{"vec"} for a named vector. See "Value". The default is \code{"vec"} unless \code{file} is not \code{NULL}.
}
\item{file}{optional; a file name to save the output if \code{type = "df"}. See \fun{write.csv}, which \code{var.name} calls. Must end in \code{.csv}.
}
\item{minimal}{whether the output should contain all variable names (i.e., all rows that appear the output of \code{bal.tab()}) or just the unique base variables. See "Details".
}
}
\value{
If \code{type = "vec"}, a character vector the the variable names both as the names and the entries.

If \code{type = "df"}, a \code{data.frame} with two columns called \code{"old"} and \code{"new"}, each with the variables as the entries.

If \code{file} is not \code{NULL}, the output will be returned invisibly.
}
\details{
The goal of the function is to make supplying new variable names to the \code{var.names} argument in \fun{love.plot} easier. Rather than manually creating a vector or \code{data.frame} with all the variable names that one desires to change, one can use \code{var.names()} to extract variable names from a \code{bak.tab} object and edit the output. Importantly, the output can be saved to a CSV file, which can be easily edited and read back into \R for use in \code{love.plot()}, as demonstrated in the Example.

When \code{minimal = TRUE}, only a minimal set of variables will be output. For example, if the variables analyzed in \code{bal.tab()} are age, race, and married, and \code{int = TRUE} in \code{bal.tab()}, many variables will appear in the output, including expansions of the factor variables, the polynomial terms, and the interactions. Rather than renaming all of these variables individually, one can rename just the three base variables, and all variables that arise from them will be accordingly renamed. Setting \code{minimal = TRUE} requests only these base variables.
}
\note{
Not all programs can properly read the Unicode characters for the polynomial terms when requested. These may appear strange in, e.g., Excel, but \R will process the characters correctly.
}

\examples{
data(lalonde, package = "cobalt")

b1 <- bal.tab(treat ~ age + race + married, data = lalonde,
             int = TRUE)
v1 <- var.names(b1, type = "vec", minimal = TRUE)
v1["age"] <- "Age (Years)"
v1["race"] <- "Race/Eth"
v1["married"] <- "Married"
love.plot(b1, var.names = v1)

\dontrun{
b2 <- bal.tab(treat ~ age + race + married + educ + nodegree +
              re74 + re75 + I(re74==0) + I(re75==0), 
              data = lalonde)
var.names(b2, file = "varnames.csv")

##Manually edit the CSV (e.g., in Excel), then save it.
v2 <- read.csv("varnames.csv")
love.plot(b2, var.names = v2)
}
}

\keyword{support functions}
