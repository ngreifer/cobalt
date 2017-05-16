\name{split.factor}
\alias{split.factor}
\title{Split Factors into Dummy Variables}
\description{
Splits factor variables into dummy (0/1) variables. This can be useful when functions do not process factor variables well or require numeric matrices to operate.
}
\usage{
split.factor(data, varname, replace = TRUE, sep = "_", drop.first = c("if2", TRUE, FALSE), drop.singleton = FALSE)
}
\arguments{
\item{data}{
A data.frame containing factor variables to be split, or a factor variable to be split.
}
\item{varname}{
If \code{data} is a data.frame, the names of the factor variables to split. If not specified, will split all factor variables in \code{data}. If \code{data} is a factor, the stem for each of the new variables to be created.
}
\item{replace}{
Whether to replace the original factor variable with the dummies (\code{TRUE}) or the append the newly created dummies to the end of the data set (\code{FALSE}).
}
\item{sep}{
A character seperating the the stem from the value of the variable for each newly created dummy.
}
\item{drop.first}{
Whether to drop the first dummy created for each factor. If \code{"if2"}, will only drop the first category if the factor has exactly two levels.
}
\item{drop.singleton}{
Whether to drop a factor variable if it only has one level.
}
}
\value{
A data.frame containing the newly created dummy variables replacing
}


\examples{

}

