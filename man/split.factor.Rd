\name{split.factor}
\alias{split.factor}
\alias{un.split.factor}
\title{Split and Unsplit Factors into Dummy Variables}
\description{
\code{split.factor()} splits factor variables into dummy (0/1) variables. This can be useful when functions do not process factor variables well or require numeric matrices to operate. \code{un.split.factor()} combines dummy variables into factor variables, undoing the operation of \code{split.factor()}.
}
\usage{
split.factor(data, var.name, replace = TRUE, sep = "_", 
             drop.level = NULL, drop.first = c(TRUE, 
             FALSE, "if2"), drop.singleton = FALSE, 
             check = TRUE)

un.split.factor(data, var.name, replace = TRUE, sep = "_", 
                dropped.level = NULL)
}
\arguments{
\item{data}{
A data.frame containing the variables to be split or unsplit. In \code{split.factor()}, can be a factor variable to be split.
}
\item{var.name}{
For \code{split.factor()}, the names of the factor variables to split. If not specified, will split all factor variables in \code{data}. If \code{data} is a factor, the stem for each of the new variables to be created. For \code{un.split.factor()}, the name of the previously split factor.
}
\item{replace}{
Whether to replace the original variable(s) with the new variable(s) (\code{TRUE}) or the append the newly created variable(s) to the end of the data set (\code{FALSE}).
}
\item{sep}{
A character seperating the the stem from the value of the variable for each dummy. For example, for \code{"race_black"}, \code{sep = "_"}.
}
\item{drop.level}{
The name of a level of \code{var.name} for which to drop the dummy variable. Only works if there is only one variable to be split.
}
\item{drop.first}{
Whether to drop the first dummy created for each factor. If \code{"if2"}, will only drop the first category if the factor has exactly two levels. The default is to always drop the first dummy (\code{TRUE}).
}
\item{drop.singleton}{
Whether to drop a factor variable if it only has one level.
}
\item{check}{
Whether to make sure the variables specified in \code{var.name} are actually factor (or character) variables. If splitting non-factor (or non-character) variables into dummies, set \code{check = FALSE}. If \code{check = FALSE} and \code{data} is a data frame, an argument to \code{var.name} must be specified.
}
\item{dropped.level}{
The value of the original factor variable whose dummy was dropped when the variable was split. If left empty and a dummy was dropped, the resulting factor will have the value \code{NA} instead of the dropped value. Can only be specified when one factor is to be unsplit.
}
}
\value{
For \code{split.factor()}, a data frame containing the original data set with the newly created dummies. For \code{un.split.factor()}. a data frame containing the original data set with the newly created factor variables.
}


\examples{
data("lalonde", package = "cobalt")

lalonde.split <- split.factor(lalonde, "race", 
                           replace = TRUE, 
                           drop.first = TRUE)
# A data set with "race_hispan" and "race_white" instead
# of "race".

lalonde.unsplit <- un.split.factor(lalonde.split, "race", 
                        replace = TRUE,
                        dropped.level = "black")

all.equal(lalonde, lalonde.unsplit) #TRUE
}

