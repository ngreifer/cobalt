\name{splitfactor}
\alias{splitfactor}
\alias{unsplitfactor}
\title{Split and Unsplit Factors into Dummy Variables}
\description{
\code{splitfactor()} splits factor variables into dummy (0/1) variables. This can be useful when functions do not process factor variables well or require numeric matrices to operate. \code{unsplitfactor()} combines dummy variables into factor variables, undoing the operation of \code{splitfactor()}.
}
\usage{
splitfactor(data, var.name, replace = TRUE, sep = "_", 
            drop.level = NULL, drop.first = TRUE, 
            drop.singleton = FALSE, 
            drop.na = TRUE, check = TRUE)

unsplitfactor(data, var.name, replace = TRUE, sep = "_", 
                dropped.level = NULL, dropped.na = TRUE)
}
\arguments{
\item{data}{
A data.frame containing the variables to be split or unsplit. In \code{splitfactor()}, can be a factor variable to be split.
}
\item{var.name}{
For \code{splitfactor()}, the names of the factor variables to split. If not specified, will split all factor variables in \code{data}. If \code{data} is a factor, the stem for each of the new variables to be created. For \code{unsplitfactor()}, the name of the previously split factor.
}
\item{replace}{
Whether to replace the original variable(s) with the new variable(s) (\code{TRUE}) or the append the newly created variable(s) to the end of the data set (\code{FALSE}).
}
\item{sep}{
A character separating the the stem from the value of the variable for each dummy. For example, for \code{"race_black"}, \code{sep = "_"}.
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
\item{drop.na}{
If \code{NA}s are present in the variable, how to handle them. If \code{TRUE}, no new dummy will be created for \code{NA} values, but all created dummies will have \code{NA} where the original variable was \code{NA}. If \code{FALSE}, \code{NA} will be treated like any other factor level, given its own column, and the other dummies will have a value of 0 where the original variable is \code{NA}.
}
\item{check}{
Whether to make sure the variables specified in \code{var.name} are actually factor (or character) variables. If splitting non-factor (or non-character) variables into dummies, set \code{check = FALSE}. If \code{check = FALSE} and \code{data} is a data frame, an argument to \code{var.name} must be specified.
}
\item{dropped.level}{
The value of the original factor variable whose dummy was dropped when the variable was split. If left empty and a dummy was dropped, the resulting factor will have the value \code{NA} instead of the dropped value. Can only be specified when one factor is to be unsplit.
}
\item{dropped.na}{
If \code{TRUE}, will assume that \code{NA}s in the variables to be unsplit correspond to \code{NA} in the unsplit factor (i.e., that \code{drop.na = TRUE} was specified in \code{split.factor()}). If \code{FALSE}, will assume there is a dummy called "var.name_stem_NA" (e.g., "x_NA") that contains 1s where the unsplit factor should be \code{NA} (i.e., that \code{drop.na = FALSE} was specified in \code{split.factor()}. If \code{NA}s are stored in a different column with the same stem, e.g., "x_miss", that name (e.g., "miss") can be entered instead.
}
}
\details{
If there are \code{NA}s in the variable to be split, the new variables created by \code{splitfactor()} will have \code{NA} where the original variable is \code{NA}.
}
\value{
For \code{splitfactor()}, a data frame containing the original data set with the newly created dummies. For \code{unsplitfactor()}. a data frame containing the original data set with the newly created factor variables.
}


\examples{
data("lalonde", package = "cobalt")

lalonde.split <- splitfactor(lalonde, "race", 
                           replace = TRUE, 
                           drop.first = TRUE)
# A data set with "race_hispan" and "race_white" instead
# of "race".

lalonde.unsplit <- unsplitfactor(lalonde.split, "race", 
                        replace = TRUE,
                        dropped.level = "black")

all.equal(lalonde, lalonde.unsplit) #TRUE
}

