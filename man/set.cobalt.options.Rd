\name{set.cobalt.options}
\alias{set.cobalt.options}
\alias{get.cobalt.options}
\title{
Set options in cobalt
}
\description{
Makes it easier to set \pkg{cobalt} options. \code{set.cobalt.options} is essentially a wrapper for \code{\link{options}} but performs several checks, and \code{get.cobalt.options} is essentially a wrapper for \code{\link{getOption}}.
}
\usage{
set.cobalt.options(..., default = FALSE)

get.cobalt.options(...)
}
\arguments{
  \item{...}{
For \code{set.cobalt.options}, \code{bal.tab} parameters and the values they should take. These should be the name of the parameter in \code{bal.tab} without \code{"cobalt_"} preceding them. See examples. If any values are \code{NULL}, the corresponding options will be set back to their defaults.

For \code{get.cobalt.options}, one or more strings containing the name of a parameter option to be retrieved. See examples. If empty, all available options and their values will be returned.
}
  \item{default}{
if \code{TRUE}, sets all \pkg{cobalt} options not named in \code{...} to their default values.
}
}
\details{
When an option is set to \code{NULL}, it is set to its default value. The defaults are not displayed but are listed on the help pages where they appear. Most options correspond to display options, which can be accessed \link[=options-display]{here}. Some others (e.g., \code{continous} and \code{binary}) are described on the \code{\link{bal.tab}} help page.
}
\seealso{
\code{\link{options}}

\code{\link{display_options}} for some arguments that can be set via options.
}
\examples{
# Set un to be TRUE to always display unadjusted 
# balance measures and set binary to "std" to 
# produce standardized mean differences for 
# binary variables.

set.cobalt.options(un = TRUE, binary = "std")

# Note: the above is equivalent to:
# options(cobalt_un = TRUE, cobalt_binary = "std")
# but performs some additional checks

get.cobalt.options("un", "binary")

# Note: the above is equivalent to:
# getOption("cobalt_un")
# getOption("cobalt_binary")

# Return all cobalt options to their defaults

set.cobalt.options(default = TRUE)

# View all available options
get.cobalt.options()
}
\keyword{support functions}