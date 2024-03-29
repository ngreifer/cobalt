% Rd macro for simplifying documentation writing

\newcommand{\fun}{\code{\link[=#1]{#1()}}}

% Because R packages need conditional use of packages in Suggests, any cross-reference to a doc in another package needs to be conditionally evaluated, too.

%\pkgfun{}{})tests whether the package is available, and, if so, produces a cross-reference to the function in the package; if not, the function name is displayed without a cross-reference. The first argument is the package, the second is the function name, e.g., \pkgfun{optmatch}{pairmatch}.
\newcommand{\pkgfun}{\ifelse{\Sexpr[results=rd,stage=render]{length(find.package("#1", quiet = TRUE, verbose = FALSE)) > 0}}{\code{\link[#1:#2]{#1::#2()}}}{\code{#1::#2()}}}
%\newcommand{\pkgfun}{\code{\link[#1:#2]{#1::#2()}}}

%E.g., \pkgfun{sandwich}{vcovCL} is the same as \code{\link[sandwich:vcovCL]{vcovCL}} if the sandwich package is installed and \code{vcovCL} if not.


%\pkgfun2{}{}{} does the same but allows the third argument to be printed, e.g., to use text that differs from the name of the function in the new package.
\newcommand{\pkgfun2}{\ifelse{\Sexpr[results=rd,stage=render]{length(find.package("#1", quiet = TRUE, verbose = FALSE)) > 0}}{\code{\link[#1:#2]{#3()}}}{\code{#3()}}}
%\newcommand{\pkgfun2}{\code{\link[#1:#2]{#3()}}}

%E.g., \pkgfun2{sandwich}{vcovCL}{meatCL} is the same as \code{\link[sandwich:vcovCL]{meatCL}} if the sandwich package is installed and \code{meatCL} if not.

