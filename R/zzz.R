### coommon method to ease the usage
#######
#' @export
#' @method print tor_obj
print.tor_obj <- function(x, ...) {
  print(tor_summarise(x, ...))
}

#######
#' @export
#' @method summary tor_obj
summary.tor_obj <- function(object, ...) {
  tor_summarise(object)
}

#######
#' @export
#' @method plot tor_obj
plot.tor_obj <- function(x, ...) {
  tor_plot(tor_obj = x,...)
}

#######

#' Display welcome message when loading the package
#'
#' Display a message when the package is being loaded.
#' @param  libname a character string giving the library directory where the package defining the namespace was found.
#' @param pkgname a character string giving the name of the package.
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "Torpor version: ", utils::packageDescription("torpor")$Version, " is loaded",
    "\n",
    "Run: vignette('example') or browseVignettes('torpor') for an example",
    "\n",
    "Find more information about the model in the companion article"
  )
}
