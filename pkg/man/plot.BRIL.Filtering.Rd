% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_outliers.R
\name{plot.BRIL.Filtering}
\alias{plot.BRIL.Filtering}
\title{Plot method for \code{BRIL.Filtering} objects}
\usage{
\method{plot}{BRIL.Filtering}(
  x,
  contents = c("p.values", "scatterplot", "dist", "hist"),
  showCenter = TRUE,
  showSelection = TRUE,
  col = "black",
  colSelection = "red",
  colCenter = "orange",
  mtextTitles = TRUE,
  mfrow,
  ...
)
}
\arguments{
\item{x}{An object of class \code{BRIL.Filtering} (see \code{\link[=filter_outliers]{filter_outliers()}})}

\item{contents}{Contents to be displayed, options are "p.values",
"scatterplot", "dist", "hist" and "all"}

\item{showCenter}{Logical value, to show the center used in the filtering}

\item{showSelection}{Logical value, to highlight the samples selected by the
filtering process}

\item{col}{Default color for non-selected samples (default: "black")}

\item{colSelection}{Color of the selected samples (default: "red")}

\item{colCenter}{Color for the center in "scatterplot" (default: "orange")}

\item{mtextTitles}{Logical value, \code{TRUE} to set smaller titles/subtitles on
top, \code{FALSE} to use the default plot title options.}

\item{mfrow}{Number of rows and columns of the figure (example:
c(4,1))}

\item{...}{Other arguments passed to or from other methods (such as pch for
the symbols, main and sub for title and subtitle, xlab, xmin, ...)}
}
\description{
Plot method for \code{BRIL.Filtering} objects
}
\details{
\emph{Red intercept lines correspond to the
selection based on the p.values exceeding the given threshold.} \cr
\emph{To display all the p.values, rerun the function \code{\link[=filter_outliers]{filter_outliers()}}
with the parameter \code{debug = TRUE}} \cr \cr

\code{contents} options:
\itemize{
\item \bold{"p.values"} provides a plot of the test p.values (in function of
the subset size)
\item \bold{"scatterplot"} displays the data in cartesian coordinates. Selected
samples are displayed in red, and the center used to compute distances as
an orange cross
\item \bold{"dist"} shows the distances of each sample to the center provided
in \code{\link[=filter_outliers]{filter_outliers()}} (in function of sample index)
\item \bold{"hist"} draws an histogram of the distances of the samples to
the center provided in \code{\link[=filter_outliers]{filter_outliers()}}
\item \bold{"all"} displays a figure with all of the options above
}
}
\examples{

# Illustrative data
XY <- rbind(
  mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
  mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
  mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
  mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
)

# Process the data
filtering <- filter_outliers(XY, median_rec(XY)$median, test = "DIP", debug = TRUE)

# Plot all default figures
plot(filtering)

# Plot P.Values and Scatterplot only
plot(filtering, contents = c("pvalues", "scatterplot"))

# Change the layout to vertical
plot(filtering, contents = c("pvalues", "scatterplot"), mfrow = c(2, 1))

# Remove title, subtitle, and axis labels
plot(filtering, contents = "scatterplot", main = "", sub = "",
  ylab = "", xlab = "")

# Other graphical options
plot(filtering, contents = "scatterplot", asp = 1,
  xlim = c(-30, 30), ylim = c(-30, 30))

plot(filtering,
  contents = "scatterplot", asp = 1, pch = 4, lwd = 2, col = "blue",
  colSelection = "green", showCenter = FALSE
)

plot(filtering, contents = "hist", main = "My Histogram",
  showSelection = FALSE, breaks = 50)

}
\seealso{
\code{\link[=filter_outliers]{filter_outliers()}}, \code{\link[=print.BRIL.Filtering]{print.BRIL.Filtering()}}, \code{\link[=bril]{bril()}},
\code{\link[=median_rec]{median_rec()}}, \code{\link[=median_mv]{median_mv()}}, \code{\link[=depth_values]{depth_values()}}
}
