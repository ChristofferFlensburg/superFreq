% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XRank.R
\name{mcri}
\alias{mcri}
\title{Turns colours into similar colours from the Murdoch Childrens Research Institute palette.}
\usage{
mcri(col = "deafult", al = 1)
}
\arguments{
\item{col}{character or numeric. Most common colours such as "green" or 'blue' are converted to MCRI. Numbers return the colour matching that number.}

\item{al}{The alpha parameter: the opaqueness. Numeric between 0 and 1.}
}
\description{
Turns colours into similar colours from the Murdoch Childrens Research Institute palette.
}
\details{
The function return a colour from the MCRI palette if match found, otherwise the input is returned unchanged. call mcri() to see the available colours.
}
\examples{
plot(1:9, rep(1,9), pch=16, cex=10, col=mcri(0:8))
text(1:9, rep(1.1,9), c('blue', 'orange', 'green', 'magenta',
                        'cyan', 'red', 'violet', 'darkblue', 'darkred'),
     col=mcri(0:8), cex=1.5, font=2)

}
