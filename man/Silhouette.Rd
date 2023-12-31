% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Silhouette.R
\name{Silhouette}
\alias{Silhouette}
\title{Silhouette Width}
\usage{
Silhouette(C, dx)
}
\arguments{
\item{C}{An integer vector specifying a k-partition of the dataset. min(C) must be 1
and max(C) must be k.}

\item{dx}{A "dist" object, which can be computed using stats::dist().}
}
\value{
A numeric matrix of class "silhouette" containing three columns
\describe{
\item{cluster}{A clustering of the dataset.}
\item{neighbor}{The clustering labels of the nearest clusters for all data points.}
\item{sil_width}{The silhouette widths of data points.}
}
}
\description{
This function computes the Silhouette Widths for all data points in the dataset.
}
\examples{
library("cluster")
dx = dist(faithful)
C = pam(dx, 2)$clustering
plot(Silhouette(C,dx))

}
\references{
Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and
validation of cluster analysis. J. Comput. Appl. Math., 20, 53–65.
}
\author{
Minh Long Nguyen \email{edelweiss611428@gmail.com}
}
