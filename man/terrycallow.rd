\name{terrycallow}
\alias{terrycallow}

\title{Apo AI Data}
\usage{data(terrycallow)}
\description{
  This data set consists 4 variables.
}
\format{
  A list with 3186 observations on 4 variables:
  \tabular{rll}{
    [ , 1] \tab slide     \tab array number \cr
    [ , 2] \tab id        \tab clone ID \cr
    [ , 3] \tab ratio     \tab \eqn{log_2\frac{R}{G}} \cr
    [ , 4] \tab intensity \tab \eqn{0.5log_2(RG)} \cr
  }
where R is background adjusted signal intensity in the Cy5 channel,
G is background adjusted signal intensity in the Cy3 channel.
}
\source{
This data set is one part of Apo AI experiment by Callow et. al.(2000).
Block one in the treatment group is chosen for this example. The whole Apo AI
data set can be found from \url{http://stat-www.berkeley.edu/users/terry/zarray/Html/apodata.html}.
}
\details{
  Apo AI dataset is used by several authors to study cDNA microarray data normalization
  and significance analysis.
}
\references{
 Matthew J. Callow, Sandrine Dudoit, Elaine L. Gong, Terence P. Speed,
 and Edward M. Rubin (2000):
 Microarray expression profiling identifies genes with altered expression in HDL-deficient mice.
  \emph{Genome Research}, Vol. 10, Issue 12, 2022-2029.
  \url{http://www.genome.org}
}
\keyword{datasets}


