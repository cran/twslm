\name{twslm}
\alias{twslm}
\alias{non.robust.twslm}
\alias{robust.twslm}
\alias{BlockByBlock}

\title{Normalization of cDNA microarray data using the two-way semi-linear model}
\description{
Normalize cDNA microarray data using the two-way semi-linear model. Two methods
are available for estimation, including a robust estimation method and the least square method.
The B-splines is used to estimate nonparametric curves in the model.
}
\usage{
twslm(sld, blk, geneid, rt, intn, df=12, degree=3, norm.only=TRUE,
    block.norm=FALSE,robust=TRUE, robust.name="Tukey",scale.constant=2.5,
    weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)

non.robust.twslm(sld, geneid, rt, intn, df=12, degree=3, norm.only=TRUE,
               tol=1e-5)

robust.twslm(sld, geneid, rt, intn, df=12, degree=3, norm.only=TRUE,
           robust=TRUE,robust.name="Tukey", scale.constant=2.5, 
           weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)

BlockByBlock(sld, blk, geneid, rt, intn, df=12, degree=3,norm.only=TRUE, 
             robust=TRUE, robust.name="Tukey", scale.constant=2.5, 
             weight.constant=4.685, ibeta=NULL,iscale=NULL,tol=1e-5)
}
\arguments{
\item{sld}{
a vector of array or slide numbers. This argument is required.
}
\item{blk}{
a vector of block numbers. This argument is required only for blockwise normalization.
}
\item{geneid}{
a vector of gene identification numbers, can be numerical numbers or gene names.
This argument is required.
}
\item{rt}{
a vector of \eqn{log_2} intensity ratio, i.e. \eqn{log_2(Cy5/Cy3)}. This argument is required.
}
\item{intn}{
a vector of average of log two total intensity, i.e. \eqn{0.5log_2(Cy5*Cy3)}. This argument is required.
}
\item{df}{
the degrees of freedom for B-spline smooth. The default is 12.
}
\item{degree}{
the order of polynomials in the B-splines. The default is 3, the cubic spline.
}
\item{norm.only}{
a logical value indicating if only normalization is carried out. The default is TRUE.
If this option is FALSE, then variance for estimated parameters of interest will be
calculated beside normalization. It will take more time for calculation if this option is FALSE.
}
\item{block.norm}{
a logical value indicating whether blockwise normalization is performed or not.
The default is FALSE, which means the default normalization is slide by slide normalization.
}
\item{robust}{
a logical value indicating if the robust procedure is incoporated in the normalization.
The default is TRUE, which means normalization is conducted using a robust method
in the two-way semilinear model. The least square method is used if this argument is FALSE.
}
\item{robust.name}{
a name for the robust procedure. The default is "Tukey", which means the
location and scale parameters are estimated iteratively with Tukey's bisquare
weight function. The other option is "Huber", which uses Huber's
weight function and the location and scale parameters are estimated iteratively.
This option works only if robust argument is TRUE.
}
\item{scale.constant}{
a constant chosen for scale estimation in the robust two-way semilinear model.
The default is 2.5.
}
\item{weight.constant}{
a constant chosen for robust location estimation.
The default is 1.345 for Huber's weight function and 4.685 for Tukey's weight function.
}
\item{ibeta}{
a vector for initalization of \eqn{\beta}. The default is NULL. The ordinary least square 
estimators for \eqn{\beta} is a good choice. Giveing this value will speed up convergence.
}
\item{iscale}{
a value for initalization of the scale parameter in the robust model. The default is NULL.
Giveing this value will speed up convergence of the algorithm.
}
\item{tol}{
a convergent criteria for iterative estimation procedure. The default is 1e-5.
}
}
\details{
Normalization is a basic step in the analysis of cDNA microarray data.
Widely used normalization method for cDNA microarray data is the Lowess normalization method
proposed by Yang et al.(2001). This method requires that at least one of the two
underlying biological assumptions, i.e. either
(i) a small fraction of genes in the experiment are differentially expressed; or
(ii) the up-regulated genes and the down-regulated genes are distributed symmetrically.
The proposed two-way semilinear model is a generalization of the semiparametric regression model.
It does not require either of the above two assumptions for normalization of cDNA microarray
data.

The proposed two-way semilinear model has the form
\deqn{y_{ij}=\phi_i(x_{ij})+\beta_j+\epsilon_{ij}.}
where \eqn{y_{ij}=log_2(Cy5/Cy3)}, \eqn{\phi_i(x_{ij})} is the intensity
dependent normalization curve for slide \eqn{i}, \eqn{x_{ij}=0.5log_2(Cy5*Cy3)},
\eqn{\beta_j} is the relative effect of gene \eqn{j}, \eqn{\epsilon_{ij}} is
the residual term, for \eqn{i=1,\ldots,n}, where \eqn{n} is the total number of slides,
\eqn{j=1,\ldots,J}, where \eqn{J} the total number of genes in the experiment.

The \code{twslm} package implements the two-way semilinear model for normalization
of cDNA microarray data. Two robust estimation procedures are implemented in the current version of \code{twslm}:
Huber's method (1981) and Tukey's method (1986). \code{twslm} can also calculate variance of estimated
parameters of interest \eqn{\beta} under the assumption of constant variance
for error terms in the model. Inference for differentially expressed genes can be carried out
based on the estmated relative gene expression levels \eqn{\hat\beta} and its variance estimator.
}
\value{
  An object of a list is returned with components:
\item{name}{
a vector of names of unique genes.
}
\item{beta}{
an estimated parameters of relative gene expression level for each gene.
}
\item{ymean}{
a mean vector for each gene after normalization. It is the arithmetic mean vector
if the least squares is used in the model. This vector will be weighted mean vector if robust methods
are used in the model.
}
\item{bvar}{
a vector of variance estimator for \eqn{\hat\beta}.
}
\item{fittedvalue}{
a vector of fitted values in the two-way semilinear model.
}
\item{bfit}{
a vector of fitted values for normalization curves.
}
\item{slide}{
a vector of slide number from the input of "twslm" function. The order is
different from the input "sld" vector.
}
\item{id}{
a vector of gene ID from the input vector "geneid" with a different order.
}
\item{ratio}{
a vector of the log two intensity ratio from the input vector "rt" with a different
order.
}
\item{intensity}{
a vector of average log two total intensity from the input vector "intn" with
a different order.
}
\item{scale}{
a scale estimator in the two-way semilinear model.
}
\item{rscale}{
a robust scale estimator if the robust method is used in the model. It is NULL for
the two-way semilinear model using ordinary least squares.
}
}
\note{ \code{twslm} is the main function to control which normalizatin method will be used.
\code{non.robust.twslm} is the function for the two-way semilinear model using the ordinary
least squares, \code{robust.twslm} is the function for robust estimation of the two-way
semilinear model, \code{BlockByBlock} is the function for blockwise normalization.
}

\references{
   Huang, J., Wang, D. & Zhang, C.H. (2005), 
   \bold{A Two-way Semi-Linear Model for Normalization and Analysis of Microarray Data}.
\emph{Journal of the American Statistical Association, 100(471):814-829}

   Wang, D., Huang, J., Xie, H., Manzella, L., Soares, M. B.,  \bold{A robust two-way semi-linear model for
 normalization of cDNA microarray data},\emph{BMC Bioinformatics 2005, 6:14}.

   Yang, Y. H., Dudoit, S., Luu, P. & Speed, T. P. (2001), \bold{Normalization for cDNA microarray}.
   In Bittner, M.L., Chen, Y., Dorsel, A.N. and Dougherty, E.R.(eds), Microarrays: Optical Technologies and
   Informatics. SPIE, Society for Optical Engineering, San Jose, CA, 4266.

   Huber, P.J. (1981), \bold{Robust Statistics}, John Wiley & Sons.

   Hampel, F. R., Ronchetti, E. M., Rousseeuw, P. J. & Stahel, W. A.(1986), \bold{Robust Statistics-The Approach
   Based on Influence Functions}, John Wiley & Sons.
}
\examples{

## Using one part of a public available dataset from Terry Speed group.
## Block one in the treatment group is chosen as an example.

data(terrycallow)
attach(terrycallow)

p<-twslm(sld=slide,geneid=id,rt=ratio,intn=intensity)

## get normalized data
normalized.data=p$ratio-p$bfit

##plot normalization curves
par(mfrow=c(3,3),cex.main=0.8,cex.lab=0.7,cex.axis=0.7,mgp=c(1,0.2,0),
         mar=c(3,3,3,1),tcl=-0.3)

for(i in 1:length(unique(p$slide))){
  plot(p$intensity[p$slide==i],p$ratio[p$slide==i],xlab="1/2log2(RG)",ylab="log2(R/G)",
     main=paste("Slide",i,sep=" "))
  ii<-order(p$intensity[p$slide==i])
  lines(p$intensity[p$slide==i][ii],p$bfit[p$slide==i][ii],col="red")
 }

detach("terrycallow")
}
\keyword{robust}
\keyword{models}
\keyword{smooth}
\keyword{nonparametric}
\author{
    Deli Wang \email{deli.wang@ccc.uab.edu} Jian Huang \email{jian@stat.uiowa.edu}
}


