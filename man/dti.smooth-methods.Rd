\name{dti.smooth-methods}
\docType{methods}
\alias{dti.smooth}
\alias{dti.smooth-methods}
\alias{dti.smooth,ANY-method}
\alias{dti.smooth,dtiData-method}
\alias{dti.smooth,dtiTensor-method}
\title{ Methods for Function `dti.smooth' in Package `dti' }
\description{
  The function provides structural adaptive smoothing for diffusion weighted image data within the context of an diffusion tensor (DTI) model.  It implements smoothing of DWI data using a structural assumption of a local (anisotropic) homogeneous  diffusion tensor model (in case a \code{"dtiData"}-object is provided).  It also implements structural adaptive smoothing of a diffusion tensor using a Riemannian metric (in case a \code{"dtiTensor"}-object is given), although we strictly recommend to use the first variant due to methodological reasons.
}
\usage{
\S4method{dti.smooth}{dtiData}(object, hmax=5, hinit=NULL, lambda=20, tau=10, rho=1, 
         graph=FALSE,slice=NULL, quant=.8, minfa=NULL, hsig=2.5, 
         lseq=NULL, method="nonlinear", rician=TRUE, 
         niter=5,result="Tensor")
 }
\arguments{
  \item{object}{ Either an object of class \code{"dtiData"} or an object of class \code{"dtiTensor"}}
  \item{hmax}{Maximal bandwidth}
  \item{hinit}{Initial bandwidth (default 1)}
  \item{lambda}{Critical parameter (default 20)}
  \item{tau}{Critical parameter for orientation scores (default 10)}
  \item{rho}{Regularization parameter for anisotropic vicinities (default 1)}
  \item{graph}{"logical": Visualize intermediate results (default FALSE)}
  \item{slice}{slice number, determines the slice used in visualization}
  \item{quant}{determines \code{minfa} as corresponding quantile of FA if \code{is.null(minfa)} }
  \item{minfa}{minimal anisotropy index (FA) to use in visualization}
  \item{hsig}{bandwidth for presmoothing of variance estimates}
  \item{lseq}{ sequence of correction factors for \code{lambda}}
  \item{method}{ Method for tensor estimation. May be \code{"linear"}, \code{"nonlinear"}}
  \item{rician}{ "logical": apply a correction for Rician bias. This is still experimental and depends on spatial independence of errors.}
  \item{niter}{ Maximum number of iterations for tensor estimates using the nonlinear model.}
  \item{result}{ Determines the created object. Alternatives are \code{"Tensor"} for create a dtiTensor-object and \code{"dtiData"}
   for a dtiData-object containing a smoothed data cube.}
}
\section{Methods}{
\describe{

  \item{object = "ANY"}{ Returns a warning. }

  \item{object = "dtiData"}{ We highly recommend to use the method \code{dti.smooth} on DWI data directly, i.e. on an object of class \code{"dtiData"}, due to methodological reasons, see Tabelow et al. (2008).  It is usually not necessary to use any other argument than \code{hmax}, which defines the maximum bandwidth of the iteration. 
  
  If \code{model=="linear"} estimates are obtained using a linearization of the tensor model. This was the estimate used in Tabelow et.al. (2008). \code{model=="nonlinear"} uses a nonlinear regression model with reparametrization that ensures the tensor to be  positive semidefinite, see  Koay et.al. (2006).  If \code{varmethod=="replicates"} the error variance is estimated from replicated gradient directions if possible, otherwise (default) an estimate is obtained from the residual sum of squares.  If \code{volseq==TRUE} the sum of location weights is fixed to \eqn{1.25^k}{1.25^k} within iteration \eqn{k}{k} (does not depend on the actual tensor). Otherwise the ellipsoid of positive location weights is determined by a bandwidth \eqn{h_k = 1.25^(k/3)}{h_k = 1.25^(k/3)}. }

}}
\value{
  An object of class \code{dtiTensor}.
}
\references{               J. Polzehl and K. Tabelow, 
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, to appear.

             K. Tabelow, H.U. Voss and J. Polzehl, 
             \emph{Modeling the orientation distribution function by 
             mixtures of angular central Gaussian distributions},
             Journal of Neuroscience Methods, to appear.

             J. Polzehl and K. Tabelow, 
             \emph{Structural adaptive smoothing in diffusion tensor imaging: 
             The R package dti}, Journal of Statistical Software, 31 (2009) pp. 1--24.
             
              K. Tabelow, J. Polzehl, V. Spokoiny and H.U. Voss. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).


              \url{ https://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\link{dtiTensor-methods}}, 
  \code{\link{dtiIndices-methods}}, 
  \code{\link{medinria}} , 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}},
  \code{\linkS4class{dtiIndices}}
} 
\keyword{ methods }
\keyword{ smooth }
