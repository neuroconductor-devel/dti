\name{dtiTensor-methods}
\docType{methods}
\alias{dtiTensor}
\alias{dtiTensor-methods}
\alias{dtiTensor,ANY-method}
\alias{dtiTensor,dtiData-method}
\title{ Methods for Function `dtiTensor' in Package `dti' }
\description{
  The method estimates, in each voxel, the diffusion tensor from the DWI data contained in an object of class \code{"dtiData"}.
}
\usage{
  \S4method{dtiTensor}{dtiData}(object, method=c( "nonlinear", "linear", "quasi-likelihood"),
          sigma = NULL, L = 1, mask=NULL, mc.cores = setCores( , reprt = FALSE))
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{method}{ Method for tensor estimation. May be \code{"linear"}, or \code{"nonlinear"}. \code{method=="quasi-likelihood"} solves the nonlinear regression problem with the
  expected value of the signal as regression function and weighting according to the signal variance.}
   \item{sigma}{(local) scale parameter of the signal's distribution.}
   \item{L}{(local) effective degrees of freedom.}
   \item{mask}{argument to specify a precomputed brain mask}
   \item{mc.cores}{Number of cores to use. Defaults to number of threads specified for openMP, see documentation of package \pkg{awsMethods}.
   Our experience suggests to use 4-6 cores if available.}
}
\value{
  An object of class \code{"dtiTensor"}.
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning. }

  \item{obj = "dtiData"}{ Estimate diffusion tensor from data in each voxel with the different options for
  the regression type and model for variance estimation.  If \code{method=="linear"} estimates are obtained
  using a linearization of the tensor model. This was the estimate used in Tabelow et.al. (2008).
  \code{method=="nonlinear"} uses a nonlinear regression model with reparametrization that ensures the
  tensor to be  positive semidefinite, see  Koay et.al. (2006).  The imlementation is based on R's internal
  C code for the
  BFGS optimization. \code{method=="quasi-likelihood"} solves the nonlinear regression problem with the
  expected value of the signal as regression function and weighting according to the signal variance.
  Tis requires additional parameters \code{sigma} and \code{L} characterizing the distribution of the signal. If \code{varmethod=="replicates"} the error variance is estimated from replicated
  gradient directions if possible, otherwise an estimate is obtained from the residual sum of squares.  If
  \code{varmodel=="global"} a homogeneous variance is assumed and estimated as the median of the local
  variance estimates.
  \code{sigma} and \code{2*L} are the scale parameter and degrees of freedom of the (local) signal distribution. \code{L} characterizes the effective number of coils. Both parameters are either scalars or arrays of the size of the images.}

}
}
\references{ J. Polzehl and K. Tabelow,
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, 44(12), 1-26 (2011).

             K. Tabelow, H.U. Voss and J. Polzehl,
             \emph{Modeling the orientation distribution function by
             mixtures of angular central Gaussian distributions},
             Journal of Neuroscience Methods, 203(1), 200-211 (2012).

             J. Polzehl and K. Tabelow,
             \emph{Structural adaptive smoothing in diffusion tensor imaging:
             The R package dti}, Journal of Statistical Software, 31(9) 1-24 (2009).

              K. Tabelow, J. Polzehl, V. Spokoiny and H.U. Voss.
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing},
              NeuroImage 39(4), 1763-1773 (2008).

              C.G. Koay, J.D. Carew, A.L. Alexander, P.J. Basser and M.E. Meyerand.
              \emph{ Investigation of Anomalous Estimates of Tensor-Derived Quantities
              in Diffusion Tensor Imaging},
              Magnetic Resonance in Medicine, 2006, 55, 930-936.

              J. Polzehl, K. Tabelow (2019). Magnetic Resonance Brain Imaging:
              Modeling and Data Analysis Using R. Springer, Use R! series. Doi:10.1007/978-3-030-29184-6.

              \url{ https://www.wias-berlin.de/projects/matheon_a3/}
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}},
  \code{\link{readDWIdata}},
  \code{\link{dtiIndices-methods}},
  \code{\link{medinria}},
  \code{\linkS4class{dtiData}},
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixtensor}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
