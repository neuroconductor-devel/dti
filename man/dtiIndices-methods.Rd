\name{dtiIndices-methods}
\docType{methods}
\alias{dtiIndices}
\alias{dtiIndices-methods}
\alias{dtiIndices,ANY-method}
\alias{dtiIndices,dtiTensor-method}
\title{ Methods for Function `dtiIndices' in Package `dti' }
\description{
  The method creates estimates of the fractional anisotropy (FA) and relative anisotropy (RA) indices, the main directions of anisotropy and several statistics used for visualization.
}
\usage{
  \S4method{dtiIndices}{dtiTensor}(object, mc.cores = setCores(,reprt=FALSE))
}
\arguments{
  \item{object}{ Object of class \code{"dtiTensor"} }
   \item{mc.cores}{Number of cores to use. Defaults to number of threads specified for openMP, see documentation of package \pkg{awsMethods}.
   Our experience suggests to use 4-6 cores if available.}
}
\value{
  An object of class \code{"dtiIndices"}.
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning. }

\item{obj = "dtiTensor"}{ Estimate tensor indices like trace, fractional and geodesic anisotropy, main diffusion direction and shape parameters. }

}
}
\references{ J. Polzehl, K. Tabelow (2019). Magnetic Resonance Brain Imaging:
Modeling and Data Analysis Using R. Springer, Use R! series. Doi:10.1007/978-3-030-29184-6.

             J. Polzehl and K. Tabelow,
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
  \code{\link{medinria}},
  \code{\link{dtiTensor-methods}},
  \code{\linkS4class{dtiTensor}},
  \code{\linkS4class{dtiIndices}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ manip }
\keyword{ models }
