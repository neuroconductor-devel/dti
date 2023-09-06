\name{medinria}
\alias{MedINRIA}
\alias{medinria}
\alias{medinria2tensor}
\alias{tensor2medinria}
\title{ Read/Write Diffusion Tensor Data from/to NIFTI File }
\description{
  Read/Write diffusion tensor data from/to NIfTI file. Interface functions to MedINRIA.
}
\usage{
  medinria2tensor(filename)
  tensor2medinria(obj, filename, xind = NULL, yind = NULL, zind = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ file name for the tensor data. }
  \item{obj}{ object of class \code{"dtiTensor"} }
  \item{xind}{ index to define a subcube in x-direction. If \code{is.null(xind)} all voxel indices are used. }
  \item{yind}{ index to define a subcube in y-direction. If \code{is.null(yind)} all voxel indices are used. }
  \item{zind}{ index to define a subcube in z-direction. If \code{is.null(zind)} all voxel indices are used. }
}
\value{
  For function \code{medinria2tensor}: object of class \code{"dtiTensor"}.
}
\references{  P. Fillard, J. Souplet and N. Toussaint 
              \emph{ Medical Image Navigation and Research Tool by INRIA (MedINRIA)}, 
              INRIA Sophia Antipolis - Research Project ASCLEPIOS 2007

              \url{https://www-sop.inria.fr/asclepios/software/MedINRIA/}
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiTensor}}, 
  \code{\link{dtiTensor-methods}}
  \code{\link{dtiIndices-methods}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ IO }
