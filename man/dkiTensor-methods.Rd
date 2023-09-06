\name{dkiTensor-methods}
\docType{methods}
\alias{dkiTensor}
\alias{dkiTensor-methods}
\alias{dkiIndices-methods}
\alias{dkiIndices}
\alias{dkiTensor,ANY-method}
\alias{dkiIndices,ANY-method}
\alias{dkiTensor,dtiData-method}
\alias{dkiIndices,dkiTensor-method}
\title{Diffusion Kurtosis Imaging (DKI)}
\description{
 These methods estimate, in each voxel, the diffusion kurtosis tensor (and the diffusion tensor) and some scalar indices.
}
\usage{
  \S4method{dkiTensor}{dtiData}(object, method=c("CLLS-QP", "CLLS-H", "ULLS", "QL", "NLR"),
                   sigma=NULL, L=1, mask=NULL,
                   mc.cores=setCores(, reprt=FALSE), verbose=FALSE)
  \S4method{dkiIndices}{dkiTensor}(object, mc.cores=setCores(, reprt=FALSE),
                   verbose=FALSE)
}
\arguments{
  \item{object}{Object of class \code{"dtiData"}}
  \item{method}{Method for tensor estimation. May be \code{"CLLS-QP"} for a qudratic
programm solution for the constrained optimization (requires package quadprog),
\code{"CLLS-H"} for a heuristic approximation described in Tabesh et al. (2011),
or \code{"ULLS"} for an unconstrained linear least squares estimation. \code{"QL"} and
\code{"NLR"} correspond to the use of unconstrained quasi-likelihood and nonlinear regression,
respectively. }
\item{sigma}{Scale parameter of intensity distribution (unprocessed). Used with \code{method="QL"}
   in the calculation of the expected intensity values.
}
\item{L}{Effective number of coils, 2*L are the degrees of freedom of the intensity
   distribution (unprocessed). The default corresponds, e.g., to the case of a SENSE reconstruction.
   Used with \code{method="QL"}
   in the calculation of the expected intensity values.
}
   \item{mask}{argument to specify a precomputed brain mask}
   \item{mc.cores}{Number of cores to use. Defaults to number of threads specified for openMP, see documentation of package \pkg{awsMethods}. Not yet fully implemented for
  these methods.}
   \item{verbose}{Verbose mode.}
}
\value{
  An object of class \code{"dkiTensor"} or \code{"dkiIndices"}.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{Returns a warning}

\item{\code{signature(object = "dtiData")}}{The method \code{"dkiTensor"} estimates the diffusion kurtosis
    model, i.e., the kurtosis tensor and the diffusion tensor.}

\item{\code{signature(object = "dkiTensor")}}{The method \code{"dkiIndices"} estimates
some scalar indices from the kurtosis tensor. The method is still experimental, some
quantities may be removed in future versions, other might be included.}

}}
\references{A. Tabesh, J.H. Jensen, B.A. Ardekani, and J.A. Helpern,
             \emph{Estimation of tensors and tensor-derived measures in diffusional
                   kurtosis imaging},
             Magnetic Resonance in Medicine, 65, 823-836 (2011).

            E.S. Hui, M.M. Cheung, L. Qi, and E.X. Wu,
            \emph{Towards better MR characterization of neural tissues using directional
            diffusion kurtosis analysis},
            Neuroimage, 42, 122-134 (2008).

            J. Polzehl, K. Tabelow (2019). Magnetic Resonance Brain Imaging:
            Modeling and Data Analysis Using R. Springer, Use R! series. Doi:10.1007/978-3-030-29184-6.

            \url{ https://www.wias-berlin.de/projects/matheon_a3/}
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}},
  \code{\link{readDWIdata}},
  \code{\linkS4class{dtiData}},
  \code{\linkS4class{dkiTensor}}
  \code{\linkS4class{dkiIndices}}
}
\keyword{methods}
\keyword{model}
