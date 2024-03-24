\name{extract-methods}
\docType{methods}
\alias{extract}
\alias{extract-methods}
\alias{extract,ANY-method}
\alias{extract,dtiData-method}
\alias{extract,dtiIndices-method}
\alias{extract,dtiTensor-method}
\alias{extract,dwiMixtensor-method}
\alias{extract,dwiQball-method}
\alias{[-methods}
\alias{[,ANY-method}
\alias{[,dtiData-method}
\alias{[,dtiIndices-method}
\alias{[,dtiTensor-method}
\alias{[,dkiIndices-method}
\alias{[,dkiTensor-method}
\alias{[,dwiMixtensor-method}
\alias{[,dwiQball-method}
\title{ Methods for Function `extract' and `[' in Package `dti' }
\description{
  The methods extract and/or compute specified statistics from object of class \code{"dtiData"}, \code{"dtiTensor"}, and \code{"dtiIndices"}. This can be restricted to a subset of voxel.
}
\usage{
  \S4method{extract}{dtiData}(x,
    what=c("data","gradient","btb","s0","sb","siq"),
    xind=TRUE, yind=TRUE, zind=TRUE)
  \S4method{extract}{dtiTensor}(x, what=c("tensor", "fa", "ga", "md", "evalues",
      "andir", "s0", "mask", "bic", "aic", "outlier"),
    xind=TRUE, yind=TRUE, zind=TRUE, mc.cores = setCores(, reprt = FALSE))
  \S4method{extract}{dwiMixtensor}(x, what=c("w0","andir", "order", "ev", "mix", "s0",
      "mask", "fa", "eorder", "bic", "aic"), xind=TRUE, yind=TRUE, zind=TRUE)
  \S4method{extract}{dtiIndices}(x, what=c("fa", "andir", "ga", "md", "bary"),
    xind=TRUE, yind=TRUE, zind=TRUE)
  \S4method{extract}{dwiQball}(x, what=c("sphcoef", "s0", "mask", "bic", "aic",
      "outlier"), xind=TRUE, yind=TRUE, zind=TRUE)
  \S4method{[}{dtiData}(x, i, j, k, drop=FALSE)
  \S4method{[}{dtiTensor}(x, i, j, k, drop=FALSE)
  \S4method{[}{dtiIndices}(x, i, j, k, drop=FALSE)
  \S4method{[}{dkiTensor}(x, i, j, k, drop=FALSE)
  \S4method{[}{dkiIndices}(x, i, j, k, drop=FALSE)
  \S4method{[}{dwiQball}(x, i, j, k, drop=FALSE)
}
\arguments{
  \item{x}{ Object of class \code{dti} }
  \item{i}{ vector of x-coordinates, defaults to whole range. }
  \item{j}{ vector of y-coordinates, defaults to whole range. }
  \item{k}{ vector of z-coordinates, defaults to whole range. }
  \item{xind}{ vector of x-coordinates, defaults to whole range. }
  \item{yind}{ vector of y-coordinates, defaults to whole range. }
  \item{zind}{ vector of z-coordinates, defaults to whole range. }
  \item{what}{ Statistic to extract. See Methods Section for details. }
  \item{drop}{ unused. }
   \item{mc.cores}{Number of cores to use. Defaults to number of threads specified for openMP,
   see documentation of package \pkg{awsMethods}. Our experience suggests to use 4-6 cores if available.}
}
\section{Methods}{
\describe{
  The generic extract function \code{"["} does what it is expected to do: it extracts parts of
  the object specified by \code{i}, \code{j}, and \code{k}.

  \item{x = "ANY"}{ Returns a warning for \code{extract}. Generic funtion for \code{"["}
  returns an object of same class with data clipped to the indices specified in arguments \code{i},
  \code{j} and \code{k}.}

  \item{x = "dtiData"}{ Extraction of squared gradient matrix \code{"btb"} or of S0 \code{"s0"},
  Sb \code{"sb"}, Si/mean(SO) \code{"siq"} or all images \code{"data"} restricted to the cube
  defined by arguments \code{i}, \code{j} and \code{k}.}

  \item{x = "dtiIndices"}{ Returns an array containing the specified statistics, i.e.
  fractional anisotropy \code{"fa"}, geodesic anisotropy \code{"ga"}, mean diffusivity \code{"md"},
  main direction of anisotropy \code{"andir"} and/or shape parameters \code{"bary"},
  as specified in argument \code{what}. Information is extracted for voxel within the
  cube defined by \code{xind}, \code{yind}, and \code{zind}. }

  \item{x = "dtiTensor"}{ Returns a list with component names corresponding to \code{what}
  containing the specified statistics, i.e. fractional anisotropy \code{"fa"},
  geodesic anisotropy \code{"ga"}, mean diffusivity \code{"md"}, eigenvalues \code{"evalues"},
  main direction of anisotropy \code{"andir"}, the tensor \code{"tensor"}
  the estimated S0 image \code{"s0"}, the values of the model selection criteia
  BIC \code{"bic"} or AIC \code{"aic"} and/or the mask used to restrict computations
  \code{"mask"}, as specified in argument \code{what}. Information is extracted
  for voxel within the  cube defined by arguments \code{xind}, \code{yind} and \code{zind}. }

  \item{x = "dwiMixtensor"}{ Returns a list with component names corresponding to
  \code{what} containing the specified statistics. Possible values for \code{what} are
  \code{"w0"} (size of isotropic department), \code{"order"} (estimated number of mixture components),
  \code{"eorder"} effective order), \code{"ev"} (eigenvalues), \code{"mix"} (mixture weights),
  \code{"andir"} (main directions of diffusion), \code{"fa"} (FA index), \code{"s0"}
  (the estimated S0 image), the values of the model selection criteia BIC \code{"bic"}
  or AIC \code{"aic"} and \code{mask} (the mask used to restrict computations).
  Information is extracted for voxel within the  cube defined by arguments \code{xind},
  \code{yind} and \code{zind}. }

  \item{x = "dwiQball"}{ Returns an array containing the specified statistics,
  the estimated coefficients with respect to the selected spherical harmonics basis
  \code{"sphcoef"}, the estimated S0 image \code{"s0"}, the values of the model selection
  criteia BIC \code{"bic"} or AIC \code{"aic"} and/or the mask used to restrict computations
  \code{"mask"}, as specified in argument \code{what}. Information is extracted for voxel
  within the  cube defined by arguments \code{xind}, \code{yind} and \code{zind}. }
}
}
\value{
  For function \code{extract} a list with components carrying the names of the options specified in
  argument \code{what}.  For \code{"["} the cutted object.
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\linkS4class{dtiData}},
  \code{\linkS4class{dtiTensor}},
  \code{\linkS4class{dtiIndices}}
  \code{\linkS4class{dwiMixtensor}},
  \code{\linkS4class{dwiQball}}
}
\keyword{ methods }
\keyword{ manip }
