%- Also NEED an '\alias' for EACH other topic documented here.
\name{setmask-methods}
\docType{methods}
\alias{setmask}
\alias{setmask-methods}
\alias{setmask,ANY-method}
\alias{setmask,dtiData-method}
\title{ Methods for Function `setmask' in Package `dti' }
\description{
   Read mask definition from NIfTI file and include it in dtiData object
}
\usage{
\S4method{setmask}{dtiData}(object, maskfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"dtiData"}}
  \item{maskfile}{NIfTI file containing mask definition. Dimension need to be 
  compatible, i.e. either equal \code{object@ddim0} or \code{object@ddim}}
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning }

  \item{obj = "dtiData"}{Set mask definition in dtiObject using information provided as NIfTI file as e.g. provided by \code{fsl_bet}. }

}
}
\value{
  The function returns an object of class \code{dtiData}.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}},
  \code{\link{readDWIdata}},
}
\keyword{ methods }
