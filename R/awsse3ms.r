# This file contains the implementation of dti.smooth() for
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.ms <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth.ms", function(object, ...) standardGeneric("dwi.smooth.ms"))

setMethod("dwi.smooth.ms",
          "dtiData",
          function(object,
                   kstar,
                   lambda = 12,
                   kappa0 = .5,
                   ncoils = 1,
                   sigma = NULL,
                   ws0 = 1,
                   level = NULL,
                   mask = NULL,
                   xind = NULL,
                   yind = NULL,
                   zind = NULL,
                   verbose = FALSE,
                   usemaxni = TRUE,
                   memrelease = TRUE){

            # make the call part of the object
            args <- sys.call(-1)
            args <- c(object@call,args)

            # we need a lot of object properties
            ddim <- object@ddim
            ngrad <- object@ngrad
            s0ind <- object@s0ind
            ns0 <- length(s0ind)
            ngrad <- ngrad - ns0
            grad <- object@gradient[,-s0ind]
            bvalues <- object@bvalue[-s0ind]
            if(is.null(mask)) mask <- object@mask
            sdcoef <- object@sdcoef
            level <- object@level
            vext <- object@voxelext[2:3]/object@voxelext[1]
            #
            sb <- object@si[,,,-s0ind]
            s0 <- object@si[,,,s0ind]
            if(is.null(kappa0)){
              #  select kappa based on variance reduction on the sphere
              warning("You need to specify  kappa0  returning unsmoothed object")
              return(object)
            }
            if(ns0>1){
              dim(s0) <- c(prod(ddim),ns0)
              s0 <- s0%*%rep(1/sqrt(ns0),ns0)
              #  make sure linear combination of s0 has same variance as original
            }
            if (is.null(mask)) mask <- getmask(object, level)$mask
            nvoxel <- sum(mask)
# handle sigma
            dsigma <- dim(sigma)
            if (length(sigma) == 1) {
              sigma <- rep(sigma, nvoxel)
              if (verbose) cat("using supplied sigma ", sigma, "\n")
            } else if (identical(dsigma, ddim[1:3])) {
              sigma <- sigma[mask]
            } else if (length(dsigma==4)&identical(dsigma[1:3], ddim[1:3])) {
              dim(sigma) <- c(prod(ddim[1:3]),dsigma[4])
              sigma <- sigma[mask,]
              if (verbose) cat("using supplied array of sigma\n")
            } else {
# estimate from first s0
              sigma <- awslsigmc(object@si[ , , , object@s0ind[1]], 12,
                ncoils = ncoils, lambda = 5, verbose = verbose, hsig = 5, mask = mask)$sigma[mask]
                if (verbose) cat("estimated sigma image from first S0\n")
            }
            #
            #  reduce data to contain only voxel from brain mask
            #
            dim(sb) <- c(prod(ddim),ngrad)
            s0 <- s0[mask]
            sb <- sb[mask,]
            z <- aws::smse3ms(sb, s0, bvalues, grad, kstar, lambda, kappa0,
                            mask, sigma, ns0=ns0, ws0=ws0, vext=vext, ncoils=ncoils,
                            verbose=verbose, usemaxni=usemaxni)
            #
            #  one s0 image only
            #
            ngrad <- ngrad+1
            si <- array(0,c(prod(ddim),ngrad))
            si[mask,1] <- z$th0
            si[mask,-1] <- z$th
            dim(si) <- c(ddim,ngrad)
            object@si <-  si
            object@gradient <- grad <- cbind(c(0,0,0),grad)
            object@bvalue <- bvalue <- c(0,object@bvalue[-object@s0ind])
            object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
            object@s0ind <- as.integer(1)
            object@replind <- replind(grad)
            object@ngrad <- as.integer(ngrad)
            object@call <- args
            attr(object,"ns0") <- ns0
            object
          }
)
