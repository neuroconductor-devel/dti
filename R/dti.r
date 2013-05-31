################################################################
#                                                              #
# Section for DTI functions                                    #
#                                                              #
################################################################

dtiTensorold <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensorold", function(object,  ...) standardGeneric("dtiTensorold"))

setMethod( "dtiTensorold", "dtiData",
           function( object, 
                     method = c( "nonlinear", "linear"),
                     mc.cores = setCores( , reprt = FALSE)) {

             ## check method! available are: 
             ##   "linear" - use linearized model (log-transformed)
             ##   "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
             method <- match.arg( method)
                          
             if(is.null(mc.cores)) mc.cores <- 1
             mc.cores <- min(mc.cores,detectCores())
             args <- sys.call(-1)
             args <- c(object@call,args)
             ngrad <- object@ngrad
             grad <- object@gradient
             ddim <- object@ddim
             nvox <- prod(ddim)
             s0ind <- object@s0ind
             ns0 <- length(s0ind)
             sdcoef <- object@sdcoef
             require(parallel)
             if(all(sdcoef[1:4]==0)) {
               cat("No parameters for model of error standard deviation found\n estimating these parameters\n You may prefer to run sdpar before calling dtiTensor\n")
               sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
             }
             z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
             #
             #  this does not scale well with openMP
             #
             si <- array(z$si,c(ngrad,ddim))
             si[si<.5] <- .5 ## to avoid ceros in log-transform
             index <- z$index
             rm(z)
             gc()
             if(method=="linear"){
               ngrad0 <- ngrad - length(s0ind)
               s0 <- si[s0ind,,,]
               si <- si[-s0ind,,,]
               if(ns0>1) {
                 dim(s0) <- c(ns0,prod(ddim))
                 s0 <- rep(1/ns0,ns0)%*%s0
                 dim(s0) <- ddim
               }
               mask <- s0 > object@level
               mask <- connect.mask(mask)
               idsi <- 1:length(dim(si))
               ttt <- -log(sweep(si,idsi[-1],s0,"/"))
               ttt[is.na(ttt)] <- 0
               ttt[(ttt == Inf)] <- 0
               ttt[(ttt == -Inf)] <- 0
               dim(ttt) <- c(ngrad0,prod(ddim))
               cat("Data transformation completed ",format(Sys.time()),"\n")
               
               btbsvd <- svd(object@btb[,-s0ind])
               solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
               D <- solvebtb%*% ttt
               cat("Diffusion tensors generated ",format(Sys.time()),"\n")
               
               res <- ttt - t(object@btb[,-s0ind]) %*% D
               rss <- res[1,]^2
               for(i in 2:ngrad0) rss <- rss + res[i,]^2
               dim(rss) <- ddim
               sigma2 <- rss/(ngrad0-6)
               D[c(1,4,6),!mask] <- 1e-6
               D[c(2,3,5),!mask] <- 0
               D <- dti3Dreg(D,mc.cores=mc.cores)
               dim(D) <- c(6,ddim)
               dim(res) <- c(ngrad0,ddim)
               cat("Variance estimates generated ",format(Sys.time()),"\n")
               th0 <- array(s0,object@ddim)
               th0[!mask] <- 0
               gc()
             } else {
               #  method == "nonlinear" 
               ngrad0 <- ngrad
               s0 <- si[s0ind,,,]
               if(ns0>1) {
                 dim(s0) <- c(ns0,nvox)
                 s0 <- rep(1/ns0,ns0)%*%s0
                 dim(s0) <- ddim
               }
               mask <- s0 > object@level
               mask <- connect.mask(mask)
               df <- sum(table(object@replind)-1)
               if(mc.cores==1){
                 cat("start nonlinear regression",format(Sys.time()),"\n")
                 z <- .Fortran("nlrdtirg",
                               as.double(si),
                               as.integer(ngrad),
                               as.integer(nvox),
                               as.logical(mask),
                               as.double(object@btb),
                               as.double(sdcoef),
                               th0=as.double(s0),
                               D=double(6*nvox),
                               as.integer(200),
                               as.double(1e-6),
                               res=double(ngrad*nvox),
                               rss=double(nvox),
                               double(ngrad),
                               PACKAGE="dti",DUP=FALSE)[c("th0","D","res","rss")]
                 res <- z$res
                 D <- z$D
                 rss <- z$rss
                 th0 <- z$th0
               } else {
                 cat("start nonlinear regression using ",mc.cores, "cores", format(Sys.time()),"\n")
                 z <- matrix(0,8+ngrad,nvox)
                 dim(si) <- c(ngrad,nvox)
                 z[,mask] <- plmatrix(si[,mask],pnlrdtirg,btb=object@btb,sdcoef=sdcoef,s0ind=s0ind,
                                      ngrad = ngrad,mc.cores=mc.cores)
                 th0 <- z[1,]
                 D <- z[2:7,]
                 rss <- z[8,]
                 res <- z[8+(1:ngrad),]
               }
               indr <- NULL
               if(any(is.na(res))){
                 dim(res) <- c(ngrad,nvox)
                 indr <- (1:nvox)[apply(is.na(res),2,any)]
                 cat("(2) NA's in res in voxel",indr,"\n")
                 res[,indr] <- 0
               }
               indD0 <- indD1 <- NULL
               if(any(is.na(D))|any(abs(D)>1e10)){
                 dim(D) <- c(6,nvox)
                 indD0 <- (1:nvox)[apply(is.na(D),2,any)]
                 #     cat("NA's in D in ", length(indD0),"voxel:",indD0,"\n")
                 D[,indD0] <- c(1,0,0,1,0,1)
                 indD1 <- (1:nvox)[apply(abs(D)>1e10,2,any)]
                 #     cat("Inf's in D in", length(indD1)," voxel:",indD1,"\n")
                 D[,indD1] <- c(1,0,0,1,0,1)
               }
               dim(th0) <- ddim
               dim(D) <- c(6,ddim)
               dim(res) <- c(ngrad,ddim)
               dim(rss) <- ddim
               dim(mask) <- c(1, ddim)
               indD <- (1:nvox)[D[2,,,, drop=FALSE]==0&D[3,,,, drop=FALSE]==0&D[5,,,, drop=FALSE]==0&mask]
               if(!is.null(indr)) indD <- unique(indD,indr)
               if(!is.null(indD0)) indD <- unique(indD,indD0)
               if(!is.null(indD1)) indD <- unique(indD,indD1)
               dim(mask) <- ddim
               # this does not work in case of 2D data: 
               if(length(indD)>0){
                 cat("Problems in ",length(indD)," voxel, using fall back strategy \n ")
                 dim(si) <- c(ngrad,nvox)
                 dim(D) <- c(6,nvox)
                 dim(res) <- c(ngrad,nvox)
                 #     D[,unique(indDna,indDInf)] <- c(1,0,0,1,0,1)
                 #     res[,unique(indDna,indDInf)] <- 0
                 #     rss[unique(indDna,indDInf)] <- 1e10
                 #     cat("NA's in D in voxel",unique(indDna,indDInf),"\n")
                 if(mc.cores==1||length(indD)<2*mc.cores){
                   for(i in indD){
                     zz <- optim(c(1,0,0,1,0,1),opttensR,method="BFGS",si=si[-s0ind,i],s0=s0[i],grad=grad[,-s0ind],sdcoef=sdcoef)
                     D[,i] <- rho2D(zz$par)
                     th0[i] <- s0[i]
                     rss[i] <- zz$value
                     res[s0ind,i] <- 0
                     res[-s0ind,i] <- tensRres(zz$par,si[-s0ind,i],s0[i],grad[,-s0ind])
                   }
                 } else {
                   zz <- pmatrix(si[,indD],pnltens,grad=grad[,-s0ind],s0ind=s0ind,
                                 sdcoef=sdcoef,mc.cores=min(mc.cores,length(indD)))
                   dim(zz) <- c(length(zz)/length(indD),length(indD))
                   if(any(is.na(zz))){
                     ind <- (1:length(indD))[apply(is.na(zz),2,any)]
                     print(zz[,ind])
                     print(si[,indD[ind]])
                   }
                   D[,indD] <- zz[1:6,]
                   th0[indD] <- zz[7,]
                   rss[indD] <- zz[8,]
                   res[s0ind,indD] <- 0
                   res[-s0ind,indD] <- zz[-(1:8),]
                 }
                 if(any(is.na(res))){
                   dim(res) <- c(ngrad,nvox)
                   indr <- (1:nvox)[apply(is.na(res),2,any)]
                   #     cat("(3) NA's in res in voxel",indr,"\n")
                   res[,indr] <- 0
                 }
                 dim(D) <- c(6,ddim)
                 dim(res) <- c(ngrad,ddim)
               }
               cat("successfully completed nonlinear regression ",format(Sys.time()),"\n")
               sigma2 <- array(0,c(1,1,1))
               rm(z)
               gc()
             }
             #
             #   get spatial correlation
             #
             if(any(is.na(res))){
               dim(res) <- c(ngrad,nvox)
               indr <- (1:nvox)[apply(is.na(res),2,any)]
               cat("NA's in res in voxel",indr,"\n")
               res[,indr] <- 0
             }
             if(any(is.na(D))|any(abs(D)>1e10)){
               dim(D) <- c(6,nvox)
               indD <- (1:nvox)[apply(is.na(D),2,any)]
               cat("NA's in D in ", length(indD),"voxel:",indD,"\n")
               D[,indD] <- c(1,0,0,1,0,1)
               mask[indD] <- FALSE
               indD <- (1:nvox)[apply(abs(D)>1e10,2,any)]
               cat("Inf's in D in", length(indD)," voxel:",indD,"\n")
               D[,indD] <- c(1,0,0,1,0,1)
               mask[indD] <- FALSE
             }
             scorr <- mcorr(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=mc.cores)
             ev <- dti3Dev(D,mask,mc.cores=mc.cores)
             dim(ev) <- c(3,ddim)   
             dim(D) <- c(6,ddim)   
             scale <- quantile(ev[3,,,][mask],.95,na.rm=TRUE)
             cat("estimated scale information",format(Sys.time()),"\n")  
             invisible(new("dtiTensor",
                           call  = args,
                           D     = D,
                           th0   = th0,
                           sigma = sigma2,
                           scorr = scorr$scorr, 
                           bw = scorr$bw, 
                           mask = mask,
                           hmax = 1,
                           gradient = object@gradient,
                           bvalue = object@bvalue,
                           btb   = object@btb,
                           ngrad = object@ngrad, # = dim(btb)[2]
                           s0ind = object@s0ind,
                           replind = object@replind,
                           ddim  = object@ddim,
                           ddim0 = object@ddim0,
                           xind  = object@xind,
                           yind  = object@yind,
                           zind  = object@zind,
                           voxelext = object@voxelext,
                           level = object@level,
                           orientation = object@orientation,
                           rotation = object@rotation,
                           source = object@source,
                           outlier = index,
                           scale = scale,
                           method = method)
             )
           })

opttensR <- function(param,si,s0,grad,sdcoef){
     z<- .Fortran("opttensR",
               as.double(param),
               as.double(si),
               as.double(s0),
               as.double(grad),
               as.integer(length(si)),
               as.double(sdcoef),
               erg=double(1),
               DUP=FALSE,
               PACKAGE="dti")$erg
     if(is.na(z)) z <- 1e12
     z
}
tensRres <- function(param,si,s0,grad){
      .Fortran("tensRres",
               as.double(param),
               as.double(si),
               as.double(s0),
               as.double(grad),
               as.integer(length(si)),
               res=double(length(si)),
               DUP=FALSE,
               PACKAGE="dti")$res
}
rho2D <- function(param){
      .Fortran("rho2D0",
               as.double(param),
               D=double(6),
               DUP=FALSE,
               PACKAGE="dti")$D
}
D2rho <- function(D){
      .Fortran("D2rho0",
               as.double(D),
               rho=double(6),
               DUP=FALSE,
               PACKAGE="dti")$rho
}
#############

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, mc.cores=setCores(,reprt=FALSE)) {
  args <- sys.call(-1)
  args <- c(object@call,args)
  ddim <- object@ddim
  n <- prod(ddim)
  z <- dtiind3D(object@D,object@mask,mc.cores=mc.cores)
  invisible(new("dtiIndices",
                call = args,
                fa = array(z$fa,ddim),
                ga = array(z$ga,ddim),
                md = array(z$md,ddim),
                andir = array(z$andir,c(3,ddim)),
                bary = array(z$bary,c(3,ddim)),
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

