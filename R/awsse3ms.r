# This file contains the implementation of dti.smooth() for 
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.ms <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth.ms", function(object, ...) standardGeneric("dwi.smooth.ms"))

setMethod("dwi.smooth.ms", "dtiData", function(object,kstar,lambda=15,kappa0=.9,ncoils=1,sigma=NULL,ws0=1,level=NULL,xind=NULL,yind=NULL,zind=NULL,verbose=FALSE,wghts=NULL){
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  level <- object@level
  vext <- object@voxelext[2:3]/object@voxelext[1]
  if(length(sigma)==1) {
     cat("using supplied sigma",sigma,"\n")
  } else {
  mask <- getmask(object,level)$mask
  sigma <- numeric(object@ngrad)
  for(i in 1:object@ngrad){
  sigma[i] <- awssigmc(object@si[,,,i],12,mask,ncoils,vext,h0=1.25,verbose=verbose)$sigma
  cat("image ",i," estimated sigma",sigma[i],"\n")
  }
 cat("quantiles of estimated sigma values",quantile(sigma),"\n")
  sigma <- median(sigma)
 cat("using median estimated sigma",sigma,"\n")
  }
#
  if(!(is.null(xind)&is.null(yind)&is.null(zind))){
  if(is.null(xind)) xind <- 1:object@ddim[1]
  if(is.null(yind)) yind <- 1:object@ddim[2]
  if(is.null(zind)) zind <- 1:object@ddim[3]
  object <- object[xind,yind,zind]
  }
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  bvalues <- object@bvalue[-s0ind]
  msstructure <- getnext3g(grad,bvalues)
  nshell <- msstructure$nbv
  sb <- object@si[,,,-s0ind]
  s0 <- object@si[,,,s0ind]
  if(is.null(kappa0)){
#  select kappa based on variance reduction on the sphere
     warning("You need to specify  kappa0  returning unsmoothed object")
     return(object)
   }
#
#  rescale so that we have Chi-distributed values
#
  sb <- sb/sigma
  s0 <- s0/sigma
  if(ns0>1){
     dim(s0) <- c(prod(ddim),ns0)
     s0 <- s0%*%rep(1/sqrt(ns0),ns0)
#  make sure linear combination of s0 has same variance as original 
     dim(s0) <- ddim
  }
#  th0 <- array(s0,c(ddim,nshell+1))
  ni0 <- array(1,ddim)
  mask <- s0>(sqrt(ns0)*level/sigma)
  gradstats <- getkappasmsh3(grad, msstructure)
#     save(gradstats,file="gradstats.rsc")
  hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
  kappa <- hseq$kappa
  cat("kappa",kappa,"\n")
  nind <- as.integer(hseq$n*1.25)
  hseqi <- hseq$h
  hseq0 <- hseq$h0 <- apply(hseqi,2,mean)
# make it nonrestrictive for the first step
  ni <- array(1,dim(sb))
  minlevel <- gamma(ncoils+0.5)/gamma(ncoils)*sqrt(2)
#  thats the mean of the central chi distribution with 2*ncoils df
  z <- list(th=sb, ni = ni, th0=s0, ni0=ni0)
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  kinit <- if(lambda<1e10) 1 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  for(k in kinit:kstar){
     gc()
     hakt <- hseqi[,k]
     hakt0 <- hseq0[k]
     thnimsh <- interpolatesphere0(z$th,z$th0,z$ni,z$ni0,msstructure,mask)
     param <- lkfullse3msh(hakt,kappa/hakt,gradstats,vext,nind) 
     param0 <- lkfulls0(hakt0,vext,nind) 
     z <- .Fortran("adsmse3c",
                as.double(sb),#y
                as.double(s0),#y0
                as.double(thnimsh$mstheta),#th
                as.double(thnimsh$msni),#ni
                as.double(thnimsh$msth0),#th0
                as.double(thnimsh$msni0),#ni0
                as.logical(mask),#mask
                as.integer(nshell+1),#ns number of shells
                as.integer(ddim[1]),#n1
                as.integer(ddim[2]),#n2
                as.integer(ddim[3]),#n3
                as.integer(ngrad),#ngrad
                as.double(lambda),#lambda
                as.double(ws0),# wghts0 rel. weight for s0 image
                as.integer(ncoils),#ncoils
                as.double(minlevel),#minlev 
                as.integer(mc.cores),#ncores
                as.integer(param$ind),#ind
                as.double(param$w),#w
                as.integer(param$n),#n
                as.integer(param0$ind),#ind0
                as.double(param0$w),#w0
                as.integer(param0$n),#n0
                th=double(prod(ddim)*ngrad),#thn
                ni=double(prod(ddim)*ngrad),#nin
                th0=double(prod(ddim)),#th0n
                ni0=double(prod(ddim)),#ni0n
                double(ngrad*mc.cores),#sw
                double(ngrad*mc.cores),#swy
                double((nshell+1)*mc.cores),#si
                double((nshell+1)*mc.cores),#thi
                double((nshell+1)*mc.cores),#nii
                DUPL=FALSE,
                PACKAGE="dti")[c("ni","th","ni0","th0")]
    dim(z$th) <- c(ddim,ngrad)
    dim(z$th0) <- c(ddim)
    dim(z$ni) <- c(ddim,ngrad)
    dim(z$ni0) <- c(ddim)
    gc()
if(verbose){
   dim(z$ni) <- c(prod(ddim),ngrad)
   cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni[mask,]),3),
  "mean of ni",signif(mean(z$ni[mask,]),3),"\n              quartiles of ni0",signif(quantile(z$ni0[mask]),3),
  "mean of ni0",signif(mean(z$ni0[mask]),3),
  " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
  } else {
      cat(".")
  }
  }
  ngrad <- ngrad+1
  si <- array(0,c(ddim,ngrad))
#
#  back to original scale
#
  si[,,,1] <-  pmax(z$th0/sqrt(ns0),minlevel)*sigma
#  go back to original s0 scale
  si[,,,-1] <- pmax(z$th,minlevel)*sigma
  object@si <-  si
  object@gradient <- grad <- cbind(c(0,0,0),grad)
  object@bvalue <- bvalue <- c(0,object@bvalue[-object@s0ind])
  object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
  object@s0ind <- as.integer(1)
  object@replind <- as.integer(1:ngrad)
  object@ngrad <- as.integer(ngrad)
  object@call <- args
  object
}
)
getkappas3 <- function(grad, trace = 0){
#
#  dist: acos(g_i%*%g_j)
#
  ngrad <- dim(grad)[2]
  kappa456 <- array(0, c(3, ngrad, ngrad))
  zbg <- betagamma(grad)
  for(i in 1:ngrad) kappa456[1,i,] <- acos(pmin(1,abs(grad[,i]%*%grad)))
  list(k456 = kappa456, bghat = zbg$bghat, dist=3)
}
getkappasmsh3 <- function(grad, msstructure, trace = 0){
ngrad <- dim(grad)[2]
nbv <- msstructure$nbv
bv <- msstructure$bv
ubv <- msstructure$ubv
bvind <- k456 <- bghat <- list(NULL)
for(i in 1:nbv){
#
#   collect information for spherical distances on each schell separately
#
ind <- (1:ngrad)[bv==ubv[i]]
z <- getkappas3(grad[,ind], trace = trace)
bvind[[i]] <- ind
k456[[i]] <- z$k456
bghat[[i]] <- z$bghat
}
list(k456 = k456, bghat = bghat, bvind = bvind, dist=3, nbv = nbv, ngrad = ngrad)
}
interpolatesphere0 <- function(theta,th0,ni,ni0,n3g,mask){
##  interpolate estimated thetas to get values on all spheres
##  n3g  generated by function  getnext3g
##  dim(theta) = c(n1,n2,n3,ngrad)
nbv <- n3g$nbv
bv <- n3g$bv
ubv <- n3g$ubv
dtheta <- dim(theta)
dth0 <- dim(th0)
ng <- dtheta[4]
mstheta <- msni <- array(0,c(nbv+1,prod(dtheta[1:3]),dtheta[4]))
dim(theta)  <- dim(ni)  <- c(prod(dtheta[1:3]),dtheta[4])
msth0 <- msni0 <- array(0,c(nbv+1,prod(dth0)))
# first component contains th0
#for(j in 1:dim(theta)[2]){
#   mstheta[1,,,,j] <- th0 
#   msni[1,,,,j] <- ni0
#}
mstheta[1,,] <- th0 
msni[1,,] <- ni0
# now the remaining shells
for(i in 1:nbv){
   for(j in 1:dim(theta)[2]){
      mstheta[i+1,mask,j] <- theta[mask,n3g$ind[i,j,]]%*%n3g$w[i,j,] 
#  correct value would be 
#  msni[1,,,,j] <- 1/(1/(ni[,n3g$ind[i,j,]])%*%(n3g$w[i,j,]^2))
#  try to be less conservative by ignorin squares in w
      msni[i+1,mask,j] <- 1/((1/ni[mask,n3g$ind[i,j,]])%*%n3g$w[i,j,])
   }
}
#  now fill vector for s0
msth0[1,] <- th0
msni0[1,] <- ni0
for(i in 1:nbv){
   indi <- (1:ng)[bv==ubv[i]]
   lindi <- length(indi)
   msth0[i+1,mask] <- theta[mask,indi]%*%rep(1/lindi,lindi)
#  correct value would be 
#  msni0[i+1,] <- 1/(1/(ni[,indi])%*%rep(1/lindi,lindi)^2)
#  try to be less conservative by ignorin squares in w
   msni0[i+1,mask] <- 1/((1/ni[mask,indi])%*%rep(1/lindi,lindi))   
}
dim(msth0) <- dim(msni0) <- c(nbv+1,dth0)
dim(mstheta) <- dim(msni) <- c(nbv+1,dtheta)
list(mstheta=mstheta,msni=msni,msth0=msth0,msni0=msni0)
}

lkfulls0 <- function(h,vext,n){
      z <- .Fortran("lkfuls0",
                    as.double(h),
                    as.double(vext),
                    ind=integer(3*n),
                    w=double(n),
                    n=as.integer(n),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("ind","w","n")]
      dim(z$ind) <- c(3,n)
list(h=h,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}

