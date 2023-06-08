.dtiOpts <- new.env(TRUE,emptyenv())

.onLoad <- function(lib, pkg){
  assign(".dtiopts", list(
               "swapx" = FALSE, "swapy" = FALSE,
               "swapz" = FALSE),
               envir = .dtiOpts)
   dti.options()
   invisible(NULL)
}

dti.options <- function(...){
  args <- list(...)
  imagepars <- get(".dtiopts", envir=.dtiOpts)
  if(length(args)>0){
  if("swapx" %in% names(args)) imagepars[["swapx"]] <- args[["swapx"]]
  if("swapy" %in% names(args)) imagepars[["swapy"]] <- args[["swapy"]]
  if("swapz" %in% names(args)) imagepars[["swapz"]] <- args[["swapz"]]
  assign(".dtiopts", imagepars, envir = .dtiOpts)
}
  invisible(imagepars)
}


sioutlier1 <- function( si, s0ind, level, mask, mc.cores = 1, verbose = TRUE){
  ##
  ##   replace si values that are larger than s0
  ##   create mask and reduce data to region covered by mask
  ##
  dsi <- dim(si)
  n <- prod(dsi[-length(dsi)])
  ng <- dsi[length(dsi)]
  ns0 <- length(s0ind)
  siind <- (1:ng)[-s0ind]
  dim(si) <- c(n,ng)
  si <- t(si)

  if (verbose) cat("outlier: ")

  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  t1 <- Sys.time()
  if(mc.cores==1||ng>250){
    z <- .Fortran(C_outlier,
                  as.double(si),
                  as.integer(n),
                  as.integer(ng),
                  as.integer(s0ind),
                  as.integer(siind),
                  as.integer(ns0),
                  si=double(n*ng),
                  index=integer(n))[c("si","index")]
    zz <- matrix(z$si,ng,n)
    index <- (1:n)[z$index]
    rm(z)
  } else {
    zz <- matrix(.Fortran(C_outlierp,
                          as.double(si),
                          as.integer(n),
                          as.integer(ng),
                          as.integer(s0ind),
                          as.integer(ns0),
                          as.integer(siind),
                          as.integer(ng-ns0),
                          si=double(n*(ng+1)),
                          as.integer(ng+1))$si,ng+1,n)
    t2 <- Sys.time()
    if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
    index <- (1:n)[as.logical(zz[ng+1,])]
    zz <- zz[1:ng,]
  }
  ng0 <- length(siind)
  s0 <- zz[s0ind,]
  si <- zz[-s0ind,]
  if(ns0>1) {
    s0 <- rep(1/ns0,ns0)%*%s0
  }
  if(is.null(mask)){
    mask <- array(s0 > level,dsi[-length(dsi)])
    mask <- connect.mask(mask)
  }
  nvox <- sum(mask)

  t2 <- Sys.time()
  if (verbose) cat( difftime( t2, t1), attr(difftime( t2, t1), "units"), "for", nvox, "voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  list(si=zz[,mask],s0=s0[mask],index=index,mask=mask)
}
sioutlier <- function( si, s0ind, mc.cores = 1, verbose = TRUE){
  dsi <- dim(si)
  n <- prod(dsi[-length(dsi)])
  ng <- dsi[length(dsi)]
  ns0 <- length(s0ind)
  siind <- (1:ng)[-s0ind]
  dim(si) <- c(n,ng)
  si <- t(si)

  if (verbose) cat("outlier: ")

  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  t1 <- Sys.time()
  if(mc.cores==1||ng>250){
    z <- .Fortran(C_outlier,
                  as.double(si),
                  as.integer(n),
                  as.integer(ng),
                  as.integer(s0ind),
                  as.integer(siind),
                  as.integer(ns0),
                  si=double(n*ng),
                  index=integer(n))[c("si","index")]
  } else {
    zz <- matrix(.Fortran(C_outlierp,
                          as.double(si),
                          as.integer(n),
                          as.integer(ng),
                          as.integer(s0ind),
                          as.integer(ns0),
                          as.integer(siind),
                          as.integer(ng-ns0),
                          si=double(n*(ng+1)),
                          as.integer(ng+1))$si,ng+1,n)
    t2 <- Sys.time()
    if (verbose) cat( difftime( t2, t1), attr(difftime( t2, t1), "units"), "for", n, "voxel\n")
    if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
    return(list(si=zz[1:ng,],index=(1:n)[as.logical(zz[ng+1,])]))
  }
  t2 <- Sys.time()
  if (verbose) cat( difftime( t2, t1), attr(difftime( t2, t1), "units"), "for", n, "voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  index <- (1:n)[z$index]
  dim(z$si) <- c(ng,dsi[-length(dsi)])
  list(si=z$si,index=index)
}

mcorr <- function(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=1){
  cat("mcorr:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(min(mc.cores,lags[1]))
  }
  t1 <- Sys.time()
  scorr <- .Fortran(C_mcorr,as.double(res),
                    as.integer(mask),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.integer(ngrad0),
                    double(prod(ddim)),
                    double(prod(ddim)),
                    scorr = double(prod(lags)),
                    as.integer(lags[1]),
                    as.integer(lags[2]),
                    as.integer(lags[3]))$scorr
  t2 <- Sys.time()
  cat(difftime(t2,t1),"\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  dim(scorr) <- lags
  scorr[is.na(scorr)] <- 0
  cat("estimated spatial correlations",format(Sys.time()),"\n")
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")
  cat("thcorr:")
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.2,.2,.2),
              upper=c(3,3,3),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0
  cat("estimated corresponding bandwidths",format(Sys.time()),"\n")
  list(scorr=scorr,bw=bw)
}

dti3Dreg <- function(D,mc.cores=1){
  nvox <- length(D)/6
  cat("dti3Dreg:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  t1 <- Sys.time()
  D <- .Fortran(C_dti3dreg,
                D=as.double(D),
                as.integer(nvox))$D
  t2 <- Sys.time()
  cat(difftime(t2,t1)," for",nvox,"voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  matrix(D,6,nvox)
}
dti3Dev <- function(D,mask,mc.cores=1){
  dimD <- dim(D)[-1]
  nvox <- prod(dimD)
  nvox0 <- sum(mask)
  dim(D) <- c(6,nvox)
  cat("dti3Dev:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  ev <- matrix(0,3,nvox)
  t1 <- Sys.time()
  ev[,mask] <- .Fortran(C_dti3dev,
                        as.double(D[,mask]),
                        as.integer(nvox0),
                        ev=double(3*nvox0))$ev
  t2 <- Sys.time()
  cat(difftime(t2,t1)," for",nvox0,"voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  dim(ev) <- c(3,dimD)
  ev
}
dti3Dand <- function(D,mask,mc.cores=1){
  dimD <- dim(D)[-1]
  nvox <- prod(dimD)
  nvox0 <- sum(mask)
  dim(D) <- c(6,nvox)
  cat("dti3Dand:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  andir <- matrix(0,3,nvox)
  t1 <- Sys.time()
  andir[,mask] <- .Fortran(C_dti3dand,
                           as.double(D[,mask]),
                           as.integer(nvox0),
                           andir=double(3*nvox0))$andir
  t2 <- Sys.time()
  cat(difftime(t2,t1)," for",nvox0,"voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  dim(andir) <- c(3,dimD)
  andir
}
dti3Dall <- function(D,mask,mc.cores=1){
  dimD <- dim(D)[-1]
  nvox <- prod(dimD)
  nvox0 <- sum(mask)
  dim(D) <- c(6,nvox)
  cat("dti3Dall:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  ev <- andir <- matrix(0,3,nvox)
  fa <- ga <- md <- numeric(nvox)
  t1 <- Sys.time()
  z <- .Fortran(C_dti3dall,
                as.double(D[,mask]),
                as.integer(nvox0),
                fa=double(nvox0),
                ga=double(nvox0),
                md=double(nvox0),
                andir=double(3*nvox0),
                ev=double(3*nvox0))[c("fa","ga","md","andir","ev")]
  t2 <- Sys.time()
  cat(difftime(t2,t1)," for",nvox0,"voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  fa[mask] <- z$fa
  ga[mask] <- z$ga
  md[mask] <- z$md
  andir[,mask] <- z$andir
  ev[,mask] <- z$ev
  list(fa=fa,ga=ga,md=md,andir=andir,ev=ev)
}
dtieigen <- function(D,mask,mc.cores=1){
  dimD <- dim(D)[-1]
  nvox <- prod(dimD)
  nvox0 <- sum(mask)
  dim(D) <- c(6,nvox)
  cat("dtieigen:")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
  }
  ev <- matrix(0,3,nvox)
  andir <- matrix(0,6,nvox)
  fa <-numeric(nvox)
  t1 <- Sys.time()
  z <- .Fortran(C_dtieigen,
                as.double(D[,mask]),
                as.integer(nvox0),
                fa=double(nvox0),
                ev=double(3*nvox0),
                andir=double(6*nvox0))[c("fa","ev","andir")]
  t2 <- Sys.time()
  cat(difftime(t2,t1)," for",nvox0,"voxel\n")
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  fa[mask] <- z$fa
  andir[,mask] <- z$andir
  ev[,mask] <- z$ev
  list(fa=fa,ev=ev,andir=andir)
}
dtiind3D <- function( D, mask, mc.cores = 1, verbose = TRUE){
  dimD <- dim(D)[-1]
  nvox <- prod(dimD)
  nvox0 <- sum(mask)
  dim(D) <- c(6,nvox)
  if (verbose) cat( "dtiind3: entering function", format( Sys.time()), "\n")
  if(mc.cores>1){
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores,reprt=verbose)
  }
  bary <- andir <- matrix(0,3,nvox)
  fa <- ga <- md <- numeric(nvox)
  t1 <- Sys.time()
  z <- .Fortran(C_dtiind3d,
                as.double(D[ , mask]),
                as.integer(nvox0),
                fa      = double(nvox0),
                ga      = double(nvox0),
                md      = double(nvox0),
                andir   = double(3*nvox0),
                bary    = double(3*nvox0))[c( "fa", "ga", "md", "andir", "bary")]
  if(mc.cores>1) setCores(mc.cores.old,reprt=FALSE)
  t2 <- Sys.time()
  if ( verbose) cat( "dtiind3: calculation took ", difftime( t2, t1), attr(difftime( t2, t1), "units"), " for", nvox0, "voxel\n")
  fa[mask] <- z$fa
  ga[mask] <- z$ga
  md[mask] <- z$md
  andir[,mask] <- z$andir
  bary[,mask] <- z$bary
  if (verbose) cat( "dtiind3: exiting function", format( Sys.time()), "\n")
  list( fa = fa, ga = ga, md = md, andir = andir, bary = bary)
}
kldist <- function(L,eta1,eta2){
  #
  #  L -number of coils
  #  eta1 - m1_1/sigma
  #  eta2 - m1_2/sigma
  n1 <- length(eta1)
  n2 <- length(eta2)
  f1 <- (2*L+eta1^2)^2/(2*L+2*eta1^2)
  c1 <- (2*L+2*eta1^2)/(2*L+eta1^2)
  f2 <- (2*L+eta2^2)^2/(2*L+2*eta2^2)
  c2 <- (2*L+2*eta2^2)/(2*L+eta2^2)
  lGf1 <- lgamma(f1/2)
  lGf2 <- lgamma(f2/2)
  flc1 <- f1/2*log(c1)
  flc2 <- f2/2*log(c2)
  psif1 <- digamma(f1/2)-outer(lGf1,lGf2,"-")-outer(flc1,flc2,"-")+outer(f1,f2,"-")/2*outer(log(c1)+psif1,rep(1,n2),"*")+
    (outer(c1,c2,"/")-1)*outer(f1,rep(1/2,n2),"*")
}


fwhm2bw <- function(hfwhm) hfwhm/sqrt(8*log(2))

replind <- function(gradient){
  #
  #  determine replications in the design that may be used for
  #  variance estimates
  #
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  replind <- numeric(ngrad)
  while(any(replind==0)){
    i <- (1:ngrad)[replind==0][1]
    ind <- (1:ngrad)[apply(abs(gradient-gradient[,i]),2,max)==0]
    replind[ind] <- i
  }
  as.integer(replind)
}



Spatialvar.gauss<-function(h,h0,d,interv=1){
  #
  #   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in
  #
  #   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
  #
  #   Spatialvar.gauss(lkern,h,h0,d)/Spatialvar.gauss(lkern,h,1e-5,d) gives the
  #   a factor for lambda to be used with bandwidth h
  #
  #
  #  interv allows for further discretization of the Gaussian Kernel, result depends on
  #  interv for small bandwidths. interv=1  is correct for kernel smoothing,
  #  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding
  #  discretisation into voxel)
  #
  h0 <- pmax(h0,1e-5)
  h <- pmax(h,1e-5)
  h<-h/2.3548*interv
  if(length(h)==1) h<-rep(h,d)
  ih<-trunc(4*h)
  ih<-pmax(1,ih)
  dx<-2*ih+1
  penl<-dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  dim(penl)<-dx
  h0<-h0/2.3548*interv
  if(length(h0)==1) h0<-rep(h0,d)
  ih<-trunc(4*h0)
  ih<-pmax(1,ih)
  dx0<-2*ih+1
  x<- ((-ih[1]):ih[1])/h0[1]
  penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
  if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
  if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
  dim(penl0)<-dx0
  penl0<-penl0/sum(penl0)
  dz<-dx+dx0-1
  z<-array(0,dz)
  if(d==1){
    for(i1 in 1:dx0) {
      ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
      ind1<-ind1[ind1<=dz][-1]
      z[-ind1]<-z[-ind1]+penl*penl0[i1]
    }
  } else if(d==2){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
    }
  } else if(d==3){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]) for(i3 in 1:dx0[3]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
      ind3<-ind3[ind3<=dz[3]][-1]
      z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
    }
  }
  sum(z^2)/sum(z)^2*interv^d
}

Varcor.gauss<-function(h){
  #
  #   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
  #
  #   in case of colored noise that was produced by smoothing with lkern and bandwidth h
  #
  h<-pmax(h/2.3548,1e-5)
  ih<-trunc(4*h)+1
  dx<-2*ih+1
  d<-length(h)
  penl <- dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl <- outer(penl,dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl <- outer(penl,outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  2*sum(penl)^2/sum(diff(penl)^2)
}


corrrisk <- function(bw,lag,data){
  z <- thcorr3D(bw,lag)
  mean((data-z)^2/outer(outer(1:lag[1],1:lag[2],"*"),1:lag(3),"*"))
}

thcorr3D <- function(bw,lag=rep(5,3)){
  g <- trunc(fwhm2bw(bw)*4)
  gw1 <- dnorm(-(g[1]):g[1],0,fwhm2bw(bw[1]))
  gw2 <- dnorm(-(g[2]):g[2],0,fwhm2bw(bw[2]))
  gw3 <- dnorm(-(g[3]):g[3],0,fwhm2bw(bw[3]))
  gwght <- outer(gw1,outer(gw2,gw3,"*"),"*")
  gwght <- gwght/sum(gwght)
  dgw <- dim(gwght)
  scorr <- .Fortran(C_thcorr,
                    as.double(gwght),
                    as.integer(dgw[1]),
                    as.integer(dgw[2]),
                    as.integer(dgw[3]),
                    scorr=double(prod(lag)),
                    as.integer(lag[1]),
                    as.integer(lag[2]),
                    as.integer(lag[3]))$scorr
  # bandwidth in FWHM in voxel units
  dim(scorr) <- lag
  scorr
}

andir2.image <- function(dtobject,slice=1,method=1,quant=0,minfa=NULL,show=TRUE,xind=NULL,yind=NULL,...){
  if(!inherits(dtobject,"dti")) stop("Not an dti-object")
  if(is.null(dtobject$anindex)) stop("No anisotropy index yet")
  #adimpro <- require(adimpro)
  anindex <- dtobject$anindex
  dimg <- dim(anindex)[1:2]
  if(is.null(xind)) xind <- 1:dimg[1]
  if(is.null(yind)) yind <- 1:dimg[2]
  if(is.null(slice)) slice <- 1
  anindex <- anindex[xind,yind,slice]
  dimg <- dim(anindex)[1:2]
  andirection <- dtobject$andirection[,xind,yind,slice]
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(3,prod(dimg))
  if(is.null(minfa)) minfa <- quantile(anindex,quant)
  if(method==1) {
    andirection[1,] <- abs(andirection[1,])
    andirection[2,] <- abs(andirection[2,])
    andirection[3,] <- abs(andirection[3,])
  } else {
    ind<-andirection[1,]<0
    andirection[,ind] <- - andirection[,ind]
    andirection[2,] <- (1+andirection[2,])/2
    andirection[3,] <- (1+andirection[3,])/2
  }
  andirection <- t(andirection)
  andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minfa)
  dim(andirection)<-c(dimg,3)
  #  if(adimpro) {
  andirection <- adimpro::make.image(andirection)
  if(show) adimpro::show.image(andirection,...)
  #  } else if(show) {
  #    dim(anindex) <- dimg
  #    image(anindex,...)
  #  }
  invisible(andirection)
}

andir.image <- function(anindex,andirection,quant=0,minfa=NULL){
  dimg <- dim(anindex)
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(3,prod(dimg))
  if(is.null(minfa)) minfa <- quantile(anindex,quant)
  andirection[1,] <- abs(andirection[1,])
  andirection[2,] <- abs(andirection[2,])
  andirection[3,] <- abs(andirection[3,])
  andirection <- t(andirection)*as.vector(anindex)*as.numeric(anindex>minfa)
  dim(andirection)<-c(dimg,3)
  adimpro::show.image(adimpro::make.image(andirection))
  invisible(NULL)
}

connect.mask <- function(mask){
  dm <- dim(mask)
  n1 <- dm[1]
  n2 <- dm[2]
  n3 <- dm[3]
  n <- n1*n2*n3
  mask1 <- .Fortran(C_lconnect,
                    as.integer(mask),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.integer((n1+1)/2),
                    as.integer((n2+1)/2),
                    as.integer((n3+1)/2),
                    integer(n),
                    integer(n),
                    integer(n),
                    mask=integer(n))$mask
  mask1 <- as.logical(mask1)
  dim(mask1) <- dm
  mask1
}

sphcoord <- function(ccoord){
  #
  #  transform cartesian into sherical coordinates
  #
  ccoord <- ccoord/sqrt(sum(ccoord^2))
  phi <- atan2(ccoord[2],ccoord[1])+2*pi*(ccoord[2]<0)
  theta <- atan2(sqrt(ccoord[2]^2+ccoord[1]^2),ccoord[3])
  c(theta,phi)
}


create.designmatrix.dti <- function(gradient) {
  dgrad <- dim(gradient)
  if (dgrad[2]==3) gradient <- t(gradient)
  dgrad <- dim(gradient)
  if (dgrad[1]!=3) stop("Not a valid gradient matrix")

  btb <- matrix(0,6,dgrad[2])
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]

  btb
}

identifyFA <- function(view,slice,xind,yind,zind){
  n1 <- switch(view,"sagittal"=length(yind),length(xind))
  n2 <- switch(view,"axial"=length(yind),length(zind))
  x <- as.vector(outer(1:n1,rep(1,n2),"*"))
  y <- as.vector(outer(rep(1,n1),1:n2,"*"))
  cat("Please use left mouse click to identify a voxel,\n terminate selection process by right mouse click\n")
  z <- identify(x,y,plot=FALSE)
  coord <- matrix(0,3,length(z))
  if(view=="sagittal"){
    coord[1,] <- slice
    coord[2,] <- yind[x[z]]
    coord[3,] <- zind[y[z]]
  } else if(view=="coronal") {
    coord[1,] <- xind[x[z]]
    coord[2,] <- slice
    coord[3,] <- zind[y[z]]
  } else {
    coord[1,] <- xind[x[z]]
    coord[2,] <- yind[y[z]]
    coord[3,] <- slice
  }
  coord
}


vcrossp <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

showFAColorScale <- function(filename = "FAcolorscale.png") {
  data("colqFA", envir = environment())
  png( filename = filename, width = 600, height = 75, bg = "white", pointsize = 12)
  par( mar = c( 2, 0.75, 0.1, 0.75))
  image( matrix( seq( 0, 1, length = 256), 256, 1), col = dti::colqFA, yaxt = "n")
  axis(1, at = seq( 0, 1, by = 0.1))
  text( 0.1, 0, "FA", pos = 4, cex = 2, font = 2, col = "white")
  dev.off()
}
hg1f1 <- function(a,b,z){
##
##  Confluent Hypergeometric 1F1 (a,b scalar, z vector)
##  rel accuracy 1e-13 for z in -1400:700 for a=-.5, .5
##  rel accuracy 2e-4 for z < -1400 for a=-.5, .5
##
   n <- length(z)
   .Fortran(C_hg1f1,
            as.double(a),
            as.double(b),
            as.double(z),
            as.integer(n),
            fz=double(n))$fz
}

unifybvals <- function(bval,dbv=51){
   nbv <- length(bval)
   nbval <- bval
   obval <- numeric(nbv)
   while(any(nbval!=obval)){
      obval <- nbval
      sbv <- sort(obval)
      dsbv <- (1:(nbv-1))[diff(sbv)<dbv]
      sbv[dsbv+1] <- sbv[dsbv]
      obv <- order(obval)
      nbval[obv] <- sbv
   }
   for(bv in unique(nbval)) nbval[nbval==bv] <- trunc(mean(bval[nbval==bv]))
   nbval
}

expanddwiobj <- function(object){
   if(!inherits(object,"dtiData")) stop("Needs dtiData object")
   ns0 <- attributes(object)$ns0
   if(is.numeric(ns0)) if(ns0>1&&object@s0ind==1){
     grad <- cbind(matrix(0,3,ns0),object@gradient[,-1])
     bvalue <- c(rep(0,ns0-1),object@bvalue)
     ngrad <- object@ngrad+ns0-1
     s0 <- object@si[,,,object@s0ind[1]]
     si <- object@si[,,,-object@s0ind]
     sinew <- array(0,c(dim(s0),ngrad))
     sinew[,,,-(1:ns0)] <- si
     sinew[,,,1:ns0] <- s0
     object@si <- sinew
     object@gradient <- grad 
     object@bvalue <- bvalue
     object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
     object@s0ind <- as.integer(1:ns0)
     object@replind <- replind(grad)
     object@ngrad <- as.integer(ngrad)
     attr(object,"ns0") <- NULL
   }
     object
}