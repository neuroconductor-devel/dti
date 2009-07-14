################################################################
#                                                              #
# Section for plot() functions                                 #
#                                                              #
################################################################

setMethod("plot", "dwi", function(x, y, ...) cat("No implementation for class dwi\n"))

##############

setMethod("plot", "dtiData", 
function(x, y,slice=1, gradient=NULL, view= "axial", show=TRUE, density=FALSE, xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), ...) {
  if(is.null(x@si)) cat("No dwi data yet")
  maxsi <- max(x@si)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if(is.null(gradient)) gradient <- x@s0ind[1]
  if(gradient<1||gradient>x@ngrad) {
    warning("gradient number out of range, show s0 image")
    gradient <- x@s0ind[1]
  }
  if(density) { 
    z <- density(x@si[xind,yind,zind,gradient])
    if(show) {
      plot(z,main="Density of S0-values") 
      lines(c(x@level,x@level),c(0,max(z$y)),col=2)
    }
    return(invisible(z))
  }
  adimpro <- require(adimpro)
  if (view == "sagittal") {
    if(slice<1||slice>x@ddim[1]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[1]%/%2
    }
    img <- x@si[slice,yind,zind,gradient]
  } else if (view == "coronal") {
    if(slice<1||slice>x@ddim[2]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[2]%/%2
    }
    img <- x@si[xind,slice,zind,gradient]
  } else {
    if(slice<1||slice>x@ddim[3]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[3]%/%2
    }
    img <- x@si[xind,yind,slice,gradient]
  }
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if(adimpro) {
    img <- make.image(65535*img/maxsi)
    if(show) show.image(img,...)
  } else if(show) {
    image(img,...)
  }
  par(oldpar)
  invisible(img)
})

##############

setMethod("plot", "dtiTensor", function(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1,what="FA", qrange=c(.01,.99),xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...) {
  if(is.null(x@D)) cat("No diffusion tensor yet")
  adimpro <- require(adimpro)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if (view == "sagittal") {
    D <- x@D[,slice,yind,zind]
    mask <- x@mask[slice,yind,zind]
    n1 <- length(yind)
    n2 <- length(zind)
  } else if (view == "coronal") {
    D <- x@D[,xind,slice,zind]
    mask <- x@mask[xind,slice,zind]
    n1 <- length(xind)
    n2 <- length(zind)
  } else {
    D <- x@D[,xind,yind,slice]
    mask <- x@mask[xind,yind,slice]
    n1 <- length(xind)
    n2 <- length(yind)
  }
  if(what=="GA"){
  z <- .Fortran("dti2Dga",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.logical(mask),
                fa=double(n1*n2),
                md=double(n1*n2),
                andir=double(3*n1*n2),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir")]
  } else {
  z <- .Fortran("dti2Dfa",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.logical(mask),
                fa=double(n1*n2),
                md=double(n1*n2),
                andir=double(3*n1*n2),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir")]
   }
   oldpar <- par(mfrow=c(3,3),mar=mar,mgp=mgp,...)
#  now draw information to graphical device
   on.exit(par(oldpar))
   img<-D[1,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dxx: mean",signif(mean(D[mask]),3),"max",signif(max(D[1,,][mask]),3)))
   img<-D[2,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxy: min",signif(min(D[2,,][mask]),3),"max",signif(max(D[2,,][mask]),3)))
   img<-D[3,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxz: min",signif(min(D[3,,][mask]),3),"max",signif(max(D[3,,][mask]),3)))
   show.image(make.image(matrix(z$fa,n1,n2)))
   if(what=="GA"){
   title(paste("Anisotropy Index  (GA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
   } else {
   title(paste("Geodesic Anisotropy (FA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
   }
   img<-D[4,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dyy: min",signif(min(D[4,,][mask]),3),"max",signif(max(D[4,,][mask]),3)))
   img<-D[5,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dyz: min",signif(min(D[5,,][mask]),3),"max",signif(max(D[5,,][mask]),3)))
   andir.image(matrix(z$fa,n1,n2),array(z$andir,c(3,n1,n2)),quant=quant,minanindex=minanindex)
   title(paste("Anisotropy directions"))
   img <- matrix(z$md,n1,n2)
   show.image(make.image(65535*img/max(img)))
   if(what=="GA"){
   title(paste("Mean log diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
   } else {
   title(paste("Mean diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
   }
   img<-D[6,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dzz: min",signif(min(D[6,,][mask]),3),"max",signif(max(D[6,,][mask]),3)))
   invisible(NULL)
})

##############

setMethod("plot", "dtiIndices", 
function(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, identify=FALSE, density=FALSE, contrast.enh=1,what="FA",xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
  if(is.null(x@fa)) cat("No anisotropy index yet")
  if(!(method %in% 1:5)) {
      warning("method out of range, reset to 1")
      method <- 1
  }
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if(density) { 
   x <- x[xind,yind,zind]
   z <- density(if(what=="FA") x@fa[x@fa>0] else x@ga[x@ga>0])
   if(show) {
      plot(z,main=paste("Density of positive",what,"-values")) 
   }
   return(invisible(z))
  }
  adimpro <- require(adimpro)
  oldpar <- par(mar=mar,mgp=mgp, ...)
#  if(what=="GA") maxga <- max(x@ga) 
#  if(what=="GA") maxga <- quantile(x@ga,0.99) 
#  resulting image needs to be rescaled 
  if (view == "sagittal") {
#    anindex <- if(what=="GA") pmin(x@ga[slice,yind,zind]/maxga, 1) else x@fa[slice,yind,zind]
    anindex <- if(what=="GA") tanh(x@ga[slice,yind,zind]) else x@fa[slice,yind,zind]
    if (method == 3) {
      andirection <- x@bary[,slice,yind,zind]
    } else {
      andirection <- x@andir[,slice,yind,zind]
    }
  } else if (view == "coronal") {
#    anindex <- if(what=="GA") pmin(x@ga[xind,slice,zind]/maxga, 1) else x@fa[xind,slice,zind]
    anindex <- if(what=="GA") tanh(x@ga[xind,slice,zind]) else x@fa[xind,slice,zind]
    if (method == 3) {
      andirection <- x@bary[,xind,slice,zind]
    } else {
      andirection <- x@andir[,xind,slice,zind]
    }
  } else {
#    anindex <- if(what=="GA") pmin(x@ga[xind,yind,slice]/maxga, 1) else x@fa[xind,yind,slice]
    anindex <- if(what=="GA") tanh(x@ga[xind,yind,slice]) else x@fa[xind,yind,slice]
    if (method == 3) {
      andirection <- x@bary[,xind,yind,slice]
    } else {
      andirection <- x@andir[,xind,yind,slice]
    }
  }
  anindex[anindex>1] <- 0
  anindex[anindex<0] <- 0
  if ((method==1) || (method==2) || (method==4)) {
    if(contrast.enh<1&&fa.contrast.enh>0) anindex <- pmin(anindex/contrast.enh,1)
    if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minanindex <- 0
    if(method==1) {
      andirection[1,,] <- abs(andirection[1,,])
      andirection[2,,] <- abs(andirection[2,,])
      andirection[3,,] <- abs(andirection[3,,])
    } else if (method==2) {
      ind<-andirection[1,,]<0
      dim(andirection) <- c(3,prod(dim(ind)))
      andirection[,ind] <- - andirection[,ind]
      andirection[2,] <- (1+andirection[2,])/2
      andirection[3,] <- (1+andirection[3,])/2
      dim(andirection) <- c(3,dim(ind))
    } else {
      andirection[1,,] <- andirection[1,,]^2
      andirection[2,,] <- andirection[2,,]^2
      andirection[3,,] <- andirection[3,,]^2
    }
    andirection <- aperm(andirection,c(2,3,1))
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      andirection <- make.image(andirection,gamma=TRUE)
      if(show) show.image(andirection,...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
        par(oldpar)
        invisible(andirection)
      }
    } else if(show) {
      dim(anindex) <- dim(andirection)[2:3]
      image(anindex,...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
        par(oldpar)
        invisible(NULL)
      }
    }
  } else if (method==3) {
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      bary <- make.image(aperm(andirection,c(2,3,1)))
      if(show) show.image(bary,...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
         par(oldpar)
         invisible(bary)
      }
    } else if(show) {
      image(andirection[1,,],...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
         par(oldpar)
         invisible(NULL)
      }
    }
  } else if (method==5) {
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      img.hsi.data <- array(0,dim=c(dim(andirection)[2:3],3))
      img.hsi.data[,,1] <- atan2(andirection[2,,],andirection[1,,])
      img.hsi.data[,,1] <- img.hsi.data[,,1] + pi*(img.hsi.data[,,1]<0)
      img.hsi.data[,,2] <- abs(acos(andirection[3,,]))
      img.hsi.data[,,3] <- anindex
      img.hsi <- make.image(img.hsi.data,gamma=TRUE,xmode="HSI")
      if(show) show.image(img.hsi,...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
         par(oldpar)
         invisible(img.hsi)
      }
    } else if(show) {
      image(andirection[1,,],...)
      if(identify){
         identify.fa(view,slice,xind,yind,zind)
      } else {
         par(oldpar)
         invisible(NULL)
      }
    }
  }
})
