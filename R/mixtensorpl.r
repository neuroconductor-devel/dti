dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor","dtiData",function(object, maxcomp=3,
          model=c("MT","MTiso","MTisoFA","MTisoEV"),
          fa=NULL, lambda=NULL, mask=NULL, reltol=1e-10, maxit=5000, ngc=1000,
          nguess=100*maxcomp^2, msc=c("BIC","AIC","AICC","none"),
          mc.cores = setCores(,reprt=FALSE)){
  #
  #  uses  S(g) = w_0 exp(-b*l_1) +\sum_{i} w_i exp(-b*l_2-b*(l_1-l_2)(g^T d_i)^2)
  #  w_i corresponds to th0* volume fraction of compartment i
  #  FA or (FA and l_2)  may be fixed
  #  Optimization method: L-BFGS-B for tensor mixture models with isotropic compartment

  ## check model
  model <- match.arg(model)
  bvalue <- object@bvalue[-object@s0ind]
  if(is.null(mask)) mask <- object@mask
# prefer mask in dtiObject over fallback based on level
  if(model=="MT"&&sd(bvalue)>.1*median(bvalue)){
    cat("b-values indicate measurements on multiple shells,\n
         model 'MT' not yet implemented\n using model 'MTiso' instead\n")
    model<-"MTiso"
  }
  if(model=="MT"){
#
#    this handles the case of MixTensor-models without isotropic compartments
#    only implemented for single shell data
#
     return(dwiMixtensorMT(object, maxcomp, mask=NULL,
                reltol=1e-10, maxit=5000,ngc=1000, nguess=100*maxcomp^2,
                msc=c("BIC","AIC","AICC","none"), mc.cores = setCores(,reprt=FALSE)))
  }
  ## check msc
  msc <- match.arg(msc)
  factr <- reltol/1e-14 ## this is 1e6
  set.seed(1)
  bvalue <- object@bvalue
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  if(5*(2+3*maxcomp)>ngrad){
    #     maxcomp <- max(1,trunc((ngrad-5)/15))
    cat("Maximal number of components reduced to", maxcomp,"due to insufficient
           number of gradient directions\n")
  }
  #
  #   which model should be used
  #
  if(model=="MTisoEV"&&(is.null(fa)||is.null(lambda))){
    cat("No eigenvalues specified with model=='MTisoEV'\n")
    if(is.null(fa)){
      cat("setting model='MTiso'\n")
      model<-"MTiso"
    } else {
      cat("setting model='MTisoFA'\n")
      model<-"MTisoFA"
    }
  }
  imodel <- switch(model,MTisoEV=0,MTisoFA=1,MTiso=2)
  alpha <- switch(model,MTisoEV=(fa^2+sqrt(fa^2*(3-2*fa^2)))/(1-fa^2),
                  MTisoFA=(fa^2+sqrt(fa^2*(3-2*fa^2)))/(1-fa^2),
                  MTiso=0)# not used
  lambda <- switch(model,MTisoEV=lambda*maxbv/(1+alpha),
                   MTisoFA=0,# will be adjusted later
                   MTiso=0)#  will not be used
  #
  #  First tensor estimates to generate eigenvalues and -vectors
  #
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object, mask=mask, mc.cores = mc.cores)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- dtieigen(tensorobj@D, tensorobj@mask, mc.cores = mc.cores)
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))*maxbv
  #
  #  rescale by bvalue to go to implemented scale
  #
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
  #
  #  prepare parameters for searching initial estimates
  #
    lambdahat <- ev[3,,,] #
    # use third ev instead of (ev[2,,,]+ev[3,,,])/2 to avoid effects from mixtures
    alphahat <- if(imodel==2) (ev[1,,,]-lambdahat)/lambdahat else alpha
    lambdahat <- median(lambdahat[!is.na(alphahat)&fa>.6])
    if(imodel==2) alphahat <- max(5,median(alphahat[!is.na(alphahat)&fa>.6]))
    fahat <- alphahat/sqrt(3+2*alphahat+alphahat^2)
    cat("Using lambda_2=",lambdahat,"fa=",fahat," and alpha=",alphahat,"in initial estimates\n")

  cat("Start search outlier detection at",format(Sys.time()),"\n")
  #
  #  replace physically meaningless S_i by mena S_0 values
  #
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
  #
  #  this does not scale well with openMP
  #
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
  #
  nvox <- prod(ddim[1:3])
  cat("means0:")
  t1 <- Sys.time()
  z <- .Fortran(C_means0,# mixtensbv.f
                  as.double(si[s0ind,,,,drop=FALSE]),
                  as.integer(nvox),
                  as.integer(ns0),
                  as.integer(object@level),
                  s0=double(nvox),
                  mask=integer(nvox),
                  PACKAGE="dti")[c("s0","mask")]
  if(is.null(mask)) mask <- array(as.logical(z$mask),ddim[1:3])
# fallback if no mask definition was provided
  s0 <- array(z$s0,ddim[1:3])
  means0 <- max(s0[mask])
  #
  #  rescale s0 for numerical reasons
  #
  si <- si/means0
  npar <- switch(model,
                 MTisoEV=1+3*(0:maxcomp),
                 MTisoFA=2+3*(0:maxcomp),
                 MTiso=3+3*(0:maxcomp))
  #
  #   compute penalty for model selection, default BIC
  #
  penIC <- switch(msc,AIC=2*npar/ngrad,BIC=log(ngrad)*npar/ngrad,
                  AICC=(1+npar/ngrad)/(1-(npar+2)/ngrad),
                  None=log(ngrad)-log(ngrad-npar))
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  avoid situations where si's are larger than s0
  #
  grad <- t(object@gradient)
  #
  #   determine initial estimates for orientations
  #
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders", envir = environment())
  polyeder <- icosa3 <- icosa3 ## WORKAROUND to make "polyeders" not global variable (for R CMD check)
  vert <- polyeder$vertices
  # remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
  #
  #  compute initial estimates (EV from grid and orientations from icosa3$vertices)
  #
  dim(si) <- c(ngrad,nvox)
    siind <- wi <- matrix(0,maxcomp+1,nvox)
    siind[1,!mask] <- -1
    krit <- numeric(nvox)
    nvoxm <- sum(mask)
    if(mc.cores<=1){
      z <-   getsiindbv(si[,mask],grad,bvalue,t(vert),alphahat,lambdahat,
                       maxcomp,maxc=maxc,nguess=nguess)
      krit[mask] <- z$krit # sqrt(sum of squared residuals) for initial estimates
      siind[,mask] <- z$siind # components 1: model order 2:
      wi[,mask] <- z$wi # weigths
      # grid index for EV 2+(1:m) index of orientations
    } else {
      mc.cores.old <- setCores(,reprt=FALSE)
      setCores(mc.cores)
      x <- si[,mask]
      nvico <- dim(vert)[2]
      dgrad <- matrix(abs(grad%*%vert),ngrad,nvico)
      dgrad <- dgrad/max(dgrad)
      dgradi <- matrix(abs(t(vert)%*%vert),nvico,nvico)
      dgradi <- dgradi/max(dgradi)
      isample <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
      nguess <- length(isample)/maxcomp
      cat("using ",nguess,"guesses for initial estimates\n")
      z <- plmatrix(x,pgetsiindbv,grad=grad,bv=bvalue,nvico=nvico,
                    dgrad=dgrad,dgradi=dgradi,isample=isample,alpha=alphahat,
                    lambda=lambdahat,maxcomp=maxcomp,maxc=maxc,nguess=nguess)
      setCores(mc.cores.old,reprt=FALSE)
      krit[mask] <- z[1,] # risk
      siind[,mask] <- z[2:(maxcomp+2),] # siind
      wi[,mask] <- z[-(1:(maxcomp+2)),] # wi

    }
  cat("Model orders for initial estimates")
  print(table(siind[1,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  #  logarithmic eigen values
  orient <- array(0,c(2,maxcomp,ddim))
  prt0 <- Sys.time()
    siind <- siind[-1,,drop=FALSE]
    if(mc.cores<=1){
      cat("Starting parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      z <- switch(imodel+1,.C(C_mixtrl0b,
                              as.integer(nvoxm),#nvoxm
                              as.integer(siind[,mask]),#siind
                              as.double(wi[,mask]),# wi
                              as.integer(ngrad),#ngrad
                              as.integer(maxcomp),#maxcomp
                              as.integer(maxit),#maxit
                              as.double(t(grad)),#grad_in
                              as.double(bvalue),#bv_in
                              lambda  = as.double(lambda),#lambda_in
                              alpha   = as.double(alpha),#alpha_in
                              as.double(factr),#factr
                              as.double(penIC),#penIC
                              as.double(krit[mask]),#best rss
                              as.double(vert),#vert
                              as.double(si[,mask]),#si_in
                              sigma2  = double(nvoxm),#sigma2_ret error variance
                              orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                              order   = integer(nvoxm),#order_ret selected order of mixture
                              mix     = double((maxcomp+1)*nvoxm),#mixture weights
                              PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
                  .C(C_mixtrl1b,
                     as.integer(nvoxm),#n1
                     as.integer(siind[,mask]),#siind
                     as.double(wi[,mask]),# wi
                     as.integer(ngrad),#ngrad
                     as.integer(maxcomp),#maxcomp
                     as.integer(maxit),#maxit
                     as.double(t(grad)),#grad_in
                     as.double(bvalue),#bv_in
                     as.double(lambdahat),#lambda_in
                     alpha = as.double(alpha),#alpha_in
                     as.double(factr),#factr
                     as.double(penIC),#penIC
                     as.double(krit[mask]),#best rss
                     as.double(vert),#vert
                     as.double(si[,mask]),#si_in
                     sigma2  = double(nvoxm),#sigma2_ret error variance
                     orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                     order   = integer(nvoxm),#order_ret selected order of mixture
                     lambda  = double(nvoxm),#lambda_ret lambda_2
                     mix     = double((maxcomp+1)*nvoxm),#mixture weights
                     PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
                  .C(C_mixtrl2b,
                     as.integer(nvoxm),#n1
                     as.integer(siind[,mask]),#siind
                     as.double(wi[,mask]),# wi
                     as.integer(ngrad),#ngrad
                     as.integer(maxcomp),#maxcomp
                     as.integer(maxit),#maxit
                     as.double(t(grad)),#grad_in
                     as.double(bvalue),#bv_in
                     as.double(lambdahat),#lambda_in
                     as.double(alphahat),#alpha_in
                     as.double(factr),#factr
                     as.double(penIC),#penIC
                     as.double(krit[mask]),#best rss
                     as.double(vert),#vert
                     as.double(si[,mask]),#si_in
                     sigma2  = double(nvoxm),#sigma2_ret error variance
                     orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                     order   = integer(nvoxm),#order_ret selected order of mixture
                     alpha   = double(nvoxm),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2
                     lambda  = double(nvoxm),#lambda_ret lambda_2
                     mix     = double((maxcomp+1)*nvoxm))[c("sigma2","orient","order","alpha","lambda","mix")])
      cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      sigma2 <-  array(0,ddim)
      sigma2[mask] <- z$sigma2
      orient <- matrix(0,2*maxcomp,nvox)
      orient[,mask] <- z$orient
      dim(orient) <- c(2, maxcomp, ddim)
      order <- array(0, ddim)
      order[mask] <- z$order
      lev <- matrix(0,2,nvox)
      lev[2,mask] <- z$lambda
      lev[1,mask] <- (z$alpha+1)*z$lambda
      dim(lev) <- c(2,ddim)
      mix <- matrix(0,maxcomp+1,nvox)
      mix[,mask] <- z$mix
      dim(mix) <- c(maxcomp+1, ddim)
    } else {
      cat("Starting parameter estimation and model selection (C-code) on",mc.cores," cores",format(Sys.time()),"\n")
      x <- matrix(0,ngrad+2+2*maxcomp,sum(mask))
      dim(si) <- c(ngrad,nvox)
      x[1:ngrad,] <- si[,mask]
      x[ngrad+1,] <- krit[mask]
      dim(siind) <- c(maxcomp,nvox)
      x[ngrad+2:(1+maxcomp),] <- siind[,mask]
      x[-(1:(ngrad+1+maxcomp)),] <- wi[,mask]
      res <- matrix(0,5+3*maxcomp,nvox)
      res[,mask] <- switch(imodel+1,
                           plmatrix(x,pmixtn0b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambda,alpha=alpha,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores),# model=0
                           plmatrix(x,pmixtn1b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambdahat,alpha=alpha,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores),# model=1
                           plmatrix(x,pmixtn2b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambdahat,alpha=alphahat,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores))# model=2
      cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      rm(x)
      gc()
      sigma2 <-  array(res[2,],ddim)
      orient <- array(res[maxcomp+5+1:(2*maxcomp),], c(2, maxcomp, ddim))
      order <- array(as.integer(res[1,]), ddim)
      lev <- array(res[3:4,], c(2,ddim))
      mix <- array(res[5+0:maxcomp,], c(maxcomp+1, ddim))
    }
    th0 <- apply(mix,2:4,sum)
    order[th0==0] <- 0
    mix <- sweep(mix[-1,,,,drop=FALSE],2:4,th0,"/")
    mix[is.na(mix)] <- 0
    if(any(mix<0)){
       cat("neg. weights:", sum(mix<0), "minimum", min(mix)," replaced by 0\n")
       mix[mix<0] <- 0
    }
    th0 <- th0*means0
    save(krit,wi,mask,z,mix,th0,file="tmp.rsc")
    method <- switch(imodel+1,"MTisoEV","MTisoFA","MTiso")
    model <- switch(imodel+1,"iso-prolate-fixedev","iso-prolate-fixedfa","iso-prolate")
  invisible(new("dwiMixtensor",
                model = model,
                call   = args,
                ev     = lev/maxbv,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = 0,
                th0    = th0,
                sigma  = sigma2*means0^2,
                scorr  = array(1,c(1,1,1)),
                bw     = c(0,0,0),
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
  )
}
)
#
#   Initial estimates
#
selisample <- function(ngrad,maxcomp,nguess,dgrad,maxc){
  saved.seed <- .Random.seed
  set.seed(1)
  isample <- matrix(sample(ngrad,maxcomp*nguess,replace=TRUE),maxcomp,nguess)
  ind <- rep(TRUE,nguess)
  if(maxcomp>1){
    ind <- .Fortran(C_selisamp,
                    as.integer(isample),
                    as.integer(nguess),
                    as.integer(maxcomp),
                    as.double(dgrad),
                    as.integer(dim(dgrad)[1]),
                    ind = integer(nguess),
                    as.double(maxc))$ind
    .Random.seed <- saved.seed
  }
  isample[,ind]
}
getsiind3 <- function(si,mask,sigma2,grad,bv,vico,th,indth,ev,fa,andir,maxcomp=3,
                      maxc=.866,nguess=100,mc.cores = setCores(,reprt=FALSE)){
  # assumes dim(grad) == c(ngrad,3)
  # assumes dim(si) == c(ngrad,n1,n2,n3)
  # SO removed
  ngrad <- dim(grad)[1]
  nvico <- dim(vico)[1]
  ddim <- dim(fa)
  nsi <- dim(si)[1]
  dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
  dgrad <- dgrad/max(dgrad)
  dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
  dgradi <- dgradi/max(dgradi)
  nth <- length(th)
  nvoxel <- prod(ddim)
  landir <- fa>.3
  landir[is.na(landir)] <- FALSE
  if(any(is.na(andir))) {
    cat(sum(is.na(andir)),"na's in andir")
    andir[is.na(andir)]<-sqrt(1/3)
  }
  if(any(is.na(landir))) {
    cat(sum(is.na(landir)),"na's in landir")
    landir[is.na(landir)]<-0
  }
  if(any(is.na(fa))) {
    cat(sum(is.na(fa)),"na's in fa")
    fa[is.na(fa)]<-0
  }
  iandir <- .Fortran(C_iandir,
                     as.double(t(vico)),
                     as.integer(nvico),
                     as.double(andir),
                     as.integer(nvoxel),
                     as.integer(landir),
                     iandir=integer(prod(ddim)))$iandir
  isample0 <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
  if(maxcomp>1) isample1 <- selisample(nvico,maxcomp-1,nguess,dgradi,maxc)
  if(maxcomp==1) isample1 <- sample(ngrad, nguess, replace = TRUE)
  #
  #  eliminate configurations with close directions
  #
  # this provides configurations of initial estimates with minimum angle between
  # directions > acos(maxc)
  nvoxel <- prod(dim(si)[-1])
  cat("using ",nguess,"guesses for initial estimates\n")
  siind <- matrix(as.integer(0),maxcomp+2,nvoxel)
  krit <- numeric(nvoxel)
  # first voxel with fa<.3
  if(sum(mask&!landir)>0){
    cat(sum(mask&!landir),"voxel with small FA\n")
    nguess <- length(isample0)/maxcomp
    if(mc.cores<=1){
      z <- .Fortran(C_getsii30,
                    as.double(si),
                    as.double(sigma2),
                    as.integer(nsi),
                    as.integer(nvoxel),
                    as.integer(maxcomp),
                    as.double(dgrad),
                    as.integer(nvico),
                    as.double(th),
                    as.integer(nth),
                    as.integer(indth),
                    double(nsi*nvico),
                    as.integer(isample0),
                    as.integer(nguess),
                    double(nsi),
                    double(nsi*(maxcomp+2)),
                    siind=integer((maxcomp+2)*nvoxel),
                    krit=double(nvoxel),
                    as.integer(maxcomp+2),
                    as.integer(mask&!landir))[c("siind","krit")]
      dim(z$siind) <- c(maxcomp+2,nvoxel)
      siind[,!landir] <- z$siind[,!landir]
      krit[!landir] <- z$krit[!landir]
    } else {
      x <- matrix(0,nsi+2,sum(mask&!landir))
      x[1:nsi,] <- matrix(si,nsi,nvoxel)[,mask&!landir]
      x[nsi+1,] <- sigma2[mask&!landir]
      x[nsi+2,] <- indth[mask&!landir]
      z <- plmatrix(x,pgetsii30,maxcomp=maxcomp,dgrad=dgrad,th=th,isample0=isample0,
                    nsi=nsi,nth=length(th),nvico=nvico,nguess=nguess,
                    mc.cores=mc.cores)
      dim(z) <- c(maxcomp+3,sum(mask&!landir))
      siind[,mask&!landir] <- z[-1,]
      krit[mask&!landir] <- z[1,]
    }
  }
  # now voxel where first tensor direction seems important
  if(sum(mask&landir)>0){
    if(maxcomp >0){
      cat(sum(mask&landir),"voxel with distinct first eigenvalue \n")
      nguess <- if(maxcomp>1) length(isample1)/(maxcomp-1) else length(isample1)
      if(mc.cores<=1){
        z <- .Fortran(C_getsii31,
                      as.double(si),
                      as.double(sigma2),
                      as.integer(nsi),
                      as.integer(nvoxel),
                      as.integer(maxcomp),
                      as.double(dgrad),
                      as.integer(nvico),
                      as.integer(iandir),
                      as.double(th),
                      as.integer(nth),
                      as.integer(indth),
                      double(nsi*nvico),
                      as.integer(isample1),
                      as.integer(nguess),
                      double(nsi),
                      double(nsi*(maxcomp+2)),
                      siind=integer((maxcomp+2)*nvoxel),
                      krit=double(nvoxel),
                      as.integer(maxcomp+2),
                      as.integer(mask&landir),
                      as.double(dgradi),
                      as.double(maxc))[c("siind","krit")]
        dim(z$siind) <- c(maxcomp+2,nvoxel)
        siind[,landir] <- z$siind[,landir]
        krit[landir] <- z$krit[landir]
      } else {
        x <- matrix(0,nsi+3,sum(mask&landir))
        x[1:nsi,] <- matrix(si,nsi,nvoxel)[,mask&landir]
        x[nsi+1,] <- sigma2[mask&landir]
        x[nsi+2,] <- indth[mask&landir]
        x[nsi+3,] <- iandir[mask&landir]
        z <- plmatrix(x,pgetsii31,maxcomp=maxcomp,dgrad=dgrad,th=th,isample1=isample1,
                      nsi=nsi,nth=length(th),nvico=nvico,nguess=nguess,
                      dgradi=dgradi,maxc=maxc,mc.cores=mc.cores)
        dim(z) <- c(maxcomp+3,sum(mask&landir))
        siind[,mask&landir] <- z[-1,]
        krit[mask&landir] <- z[1,]
      }
    }
  }
  list(siind=array(siind,c(maxcomp+2,dim(si)[-1])),
       krit=array(krit,dim(si)[-1]))
}
#
#   create initial estimates for directions
#

getsiindbv <- function(si,grad,bv,vico,alpha,lambda,maxcomp=3,maxc=.866,nguess=100){
  # assumes dim(grad) == c(ngrad,3)
  # assumes dim(si) == c(ngrad,n1,n2,n3)
  # SO removed
  ngrad <- dim(grad)[1]
  nvico <- dim(vico)[1]
  nsi <- dim(si)[1]
  dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
  dgrad <- dgrad/max(dgrad)
  dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
  dgradi <- dgradi/max(dgradi)
  isample <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
  nguess <- length(isample)/maxcomp
  #
  #  eliminate configurations with close directions
  #
  # this provides configurations of initial estimates with minimum angle between
  # directions > acos(maxc)
  nvoxel <- prod(dim(si)[-1])
  cat("using ",nguess,"guesses for initial estimates\n")
  z <- .Fortran(C_getsiibv,
                as.double(si),
                as.integer(nsi),
                as.integer(nvoxel),
                as.integer(maxcomp),
                as.double(dgrad),
                as.double(bv),
                as.integer(nvico),
                as.double(alpha),
                as.double(lambda),
                double(nsi*nvico),
                as.integer(isample),
                as.integer(nguess),
                double(nsi),
                double(nsi),#z0
                double(nsi*(maxcomp+1)),
                siind=integer((maxcomp+1)*nvoxel),
                wi=double((maxcomp+1)*nvoxel),
                krit=double(nvoxel),
                as.integer(maxcomp+1),
                PACKAGE="dti")[c("siind","krit","wi")]
  dim(z$siind) <- dim(z$wi) <- c(maxcomp+1,nvoxel)
  siind <- z$siind
  wi <- z$wi
  krit <- z$krit
  # now voxel where first tensor direction seems important
  list(siind=array(siind,c(maxcomp+1,dim(si)[-1])),
       wi=array(wi,c(maxcomp+1,dim(si)[-1])),
       krit=array(krit,dim(si)[-1]))
}

dwiMixtensorMT <- function(object, maxcomp=3, mask=NULL,
              reltol=1e-10, maxit=5000,ngc=1000, nguess=100*maxcomp^2,
              msc=c("BIC","AIC","AICC","none"), mc.cores = setCores(,reprt=FALSE)){
  #
  #  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
  #
  #  Optimization methods:
  #     BFGS for tensor mixture models without isotropic compartment
  #     L-BFGS-B for tensor mixture models with isotropic compartment
  #
  ## check msc
  msc <- match.arg(msc)
  factr <- reltol/1e-14 ## this is 1e6
  set.seed(1)
  bvalue <- object@bvalue[-object@s0ind]
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
    cat("Maximal number of components reduced to", maxcomp,"due to insufficient
           number of gradient directions\n")
  }
  imodel <- 3
  alpha <- 0# not used
  lambda <- 0# will not be used
  #
  #  First tensor estimates to generate eigenvalues and -vectors
  #
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object, mask=mask, mc.cores = mc.cores)
  mask <- tensorobj@mask
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- dtieigen(tensorobj@D, tensorobj@mask, mc.cores = mc.cores)
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))*maxbv
  #
  #  rescale by bvalue to go to implemented scale
  #
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
  #
  #  prepare parameters for searching initial estimates
  #
    nth <- 5
    th <- (ev[1,,,] - (ev[2,,,]+ev[3,,,])/2)
    falevel <- min(quantile(fa[fa>0],.75),.4)
    cat("falevel",falevel,"\n")
    qth <- unique(quantile(th[fa>=falevel&fa<.95],seq(.8,.99,length=nth)))
    nth <- length(qth)
    if(nth>1){
      indth <- cut(th,qth,labels=FALSE)
      indth[th<=qth[1]] <- 1
      indth[th>=qth[nth]] <- nth
      th <- qth
    } else {
      indth <- rep(1,nvox)
      th <- qth
    }
    cat("using th:::",th,"\n")
  cat("Start search outlier detection at",format(Sys.time()),"\n")
  #
  #  replace physically meaningless S_i by mena S_0 values
  #
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
  #
  #  this does not scale well with openMP
  #
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
  #
  nvox <- prod(ddim[1:3])
  cat("sweeps0:")
  t1 <- Sys.time()
  if(mc.cores==1||ngrad0>250){
    z <- .Fortran(C_sweeps0,# mixtens.f
                  as.double(si[-s0ind,,,,drop=FALSE]),
                  as.double(si[s0ind,,,,drop=FALSE]),
                  as.integer(nvox),
                  as.integer(ns0),
                  as.integer(ngrad0),
                  as.integer(object@level),
                  siq=double(nvox*ngrad0),
                  s0=double(nvox),
                  vsi=double(nvox),
                  mask=integer(nvox))[c("siq","s0","vsi","mask")]
    z$mask <- as.logical(z$mask)
    t2 <- Sys.time()
    cat(difftime(t2,t1),"for",nvox,"voxel\n")
    s0 <- array(z$s0,ddim[1:3])
    siq <- array(z$siq,c(ngrad0,ddim[1:3]))
    #
    #  siq is permutated c(4,1:3)
    #
    sigma2 <- array(z$vsi,ddim[1:3])
    mask <- array(z$mask&mask,ddim[1:3])
  } else {
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
    z <- matrix(.Fortran(C_sweeps0p,# mixtens.f
                         as.double(si[-s0ind,,,,drop=FALSE]),
                         as.double(si[s0ind,,,,drop=FALSE]),
                         as.integer(nvox),
                         as.integer(ns0),
                         as.integer(ngrad0),
                         as.integer(object@level),
                         siq=double(nvox*(ngrad0+3)),
                         as.integer(ngrad0+3))$siq,ngrad0+3,nvox)
    t2 <- Sys.time()
    cat(difftime(t2,t1),"for",nvox,"voxel\n")
    setCores(mc.cores.old,reprt=FALSE)
    s0 <- array(z[ngrad0+1,],ddim[1:3])
    siq <- array(z[1:ngrad0,],c(ngrad0,ddim[1:3]))
    #
    #  siq is permutated c(4,1:3)
    #
    sigma2 <- array(z[ngrad0+2,],ddim[1:3])
    mask <- array(as.logical(z[ngrad0+3,])&mask,ddim[1:3])
  }
  rm(si)
  rm(z)
  gc()
  npar <- 1+3*(0:maxcomp)
  #
  #   compute penalty for model selection, default BIC
  #
  penIC <- switch(msc,AIC=2*npar/ngrad0,BIC=log(ngrad0)*npar/ngrad0,
                  AICC=(1+npar/ngrad0)/(1-(npar+2)/ngrad0)-1,
                  None=log(ngrad0)-log(ngrad0-npar))
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  avoid situations where si's are larger than s0
  #
  grad <- t(object@gradient[,-s0ind])
  #
  #   determine initial estimates for orientations
  #
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders", envir = environment())
  polyeder <- icosa3 <- icosa3 ## WORKAROUND to make "polyeders" not global variable (for R CMD check)
  vert <- polyeder$vertices
  # remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
  #
  #  compute initial estimates (EV from grid and orientations from icosa3$vertices)
  #
  dim(siq) <- c(ngrad0,nvox)
    siind <-  getsiind3(siq,mask,sigma2,grad,bvalue,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess,mc.cores=mc.cores)
    krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
    siind <- siind$siind # components 1: model order 2:
    # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  #  logarithmic eigen values
  orient <- array(0,c(2,maxcomp,ddim))
  prt0 <- Sys.time()
  #
  #   loop over voxel in volume
  #
    if(mc.cores<=1){
      cat("Starting parameter estimation and model selection",format(Sys.time()),"\n")
      #  dim(siind) <- c(2+maxcomp,nvox)
      nvoxm <- sum(mask)
      z <- .C(C_mixture,
              as.integer(nvoxm),
              as.integer(siind[,mask]),
              as.integer(ngrad0),
              as.integer(maxcomp),
              as.integer(maxit),
              as.double(100),# penalty for unfeasible parameters
              as.double(t(grad)),
              as.double(reltol),
              as.double(th),
              as.double(penIC),
              as.double(sigma2[mask]),
              as.double(vert),
              as.double(siq[,mask]),
              sigma2  = double(nvoxm),# error variance
              orient  = double(2*maxcomp*nvoxm), # phi/theta for all mixture tensors
              order   = integer(nvoxm),   # selected order of mixture
              lev     = double(2*nvoxm),         # logarithmic eigenvalues
              mix     = double(maxcomp*nvoxm))[c("sigma2","orient","order","lev","mix")]
      cat("End parameter estimation and model selection ",format(Sys.time()),"\n")
      sigma2 <-  array(0,ddim)
      sigma2[mask] <- z$sigma2
      orient <- matrix(0,2*maxcomp,nvox)
      orient[,mask] <- z$orient
      dim(orient) <- c(2, maxcomp, ddim)
      order <- array(0, ddim)
      order[mask] <- z$order
      lev <- matrix(0,2,nvox)
      lev[,mask] <- z$lev
      dim(lev) <- c(2,ddim)
      mix <- matrix(0,maxcomp,nvox)
      mix[,mask] <- z$mix
      dim(mix) <- c(maxcomp, ddim)
    } else {
      cat("Starting parameter estimation and model selection on",mc.cores," cores",format(Sys.time()),"\n")
      x <- matrix(0,ngrad0+3+maxcomp,sum(mask))
      dim(siq) <- c(ngrad0,nvox)
      x[1:ngrad0,] <- siq[,mask]
      x[ngrad0+1,] <- sigma2[mask]
      dim(siind) <- c(2+maxcomp,nvox)
      x[ngrad0+2:(3+maxcomp),] <- siind[,mask]
      res <- matrix(0,4+3*maxcomp,nvox)
      res[,mask] <- plmatrix(x,pmixtens,
                             ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                             pen=100,grad=grad,reltol=reltol,th=th,
                             penIC=penIC,vert=vert,
                             mc.cores=mc.cores)
      cat("End parameter estimation and model selection ",format(Sys.time()),"\n")
      rm(x)
      gc()
      sigma2 <-  array(res[2,],ddim)
      orient <- array(res[maxcomp+4+1:(2*maxcomp),], c(2, maxcomp, ddim))
      order <- array(as.integer(res[1,]), ddim)
      lev <- array(res[3:4,], c(2,ddim))
      lev[1,,,] <- lev[1,,,] + lev[2,,,] ## lev[1,...] contained difference before
      mix <- array(res[4+(1:maxcomp),], c(maxcomp, ddim))
    }
    method <- "MT"
    model <- "prolate"
   invisible(new("dwiMixtensor",
                model = model,
                call   = args,
                ev     = lev/maxbv,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = 0,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)),
                bw     = c(0,0,0),
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
  )
}
