//  L-BFGS-B  version
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

int ngradcc = 0;
int iibv = 0;
double* si_init, *bv, *gradients;
double alpha, lambda;

extern void F77_NAME(rskmixb2)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* result);

extern void F77_NAME(rskmixb1)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* alpha, double* result);

extern void F77_NAME(rskmixb0)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* lambda,
             double* alpha, double* result);

extern void F77_NAME(drskmb2)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* result);

extern void F77_NAME(drskmb1)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* alpha, double* result);

extern void F77_NAME(drskmb0)(double* param, int* npar, double* si,
             double* grad, double* bv, int* ng, double* lambda,
             double* alpha, double* result);

extern void paroforient(double *dir, double *angles);

typedef struct
{
  int ngrad;
  double* si;
  double* grad;
  double* bv;
  int fnscale;
} optimex2;

typedef struct
{
  int ngrad;
  double* si;
  double* grad;
  double* bv;
  double alpha;
  int fnscale;
} optimex1;

typedef struct
{
  int ngrad;
  double* si;
  double* grad;
  double* bv;
  double lambda;
  double alpha;
  int fnscale;
} optimex0;

typedef struct{
  int order;
  double lambda;
  double alpha;
  double w0;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskmb2_ret;

typedef struct{
  int order;
  double lambda;
  double w0;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskmb1_ret;

typedef struct{
  int order;
  double w0;
  double* mix;
  double* orient;
  double* param;
  double value;
} mfunrskmb0_ret;





double rskmixb2(int param_length, double *param, void* ex){

  double result = 0;

  double* si = calloc(ngradcc, sizeof(double));
  int i;
  for(i=0;i<ngradcc;i++){
    si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(rskmixb2)(param, &param_length, si, gradients, bv, &ngradcc, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  free(si);

  return result;
}

double rskmixb1(int param_length, double *param, void* ex){

  double result = 0;

  double* si = calloc(ngradcc, sizeof(double));
  int i;
  for(i=0;i<ngradcc;i++){
    si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(rskmixb1)(param, &param_length, si, gradients, bv, &ngradcc, &alpha, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  free(si);

  return result;
}

double rskmixb0(int param_length, double *param, void* ex){

  double result = 0;

  double* si = calloc(ngradcc, sizeof(double));
  int i;
  for(i=0;i<ngradcc;i++){
    si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(rskmixb0)(param, &param_length, si, gradients, bv, &ngradcc, &lambda, &alpha, &result);

  //check for infinity
//  if(result == R_PosInf || result == R_NegInf){
//    return 0;
//  }

  free(si);

  return result;
}

void drskmb2(int param_length, double* param, double* result, void* ex){

  int i;
  double* si = calloc(ngradcc, sizeof(double));
  for(i=0;i<ngradcc;i++){
   si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(drskmb2)(param, &param_length, si, gradients, bv, &ngradcc, result);

  free(si);
}

void drskmb1(int param_length, double* param, double* result, void* ex){

  int i;
  double* si = calloc(ngradcc, sizeof(double));
  for(i=0;i<ngradcc;i++){
   si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(drskmb1)(param, &param_length, si, gradients, bv, &ngradcc, &alpha, result);

  free(si);
}

void drskmb0(int param_length, double* param, double* result, void* ex){

  int i;
  double* si = calloc(ngradcc, sizeof(double));
  for(i=0;i<ngradcc;i++){
   si[i] = si_init[i+iibv*ngradcc];
  }

  F77_CALL(drskmb0)(param, &param_length, si, gradients, bv, &ngradcc, &lambda, &alpha, result);

  free(si);
}

mfunrskmb2_ret getparam2b(int param_length, double* param, double fmin){

  mfunrskmb2_ret ret_val;
  int i,j;

  int c_ord = (param_length-3)/3;

  double  w0, lambda, alpha;
  double* w_tmp = calloc(param_length, sizeof(double));
  double* param_work = calloc(param_length, sizeof(double));
  int* o = calloc(c_ord, sizeof(int));

   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
 //  calculate weights
  double* mix = (double*) R_alloc((long unsigned int) c_ord, sizeof(double));
  for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i];
      o[i] = i;
   }
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
//  sorted index in o
   double* orient = (double*) R_alloc(2*(long unsigned int) c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[1+2*i];
   }
  lambda = param[3*c_ord+1];
  alpha = param[3*c_ord+2];
  w0 = param[3*c_ord];
  ret_val.order = c_ord;
  ret_val.lambda = lambda;
  ret_val.alpha = alpha;
  ret_val.w0 = w0;
  ret_val.mix = mix;
  ret_val.orient = orient;
  ret_val.param = param;
  ret_val.value = fmin;
//  Rprintf("getparam2: npar %i, fmin= %f ", param_length, fmin);
//  Rprintf("\n");

  free(o);
  free(param_work);
  free(w_tmp);

  return ret_val;
}

mfunrskmb1_ret getparam1b(int param_length, double* param, double fmin){

  mfunrskmb1_ret ret_val;
  int i,j;

  int c_ord = (param_length-2)/3;

  double* w_tmp = calloc(param_length, sizeof(double));
  double* param_work = calloc(param_length, sizeof(double));
  int* o = calloc(c_ord, sizeof(int));

   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
//  calculate weights
   double* mix = (double*) R_alloc((long unsigned int) c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i];
      o[i] = i;
   }
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
//  sorted index in o
   double* orient = (double*) R_alloc(2*(long unsigned int) c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[2*i+1];
   }
  ret_val.w0 = param[3*c_ord];
  ret_val.lambda = param[3*c_ord+1];
  ret_val.order = c_ord;
  ret_val.param = param;
  ret_val.mix = mix;
  ret_val.orient = orient;
  ret_val.value = fmin;

  free(o);
  free(param_work);
  free(w_tmp);

  return ret_val;
}

mfunrskmb0_ret getparam0b(int param_length, double* param, double fmin){

  mfunrskmb0_ret ret_val;
  int i,j;

  int c_ord = (param_length-1)/3;

  double* w_tmp = calloc(param_length, sizeof(double));
  double* param_work = calloc(param_length, sizeof(double));
  int* o = calloc(c_ord, sizeof(int));

   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }
//  calculate weights
   double* mix = (double*) R_alloc((long unsigned int) c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      mix[i] = param[3*i];
      o[i] = i;
   }
   revsort(mix, o, c_ord);
//  this sorts in decreasing order, now rearrange parameters
   double* orient = (double*) R_alloc(2*(long unsigned int) c_ord, sizeof(double));
   for( i = 0; i < c_ord; i++){
      w_tmp[i] = param[3*o[i]];
      orient[2*i] = param[3*o[i]+1];
      orient[2*i+1] = param[3*o[i]+2];
   }
   for(j = 0; j < c_ord; j++){
   while (orient[2*j] < 0) {
      orient[2*j] = orient[2*j] + M_PI;
   }
   while (orient[2*j] > M_PI) {
      orient[2*j] = orient[2*j] - M_PI;
   }
   while (orient[1+2*j] < 0) {
      orient[1+2*j] = orient[1+2*j] + 2 * M_PI;
   }
   while (orient[1+2*j] > 2 * M_PI) {
      orient[1+2*j] = orient[1+2*j] - 2 * M_PI;
   }
  }
   for( i = 0; i < c_ord; i++){
      param[3*i] = w_tmp[i];
      param[3*i+1] = orient[2*i];
      param[3*i+2] = orient[2*i+1];
   }
  ret_val.w0 = param[3*c_ord];
  ret_val.order = c_ord;
  ret_val.mix = mix;
  ret_val.orient = orient;
  ret_val.param = param;
  ret_val.value = fmin;

  free(o);
  free(param_work);
  free(w_tmp);

  return ret_val;
}

void mixtrl2b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp, int* maxit,
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,
          double* factr, double* penIC, double* sigma2, double* vert,
          double* si_in, double* sigma2_ret, double* orient_ret, int* order_ret,
          double* alpha_ret, double* lambda_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, general EV
// optmethod: L-BFGS-B
// n1 - number of voxel
// siind - initial parameters for orientations
// wi - compartment fraction * th0
// ngrad - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// si_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// alpha_ret  - estimated alpha values (n1)
// lambda_ret - estimated lambda values (n1)
// mix_ret - estimated mixture coefficients (maxcomp+1,n1)

   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init;
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0, smallw = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0,
          *param_tmp = 0, *si = 0, *lower = 0, *upper = 0;
   char msg[60];

   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim
    //Setting global variables
   si_init = si_in; gradients = grad_in, ngradcc = *ngrad; bv = bv_in;
   alpha = *alpha_in;
   lambda = *lambda_in;

   //Gradient vectors corresponding to minima in spherical coordinates

   param_length_init=3*maxcompc+3;

   param = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_work = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_last = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   lower = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   upper = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc((long unsigned int) param_length_init, sizeof(int));
   si = (double *) R_alloc((long unsigned int) ngradcc, sizeof(double));

   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }
   for(iibv = 0; iibv < *n1; iibv++){
       for(i=0; i<param_length_init; i++){
          lower[i] = R_NegInf;
          upper[i] = R_PosInf;
	        nbd[i] = 0;
       }
      param_length = param_length_init;
      sigma2_ret[iibv] = sigma2[iibv];

         ord=maxcompc+1;

         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*iibv]-1;

            if(iv==-1) iv = j; // this should never happen

            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];

            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*iibv] = angles[0];
            orient_ret[1+2*j+2*maxcompc*iibv] = angles[1];

            param[3*j] = wi[j+1+(maxcompc+1)*iibv];
//initial value for w(j) corresponds to volume w(j)/(1+sum_k w(j)), i.e. a 17% isotropic compartment
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 0;
            upper[3*j] = 1;
	          nbd[3*j] = 2;
         }
// lower values for lambda and alpha to get lambda_1>=lambda_2>=0
         lower[3*maxcompc] = 0;//w0
         upper[3*maxcompc] = 1;//w0
         lower[3*maxcompc+1] = 1e-4;//lambda
         upper[3*maxcompc+1] = 10;
         lower[3*maxcompc+2] = 0.5;//alpha
         upper[3*maxcompc+2] = 10;
	       nbd[3*maxcompc+1] = 2;
	       nbd[3*maxcompc+2] = 2;
	       nbd[3*maxcompc] = 2;
         // initialize EV-parameter and w0
	       param[3*maxcompc] = wi[(maxcompc+1)*iibv];
         param[3*maxcompc+1] = lambda;//lambda_2
         param[3*maxcompc+2] = alpha;
//alpha; lambda_1=(1+alpha)*lambda; FA=alpha/sqrt((1+alpha)^2+2)
//alpha=1.7 corresponds to a FA of 0.7

         sigma_init = sigma2[iibv];
         krit = log(sigma_init) + penIC[0];


         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
         }


         // use AIC/ngrad, BIC/ngrad or AICC/ngrad respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
// reorder such that parameters for smallest compartment are last
	            	  smallw = param_work[3*k+3];
		              for (l = 0; l < k; l++){
		                 if(param_work[3*l] < smallw){
		                    smallw = param_work[3*l];
		                  	param_work[3*l] = param_work[3*k];
		                  	param_work[3*k] = smallw;
		                    ttt = param_work[3*l+1];
		                  	param_work[3*l+1] = param_work[3*k+1];
		                  	param_work[3*k+1] = ttt;
		                    ttt = param_work[3*l+2];
		                  	param_work[3*l+2] = param_work[3*k+2];
		                  	param_work[3*k+2] = ttt;
		                 }
	            	  }
                  for (l = 0; l < param_length-3; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_tmp[3*k] = param_work[3*k+3];
                  param_tmp[3*k+1] = param_work[3*k+4];
                  param_tmp[3*k+2] = param_work[3*k+5];
                  param_length=3*k+3;
                  lower[3*k] = 0.;
                  upper[3*k] = 1.; // si are standardized by maxs0
                  lower[3*k+1] = 0.01;//lambda
                  upper[3*k+1] = 10;//lambda
                  lower[3*k+2] = 0.5;//alpha
                  upper[3*k+2] = 20;//alpha
 	                nbd[3*k] = 2;
	                nbd[3*k+1] = 2;
 	                nbd[3*k+2] = 2;
                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }

               optimex2 myoptimpar;

               //calculation and estimation of the best fitting model

               // initialize ret_val to avoid warnings
               mfunrskmb2_ret ret_val;
               ret_val.order = 0;
               ret_val.w0 = 0;
               ret_val.lambda = 0;
               ret_val.alpha = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;

               for(l=0;l<ngradcc;l++){
                  si[l] = si_init[l+iibv*ngradcc];
               }

               lbfgsb(param_length, lmm, param_work, lower, upper, nbd,
                      &Fmin, rskmixb2, drskmb2, &fail, &myoptimpar, *factr, pgtol,
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);

               ret_val = getparam2b(param_length, param_work, Fmin);

               value = ret_val.value;
               ord = ret_val.order;

               if(ord < k)
               {
                  ttt = krit;

                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradcc;

                  ttt = log(si2new) + penIC[ord];

                  for(l = 0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }


               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){
                  krit = ttt;

                  order_ret[iibv] = ret_val.order;

                  lambda_ret[iibv] = ret_val.lambda;
                  alpha_ret[iibv] = ret_val.alpha;
                  mix_ret[(maxcompc+1)*iibv] = ret_val.w0;
                  for(l = 0; l < ord; l++){
                     mix_ret[l+1+(maxcompc+1)*iibv] = ret_val.mix[l];
                     orient_ret[2*l+2*maxcompc*iibv] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*iibv] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+1+(maxcompc+1)*iibv] = 0;
                  }

                  sigma2_ret[iibv] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end iibv
}
void mixtrl1b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp, int* maxit,
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,
          double* factr, double* penIC, double* sigma2, double* vert,
          double* si_in, double* sigma2_ret, double* orient_ret,
          int* order_ret, double* lambda_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, fixed FA
// optmethod: L-BFGS-B
// n1 - number of voxel
// siind - initial parameters for orientations
// ngrad - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// alpha_in - (lambda1-lambda_2)/lambda_2
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// si_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// lambda_ret - estimated lambda values (n1)
// mix_ret - estimated mixture coefficients (maxcomp+1,n1)

   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init;
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0,
          *param_tmp = 0, *si = 0, *lower = 0, *upper = 0;
   char msg[60];

   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim

   //Setting global variables
   si_init = si_in; gradients = grad_in, ngradcc = *ngrad; bv = bv_in;
   alpha = *alpha_in;
   lambda = *lambda_in;


   //Gradient vectors corresponding to minima in spherical coordinates

   param_length_init=3*maxcompc+2;

   param = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_work = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_last = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   lower = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   upper = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc((long unsigned int) param_length_init, sizeof(int));
   si = (double *) R_alloc((long unsigned int) ngradcc, sizeof(double));

   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }

   for(iibv = 0; iibv < *n1; iibv++){
      for(i=0; i<param_length_init; i++){
         lower[i] = R_NegInf;
         upper[i] = R_PosInf;
	       nbd[i] = 0;
      }
      param_length = param_length_init;
      sigma2_ret[iibv] = sigma2[iibv];

         ord=maxcompc+1;

         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*iibv]-1;

            if(iv==-1) iv = j; // this should never happen

            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];

            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*iibv] = angles[0];
            orient_ret[1+2*j+2*maxcompc*iibv] = angles[1];

            param[3*j] = wi[j+1+(maxcompc+1)*iibv];
//initial value for w(j) corresponds to volume w(j)/(1+sum_k w(j)), i.e. a 17% isotropic compartment
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 0;
            upper[3*j] = 1;
	          nbd[3*j] = 2;
         }
         lower[3*maxcompc] = 0;//w0
         upper[3*maxcompc] = 1;//w0
         lower[3*maxcompc+1] = 1e-4;//lambda
         upper[3*maxcompc+1] = 10;
         nbd[3*maxcompc] = 2;
	       nbd[3*maxcompc+1] = 2;
// lower value for alpha to get lambda_1>=lambda_2
         // initialize EV-parameter
         param[3*maxcompc] = wi[(maxcompc+1)*iibv];
         param[3*maxcompc+1] = lambda;//lambda_2
//alpha; lambda_1=(1+alpha)*lambda; FA=alpha/sqrt((1+alpha)^2+2)
//alpha=1.7 corresponds to a FA of 0.7

         sigma_init = sigma2[iibv];
         krit = log(sigma_init) + penIC[0];


         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
         }


         // use AIC/ngrad, BIC/ngrad or AICC/ngrad respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length-1; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_tmp[3*k] = param_work[3*k+3];
                  param_tmp[3*k+1] = param_work[3*k+4];
                  param_length=3*k+2;
                  lower[3*k] = 0.;
                  upper[3*k] = 1.; // si are standardized by maxs0
                  lower[3*k+1] = 0.01;//lambda
                  upper[3*k+1] = 10;//lambda
                  nbd[3*k] = 2;
                  nbd[3*k+1] = 2;

                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }

               optimex1 myoptimpar;

               //calculation and estimation of the best fitting model

               // initialize ret_val to avoid warnings
               mfunrskmb1_ret ret_val;
               ret_val.w0 = 0;
               ret_val.order = 0;
               ret_val.lambda = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;

               for(l=0;l<ngradcc;l++){
                  si[l] = si_init[l+iibv*ngradcc];
               }

               lbfgsb(param_length, lmm, param_work, lower, upper, nbd,
                      &Fmin, rskmixb1, drskmb1, &fail, &myoptimpar, *factr, pgtol,
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);

               ret_val = getparam1b(param_length, param_work, Fmin);

               value = ret_val.value;
               ord = ret_val.order;

               if(ord < k)
               {
                  ttt = krit;

                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradcc;

                  ttt = log(si2new) + penIC[ord];

                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }


               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){
                  krit = ttt;

                  order_ret[iibv] = ret_val.order;

                  lambda_ret[iibv] = ret_val.lambda;
                  mix_ret[(maxcompc+1)*iibv] = ret_val.w0;

                  for(l = 0; l < ord; l++){
                     mix_ret[l+1+(maxcompc+1)*iibv] = ret_val.mix[l];
                     orient_ret[0+2*l+2*maxcompc*iibv] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*iibv] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+(maxcompc+1)*iibv] = 0;
                  }

                  sigma2_ret[iibv] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end iibv
}
void mixtrl0b( int* n1, int* siind, double* wi, int* ngrad, int* maxcomp, int* maxit,
          double* grad_in, double* bv_in, double* lambda_in, double* alpha_in,
          double* factr, double* penIC, double* sigma2, double* vert,
          double* si_in, double* sigma2_ret, double* orient_ret,
          int* order_ret, double* mix_ret){
// mixtensor: prolate tensors with isotropic compartment, fixed eigenvalues
// optmethod: L-BFGS-B
// n1 - number of voxel
// siind - initial parameters for orientations
// ngrad - number of gradient directions
// maxcomp - maximum number of mixture components
// maxit - maxium number of iterations in L-BFGS-B
// grad_in - array of gradients (3,ngrad)
// bv_in   - vector of b-values
// lambda_in - lambda_2
// alpha_in - (lambda1-lambda_2)/lambda_2
// factr  - control-parameter in  L-BFGS-B (default 1e7)
// penIC   - penalty for model complexity
// sigma2  - vector of  residual variances for best model with initial estimates
// vert    - vertices in polyeder used for initioal orientations
// si_in  - array of si/s0 (ngrad,n1)
// sigma2_ret - estimated residual variances
// orient_ret - estimated orientations (2,maxcomp,n1)
// order_ret - estimated model order (n1)
// mix_ret - estimated mixture coefficients (maxcomp,n1)

   int maxcompc = *maxcomp;
   int trace = 0, nREPORT = 1, lmm = 5;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init;
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, pgtol = 0;
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0,
          *param_tmp = 0, *si = 0, *lower = 0, *upper = 0;
   char msg[60];

   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim

   //Setting global variables
   alpha = *alpha_in; lambda = *lambda_in;;
   si_init = si_in; gradients = grad_in, ngradcc = *ngrad; bv = bv_in;


   //Gradient vectors corresponding to minima in spherical coordinates

   param_length_init=3*maxcompc + 1;

   param = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_work = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_last = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   lower = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   upper = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   int *nbd = (int *) R_alloc((long unsigned int) param_length_init, sizeof(int));
   si = (double *) R_alloc((long unsigned int) ngradcc, sizeof(double));

   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }

   for(iibv = 0; iibv < *n1; iibv++){
       for(i=0; i<param_length_init; i++){
          lower[i] = R_NegInf;
          upper[i] = R_PosInf;
	        nbd[i] = 0;
       }
      param_length = param_length_init;
      sigma2_ret[iibv] = sigma2[iibv];

         ord=maxcompc+1;

         for (j = 0; j < maxcompc; j++){
            iv = siind[j+maxcompc*iibv]-1;

            if(iv==-1) iv = j; // this should never happen

            dir[0]=vert[iv*3];  //vert dim(max(siind),3)
            dir[1]=vert[1+iv*3];
            dir[2]=vert[2+iv*3];

            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*iibv] = angles[0];
            orient_ret[1+2*j+2*maxcompc*iibv] = angles[1];

            param[3*j] = wi[j+1+(maxcompc+1)*iibv];
//initial value for w(j) corresponds to volume w(j)/(1+sum_k w(j)), i.e. a 17% isotropic compartment
            param[3*j+1] = angles[0];
            param[3*j+2] = angles[1];
// set values for lower corresponding to weights to 0
            lower[3*j] = 0;
            upper[3*j] = 1;
	          nbd[3*j] = 2;
         }
// lower value for alpha to get lambda_1>=lambda_2
         lower[3*maxcompc] = 0;//w0
         upper[3*maxcompc] = 1;//w0
         nbd[3*maxcompc] = 2;
         // initialize w0
         param[3*maxcompc] = wi[(maxcompc+1)*iibv];

         sigma_init = sigma2[iibv];
         krit = log(sigma_init) + penIC[0];

//          Rprintf("param from orient:\n");

         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
//             Rprintf(" %f ", param[l]);
         }

//          Rprintf("\n");

         // use AIC/ngrad, BIC/ngrad or AICC/ngrad respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length-1; l++){
                     param_tmp[l] = param_work[l];
                  }
                  param_tmp[3*k] = param_work[3*k+3];
                  param_length=3*k+1;
                  lower[3*k] = 0.;
                  upper[3*k] = 1.; // si are standardized by maxs0
                  nbd[3*k] = 2;

                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
                  }
               }


               //neuer check wenn elses drin sind

               optimex0 myoptimpar;

               //calculation and estimation of the best fitting model

               // initialize ret_val to avoid warnings
               mfunrskmb0_ret ret_val;
               ret_val.order = 0;
               ret_val.w0 = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;

               for(l=0;l<ngradcc;l++){
                  si[l] = si_init[l+iibv*ngradcc];
               }

               lbfgsb(param_length, lmm, param_work, lower, upper, nbd,
                      &Fmin, rskmixb0, drskmb0, &fail, &myoptimpar, *factr, pgtol,
                      &fncount, &grcount, *maxit, msg, trace,  nREPORT);

               ret_val = getparam0b(param_length, param_work, Fmin);

               value = ret_val.value;
               ord = ret_val.order;

               if(ord < k)
               {
                  ttt = krit;

                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngradcc;

                  ttt = log(si2new) + penIC[ord];

                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }


               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){
                  krit = ttt;

                  order_ret[iibv] = ret_val.order;

                  mix_ret[(maxcompc+1)*iibv] = ret_val.w0;
                  for(l = 0; l < ord; l++){
                     mix_ret[l+1+(maxcompc+1)*iibv] = ret_val.mix[l];
                     orient_ret[2*l+2*maxcompc*iibv] = ret_val.orient[l*2];
                     orient_ret[1+2*l+2*maxcompc*iibv] = ret_val.orient[1+l*2];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+(maxcompc+1)*iibv] = 0;
                  }

                  sigma2_ret[iibv] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end iibv
}
