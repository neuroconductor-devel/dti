#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

int ngrad0c = 0;
int i1 = 0;
double* siq_init, *penc, *grad, pen = 100;

extern void F77_NAME(mfunpl0)(double* param, double* siq, double* grad, int* m,
             int* lpar, int* ngrad, double* pen, double* tmp,
             double* tmp2, double* result);

extern void F77_NAME(mfunpl0g)(double* param, double* siq, double* grad,
            int* m, int* lpar, int* ngrad,
            double* tmp1,     //(ngrad*m),# z(n,m)
            double* tmp2,//(m*m),# v(m,m)
            double* tmp3,//(ngrad),# w(n) need w((m+1):n) for solver dgelsy
            double* tmp4,//(ngrad*m),# dkgj(n,m)
            double* tmp5,//(ngrad*m),# dkgj2(n,m)
            double* tmp6,//(ngrad*m),# ddkdphig(n,m)
            double* tmp7,//(ngrad*m),# ddkdetag(n,m)
            double* tmp8,//(m*m),# dvdth(m,m)
            double* tmp9,//(m*m*m),# dvdphi(m,m,m)
            double* tmp10,//(m*m*m),# dvdeta(m,m,m)
            double* tmp11,//(ngrad*m*3),# dzdpars(n,m,3)
            double* tmp12,//(m*lpar),# dwdpars(m,lpar)
            double* tmp13,//(m*lpar),# dwdpars2(m,lpar)
            double* tmp14,//(ngrad*m),# zs(n,m)
            double* tmp15,//(ngrad*m),# work1(n,m)
            double* tmp16,//(ngrad*m),# work2(n,m)
            double* tmp17,//(ngrad),# scopy(n)
            double* pen,
            double* result //result dfdpar double(lpar)
         );


extern void F77_NAME(mfunpl0w)(double* param, double* w, double* siq,
                double* grad, int* m, int* param_length,
                int* ngrad, double* tmp1, double* result);


typedef struct
{
  int ngrad;
  double* siq;
  double* grad;
  double* pen;
  int fnscale;
} optimex;


typedef struct{
  int order;
  double lev_par;
  double lev_sum_w;
  double* mix;
  double* orient;
  double* param;
  double value;
  double* w;
} mfunplwghts0_ret;

/* C++ code
public static T[] SubArrayDeepClone<T>(this T[] data, int index, int length){
  // ggf. fuer schoenes array kopieren (https://stackoverflow.com/questions/943635/c-arrays-getting-a-sub-array-from-an-existing-array)
    T[] arrCopy = new T[length];
    Array.Copy(data, index, arrCopy, 0, length);
    using (MemoryStream ms = new MemoryStream())
    {
        var bf = new BinaryFormatter();
        bf.Serialize(ms, arrCopy);
        ms.Position = 0;
        return (T[])bf.Deserialize(ms);
    }
    }*/


static R_INLINE int get_ind2d(int x, int y, int dim){
  return x + y*dim;
}


int compare_doubles (const void *a, const void *b){
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da < *db) - (*da > *db);
}


int contains(double* array, int length, double value){
  int i;
  for(i = 0; i < length; i++){
    if(value == array[i]){
      return 1;
    }
  }
  return 0;
}

int contains_int(int* array, int length, int value){
  int i;
  for(i = 0; i < length; i++){
    if(value == array[i]){
      return 1;
    }
  }
  return 0;
}



void paroforient(double *dir, double *angles){
  angles[0] = acos(dir[2]);        // angles[0]=theta; angles[1]=phi
  double sth = sin(angles[0]);
  angles[1] = 0;

  if (sth < 1e-8){
    angles[0] = 0;
  }else{
    double z = dir[0]/sth;
    if(fabs(z)>=1e0){
      angles[1] = z < 0 ? 0 : M_PI;
    }else{
      angles[1] = acos(z) * sign(dir[1]);
    }
    if(angles[1]<0) angles[1]=angles[1]+2*M_PI;
  }

  return;
}


double mfunpl0(int param_length, double *param, void* ex){

  int m = (param_length-1)/2;
  double* z = Calloc(ngrad0c*m, double);
  double* w = Calloc(ngrad0c, double);
  double result = 0;

  double* siq = Calloc(ngrad0c,double);
  int i;
  for(i=0;i<ngrad0c;i++){
    siq[i] = siq_init[i+i1*ngrad0c];
  }

  F77_CALL(mfunpl0)(param, siq, grad, &m, &param_length, &ngrad0c,
          &pen, z, w, &result);


  //check for infinity
  if(result == R_PosInf || result == R_NegInf){
    return 0;
  }

  Free(z);
  Free(w);
  Free(siq);

  return result;
}




void gmfunpl0(int param_length, double* param, double* result, void* ex){

  int m = (param_length-1)/2;

  double* z = Calloc(ngrad0c*m, double);
  double* v = Calloc(m*m, double);
  double* w = Calloc(ngrad0c, double);
  double* dkgj = Calloc(ngrad0c*m, double);
  double* dkgj2 = Calloc(ngrad0c*m, double);
  double* ddkdphig = Calloc(ngrad0c*m, double);
  double* ddkgetag = Calloc(ngrad0c*m, double);
  double* dvdth = Calloc(m*m, double);
  double* dvdphi = Calloc(m*m*m, double);
  double* dvdeta = Calloc(m*m*m, double);
  double* dzdpars = Calloc(ngrad0c*m*3, double);
  double* dwdpars = Calloc(m*param_length, double);
  double* dwdpars2 = Calloc(m*param_length, double);
  double* zs = Calloc(ngrad0c*m, double);
  double* work1 = Calloc(ngrad0c*m, double);
  double* work2 = Calloc(ngrad0c*m, double);
  double* scopy = Calloc(ngrad0c, double);

  int i;
  double* siq = Calloc(ngrad0c, double);
  for(i=0;i<ngrad0c;i++){
   siq[i] = siq_init[get_ind2d(i,i1,ngrad0c)];
  }

  F77_CALL(mfunpl0g)(param, siq, grad, &m, &param_length, &ngrad0c,
           z, v, w, dkgj, dkgj2, ddkdphig, ddkgetag,
           dvdth, dvdphi, dvdeta, dzdpars, dwdpars, dwdpars2, zs, work1, work2, scopy,
           &pen, result);
  Free(siq);
  Free(z);
  Free(v);
  Free(w);
  Free(dkgj);
  Free(dkgj2);
  Free(ddkdphig);
  Free(ddkgetag);
  Free(dvdth);
  Free(dvdphi);
  Free(dvdeta);
  Free(dzdpars);
  Free(dwdpars);
  Free(dwdpars2);
  Free(zs);
  Free(work1);
  Free(work2);
  Free(scopy);
}




/**
 *  Fall mixtensor, BFGS/CG
 *
 */
mfunplwghts0_ret mfunplwghts0(int param_length, double* param, double* siq){

  mfunplwghts0_ret ret_val;
  int i,j;

  int m = (param_length-1)/2, flag = 0, ord = 0;

  if (param[0] < 0){
    param[0] = 0;
  }

  double pen = 100;
  double tmp2 = 0;
  double* param_work = Calloc(param_length, double);
  double* tmp = Calloc(ngrad0c*m, double);
  double* w = Calloc(ngrad0c, double);
  double* w_red = Calloc(m, double);
  double* w_red2 = Calloc(m, double);

  //copy param
   for( i = 0; i < param_length; i++){
      param_work[i] = param[i];
   }

  //initialize w
  for(i=0;i<ngrad0c;i++){
      w[i] = 0;
  }

  F77_CALL(mfunpl0)(param_work, siq, grad, &m, &param_length, &ngrad0c,
          &pen, tmp, w, &tmp2);

  tmp2 = 0;
//  Free(tmp);
//  tmp = Calloc(ngrad0c*m, double);

  for(i = 0; i < m; i++){
    if(w[i] > 0){
      w_red[i] = w[i];
    }else{
      w_red[i] = 0;
    }
  }

  for(i = 0; i < m; i++){
    if(w_red[i] == R_PosInf || w_red[i] == R_NegInf){
      flag = 1;
    }
  }


  if(flag == 0){
    F77_CALL(mfunpl0w)(param_work, w_red, siq, grad, &m, &param_length, &ngrad0c, tmp, &tmp2);
  }else{
    tmp2 = pow(10,20);
  }


  for(i = 0; i < m; i++){
   w_red2[i] = w_red[i];
  }

  qsort(w_red, (long unsigned int) m, sizeof(double), compare_doubles);
  int* o = Calloc(m, int);
  for(i = 0; i<m; i++){
   o[i] = -1;
  }

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      if(w_red2[j] == w_red[i] && contains_int(o,m,j) == 0){
      o[i] = j;
      j = m;
      }
    }
  }

  double sum_w = 0;
  for (i = 0; i < m; i++){
    if(w_red[i] > 0){
      ord++;
      sum_w += w_red[i];
    }
  }


  /** correct order:
    * if order == 0, C for-loops wouldn't start
    * and arrays wouldn't be allocated
    * R solves this with o <- o[1:ord] = o[1:0] = o[1]
    * so c_ord will be 1 iff ord == 0, else c_ord = ord
    **/
  int c_ord = ord;
  if (c_ord == 0) c_ord = 1;

  // here ord will not be replaced, because zero-check is important
  double* mix = (double*) R_alloc((long unsigned int) ord, sizeof(double));

  if(ord > 0){
   for(i = 0; i < ord; i++){
      mix[i] = w_red[i]/sum_w;
   }
  }else{
   mix = NULL;
  }

  double* or = (double*) R_alloc(2*(long unsigned int) c_ord, sizeof(double));
  for(i = 0; i < c_ord; i++){
   or[get_ind2d(0,i,2)] = param[2*o[i]+1];
   or[get_ind2d(1,i,2)] = param[2*o[i]+2];
  }

  for(j = 0; j < c_ord; j++){
   while (or[get_ind2d(0,j,2)] < 0) {
      or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] + M_PI;
   }
   while (or[get_ind2d(0,j,2)] > M_PI) {
      or[get_ind2d(0,j,2)] = or[get_ind2d(0,j,2)] - M_PI;
   }
   while (or[get_ind2d(1,j,2)] < 0) {
      or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] + 2 * M_PI;
   }
   while (or[get_ind2d(1,j,2)] > 2 * M_PI) {
      or[get_ind2d(1,j,2)] = or[get_ind2d(1,j,2)] - 2 * M_PI;
   }
  }

  for(i = 0; i < c_ord; i++){       //par[0] stays
   param[2*i+1] = or[2*i];
   param[2*i+2] = or[2*i+1];

  }

  ret_val.lev_par = param[0];  //lev[1]
  ret_val.lev_sum_w = -log(sum_w);   //lev[2]
  ret_val.order = ord;
  ret_val.param = param;
  ret_val.value = tmp2;
  ret_val.orient = or;
  ret_val.mix = mix;

  Free(tmp);
  Free(param_work);
  Free(w);
  Free(w_red);
  Free(w_red2);
  Free(o);

  return ret_val;
}

void mixture( int* n1, int* siind, int* ngrad0, int* maxcomp, int* maxit,
          double* pen, double* grad_in, double* reltol,
          double* th, double* penIC, double* sigma2, double* vert,
          double* siq_in,
          double* sigma2_ret, double* orient_ret,
          int* order_ret, double* lev_ret, double* mix_ret){


   // mixtensor optmethod: BFGS
   //
   // i1 and n1,  current/total respectively
   // ngrad : # gradients
   // maxcomp : maximal order of mix ???

   int maxcompc = *maxcomp, *mask;
   int trace = 0, nREPORT = 1;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int ord, iv, i, j, k, l, param_length, param_length_init, tmp_int;
   double krit, sigma_init, si2new = 0, value = 0;
   double ttt = 0, abstol = R_NegInf;
//   double alpha = 1, beta = 0.5, gamma = 2; //Optim parameters for Nelder Mead
   double angles[2], dir[3];
   double *param = 0, *param_work = 0, *param_last = 0,
   *param_tmp = 0, *siq = 0;//, *w = 0;

   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim

   //Setting global variables
   siq_init = siq_in; penc = pen; grad = grad_in, ngrad0c = *ngrad0;

   // not yet initialized
//   w = (double*) R_alloc(maxcompc, sizeof(double));

   //Gradient vectors corresponding to minima in spherical coordinates

   param_length_init=2*maxcompc+1;

   param = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_work = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_last = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   param_tmp = (double *) R_alloc((long unsigned int) param_length_init, sizeof(double));
   mask = (int *) R_alloc((long unsigned int) param_length_init, sizeof(int));
   for (i = 0; i < param_length_init; i++) mask[i] = 1;

   siq = (double *) R_alloc((long unsigned int) ngrad0c, sizeof(double));

   //initialize param
   for(i=0; i<param_length_init; i++){
       param[i] = 0;
   }

   for(i1 = 0; i1 < *n1; i1++){
      param_length = param_length_init;
      sigma2_ret[i1] = sigma2[i1];

         ord=maxcompc+1;

         for (j = 0; j < maxcompc; j++){
            iv = siind[(j+2)+(maxcompc+2)*i1]-1;

            if(iv==-1) iv = j; // this should never happen

            dir[0]=vert[get_ind2d(0,iv,3)];  //vert dim(max(siind),3)
            dir[1]=vert[get_ind2d(1,iv,3)];
            dir[2]=vert[get_ind2d(2,iv,3)];

            paroforient(dir, angles);
            orient_ret[0+2*j+2*maxcompc*i1] = angles[0];
            orient_ret[1+2*j+2*maxcompc*i1] = angles[1];

            param[2*j+1] = angles[0];
            param[2*j+2] = angles[1];
         }

         // initialize EV-parameter
         tmp_int = siind[1+(maxcompc+2)*i1];
         if(tmp_int > 0){
            param[0] = th[tmp_int-1];
         }else{
            param[0] = 0.001;
         }

         sigma_init = sigma2[i1];
         krit = log(sigma_init) + penIC[0];

//          Rprintf("param from orient:\n");

         for(l = 0; l < param_length; l++){
            param_work[l] = param[l];
            param_last[l] = param[l];
//             Rprintf(" %f ", param[l]);
         }

//          Rprintf("\n");

         // use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
         for(k = maxcompc; k > 0; k--){
            if(k < ord){
               if(k != maxcompc){
                  for (l = 0; l < param_length; l++){
                     param_tmp[l] = param_work[l];
                  }
                        param_length=2*k+1;

                  for(l= 0; l < param_length; l++){
                     param_work[l] = param_tmp[l];
                     param_last[l] = param_tmp[l];
//                   Rprintf(" %f ", param_work[l]);
                  }
//                   Rprintf("\n");
               }


               //neuer check wenn elses drin sind

               optimex myoptimpar;

               //calculation and estimation of the best fitting model

               // initialize ret_val to avoid warnings
               mfunplwghts0_ret ret_val;
               ret_val.order = 0;
               ret_val.lev_par = 0;
               ret_val.lev_sum_w = 0;
               ret_val.mix = NULL;
               ret_val.orient = NULL;
               ret_val.param = NULL;
               ret_val.value = 0;
               ret_val.w = NULL;

               for(l=0;l<ngrad0c;l++){
                  siq[l] = siq_init[get_ind2d(l,i1,ngrad0c)];
               }

               vmmin(param_length, param_work, &Fmin, mfunpl0, gmfunpl0,
                     *maxit, trace, mask, abstol, *reltol, nREPORT,
                     &myoptimpar, &fncount, &grcount, &fail);

               ret_val = mfunplwghts0(param_length, param_work, siq);
 //            Rprintf("ret_val %f ", ret_val.value);
 //            Rprintf("\n");

               value = ret_val.value;

               ord = ret_val.order;

               if( (ret_val.lev_par < 0) || (ret_val.lev_sum_w < 0) || (ord < k) )
               {
                  ttt = krit;

                  for(l =0; l < param_length; l++){
                     param_work[l] = param_last[l];
                  }
               }
               else
               {
                  si2new = value/ngrad0c;

                  ttt = log(si2new) + penIC[ord];

                  for(l =0; l < param_length; l++){
                     param_work[l] = ret_val.param[l];
                  }
               }


               //use directions corresponding to largest weights as initial directions
               if(ttt < krit){
                  krit = ttt;

                  order_ret[i1] = ret_val.order;

                  lev_ret[0+2*i1] = ret_val.lev_par;
                  lev_ret[1+2*i1] = ret_val.lev_sum_w;

                  for(l = 0; l < ord; l++){
                     mix_ret[l+maxcompc*i1] = ret_val.mix[l];
                     orient_ret[0+2*l+2*maxcompc*i1] = ret_val.orient[get_ind2d(0,l,2)];
                     orient_ret[1+2*l+2*maxcompc*i1] = ret_val.orient[get_ind2d(1,l,2)];
                  }
                  for(l = ord; l < maxcompc; l++){
                     mix_ret[l+maxcompc*i1] = 0;
                  }

                  sigma2_ret[i1] = si2new;
               } // end if(ttt < krit)
            } // end if(k < ord)
         } // end model order
      R_CheckUserInterrupt();
   } // end i1
}
