#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

//sumweight over j
SEXP w2(SEXP w,SEXP id,SEXP gno ){
  int i,I=length(w),J=INTEGER(gno)[0];
  SEXP sw;
  PROTECT(sw=allocVector(REALSXP,J));
  for(i=0;i<J;i++) REAL(sw)[i]=0.0;
  for(i=0;i<I;i++) REAL(sw)[INTEGER(id)[i]-1]+=REAL(w)[INTEGER(id)[i]-1]
             *REAL(w)[INTEGER(id)[i]-1];
  UNPROTECT(1);

  return(sw);
}

//calculate mean of y by unique geneid;
void ymean(double *y,int *gno,int *gid,int *obs,int *nj,double *ymean ){

int j,I=*obs,J=*gno;
for(j=0;j<J;j++) ymean[j]=0.0;     //initialize ymean;
for(j=0;j<I;j++) ymean[gid[j]-1] +=y[j];  //sum ratio over j=1-J;
for(j=0;j<J;j++) ymean[j] /=(double)nj[j];        //take mean;

}

SEXP Wymean(SEXP y, SEXP id,SEXP gno,SEXP sumw ){

int j,J=INTEGER(gno)[0],I=length(y);
SEXP mean;

 PROTECT(mean=allocVector(REALSXP,J));

 for(j=0;j<J;j++) REAL(mean)[j]=0.0;    //initialize mean;
 for(j=0;j<I;j++) REAL(mean)[INTEGER(id)[j]-1]+=1.0/
    (REAL(sumw)[INTEGER(id)[j]-1]==0.0?1.0:REAL(sumw)[INTEGER(id)[j]-1])*REAL(y)[j];

 UNPROTECT(1);

 return(mean);
}

//residual function;
SEXP Yresid(SEXP y,SEXP beta,SEXP id ){

 int i,I=length(y);
 SEXP r;

 PROTECT(r=allocVector(REALSXP,I));
 for(i=0;i<I;i++) REAL(r)[i]=0.0;

 for(i=0;i<I;i++) REAL(r)[i]=REAL(y)[i]-REAL(beta)[INTEGER(id)[i]-1];
 UNPROTECT(1);

 return(r);
}

//Anova mean,balanced design;
SEXP Anovamean(SEXP y,SEXP id,SEXP gno,SEXP firstg ){
 int i,I=INTEGER(gno)[0],fg=INTEGER(firstg)[0];
 double ysum=0.0,meanf=0.0;
 SEXP mean;

 PROTECT(mean=allocVector(REALSXP,I));
 for(i=0;i<I;i++) REAL(mean)[i]=0.0;
 for(i=0;i<length(y);i++) {ysum+=REAL(y)[i]; REAL(mean)[INTEGER(id)[i]-1]+=REAL(y)[i];}
 for(i=1;i<I;i++) {REAL(mean)[i]=REAL(mean)[i]/(fg*1.0)-ysum/(fg*I*1.0); meanf+=REAL(mean)[i];}  //calculate beta's;
 REAL(mean)[0]=-meanf;

 UNPROTECT(1);

 return(mean);
}

//Unbalanced design of Anova mean;
SEXP UAnovamean(SEXP y, SEXP nj, SEXP gno,SEXP id ){
 int i,I=length(y),j,J=INTEGER(gno)[0];
 double sumbeta=0.0,sym=0.0,N=0.0;
 SEXP beta;

 PROTECT(beta=allocVector(REALSXP,J));

 for(i=0;i<J;i++) {
        REAL(beta)[i]=0.0;
        N+=(REAL(nj)[i]==0.0)?0.0:(1.0/REAL(nj)[i]);
 }

 for(i=0;i<I;i++){
   sym+=((REAL(nj)[INTEGER(id)[i]-1]==0.0)?0.0:(1.0/REAL(nj)[INTEGER(id)[i]-1]))*REAL(y)[i];  //sum over all genes;
   REAL(beta)[INTEGER(id)[i]-1]+=(REAL(nj)[INTEGER(id)[i]-1]==0.0? 0.0:(1.0/REAL(nj)[INTEGER(id)[i]-1]))*REAL(y)[i];
  }
 for(j=1;j<J;j++) {
   REAL(beta)[j]-=(REAL(nj)[j]==0.0)?0.0:1.0/REAL(nj)[j]*sym/N;
   sumbeta+=REAL(beta)[j];
   }
 REAL(beta)[0]=-sumbeta;

 UNPROTECT(1);

 return(beta);
}

//Anova b mean function for unbalance design;
void UAnovabmean(double *b,int *nr,int *nc,int *gid,int *gno,double *nj,
                double *bmean ){

int i,j,k,J=*gno,K=*nc,row=*nr;
double b1=0.0,sym=0.0,N=0.0;

  for(j=0;j<J;j++){
    N+=(nj[j]==0.0)?0.0:(1.0/nj[j]);     //calculate N
   for(k=0;k<K;k++) bmean[j+k*J]=0.0;  //initialize bmean;
   }
  for(k=0;k<K;k++){
     for(i=0;i<row;i++){
       bmean[gid[i]-1+k*J]+=(nj[gid[i]-1]==0.0?0.0:1.0/nj[gid[i]-1])*b[i+k*row]; //calculate mean of genes
       sym+=bmean[gid[i]-1+k*J];
      }
      bmean[k*J]=sym;
      sym=0.0;
  }

  for(k=0;k<K;k++){
    for(j=1;j<J;j++){
      bmean[j+k*J]-=(nj[j]==0.0)?0.0:(1.0/nj[j]*b[k*J]/N);  //sum weighted b over j for all ks;
      b1+=bmean[j+k*J];
    }
   bmean[k*J]=-b1;
   b1=0.0;
  }

}

//Bresid
SEXP Bresid(SEXP b,SEXP bmean,SEXP id ){
  int i,I=nrows(b),j,J=nrows(bmean),K=ncols(b);
  SEXP r;
  PROTECT(r=allocMatrix(REALSXP,I,K));

  for(i=0;i<K;i++)
  for(j=0;j<I;j++) REAL(r)[j+i*I]=0.0;

  for(i=0;i<K;i++)
    for(j=0;j<I;j++) REAL(r)[j+i*I]=REAL(b)[j+i*I]-REAL(bmean)[INTEGER(id)[j]-1+i*J];  //bresid;

  UNPROTECT(1);

  return(r);
}

//weighted y residual;
SEXP Wyresid(SEXP y,SEXP id,SEXP idno,SEXP w ){
 int i,I=length(y),J=INTEGER(idno)[0];
 SEXP mean,sw,r;

 PROTECT(mean=allocVector(REALSXP,J));
 PROTECT(sw=allocVector(REALSXP,J));
 PROTECT(r=allocVector(REALSXP,I));

 for(i=0;i<I;i++) REAL(r)[i]=0.0;
 for(i=0;i<J;i++) {REAL(mean)[i]=0.0;REAL(sw)[i]=0.0;}

 for(i=0;i<I;i++) {
     REAL(mean)[INTEGER(id)[i]-1]+=REAL(w)[i]*REAL(y)[i];
     REAL(sw)[INTEGER(id)[i]-1]+=REAL(w)[i];
  }
 for(i=0;i<I;i++) REAL(r)[i]=REAL(y)[i]*sqrt(REAL(w)[i])
           -sqrt(REAL(w)[i])*(REAL(mean)[INTEGER(id)[i]-1]/REAL(sw)[INTEGER(id)[i]-1]);

 UNPROTECT(3);

 return (r);
}

//determine Huber's robust weight;
void Huber(double *resid,double *weight,int *obs,double *sigma,double *k ){

  int i,I=*obs;
  double s=*sigma,K=*k;
  for(i=0;i<I;i++) weight[i]=0.0;

 for(i=0;i<I;i++)
     weight[i]=fabs(resid[i])>K*s ?(K/(fabs(resid[i]/s))):1.0;
}

//Huber scale calculation;
SEXP Huberscale(SEXP r,SEXP p,SEXP initsigma,SEXP hconst ){

 int i;
 double beta,d=REAL(hconst)[0],tmp=0.0;
 SEXP newsigma;

 PROTECT(newsigma=allocVector(REALSXP,1));
 REAL(newsigma)[0]=0.0;

 beta=pchisq(d*d,3.0,1,0)+d*d*(1.0-pchisq(d*d,1.0,1,0));  //constant in huber's scale function;
 for(i=0;i<length(r);i++){
        if(REAL(r)[i]<-d*REAL(initsigma)[0]) tmp=-d*REAL(initsigma)[0];
  else if(REAL(r)[i]>d*REAL(initsigma)[0]) tmp=d*REAL(initsigma)[0];
  else tmp=REAL(r)[i];

  REAL(newsigma)[0]+=tmp*tmp;
 }

 REAL(newsigma)[0]/=((length(r)-INTEGER(p)[0])*beta);

 UNPROTECT(1);

 return(newsigma);
}
//Huber weight function; SOME PROBLEM OF THIS FUNCTION BECAUSE NO SCALE PARAMETERS;
SEXP HuberWeight(SEXP r ){
 int i, I=length(r);
 SEXP w;

 PROTECT(w=allocVector(REALSXP,I));
 for(i=0;i<I;i++) REAL(w)[i]=0.0;

 for(i=0;i<I;i++)
  REAL(w)[i]=(fabs(REAL(r)[i]))>2.5 ?
    (2.5/(fabs(REAL(r)[i]))):1.0;

 UNPROTECT(1);

 return(w);
}

//Tukeyscale function; tconst use 2.5
SEXP Tukeyscale(SEXP r,SEXP pp,SEXP initsigma,SEXP tconst ){

 int i,N=length(r);
 double beta=0.0,K=REAL(tconst)[0],w=0.0,a=0.0,b=0.0,c=0.0,d=0.0,p=0.0,q=0.0;
 SEXP newsigma;

 PROTECT(newsigma=allocVector(REALSXP,1));
 REAL(newsigma)[0]=0.0;

 beta=3.0/(K*K)*pchisq(K*K,3.0,1,0)-9.0*pchisq(K*K,5.0,1,0)/(K*K*K*K)
     +15.0*pchisq(K*K,7.0,1,0)/(K*K*K*K*K*K)
     +2.0*(1.0-pnorm(K,0.0,1.0,1,0));       //constant in Tukey's scale function;

 for(i=0;i<N;i++){
  if(fabs(REAL(r)[i])>K*REAL(initsigma)[0]) w+=1.0;
  else {
         a+=pow(REAL(r)[i]/K,6.0);
         b+=-3.0*pow(REAL(r)[i]/K,4.0);
         c+=3.0*pow(REAL(r)[i]/K,2.0);
       }
 }

 d=w-(N-INTEGER(pp)[0])*1.0*beta;
 p=c/a-1.0/3.0*pow(b/a,2.0);
 q=d/a-b*c/(3.0*a*a)+2.0*b*b*b/(27.0*a*a*a);

 REAL(newsigma)[0]=1.0/sqrt(pow(-q/2.0+sqrt(q*q/4.0+p*p*p/27.0),1.0/3.0)-pow(q/2.0+sqrt(q*q/4.0+p*p*p/27.0),1.0/3.0)-b/(3.0*a));
 UNPROTECT(1);

 return(newsigma);
}
//Based on Cardan formula

//determine Tukey's robust weight
SEXP TukeyWeight(SEXP resid,SEXP m,SEXP k ){
  int i,I=length(resid);
  SEXP weight;

  PROTECT(weight=allocVector(REALSXP,I));

  for(i=0;i<I;i++) REAL(weight)[i]=0.0;

  for(i=0;i<I;i++){
    REAL(weight)[i]=fabs(REAL(resid)[i])<(REAL(k)[0]*REAL(m)[0]) ?
      pow(1.0-pow(REAL(resid)[i]/(REAL(k)[0]*REAL(m)[0]),2.0),2.0):0.0;

  }
   UNPROTECT(1);
 return(weight);
}

//calculate residual of y;
SEXP resid(SEXP y, SEXP gid, SEXP ymean ){

int i,I=length(y);
SEXP r;

PROTECT(r=allocVector(REALSXP,I));

for(i=0;i<I;i++) REAL(r)[i]=0.0;   //initialize yresid;
for(i=0;i<I;i++) REAL(r)[i]=REAL(y)[i]-REAL(ymean)[INTEGER(gid)[i]-1];

UNPROTECT(1);

return(r);
}

//fitted function;
SEXP fittedvalue(SEXP y,SEXP beta,SEXP id ){

 int i,I=length(y);
 SEXP fit;

 PROTECT(fit=allocVector(REALSXP,I));
 for(i=0;i<I;i++) REAL(fit)[i]=0.0;

 for(i=0;i<I;i++)
   REAL(fit)[i]=REAL(y)[i]+REAL(beta)[INTEGER(id)[i]-1];

 UNPROTECT(1);

 return(fit);
}

//score calculation;
SEXP score(SEXP D, SEXP gammae,SEXP r ){
 int i,I=length(r),k,K=ncols(D);
 SEXP s;

 PROTECT(s=allocVector(REALSXP,K));
 for(k=0;k<K;k++) REAL(s)[k]=0.0;
 for(k=0;k<K;k++)
  for(i=0;i<I;i++)
   REAL(s)[k]+=REAL(r)[i]*REAL(r)[i]*REAL(D)[i+I*k]*REAL(gammae)[i]-REAL(D)[i+I*k];

  UNPROTECT(1);

  return(s);
}
//-Hessian matrix for Newton method
SEXP Hessian(SEXP D, SEXP gammae, SEXP r ){
 int i,j,I=length(r),k,K=ncols(D);
 SEXP h;

 PROTECT(h=allocMatrix(REALSXP,K,K));
 for(k=0;k<K;k++)
   for(i=0;i<K;i++) REAL(h)[i+k*K]=0.0;
//caculate upper trangle of Hessian matrix
 for(k=0;k<K;k++)
   for(j=0;j<k+1;j++)
     for(i=0;i<I;i++)
       REAL(h)[j+k*K]+=REAL(r)[i]*REAL(r)[i]*REAL(gammae)[i]*REAL(D)[i+k*I]*REAL(D)[i+j*I];
//caculate lower trangle of Hessian matrix
 for(k=0;k<K;k++)
   for(j=k+1;j<K;j++) REAL(h)[j+k*K]=REAL(h)[k+j*K];

 UNPROTECT(1);

 return(h);
}

//score calculation for robust regression;
SEXP scoreR(SEXP D, SEXP gammae,SEXP r,SEXP robustname ){
 int i,I=length(r),k,K=ncols(D);
 double tmp=0.0;
 SEXP s;

 PROTECT(s=allocVector(REALSXP,K));
 for(k=0;k<K;k++) REAL(s)[k]=0.0;
 for(k=0;k<K;k++){
  for(i=0;i<I;i++){
     if(INTEGER(robustname)[0]==1){
       if(REAL(r)[i]*sqrt(REAL(gammae)[i])>2.5) tmp=2.5;
       else if(REAL(r)[i]*sqrt(REAL(gammae)[i])< -2.5) tmp=-2.5;
       else tmp=REAL(r)[i]*sqrt(REAL(gammae)[i]);
      }
     else if(INTEGER(robustname)[0]==2){
        if(fabs(REAL(r)[i]*sqrt(REAL(gammae)[i]))>6.0) tmp=0.0;
        else tmp=REAL(r)[i]*sqrt(REAL(gammae)[i])*
       pow(1.0-pow(REAL(r)[i]*sqrt(REAL(gammae)[i])/6.0,2.0),2.0);
      }
     REAL(s)[k]+=REAL(r)[i]*sqrt(REAL(gammae)[i])*REAL(D)[i+I*k]*tmp-REAL(D)[i+I*k];
  }
 }
 UNPROTECT(1);

  return(s);
}
//-Hessian matrix for Newton method for Robust regression
SEXP HessianR(SEXP D, SEXP gammae, SEXP r,SEXP robustname ){
 int i,j,I=length(r),k,K=ncols(D);
 double tmp=0.0,tmp1=0.0;
 SEXP h;

 PROTECT(h=allocMatrix(REALSXP,K,K));
 for(k=0;k<K;k++)
   for(i=0;i<K;i++) REAL(h)[i+k*K]=0.0;
//caculate upper trangle of Hessian matrix
 for(k=0;k<K;k++){
   for(j=0;j<k+1;j++){
     for(i=0;i<I;i++){
       if(INTEGER(robustname)[0]==1){
         if(REAL(r)[i]*sqrt(REAL(gammae)[i])>2.5){ tmp=2.5;tmp1=0.0;}
         else if(REAL(r)[i]*sqrt(REAL(gammae)[i])< -2.5){ tmp=-2.5;tmp1=0.0;}
         else{ tmp=REAL(r)[i]*sqrt(REAL(gammae)[i]); tmp1=1.0;}
        }
       else if(INTEGER(robustname)[0]==2){
          if(fabs(REAL(r)[i]*sqrt(REAL(gammae)[i]))>6.0){ tmp=0.0;tmp1=0.0;}
          else {tmp=REAL(r)[i]*sqrt(REAL(gammae)[i])*pow(1.0-pow(REAL(r)[i]*sqrt(REAL(gammae)[i])/(6.0),2.0),2.0);
             tmp1=1.0-6.0*pow(REAL(r)[i]*sqrt(REAL(gammae)[i])/(6.0),2.0)+
               5.0*pow(REAL(r)[i]*sqrt(REAL(gammae)[i])/(6.0),4.0);
                }
    }
       REAL(h)[j+k*K]+=(REAL(r)[i]*sqrt(REAL(gammae)[i])*REAL(D)[i+k*I]*REAL(D)[i+j*I]*tmp+
                         pow(REAL(r)[i]*sqrt(REAL(gammae)[i]),2.0)*REAL(D)[i+k*I]*REAL(D)[i+j*I]*tmp1);
     }
   }
 }
//caculate lower trangle of Hessian matrix
 for(k=0;k<K;k++)
   for(j=k+1;j<K;j++) REAL(h)[j+k*K]=REAL(h)[k+j*K];

 UNPROTECT(1);

 return(h);
}


//varaiance of beta;
SEXP betavar(SEXP bf,SEXP b,SEXP nj,SEXP gno ){

 int i,I=nrows(bf),j,K=ncols(bf),k,J=INTEGER(gno)[0];
 double tmp=0.0,tmp1=0.0,NN=0.0;
 SEXP varbeta;

 PROTECT(varbeta=allocVector(REALSXP,J));
 for(i=0;i<J;i++) {REAL(varbeta)[i]=0.0;NN+=(REAL(nj)[i]==0.0?0.0:1.0/REAL(nj)[i]);}

 //calculate diagno components: BF%*%B part;
 for(j=0;j<I;j++)
   for(k=0;k<K;k++) REAL(varbeta)[j+1]+=REAL(bf)[j+k*I]*REAL(b)[k+j*K];

 //calculate diagnal components: Z'Z part;
 for(i=1;i<J;i++){
   for(j=0;j<J;j++) REAL(varbeta)[i]+=(REAL(nj)[j]==0.0?0.0:1.0/REAL(nj)[j])*
            (REAL(nj)[i]==0.0?0.0:1.0/REAL(nj)[i])/NN;
   REAL(varbeta)[i]-=(REAL(nj)[i]==0.0)?0.0:
                      1.0/(REAL(nj)[i]*REAL(nj)[i]*NN);
   tmp+=REAL(varbeta)[i];
   }
  //calculate covariance of beta1;
  for(j=1;j<I;j++){
   for(i=0;i<j;i++){
     tmp1+=-(REAL(nj)[i]==0.0?0.0:1.0/REAL(nj)[i])*(REAL(nj)[j]==0.0?0.0:1.0/REAL(nj)[j])/NN;
     for(k=0;k<K;k++) tmp1+=REAL(bf)[i+k*I]*REAL(b)[k+j*K];
    //for balanced design tmp1+=-1.0/(INTEGER(n1)[0]*J*1.0);
    }
  }

 REAL(varbeta)[0]=tmp+2.0*tmp1;  //variance of beta1;
 UNPROTECT(1);

 return(varbeta);
 }
//variance calculation for Terry speed data control-treatment;
SEXP TSbetavar(SEXP bf,SEXP b,SEXP nj,SEXP gno ){

 int i,I=nrows(bf),j,K=ncols(bf),k,J=INTEGER(gno)[0];
 int halfJ=J/2;
 double tmp=0.0,tmp1=0.0,N=0.0;
 SEXP varbeta,covbeta,var;

 PROTECT(varbeta=allocVector(REALSXP,J));
 PROTECT(covbeta=allocVector(REALSXP,halfJ));

 for(i=0;i<J;i++) {
    REAL(varbeta)[i]=0.0;
    if(i<J/2) REAL(covbeta)[i]=0.0;
   }
 for(i=0;i<J;i++) N+=1.0/REAL(nj)[i];   //calculate N

 //calculate diagno components: BF%*%B part;
 for(j=0;j<I;j++)
   for(k=0;k<K;k++) REAL(varbeta)[j+1]+=REAL(bf)[j+k*I]*REAL(b)[k+j*K];

 //calculate diagnal components: Z'Z part;
 for(i=1;i<J;i++){
   for(j=0;j<J;j++) REAL(varbeta)[i]+=1.0/(REAL(nj)[j]*REAL(nj)[i]*N);
   REAL(varbeta)[i]-=1.0/(REAL(nj)[i]*REAL(nj)[i]*N);
   tmp+=REAL(varbeta)[i];
   }
  //calculate covariance of beta1;
  for(j=1;j<I;j++){
   for(i=0;i<j;i++){
     tmp1+=-1.0/(REAL(nj)[i]*REAL(nj)[j]*N);
     for(k=0;k<K;k++) tmp1+=REAL(bf)[i+k*I]*REAL(b)[k+j*K];
    //for balanced design tmp1+=-1.0/(INTEGER(n1)[0]*J*1.0);
    }
  }
  REAL(varbeta)[0]=tmp+2.0*tmp1;  //variance of beta1;

 //calculate covariance between betaj and betaj*,j=2,...J;
  tmp=0.0;tmp1=0.0;
  for(j=1;j<halfJ;j++){
//    for(i=halfJ+1;i<J;i++){
      REAL(covbeta)[j]+=-1.0/(REAL(nj)[j]*REAL(nj)[j+halfJ]*N);
     for(k=0;k<K;k++) REAL(covbeta)[j]+=REAL(bf)[j-1+k*I]*REAL(b)[k+(halfJ+j-1)*K];
  }
 //calculate covariance of beta1 and beta1*;
 for(i=0;i<halfJ-1;i++){        //J-1 terms
  tmp+=-1.0/(REAL(nj)[i+1]*REAL(nj)[halfJ]*N);
  for(k=0;k<K;k++) tmp1+=REAL(bf)[i+k*I]*REAL(b)[k+(halfJ-1)*K];
 }
 //Rprintf("tmp=%f\ttmp1=%f\n",tmp,tmp1);
 for(i=halfJ;i<J-1;i++){                               //J-1 terms
  tmp+=-1.0/(REAL(nj)[i+1]*REAL(nj)[halfJ]*N);
  for(k=0;k<K;k++) tmp1+=REAL(bf)[i+k*I]*REAL(b)[k+(halfJ-1)*K];
 }
 //Rprintf("tmp=%f\ttmp1=%f\n",tmp,tmp1);
 REAL(covbeta)[0]+=tmp-tmp1-REAL(varbeta)[halfJ];  //covariance between beta1 and beta1*;

 PROTECT(var=allocVector(VECSXP,2));   //output list of variance and covariance;
 SET_VECTOR_ELT(var,0,varbeta);       //first element is variance;
 SET_VECTOR_ELT(var,1,covbeta);       //second element is covariance;

 UNPROTECT(3);

 return(var);
 }

//robust varaiance of beta;
SEXP Rbetavar(SEXP bf,SEXP b,SEXP w ){
 int i,J=nrows(bf),j,K=ncols(bf),k;
 SEXP varbeta;

 PROTECT(varbeta=allocVector(REALSXP,J));
 for(i=0;i<J;i++) REAL(varbeta)[i]=0.0;

 for(j=0;j<J;j++)
   for(k=0;k<K;k++) REAL(varbeta)[j]+=REAL(bf)[j+k*J]*REAL(b)[k+j*K];
 for(i=0;i<J;i++) REAL(varbeta)[i]+=1.0/REAL(w)[i];

 UNPROTECT(1);

 return(varbeta);
 }


void ymean1(double *y,int *gno,int *gid,int *obs,int *nj,double *ymean,double *yresid ){

int i,j,J=*gno,I=*obs;
for(j=0;j<J;j++) ymean[j]=0.0;     //initialize ymean;
for(i=0;i<I;i++) yresid[i]=0.0;   //initialize yresid;

for(i=0;i<I;i++) ymean[gid[i]-1] +=y[i];  //sum ratio over j=1-J;
for(j=0;j<J;j++) ymean[j] /=nj[j];        //take mean;
for(i=0;i<I;i++) yresid[i]=y[i]-ymean[gid[i]-1];
}

void bmean(double *b,int *nr,int *nc,int *gid,int *nj,int *gno,double *bmean,double *bresid){

int i,j,k,J=*gno,K=*nc,row=*nr;

for(j=0;j<J;j++)
 for(k=0;k<K;k++) bmean[j+k*J]=0.0;  //initialize bmean;
for(i=0;i<row;i++)
 for(k=0;k<K;k++) bresid[i+k*row]=0.0;   //initialize bresid;

for(i=0;i<row;i++)
 for(k=0;k<K;k++) bmean[gid[i]-1+k*J]+=b[i+k*row];  //sum b over j;
for(j=0;j<J;j++)
 for(k=0;k<K;k++) bmean[j+k*J]/=nj[j];             //take mean;

for(i=0;i<row;i++)
 for(k=0;k<K;k++) bresid[i+k*row]=b[i+k*row]-bmean[gid[i]-1+k*J];  //bresid;

}

//Anova b mean function;
void Anovabmean(double *b,int *nr,int *nc,int *gid,int *gno,int *fgno,double *bmean,double *bresid
 ){

int i,j,k,J=*gno,K=*nc,row=*nr;
double ysum,msum;

for(j=0;j<J;j++)
 for(k=0;k<K;k++) bmean[j+k*J]=0.0;  //initialize bmean;
for(i=0;i<row;i++)
 for(k=0;k<K;k++) bresid[i+k*row]=0.0;   //initialize bresid;

for(k=0;k<K;k++){
  ysum=0.0;msum=0.0;
  for(i=0;i<row;i++){
    bmean[gid[i]-1+k*J]+=b[i+k*row];  //sum b over j;
    ysum+=b[i+k*row];
  }
  for(j=1;j<J;j++) {bmean[j+k*J]=bmean[j+k*J]/(*fgno)-ysum/(*fgno*J);msum+=bmean[j+k*J];}
  bmean[k*J]=-msum;        //the first gene mean;
 }

for(i=0;i<row;i++)
 for(k=0;k<K;k++) bresid[i+k*row]=b[i+k*row]-bmean[gid[i]-1+k*J];  //bresid;

}

//Weight B mean;
void Wbmean(double *b,int *nr,int *nc,int *gid, double *w,int *gno,double *sweight,
double *bmean,double *bresid ){

int i,j,k,J=*gno,K=*nc,row=*nr;

for(j=0;j<J;j++){
 sweight[j]=0.0;                      //initialize sweight;
 for(k=0;k<K;k++) bmean[j+k*J]=0.0;  //initialize bmean;
}

for(i=0;i<row;i++) for(k=0;k<K;k++) bresid[i+k*row]=0.0;   //initialize bresid;

for(i=0;i<row;i++) sweight[gid[i]-1]+=w[i];   //sum weight over j;

for(i=0;i<row;i++)
 for(k=0;k<K;k++) bmean[gid[i]-1+k*J]+=w[i]*b[i+k*row];  //sum weighted b over j;

for(j=0;j<J;j++)
 for(k=0;k<K;k++) bmean[j+k*J]/=sweight[j];             //calculate weighted mean;

for(i=0;i<row;i++)
 for(k=0;k<K;k++) bresid[i+k*row]=sqrt(w[i])*(b[i+k*row]-bmean[gid[i]-1+k*J]);  //bresid;

}
//sum the weight for each gene;

SEXP sumweight(SEXP w,SEXP id, SEXP gno ){

  int i,J=INTEGER(gno)[0];
  SEXP sw;

  PROTECT(sw=allocVector(REALSXP,J));
  for(i=0;i<J;i++) REAL(sw)[i]=0.0;

  for(i=0;i<length(w);i++) REAL(sw)[INTEGER(id)[i]-1]+=REAL(w)[i];

  UNPROTECT(1);

 return(sw);
}

//Asymtotic Robust variance according to Huber's formular;
SEXP Robustvar(SEXP r,SEXP p,SEXP d,SEXP sigma,SEXP rn ){
 int i,I=length(r);
 double K=0.0,m=0.0,var=0.0,psi2=0.0,tmpk=0.0,s=REAL(sigma)[0];
 SEXP ss;

 PROTECT(ss=allocVector(REALSXP,1));
 REAL(ss)[0]=0.0;

 if(INTEGER(rn)[0]==1) {
    for(i=0;i<I;i++){                     //Huber's
       m+=fabs(REAL(r)[i]/s)>REAL(d)[0]?0.0:1.0;
       psi2+=fabs((REAL(r)[i]/s)>REAL(d)[0])?(REAL(d)[0]*REAL(d)[0]):(REAL(r)[i]*REAL(r)[i]/(s*s));//sum psi square
     }
     m/=length(r);                                   //E(psi');
     tmpk=psi2/((length(r)-INTEGER(p)[0])*m*m);
    for(i=0;i<I;i++) var+=fabs((REAL(r)[i]/s)>REAL(d)[0])?(m*m):(pow(1.0-m,2.0));
  }
 else if(INTEGER(rn)[0]==2){                 //Tukey's
   for(i=0;i<I;i++){
      m+=fabs(REAL(r)[i]/s)>REAL(d)[0]?0.0:(1.0-6.0*pow(REAL(r)[i]/(s*INTEGER(d)[0]),2.0)
                                      +5.0*pow(REAL(r)[i]/(s*INTEGER(d)[0]),4.0));
      psi2+=fabs((REAL(r)[i]/s)>REAL(d)[0])?
         0.0:(pow(REAL(r)[i]/s*pow(1.0-pow(REAL(r)[i]/(s*INTEGER(d)[0]),2.0),2.0),2.0)); //sum psi square
    }
     m/=length(r);                                   //E(psi');
     tmpk=psi2/((length(r)-INTEGER(p)[0])*m*m);

    for(i=0;i<I;i++) var+=fabs((REAL(r)[i]/s)>REAL(d)[0])?(m*m):(pow(1.0-6.0*pow(REAL(r)[i]/(s*INTEGER(d)[0]),2.0)
                                      +5.0*pow(REAL(r)[i]/(s*INTEGER(d)[0]),4.0)-m,2.0));

 }

 K=1.0+(INTEGER(p)[0]/length(r))*(var/(length(r)*m*m));               //K
 REAL(ss)[0]=K*K*tmpk;

 UNPROTECT(1);

 return(ss);
}

//Robust objective function;
SEXP Robustobj(SEXP r,SEXP sigma, SEXP cc,SEXP rn ){

  int i,I=length(r);
  SEXP ss;

  PROTECT(ss=allocVector(REALSXP,1));
  REAL(ss)[0]=0.0;

  if(INTEGER(rn)[0]==1){
    for(i=0;i<I;i++) REAL(ss)[0]+=fabs(REAL(r)[i])>REAL(cc)[0]*REAL(sigma)[0]?
                    (REAL(cc)[0]*fabs(REAL(r)[i]/REAL(sigma)[0])-pow(REAL(cc)[0],2.0)/2.0)
                   :(0.5*pow(REAL(r)[i]/REAL(sigma)[0],2.0));
   }
  else if(INTEGER(rn)[0]==2){
    for(i=0;i<I;i++) REAL(ss)[0]+=fabs(REAL(r)[i])>REAL(cc)[0]*REAL(sigma)[0]?
                 //(REAL(cc)[0]*REAL(cc)[0]/6.0):
                 1.0:(pow(1.0-pow(1.0-pow(REAL(r)[i]/(REAL(sigma)[0]*REAL(cc)[0]),2.0),2.0),3.0));
  }
  REAL(ss)[0]=sqrt(REAL(ss)[0]);

  UNPROTECT(1);

  return(ss);
}
