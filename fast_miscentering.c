#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define TOL 1e-4
#define key 4 //Rule for QAG
#define ws_size 8000

typedef struct integrand_params{
  gsl_spline*Sspl;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  double lxlo;double lxhi;
  double Rm;double Rp;
  double*R;double*Sigma;
  int NR;int N;
  double alpha,A; //Power law variables
  double *u; //Critical points for angular integral
}integrand_params;


//Sigma(R)
double get_Sigma(double u,double x,double Rm,double Rp,
		 double*R,double*Sigma,int NR,double alpha,double A,
		 gsl_spline*Sspl,gsl_interp_accel*acc){
  double prod = 2*Rm*Rm*x;
  double arg = sqrt(Rp*Rp+prod+2*Rp*sqrt(prod)*u);
  if (arg < R[0]){// Assume power law at extremes
    return A*pow(arg,alpha);
  }else if (arg > R[NR-1]){
    return 0.0;//Sigma(R) is zero or negative at large R
  }
  return gsl_spline_eval(Sspl,arg,acc);
}

double calc_Sigma_ang(double x,double Rm,double Rp,
		      double*R,double*Sigma,int NR,int N,double*u,
		      double alpha, double A,
		      gsl_spline*Sspl,gsl_interp_accel*acc){
  double sum = 0;
  int i;
  for(i=0;i<N;i++)
    sum += get_Sigma(u[i],x,Rm,Rp,R,Sigma,NR,alpha,A,Sspl,acc);
  return sum;
}

double integrand(double lx,void*pars){
  double x = exp(lx);
  integrand_params*params = (integrand_params*)pars;  
  double Rm=params->Rm,Rp=params->Rp;
  double*R=params->R;
  double*Sigma=params->Sigma;
  gsl_spline*Sspl=params->Sspl;
  gsl_interp_accel*acc=params->acc;
  int NR=params->NR,N=params->N;
  double alpha=params->alpha,A=params->A;
  double*u=params->u;
  double ret = x*exp(-x)*calc_Sigma_ang(x,Rm,Rp,R,Sigma,NR,N,u,alpha,A,Sspl,acc);
  return ret;
}

double calc_Sigma_misc_at_R(double Rp,integrand_params*params){
  params->Rp=Rp;
  //double lxlo=params->lxlo;
  //double lxhi=params->lxhi;
  gsl_integration_workspace*ws=params->workspace;
  gsl_function F;
  F.function=&integrand;
  F.params=params;
  double result,abserr;
  int status=0;
  status|=gsl_integration_qag(&F,-10,+20,TOL,TOL/10.,ws_size,key,ws,&result,&abserr);
  //status=gsl_integration_qag(&F,lxlo-3,lxhi,TOL,TOL/10.,ws_size,key,ws,&result,&abserr);
  //printf("res=%e\terr=%e\n",result,abserr);
  if (status)
    return -1e9;
  return result;
}

int calc_Sigma_misc(double Rm,double*R,double*Sigma,
		      double*Sigma_misc,int NR,int N){
  int i;

  gsl_spline*Sspl=gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(Sspl,R,Sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace
    =gsl_integration_workspace_alloc(ws_size);
  double lxlo=log(R[0]*R[0]/(2.*Rm*Rm));
  double lxhi=log(R[NR-1]*R[NR-1]/(2.*Rm*Rm));

  integrand_params*params=malloc(sizeof(integrand_params));
  params->Sspl=Sspl;
  params->acc=acc;
  params->workspace=workspace;
  params->lxlo=lxlo,params->lxhi=lxhi;
  params->Rm=Rm;
  params->R=R,params->Sigma=Sigma;
  params->NR=NR,params->N=N;

  //Calculate the power law at small radii
  params->alpha = log(Sigma[1]/Sigma[0])/log(R[1]/R[0]);
  params->A = Sigma[0]/pow(R[0],params->alpha);

  double*u=malloc(sizeof(double)*N);
  double w = PI/N;
  for (i=0;i<N;i++)
    u[i] = cos((i+0.5)*w);
  params->u=u;

  //#pragma omp parallel shared(Rm,R,Sigma,Sigma_misc,NR,N)
  //#pragma omp for
  for(i=0;i<NR;i++){
    //printf("%d\t%e\t",i,R[i]);
    Sigma_misc[i]=calc_Sigma_misc_at_R(R[i],params)/N;
  }
  gsl_spline_free(Sspl),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}
