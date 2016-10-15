#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define TOL 1e-4
#define workspace_size 8000

typedef struct integrand_params{
  gsl_spline*Sspl;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  double Rl;
  double lrmin;double lrmax;
  double rmin;double rmax;
  double Rmis;double invRmis2;
  int NR;int N;
  double*R;double*Sigma;
}integrand_params;

//Sigma(R)
double get_Sigma(double u,double Rl,double Rp,double*R,
		 double*Sigma,int NR,gsl_spline*Sspl,gsl_interp_accel*acc){
  double arg = sqrt(Rl*Rl+Rp*Rp+2*Rl*Rp*u);
  double Rmin = R[0];
  double Rmax = R[NR-1];
  double alpha,A;
  if (arg < Rmin){
    alpha = log(Sigma[1]/Sigma[0])/log(R[1]/R[0]);
    A = Sigma[0]/pow(R[0],alpha);
    return A*pow(arg,alpha);
  }else if (arg > Rmax){
    alpha = log(Sigma[NR-1]/Sigma[NR-2])/log(R[NR-1]/R[NR-2]);
    A = Sigma[NR-1]/pow(R[NR-1],alpha);
    return A*pow(arg,alpha);
  }// Assume power laws at extremes
  return gsl_spline_eval(Sspl,arg,acc);
}

double calc_Sigma_angular_at_R(double Rp,double Rl,double*R,double*Sigma,
			       int NR,int N,gsl_spline*Sspl,gsl_interp_accel*acc){
  double x,f;
  double w = PI/N;

  double sum = 0;
  int i;
  for(i=1;i<=N;i++){
    x = cos((i-0.5)*w);
    f = get_Sigma(x,Rl,Rp,R,Sigma,NR,Sspl,acc);
    sum += f;
  }
  return w*sum/PI; // properly normalized
}

int calc_Sigma_angular(double Rp,double*R,double*Sigma,double*Sigma_angular,int NR,int N){
  int i;
  gsl_spline*Sspl = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(Sspl,R,Sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

#pragma omp parallel shared(Rp,R,Sigma,NR,N)
#pragma omp for
    for(i=0;i<NR;i++)
      Sigma_angular[i] = calc_Sigma_angular_at_R(Rp,R[i],R,Sigma,NR,N,Sspl,acc);
  
    return 0;
}


double P_mis(double Rp,double invRmis2){
  //invRmis2 = 1./Rmis/Rmis/2.;
  return Rp*2*invRmis2*exp(-Rp*Rp*invRmis2);
}

double integrand(double lRp,void*params){
  double Rp = exp(lRp);
  integrand_params*pars = (integrand_params*)params;
  double Rl = pars->Rl;
  double invRmis2 = pars->invRmis2;
  double*R = pars->R;
  double*Sigma = pars->Sigma;
  gsl_spline*Sspl = pars->Sspl;
  gsl_interp_accel*acc = pars->acc;
  int NR = pars->NR;
  int N = pars->N;

  double SA = calc_Sigma_angular_at_R(Rp,Rl,R,Sigma,NR,N,Sspl,acc);
  return Rp*Rp*2*invRmis2*exp(-Rp*Rp*invRmis2)*SA;
}

double calc_Sigma_mis_at_R(double Rl,integrand_params*params){
  params->Rl = Rl;
  double lrmin = params->lrmin,lrmax = params->lrmax;
  gsl_integration_workspace*workspace=params->workspace;
  gsl_function F;
  F.function = &integrand;
  F.params = params;

  double result,abserr;
  int status = 0;
  status = gsl_integration_qag(&F,lrmin,lrmax,TOL,TOL/10.,workspace_size,6,workspace,&result,&abserr);
  return result;
}

double calc_Sigma_miscentered(double Rmis,double*R,double*Sigma,
			      double*Sigma_miscentered,int NR,int N){
  gsl_spline*Sspl = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(Sspl,R,Sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->Sspl=Sspl;
  params->workspace=workspace;
  params->lrmin=log(R[0]),params->lrmax=log(R[NR-1]);
  params->rmin=R[0],params->rmax=R[NR-1];
  params->Rmis=Rmis;
  params->invRmis2=1./(Rmis*Rmis*2.);
  params->NR=NR,params->N=N;
  params->R=R,params->Sigma=Sigma;

  int i;
  for (i = 0; i < NR; i++){
    params->Rl = R[i];
    Sigma_miscentered[i] = calc_Sigma_mis_at_R(R[i],params);
    printf("%e %e\n",R[i],Sigma_miscentered[i]);
  }

  gsl_spline_free(Sspl),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);

  return 0;
}
