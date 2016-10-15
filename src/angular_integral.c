#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979

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

double calc_Sigma_angular_at_R(double Rp,double Rl,double*R,double*Sigma,int NR,int N){
  double x,f;
  double w = PI/N;

  gsl_spline*Sspl = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(Sspl,R,Sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    x = cos((i-0.5)*w);
    f = get_Sigma(x,Rl,Rp,R,Sigma,NR,Sspl,acc);
    sum += f;
  }
  return w*sum/PI; // properly normalized
}

int calc_Sigma_angular(double Rp,double*R,double*Sigma,double*Sigma_angular,int NR,int N){
  int i;
  
#pragma omp parallel shared(Rp,R,Sigma,NR,N)
#pragma omp for
    for(i=0;i<NR;i++)
      Sigma_angular[i] = calc_Sigma_angular_at_R(Rp,R[i],R,Sigma,NR,N);
  
    return 0;
}

double Gaussian(double Rp,double Rmis){
  double invRmis2 = 1./Rmis/Rmis;
  return Rp*invRmis2*exp(-Rp*Rp*invRmis2);
}
