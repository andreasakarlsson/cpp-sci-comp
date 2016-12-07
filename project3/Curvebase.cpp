#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

#include "Curvebase.h"

double Curvebase::integrate(double (Curvebase::*funcp)(double), double a, double b, double tol1){
  double Int1 = 0.0;
  double Int2 = 0.0;
  double errest = 10;
  int n = 1;

  while (errest > 15*tol1){
    std::vector<double> lims (n+1);
    double h = (b-a)/(n);
    for(int i = 0; i<n+1; i++) {lims[i] = a+i*h;}
    for(int j = 0; j<n; j++) {Int1 += I(funcp, lims[j],lims[j+1]);}
    errest = abs(Int1-Int2);
    n *= 2;  // update number of intervals used
    Int2 = Int1;
    Int1 = 0;
  }
  return Int2;
}
double Curvebase::I(double (Curvebase::*funcp)(double),double a, double b){
  return (b-a)/6*((this->*funcp)(a)+4*(this->*funcp)((a+b)/2)+(this->*funcp)(b));
}

double Curvebase::newton(double (Curvebase::*fp1)(double,double), double (Curvebase::*dfp1)(double), double s, double initG, double tol1){
  double x1;
  double x = initG;
  int it = 0; int maxit = 1000;
  double err = tol1 + 1;
  while( err > tol1 && it < maxit){
    x1 = x-(this->*fp1)(x,s)/(this->*dfp1)(x);
    err = abs(x1-x);
    x = x1;
    it++;
  }
  if(err >= tol1) {std::cout << s  << initG << " - No convergence Newton. Res = " << x << x1 << std::endl;}
  return x;
}


double Curvebase::f(double p){
  return sqrt(dxp(p)*dxp(p)+dyp(p)*dyp(p)); // The Integrand of the arc length integral.
}
double Curvebase::fp(const double p, const double s){
  return integrate(l,a_,p,tol/100)-s*lb; // Equation to determine p from s.
}

Curvebase::Curvebase() {}

double (Curvebase::*l)(double) = &Curvebase::f; // class member function pointer.
double (Curvebase::*fpP)(double,double) = &Curvebase::fp;

double Curvebase::x(double s){
  if((s<0)||(s>1)){ std::cout << "Got invalid s" << std::endl;}
  double p = newton(fpP,l,s,0.5*(a_+b_),tol);
  return xp(p);
}
double Curvebase::y(double s){
  if((s<0)||(s>1)){ std::cout << "Got invalid s" << std::endl;}
  double p = newton(fpP,l,s,0.5*(a_+b_),tol);
  return yp(p);
}

Curvebase::~Curvebase(){}
