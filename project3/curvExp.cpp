#include <iostream>
#include <cmath>
#include "curvExp.h"

curvExp::curvExp(double a1, double b1) : Curvebase() {
  if ( b1 < a1 ) {
    std::cout << "\n changed a & b order" << std::endl;
    a_ = b1;
    b_ = a1;
  } else {
    a_ = a1;
    b_ = b1;
  }
  lb = Curvebase::integrate(l,a_,b_,tol/100); // total length of curve.
}

curvExp::~curvExp(){ }

double curvExp::xp(double p){
  return p;
}

double curvExp::yp(double p){
  if(p<-3)
    return 0.5/(1+exp(-3*(p+6)));
  else
    return 0.5/(1+exp(3*p));
}

double curvExp::dxp(double p){
  return 1.0;
}

double curvExp::dyp(double p){
  if(p<-3){
    return (3/2)*exp(-3*(p+6))/((1+exp(-3*(p+6)))*(1+exp(-3*(p+6))));}
  else
    return -(3/2)*exp(3*p)/((1+exp(3*p))*(1+exp(3*p)));
}
