#include <iostream>
#include "Curvebase.h"
#include "curvStraight.h"

curvStraight::curvStraight(double a1, double b1, int orientation, double SecondDim) : Curvebase::Curvebase() {
  if ( b1 < a1 ) {
    std::cout << "\n changed a & b order" << std::endl;
    a_ = b1;
    b_ = a1;
  } else {
    a_ = a1;
    b_ = b1;
  }
  o_ = orientation;
  Sdim_ = SecondDim;
  lb = Curvebase::integrate(l,a_,b_,tol/100);  // total length of curve.
}

curvStraight::~curvStraight(){}


double curvStraight::xp(double p){
  if(o_ == 0) return p;
  else return Sdim_;
}
double curvStraight::yp(double p){
  if(o_ == 0) return Sdim_;
  else return p;
}
double curvStraight::dxp(double p){
  if(o_ == 0) return 1.0;
  else return 0.0;
}
double curvStraight::dyp(double p){
  if(o_ == 0) return 0.0;
  else return 1.0;
}

int o_;
double Sdim_;
