#ifndef CURVSTRAIGHT_H
#define CURVSTRAIGHT_H

#include "Curvebase.h"

class curvStraight : public Curvebase::Curvebase {

public:
  curvStraight(double a1, double b1, int orientation, double SecondDim);
  ~curvStraight();


private:
  double xp(double p);
  double yp(double p);
  double dxp(double p);
  double dyp(double p);
  int o_;
  double Sdim_;

};

#endif
