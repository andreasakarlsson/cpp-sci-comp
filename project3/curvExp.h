#ifndef CURVEXP_H
#define CURVEXP_H

#include "Curvebase.h"

class curvExp : public Curvebase::Curvebase {

public:
  curvExp(double a1, double b1);
  ~curvExp();

private:

  double xp(double p);
  double yp(double p);
  double dxp(double p);
  double dyp(double p);

};
#endif
