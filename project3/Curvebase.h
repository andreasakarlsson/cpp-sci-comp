#ifndef CURVEBASE_H
#define CURVEBASE_H

class Curvebase{
 public:
  Curvebase();
  double x(double s);
  double y(double s);
  ~Curvebase();

 protected:
  double integrate(double (Curvebase::*funcp)(double), double a, double b, double tol1);
  double I(double (Curvebase::*funcp)(double),double a, double b);
  double newton(double (Curvebase::*fp1)(double,double), double (Curvebase::*dfp1)(double), double s, double initG, double tol1);

  virtual double xp(double p) = 0;
  virtual double yp(double p) = 0;
  virtual double dxp(double p) = 0;
  virtual double dyp(double p) = 0;

  double f(double p);
  double fp(const double p, const double s);

  double a_;
  double b_;
  double lb; // Length of curve. Initialization is done in instantiated class.

  double (Curvebase::*l)(double)=&Curvebase::f; // class member function pointer.
  double (Curvebase::*fpP)(double,double)=&Curvebase::fp;
  double const tol = 1e-6; // tolerance.

};

#endif
