#ifndef DOMAIN_H
#define DOMAIN_H

#include "Curvebase.h"

class Domain {

 public:
  // Domain();
  Domain(Curvebase&,Curvebase&,Curvebase&,Curvebase&);
  virtual ~Domain();
  bool check_consistency(double epsilon = 1e-3);
  void generate_grid(int n, int m, bool stretch = false);
  int n_;
  int m_;
  double h1;
  double h2;

  // The eight corner values (x & y). Calculated once
  // before the loop instead of for every iteration to
  // increase speed.

  double s0x0;
  double s1x0;
  double s3x1;
  double s2x1;

  double s0y0;
  double s1y0;
  double s3y1;
  double s2y1;


  // Values for the sides are calculated before the loop
  // this way we do m+1 calculations instead of
  // (m+1)*(m+1).
  double x_i_s0;
  double x_i_s2;
  double y_i_s0;
  double y_i_s2;

  double x_j_s3;
  double x_j_s1;
  double y_j_s3;
  double y_j_s1;

  void save2file(const char* fname = "outfile.bin");

  int xsize();
  int ysize();
  bool grid_valid();

  double* xvector();
  double* yvector();

  private:
  Curvebase* sides[4];

  inline double phi1(double w);
  inline double phi2(double w);
  
  double *x_,*y_;
};
#endif
