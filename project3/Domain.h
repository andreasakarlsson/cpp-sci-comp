#ifndef DOMAIN_H
#define DOMAIN_H


#include "Curvebase.h"

class Domain {

 public:
  // Domain();
  Domain(std::shared_ptr<Curvebase> s1, std::shared_ptr<Curvebase> s2, std::shared_ptr<Curvebase> s3, std::shared_ptr<Curvebase> s4);
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

  // Get size of grid
  inline int xsize(){ return n_; }
  inline int ysize(){ return m_; }

  inline bool grid_valid() { return m_ != 0; }

  // access to x values
  inline double x(int i){
      if (i < 0 || i > n_) {
          std::cout << "Error: Gird index out of bounds" << std::endl; // print error message
          throw std::exception();
      }
      return x_[0+i*(m_+1)];
  }

  // access to y values
  inline double y(int j){
      if (j < 0 || j > m_) {
        std::cout << "Error: Gird index out of bounds" << std::endl; // print error message
        throw std::exception();
      }
      return y_[j+0*(m_+1)];
  }

  // access to x_ and y_ arrays
  inline double* xvector(){ return x_; }
  inline double* yvector(){ return y_; }


  private:
  
  std::shared_ptr<Curvebase> sides[4];
  //Curvebase* sides[4];

  inline double phi1(double w){ return 1.0 - w; }
  inline double phi2(double w){ return w; }
  
  double *x_,*y_;

};
#endif
