#include <cmath>
#include <iostream>
#include <cstdio>
#include <memory>   // for shared_ptr
#include "Curvebase.h"
#include "Domain.h"
using namespace std;


Domain::Domain (shared_ptr<Curvebase> s1, shared_ptr<Curvebase> s2, shared_ptr<Curvebase> s3, shared_ptr<Curvebase> s4){

  sides[0] = s1;
  sides[1] = s2;
  sides[2] = s3;
  sides[3] = s4;

  if(!check_consistency(1e-5)){
    sides[0] = sides[1] = sides[2] = sides[3] = nullptr;
    std::cout << " Not consistent. All curves set to zero.\n" << std::endl;
  }
}


bool Domain::check_consistency(double epsilon){


  // Check that the end of a curve is connected to the start of the next.
  if (abs(sides[0]->x(1) - sides[1]->x(0))>epsilon || abs(sides[0]->y(1) - sides[1]->y(0))>epsilon ){
    std::cout << "\n Curve 0 and 1 not connected \n" << std::endl;
    return false;
  }
  if ( abs(sides[1]->x(1) - sides[2]->x(1))>epsilon || abs(sides[1]->y(1) - sides[2]->y(1))>epsilon ){
    std::cout << "\n Curve 1 and 2  not connected \n" << std::endl;
    return false;
  }
  if ( abs(sides[2]->x(0) - sides[3]->x(1))>epsilon || abs(sides[2]->y(0) - sides[3]->y(1))>epsilon ){
    std::cout << "\n Curve 2 and 3 not connected \n" << std::endl;
    return false;
  }
  if ( abs(sides[3]->x(0) - sides[0]->x(0))>epsilon || abs(sides[3]->y(0) - sides[0]->y(0))>epsilon ){
    std::cout << "\n Curve 3 and 0  not connected \n" << std::endl;
    return false;
  }

  return true;
}

void Domain::generate_grid(int n, int m, bool stretch){
  if(sides[0]==nullptr){
    std::cout << " Grid point generation aborted due to nullptr sides.\n" << std::endl;
    return;
  }
  if((n<1)||(m<1)) {m=20;n=20; std::cout << "set m and n to 20" << std::endl;}
  // if(n_ != 0) { // if n_0 non-zero => grid allready exist. Must be deleted.
  //   delete [] x_;
  //   delete [] y_;
  // }
  n_ = n;
  m_ = m;
  x_ = new double[(m_+1)*(n_+1)];
  y_ = new double[(m_+1)*(n_+1)];
  double h1 = 1.0/n_;
  double h2 = 1.0/m_;


  // The eight corner values (x & y). Calculated once
  // before the loop instead of for every iteration to
  // increase speed.
  double s0x0 = sides[0]->x(0);
  double s1x0 = sides[1]->x(0);
  double s3x1 = sides[3]->x(1);
  double s2x1 = sides[2]->x(1);

  double s0y0 = sides[0]->y(0);
  double s1y0 = sides[1]->y(0);
  double s3y1 = sides[3]->y(1);
  double s2y1 = sides[2]->y(1);


  // Values for the sides are calculated before the loop
  // this way we do m+1 calculations instead of
  // (m+1)*(m+1).
  double x_i_s0[(n_+1)];
  double x_i_s2[(n_+1)];
  double y_i_s0[(n_+1)];
  double y_i_s2[(n_+1)];

  for(int i = 0; i<= n_; i++){
    x_i_s0[i] = sides[0]->x(i*h1);
    x_i_s2[i] = sides[2]->x(i*h1);
    y_i_s0[i] = sides[0]->y(i*h1);
    y_i_s2[i] = sides[2]->y(i*h1);
  }

  double x_j_s3[(m_+1)];
  double x_j_s1[(m_+1)];
  double y_j_s3[(m_+1)];
  double y_j_s1[(m_+1)];

  for(int j = 0; j<= m_; j++){
    double ss; // ss is a stretched version of j*h2
    if (stretch) ss = (exp(1.5*j*h2)-1)/(exp(1.5)-1); else ss = j*h2;
    x_j_s3[j] = sides[3]->x(ss);
    x_j_s1[j] = sides[1]->x(ss);
    y_j_s3[j] = sides[3]->y(ss);
    y_j_s1[j] = sides[1]->y(ss);
  }


  for(int i = 0; i <= n_; i++){
    for(int j = 0; j <= m_; j++){
      x_[j+i*(m_+1)] = phi1(i*h1)*x_j_s3[j]
  	+ phi2(i*h1)*x_j_s1[j]
  	+ phi1(j*h2)*x_i_s0[i]
  	+ phi2(j*h2)*x_i_s2[i]
  	- phi1(i*h1)*phi1(j*h2)*s0x0
  	- phi2(i*h1)*phi1(j*h2)*s1x0
  	- phi1(i*h1)*phi2(j*h2)*s3x1
  	- phi2(i*h1)*phi2(j*h2)*s2x1;


      y_[j+i*(m_+1)] = phi1(i*h1)*y_j_s3[j]
  	+ phi2(i*h1)*y_j_s1[j]
  	+ phi1(j*h2)*y_i_s0[i]
  	+ phi2(j*h2)*y_i_s2[i]
  	- phi1(i*h1)*phi1(j*h2)*s0y0
  	- phi2(i*h1)*phi1(j*h2)*s1y0
  	- phi1(i*h1)*phi2(j*h2)*s3y1
  	- phi2(i*h1)*phi2(j*h2)*s2y1;

    }
  }
}

// Saves the x & y-values as a single array in a binary file.
void Domain::save2file(const char* fname){
  // Create new array that contains x_ and y_
  int sizeV = (m_+1)*(n_+1);
  double * result = new double[sizeV + sizeV];
  std::copy(x_, x_ + sizeV, result);
  std::copy(y_, y_ + sizeV, result + sizeV);

  FILE *fil;
  fil = fopen(fname,"wb");
  fwrite(result,sizeof(double),2*sizeV,fil);
  fclose(fil);
}

int Domain::xsize(){ return n_; }
int Domain::ysize(){ return m_; }
bool Domain::grid_valid() { return m_ != 0; }

double* Domain::xvector(){ return x_; }
double* Domain::yvector(){ return y_; }


inline double Domain::phi1(double w){ return 1.0 - w; }
inline double Domain::phi2(double w){ return w; }

Domain::~Domain(){
  if (n_ > 0){
    delete [] x_;
    delete [] y_;
  }
}

int m_ = 0; // initialized to zero to enable class to check if
int n_ = 0; // a grid generation been done.