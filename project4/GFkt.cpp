#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>		// fopen etc
#include <stdlib.h>		// abort()
#include <memory>		// for shared_ptr

#include "GFkt.h"
#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"


using namespace std;


GFkt::GFkt(shared_ptr<Domain> grid_) : u(grid_->xsize()+1,grid_->ysize()+1,0.0), grid(grid_) {}

// copy constructor
GFkt::GFkt(const GFkt& U) : u(U.u), grid(U.grid) {}

// copy assignment constructor
GFkt& GFkt::operator=(const GFkt& U) {
	return *this;
}

// Addition
GFkt GFkt::operator+(const GFkt& U) const {
	if(grid == U.grid) { // must be same grid to allow addition
		GFkt tmp(grid);
		tmp.u = u + U.u; 
		return tmp;
	}
	else cout << "error" << endl; // give error somehow
	abort();
}

// Subtraction
GFkt GFkt::operator-(const GFkt& U) const {
	if(grid == U.grid) { // must be same grid to allow subtraction
		GFkt tmp(grid);
		tmp.u = u - U.u; 
		return tmp;
	}
	else cout << "error" << endl; // give error somehow
	abort();
}

// Pointwise multiplication
GFkt GFkt::operator*(const GFkt& U) const {
	if(grid == U.grid) { // must be same grid to allow pointwise multiplication
		GFkt tmp(grid);
		for (int j = 0; j <= grid->ysize(); j++)
			for (int i = 0; i <= grid->xsize(); i++)
				tmp.u(i,j) = u(i,j)*U.u(i,j);
		return tmp;
	}
	else cout << "error" << endl; // give error somehow
	abort();
}

// Multiplication by scalar
GFkt GFkt::operator*(double d) const {
	GFkt tmp(grid);
	for (int j = 0; j <= grid->ysize(); j++)
		for (int i = 0; i <= grid->xsize(); i++)
			tmp.u(i,j) = u(i,j)*d;
	return tmp;
}

// Sets the MatrixNEW with a function fp
void GFkt::setfunction(double (*fp)(double,double)){
	double* x_ = grid->xvector();	// get real x & y values
	double* y_ = grid->yvector();

	for (int j = 0; j <= grid->ysize(); j++)
		for (int i = 0; i <= grid->xsize(); i++)
			u(i,j) = fp(x_[j+i*(grid->ysize()+1)], y_[j+i*(grid->ysize()+1)]);  // TODO: set real x and y values from grid
}

// Saves MatrixNEW to binary file
void GFkt::save2file(const char* fname){
		FILE *fil;
		fil = fopen(fname,"wb");
		fwrite(u.getMatrix(),sizeof(double),(grid->xsize()+1)*(grid->ysize()+1),fil);
		fclose(fil);
}

// D0x Partial derivative   (Border not fixed)
void GFkt::Dx(const char* fname){
	int n_ = grid->xsize(); int m_ = grid->ysize(); 
	MatrixNEW Dx(n_+1,m_+1,0.0);  // size 50 x 20
	double* x_ = grid->xvector();
	for(int i=1; i<=n_-1; i++){   // from i=1 to i=48 i.e. borders = 0.0
		for(int j=0; j<=m_; j++){ // from j=0 to i=49 i.e. also borders
			Dx(i,j) = (u(i+1,j) - u(i-1,j))/(x_[j+(i+1)*(m_+1)]-x_[j+(i-1)*(m_+1)]); // correct??
		}
	}
	for(int j=0; j<=m_; j++){   // fix borders (first and last column)
		// one sides expressions
		Dx(0,j) =  (u(1,j) - u(0,j))/(x_[j+(1)*(m_+1)]-x_[j+(0)*(m_+1)]); // D+
		Dx(n_,j) = (u(n_,j) - u(n_-1,j))/(x_[j+(n_)*(m_+1)]-x_[j+(n_-1)*(m_+1)]); // D-
	}
	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(Dx.getMatrix(),sizeof(double),(n_+1)*(m_+1),fil);
	fclose(fil);
}

// D0y Partial derivative   (Border not fixed)
void GFkt::Dy(const char* fname){
	int n_ = grid->xsize(); int m_ = grid->ysize();
	MatrixNEW Dy(n_+1,m_+1,0.0);  // size 50 x 20
	double* y_ = grid->yvector();
	for(int i=0; i<=n_; i++){   // from i=0 to i=49 i.e. also borders
		for(int j=1; j<=m_-1; j++){ // from i=1 to i=48 i.e. borders = 0.0
			Dy(i,j) = (u(i,j+1) - u(i,j-1))/(y_[(j+1)+i*(m_+1)]-y_[(j-1)+i*(m_+1)]); // correct??
		}
	}
	for(int i=0; i<=n_; i++){   // fix borders (first and last row)
		// one sides expressions
		Dy(i,0) =  (u(i,1) - u(i,0))/(y_[1+(i)*(m_+1)]-y_[0+(i)*(m_+1)]); // D+
		Dy(i,m_) = (u(i,m_) - u(i,m_-1))/(y_[m_+(i)*(m_+1)]-y_[(m_-1)+(i)*(m_+1)]); // D-
	}
	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(Dy.getMatrix(),sizeof(double),(n_+1)*(m_+1),fil);
	fclose(fil);
}

// D0y Partial derivative   (Border not fixed)
void GFkt::Laplacian(const char* fname){
	MatrixNEW L(grid->xsize()+1,grid->ysize()+1,0.0);  // size 50 x 20
	//double* x_ = grid->xvector();   // not used now, but should be used in real solution
	//double* y_ = grid->yvector();   
	for(int i=1; i<=grid->xsize()-1; i++){   // from i=1 to i=48 i.e. borders = 0.0
		for(int j=1; j<=grid->ysize()-1; j++){
			L(i,j) = (-4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)); // just to try.
		}
	}
	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(L.getMatrix(),sizeof(double),(grid->xsize()+1)*(grid->ysize()+1),fil);
	fclose(fil);
}