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


GFkt::GFkt(shared_ptr<Domain> grid_) : u(grid_->xsize()+1,grid_->ysize()+1,0.0), 
	ux(grid_->xsize()+1,grid_->ysize()+1,0.0), 
	uy(grid_->xsize()+1,grid_->ysize()+1,0.0),
	uxx(grid_->xsize()+1,grid_->ysize()+1,0.0),grid(grid_) {}

// copy constructor
GFkt::GFkt(const GFkt& U) : u(U.u), ux(U.ux), uy(U.uy), uxx(U.uxx), grid(U.grid) {}

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
	int n_ = grid->xsize(); int m_ = grid->ysize(); 
	double* x_ = grid->xvector();	// get real x & y values
	double* y_ = grid->yvector();

	for (int j = 0; j <= m_; j++)
		for (int i = 0; i <= n_; i++)
			u(i,j) = fp(x_[j+i*(m_+1)], y_[j+i*(m_+1)]);
}

// Saves MatrixNEW to binary file
void GFkt::save2file(const char* fname){
		FILE *fil;
		fil = fopen(fname,"wb");
		fwrite(u.getMatrix(),sizeof(double),(grid->xsize()+1)*(grid->ysize()+1),fil);
		fclose(fil);
}

// D0x Partial derivative
void GFkt::Dx(const char* fname){
	int n_ = grid->xsize(); int m_ = grid->ysize(); 
	//MatrixNEW Dx(n_+1,m_+1,0.0);  // size 50 x 20
	double* x_ = grid->xvector();
	for(int i=1; i<=n_-1; i++){   // from i=1 to i=48 i.e. borders = 0.0
		for(int j=0; j<=m_; j++){ // from j=0 to i=49 i.e. also borders
			ux(i,j) = (u(i+1,j) - u(i-1,j))/(x_[j+(i+1)*(m_+1)]-x_[j+(i-1)*(m_+1)]); 
		}
	}
	for(int j=0; j<=m_; j++){   // fixes borders (first and last column)
		// one sides expressions // Better expressions should be use for improved accuracy
		ux(0,j) =  (u(1,j) - u(0,j))/(x_[j+(1)*(m_+1)]-x_[j+(0)*(m_+1)]); // D+
		ux(n_,j) = (u(n_,j) - u(n_-1,j))/(x_[j+(n_)*(m_+1)]-x_[j+(n_-1)*(m_+1)]); // D-
	}
	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(ux.getMatrix(),sizeof(double),(n_+1)*(m_+1),fil);
	fclose(fil);
}

// D0y Partial derivative
void GFkt::Dy(const char* fname){
	int n_ = grid->xsize(); int m_ = grid->ysize();
	//MatrixNEW Dy(n_+1,m_+1,0.0);  // size 50 x 20
	double* y_ = grid->yvector();
	for(int i=0; i<=n_; i++){   // from i=0 to i=49 i.e. also borders
		for(int j=1; j<=m_-1; j++){ // from i=1 to i=48 i.e. borders = 0.0
			uy(i,j) = (u(i,j+1) - u(i,j-1))/(y_[(j+1)+i*(m_+1)]-y_[(j-1)+i*(m_+1)]);
		}
	}
	for(int i=0; i<=n_; i++){   // fixes borders (first and last row)
		// one sides expressions // Better expressions should be use for improved accuracy
		uy(i,0) =  (u(i,1) - u(i,0))/(y_[1+(i)*(m_+1)]-y_[0+(i)*(m_+1)]); // D+
		uy(i,m_) = (u(i,m_) - u(i,m_-1))/(y_[m_+(i)*(m_+1)]-y_[(m_-1)+(i)*(m_+1)]); // D-
	}
	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(uy.getMatrix(),sizeof(double),(n_+1)*(m_+1),fil);
	fclose(fil);
}

// Laplacian    
void GFkt::Laplacian(const char* fname){
	int n_ = grid->xsize(); int m_ = grid->ysize();
	//MatrixNEW L(n_+1,m_+1,0.0);  // size 50 x 20
	double* x_ = grid->xvector();   // not used now, but should be used in real solution
	//double* y_ = grid->yvector();   
	for(int i=2; i<=n_-1; i++){   // from i=2 to i=48 i.e. borders = 0.0
		for(int j=0; j<=m_; j++){
			//L(i,j) = (-4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)); // just to try.
			uxx(i,j) = (ux(i+1,j)-ux(i-1,j))/(x_[j+(i+1)*(m_+1)]-x_[j+(i-1)*(m_+1)]);
		}
	}

	  
	for(int j=0; j<=m_; j++){   // one sided expression supposed to fix first, second and last column.

		for(int i=0; i<=1; i++){ 
			double x0 = x_[j+(i+0)*(m_+1)];
			double x1 = x_[j+(i+1)*(m_+1)];
			double x2 = x_[j+(i+2)*(m_+1)];
			double u0 = u(i+0,j);
			double u1 = u(i+1,j);
			double u2 = u(i+2,j);

			uxx(i,j) =  2*(-u0*(x2-x1)+u1*(x2-x0)-u2*(x1-x0))/((x2-x0)*(x1-x0)*(x1-x2));
		}
		for(int i=n_; i<= n_; i++){
			double xm = x_[j+(i+0)*(m_+1)];
			double xm1 = x_[j+(i-1)*(m_+1)];
			double xm2 = x_[j+(i-2)*(m_+1)];
			double um = u(i-0,j);
			double um1 = u(i-1,j);
			double um2 = u(i-2,j);

			uxx(i,j) =  2*(-um*(xm1-xm2)+um1*(xm-xm2)-um2*(xm-xm1))/((xm-xm2)*(xm-xm1)*(xm2-xm1));
		}
	}

	FILE *fil;
	fil = fopen(fname,"wb");
	fwrite(uxx.getMatrix(),sizeof(double),(n_+1)*(m_+1),fil);
	fclose(fil);
}