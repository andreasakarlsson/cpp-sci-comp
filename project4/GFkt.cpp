#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>		// fopen etc
#include <stdlib.h>		// abort()
#include <memory>		// for shared_ptr
#include <exception>

#include "GFkt.h"
#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"


using namespace std;


GFkt::GFkt(shared_ptr<Domain> grid_) : u(grid_->xsize()+1,grid_->ysize()+1,0.0), grid(grid_) {}

// copy-constructor
GFkt::GFkt(const GFkt& U) : u(U.u), grid(U.grid) {}

// copy-assignment constructor
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
	else {
		cout << "error" << endl; // give error somehow
		throw std::exception();
	}
}

// Subtraction
GFkt GFkt::operator-(const GFkt& U) const {
	if(grid == U.grid) { // must be same grid to allow subtraction
		GFkt tmp(grid);
		tmp.u = u - U.u; 
		return tmp;
	}
	else {
		cout << "error: not same grid" << endl; // give error somehow
		throw std::exception();
	}
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
	else {
		cout << "error: not same grid" << endl; // give error somehow
		throw std::exception();
	}
}

// Multiplication by scalar from right
GFkt GFkt::operator*(double d) const {
	GFkt tmp(grid);
	for (int j = 0; j <= grid->ysize(); j++)
		for (int i = 0; i <= grid->xsize(); i++)
			tmp.u(i,j) = u(i,j)*d;
	return tmp;
}

// Multiplication by scalar from left
GFkt operator*(double d, GFkt& G){
	return G*d;
}

// Sets the MatrixNEW u with a function fp
void GFkt::setfunction(double (*fp)(double,double)){
	int n_ = grid->xsize(); int m_ = grid->ysize(); 
	double* x_ = grid->xvector();	// get real x & y values
	double* y_ = grid->yvector();

	for (int j = 0; j <= m_; j++)
		for (int i = 0; i <= n_; i++)
			u(i,j) = fp(x_[j+i*(m_+1)], y_[j+i*(m_+1)]);
}

// Saves MatrixNew u (i.e. function values) to a binary file.
void GFkt::save2file(const char* fname){
		FILE *fil;
		fil = fopen(fname,"wb");
		fwrite(u.getMatrix(),sizeof(double),(grid->xsize()+1)*(grid->ysize()+1),fil);
		fclose(fil);
}

// D0x Partial derivative
GFkt GFkt::Dx(){
	int n_ = grid->xsize(); int m_ = grid->ysize(); 
	MatrixNEW ux(n_+1,m_+1,0.0);
	double* x_ = grid->xvector();
	for(int i=1; i<=n_-1; i++){   // from i=1 to i=n_-1 i.e. borders = 0.0
		for(int j=0; j<=m_; j++){ // from j=0 to i=m_ i.e. also borders
			ux(i,j) = (u(i+1,j) - u(i-1,j))/(x_[j+(i+1)*(m_+1)]-x_[j+(i-1)*(m_+1)]); 
		}
	}

	// One-sided expression for borders
	for(int j=0; j<=m_; j++){  
		double h1 = x_[j+(1)*(m_+1)] - x_[j+(0)*(m_+1)];
		double h2 = x_[j+(2)*(m_+1)] - x_[j+(0)*(m_+1)];
		ux(0,j) = (-u(0,j)*(h2*h2-h1*h1) + h2*h2*u(1,j) - h1*h1*u(2,j))/(h1*h2*(h2-h1));

		double hn1 = x_[j+(n_)*(m_+1)] - x_[j+(n_-1)*(m_+1)];
		double hn2 = x_[j+(n_)*(m_+1)] - x_[j+(n_-2)*(m_+1)];
		ux(n_,j) = (-u(n_,j)*(hn2*hn2-hn1*hn1) + hn2*hn2*u(n_-1,j) - hn1*hn1*u(n_-2,j))/(hn1*hn2*(hn1-hn2));
	}

	GFkt tmp(grid);
	tmp.u = ux;
	return tmp;
}

// D0y Partial derivative
GFkt GFkt::Dy(){
	int n_ = grid->xsize(); int m_ = grid->ysize();
	MatrixNEW uy(n_+1,m_+1,0.0);
	double* y_ = grid->yvector();
	for(int i=0; i<=n_; i++){   // from i=0 to i=n_ i.e. also borders
		for(int j=1; j<=m_-1; j++){ // from i=1 to i=m_-1 i.e. borders = 0.0
			uy(i,j) = (u(i,j+1) - u(i,j-1))/(y_[(j+1)+i*(m_+1)]-y_[(j-1)+i*(m_+1)]);
		}
	}

	// One-sided expression for borders
	for(int i=0; i<=n_; i++){  
	double h1 = y_[1+i*(m_+1)] - y_[0+i*(m_+1)];
	double h2 = y_[2+i*(m_+1)] - y_[0+i*(m_+1)];
	uy(i,0) = (-u(i,0)*(h2*h2-h1*h1) + h2*h2*u(i,1) - h1*h1*u(i,2))/(h1*h2*(h2-h1));

	double hn1 = y_[m_+i*(m_+1)] - y_[(m_-1)+i*(m_+1)];
	double hn2 = y_[m_+i*(m_+1)] - y_[(m_-2)+i*(m_+1)];
	uy(i,m_) = (-u(i,m_)*(hn2*hn2-hn1*hn1) + hn2*hn2*u(i,m_-1) - hn1*hn1*u(i,m_-2))/(hn1*hn2*(hn1-hn2));
	}


	GFkt tmp(grid);
	tmp.u = uy;
	return tmp;
}
    
GFkt GFkt::Laplacian(GFkt ux_, GFkt uy_){
	MatrixNEW ux = ux_.u;
	MatrixNEW uy = uy_.u;

	int n_ = grid->xsize(); int m_ = grid->ysize();
	MatrixNEW uxx(n_+1,m_+1,0.0);  
	MatrixNEW uyy(n_+1,m_+1,0.0);
	double* x_ = grid->xvector();
	double* y_ = grid->yvector();   

	for(int i=1; i<=n_-1; i++){   // from i=1 to i=n_-1 i.e. borders = 0.0
		for(int j=0; j<=m_; j++){ // no borders.
			uxx(i,j) = (ux(i+1,j)-ux(i-1,j))/(x_[j+(i+1)*(m_+1)]-x_[j+(i-1)*(m_+1)]);
		}
	}
	for(int i=0; i<=n_; i++){   // from i=0 to i=n_ i.e. also borders
		for(int j=1; j<=m_-1; j++){ // from j=1 to j=m_-1 i.e. border = 0.0
			uyy(i,j) = (uy(i,j+1)-uy(i,j-1))/(y_[(j+1)+i*(m_+1)]-y_[(j-1)+i*(m_+1)]);
		}
	}


	// One-sided expressions for borders (4-stencil)

		double u0;
		double u1;
		double u2;
		double u3;
		double h1;
		double h2;
		double h3;


	for(int j=0; j<=m_; j++){   // one sided expression supposed to fix first and last column (for uxx).

		int i = 0;
		u0 = u(i+0,j);
		u1 = u(i+1,j);
		u2 = u(i+2,j);
		u3 = u(i+3,j);
		h1 = x_[j+(1)*(m_+1)] - x_[j+(0)*(m_+1)];
		h2 = x_[j+(2)*(m_+1)] - x_[j+(0)*(m_+1)];
		h3 = x_[j+(3)*(m_+1)] - x_[j+(0)*(m_+1)];

		uxx(i,j) = 2*(-u0*((h3*h3*h3 -h1*h1*h1)*(h2*h3*h3*h3 - h3*h2*h2*h2)-(h3*h3*h3 - h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1)) 
		+ (h2*h3*h3*h3 - h3*h2*h2*h2)*(h3*h3*h3*u1 - h1*h1*h1*u3) - (h1*h3*h3*h3 - h3*h1*h1*h1)*(h3*h3*h3*u2 - h2*h2*h2*u3) )
		/((h1*h1*h3*h3*h3 - h3*h3*h1*h1*h1) * (h2*h3*h3*h3 - h3*h2*h2*h2) - (h2*h2*h3*h3*h3 - h3*h3*h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1));

		i = n_;
		u0 = u(i-0,j);
		u1 = u(i-1,j);
		u2 = u(i-2,j);
		u3 = u(i-3,j);
		h1 = x_[j+(i-1)*(m_+1)] - x_[j+(i-0)*(m_+1)];
		h2 = x_[j+(i-2)*(m_+1)] - x_[j+(i-0)*(m_+1)];
		h3 = x_[j+(i-3)*(m_+1)] - x_[j+(i-0)*(m_+1)];

		uxx(i,j) = 2*(-u0*((h3*h3*h3 -h1*h1*h1)*(h2*h3*h3*h3 - h3*h2*h2*h2)-(h3*h3*h3 - h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1)) 
		+ (h2*h3*h3*h3 - h3*h2*h2*h2)*(h3*h3*h3*u1 - h1*h1*h1*u3) - (h1*h3*h3*h3 - h3*h1*h1*h1)*(h3*h3*h3*u2 - h2*h2*h2*u3) )
		/((h1*h1*h3*h3*h3 - h3*h3*h1*h1*h1) * (h2*h3*h3*h3 - h3*h2*h2*h2) - (h2*h2*h3*h3*h3 - h3*h3*h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1));
		
	}

	for(int i=0; i<=n_; i++){   // one sided expression supposed to fix first and last row (for uyy).

		int j = 0;
		u0 = u(i,j+0);
		u1 = u(i,j+1);
		u2 = u(i,j+2);
		u3 = u(i,j+3);
		h1 = y_[1+i*(m_+1)] - y_[0+i*(m_+1)];
		h2 = y_[2+i*(m_+1)] - y_[0+i*(m_+1)];
		h3 = y_[3+i*(m_+1)] - y_[0+i*(m_+1)];

		uyy(i,j) = 2*(-u0*((h3*h3*h3 -h1*h1*h1)*(h2*h3*h3*h3 - h3*h2*h2*h2)-(h3*h3*h3 - h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1)) 
		+ (h2*h3*h3*h3 - h3*h2*h2*h2)*(h3*h3*h3*u1 - h1*h1*h1*u3) - (h1*h3*h3*h3 - h3*h1*h1*h1)*(h3*h3*h3*u2 - h2*h2*h2*u3) )
		/((h1*h1*h3*h3*h3 - h3*h3*h1*h1*h1) * (h2*h3*h3*h3 - h3*h2*h2*h2) - (h2*h2*h3*h3*h3 - h3*h3*h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1));
		
		j = m_;
		u0 = u(i,j-0);
		u1 = u(i,j-1);
		u2 = u(i,j-2);
		u3 = u(i,j-3);
		h1 = y_[(j-1)+i*(m_+1)] - y_[(j-0)+i*(m_+1)];
		h2 = y_[(j-2)+i*(m_+1)] - y_[(j-0)+i*(m_+1)];
		h3 = y_[(j-3)+i*(m_+1)] - y_[(j-0)+i*(m_+1)];

		uyy(i,j) = 2*(-u0*((h3*h3*h3 -h1*h1*h1)*(h2*h3*h3*h3 - h3*h2*h2*h2)-(h3*h3*h3 - h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1)) 
		+ (h2*h3*h3*h3 - h3*h2*h2*h2)*(h3*h3*h3*u1 - h1*h1*h1*u3) - (h1*h3*h3*h3 - h3*h1*h1*h1)*(h3*h3*h3*u2 - h2*h2*h2*u3) )
		/((h1*h1*h3*h3*h3 - h3*h3*h1*h1*h1) * (h2*h3*h3*h3 - h3*h2*h2*h2) - (h2*h2*h3*h3*h3 - h3*h3*h2*h2*h2)*(h1*h3*h3*h3 - h3*h1*h1*h1));

	}

	GFkt tmp(grid);
	tmp.u = uxx + uyy; // Laplacian 
	return tmp;
}