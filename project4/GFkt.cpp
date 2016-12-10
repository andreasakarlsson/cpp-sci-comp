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

GFkt::GFkt(const GFkt& U) : u(U.u), grid(U.grid) {}

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
	for (int j = 0; j <= grid->ysize(); j++)
		for (int i = 0; i <= grid->xsize(); i++)
			u(i,j) = fp(i,0.05*j);  // TODO: set real x and y values from grid
}

// Saves MatrixNEW to binary file
void GFkt::save2file(const char* fname){
		FILE *fil;
		fil = fopen(fname,"wb");
		fwrite(u.getMatrix(),sizeof(double),(grid->xsize()+1)*(grid->ysize()+1),fil);
		fclose(fil);
}