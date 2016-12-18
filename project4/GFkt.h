#ifndef GFKT_H
#define GFKT_H


#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"

using namespace std;

class GFkt{

private:

	MatrixNEW u;
	MatrixNEW ux;
	MatrixNEW uy;
	MatrixNEW L;
	shared_ptr<Domain> grid;


public:

	GFkt(shared_ptr<Domain> grid_);
	GFkt(const GFkt& U);

	GFkt& operator=(const GFkt& U);
	GFkt operator+(const GFkt& U) const;
	GFkt operator-(const GFkt& U) const;
	GFkt operator*(const GFkt& U) const;
	GFkt operator*(double d) const;
	void setfunction(double (*fp)(double,double));
	void save2file(const char* fname);
	void Dx(const char* fname);
	void Dy(const char* fname);
	void Laplacian(const char* fname);

};

#endif