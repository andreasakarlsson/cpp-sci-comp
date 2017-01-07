#ifndef GFKT_H
#define GFKT_H


#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"

class GFkt{

private:

	MatrixNEW u;
	std::shared_ptr<Domain> grid;


public:

	GFkt(std::shared_ptr<Domain> grid_);
	GFkt(const GFkt& U);

	GFkt& operator=(const GFkt& U);
	GFkt operator+(const GFkt& U) const;
	GFkt operator-(const GFkt& U) const;
	GFkt operator*(const GFkt& U) const;
	GFkt operator*(double d) const;
	void setfunction(double (*fp)(double,double));
	void save2file(const char* fname);
	GFkt Dx();
	GFkt Dy();
	GFkt Laplacian(GFkt ux_, GFkt uy_);

};

	GFkt operator*(double d, GFkt& G);

#endif