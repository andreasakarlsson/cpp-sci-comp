#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <memory>		// for shared_ptr

#ifdef _WIN32
#include <conio.h> // loads getch();
#endif

#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"

using namespace std;


class GFkt{
	
private:

	MatrixNEW u;
	shared_ptr<Domain> grid;

public:

	GFkt(shared_ptr<Domain> grid_) : u(grid_->xsize()+1,grid_->ysize()+1,0.0), grid(grid_) {}

	GFkt(const GFkt& U) : u(U.u), grid(U.grid) {}

	GFkt& operator=(const GFkt& U) {
		return *this;
	}

	GFkt operator+(const GFkt& U) const {
		if(grid == U.grid) { // must be same grid to allow addition
			GFkt tmp(grid);
			tmp.u = u + U.u; 
			return tmp;
		}
		else cout << "error" << endl; // give error somehow
	}

	GFkt operator*(const GFkt& U) const {
		if(grid == U.grid) { // must be same grid to allow pointwise multiplication
			GFkt tmp(grid);
			for (int j = 0; j <= grid->ysize(); j++)
				for (int i = 0; i <= grid->xsize(); i++)
					tmp.u(i,j) = u(i,j)*U.u(i,j);
			return tmp;
		}
		else cout << "error" << endl; // give error somehow
	}

};



int main() {


	curvStraight A(-10,5,0,0); // For initial test use lower boundry that is straight
	curvStraight B(0,3,1,5);
	curvStraight C(-10,5,0,3);
	curvStraight D(0,3,1,-10);

	shared_ptr<Domain> Grid = make_shared<Domain>(A,B,C,D);
	Grid->generate_grid(49,19,false);
	Grid->save2file("task1.bin");

	GFkt GA(Grid);




// #ifdef _WIN32
// 	cout << "\n Press any key to quit..." << endl;
// 	getch();
// #endif
	return 0;
}
