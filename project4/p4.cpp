#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>		// fopen etc
#include <stdlib.h>		// abort()
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <memory>		// for shared_ptr
#include <vector>

#include "GFkt.h"
#include "MatrixNEW.h"
#include "../project3/Domain.h"
#include "../project3/Curvebase.h"
#include "../project3/curvStraight.h"

using namespace std;


double u(double x, double y){
	return sin(pow(x/10.0,2))*cos(x/10.0)+y;
}

int main() {

	// create curves
	std::shared_ptr<curvStraight> A =  make_shared<curvStraight>(-10,5,0,0);
	std::shared_ptr<curvStraight> B =  make_shared<curvStraight>(0,3,1,5);
	std::shared_ptr<curvStraight> C =  make_shared<curvStraight>(-10,5,0,3);
	std::shared_ptr<curvStraight> D =  make_shared<curvStraight>(0,3,1,-10);

	// create domain
	shared_ptr<Domain> Grid = make_shared<Domain>(A,B,C,D);
	Grid->generate_grid(49,19,false);
	Grid->save2file("task1.bin");

	// create GFkt object
	GFkt GA(Grid);
	double (*fp)(double,double) = u;  // function pointer to u()
  	GA.setfunction(fp);
  	GA.save2file("task3-1.bin");


	clock_t t;
	t = clock(); // start timing

  	// create and save discrete partial derivatives
	GFkt GAux = GA.Dx();
	GAux.save2file("task3-2.bin");

	GFkt GAuy = GA.Dy();
	GAuy.save2file("task3-3.bin");

	GFkt GAL = GA.Laplacian(GAux, GAuy);
	GAL.save2file("task3-4.bin");

	t = clock() - t;
  	printf ("\n It took %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

	return 0;
}