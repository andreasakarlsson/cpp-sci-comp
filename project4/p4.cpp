#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>		// fopen etc
#include <stdlib.h>		// abort()
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <memory>		// for shared_ptr

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


	curvStraight A(-10,5,0,0); // For initial test use lower boundry that is straight
	curvStraight B(0,3,1,5);
	curvStraight C(-10,5,0,3);
	curvStraight D(0,3,1,-10);

	shared_ptr<Domain> Grid = make_shared<Domain>(A,B,C,D);
	Grid->generate_grid(49,19,false);
	Grid->save2file("task1.bin");

	GFkt GA(Grid);


	double (*fp)(double,double) = u;  // function pointer to u()

  	cout << u(1,0) << endl;	
  	cout << u(9,0) << endl;
  	cout << u(16,0) << endl;

  	GA.setfunction(fp);

  	GA.save2file("task3-1.bin");

  	GA.Dx("task3-2.bin");
  	GA.Dy("task3-3.bin");

  	GA.Laplacian("task3-4.bin");

  	/* SHOW taks3-2.bin with MATLAB:

	fid = fopen('task1.bin','r');
	c = fread(fid,'double');
	fclose(fid);

	x = c(1:length(c)/2);
	y = c(length(c)/2+1:end);

	figure(101)
	plot(x,y,'.')
	axis([-12 7 -1 4])


	fid = fopen('task3-2.bin','r');
	u = fread(fid,'double');
	fclose(fid);
	A = vec2mat(u,50);
	figure(102)
	imagesc(A)
	figure(103)
	X = vec2mat(x,20)';
	Y = vec2mat(y,20)';
	Y(end:-1:1) = Y;
	surf(X,Y,A)
	figure(104)
	plot(X(end-1,:),A(end-1,:))

	*/

	return 0;
}