#include <iostream>
#include <cmath>
#include <cstdio>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "curvStraight.h"
#include "curvExp.h"
#include "Domain.h"

#ifdef _WIN32
#include <conio.h> // loads getch();
#endif


using namespace std;

int main() {

	curvExp E(-10,5); // lower boundary

	curvStraight A(-10,5,0,0); // For initial test use lower boundry that is straight
	curvStraight B(0,3,1,5);
	curvStraight C(-10,5,0,3);
	curvStraight D(0,3,1,-10);

	// PRINT OUT TO SHOW THAT THE CURVE OBJECTS GIVE CORRECT VALUES

	cout << "\n Curve E:" << endl;

	cout << " (" << E.x(0) << ", " << E.y(0) << ")" << endl;
	cout << " (" << E.x(0.25) << ", " << E.y(0.25) << ")" << endl;
	cout << " (" << E.x(0.5) << ", " << E.y(0.5) << ")" << endl;
	cout << " (" << E.x(0.75) << ", " << E.y(0.75) << ")" << endl;
	cout << " (" << E.x(1) << ", " << E.y(1) << ")" << endl;

	cout << "\n Curve A:" << endl;

	cout << " (" << A.x(0) << ", " << A.y(0) << ")" << endl;
	cout << " (" << A.x(0.25) << ", " << A.y(0.25) << ")" << endl;
	cout << " (" << A.x(0.5) << ", " << A.y(0.5) << ")" << endl;
	cout << " (" << A.x(0.75) << ", " << A.y(0.75) << ")" << endl;
	cout << " (" << A.x(1) << ", " << A.y(1) << ")" << endl;


	cout << "\n Curve B:" << endl;

	cout << " (" << B.x(0) << ", " << B.y(0) << ")" << endl;
	cout << " (" << B.x(0.25) << ", " << B.y(0.25) << ")" << endl;
	cout << " (" << B.x(0.5) << ", " << B.y(0.5) << ")" << endl;
	cout << " (" << B.x(0.75) << ", " << B.y(0.75) << ")" << endl;
	cout << " (" << B.x(1) << ", " << B.y(1) << ")" << endl;


	cout << "\n Curve C:" << endl;

	cout << " (" << C.x(0) << ", " << C.y(0) << ")" << endl;
	cout << " (" << C.x(0.25) << ", " << C.y(0.25) << ")" << endl;
	cout << " (" << C.x(0.5) << ", " << C.y(0.5) << ")" << endl;
	cout << " (" << C.x(0.75) << ", " << C.y(0.75) << ")" << endl;
	cout << " (" << C.x(1) << ", " << C.y(1) << ")" << endl;


	cout << "\n Durve D:" << endl;

	cout << " (" << D.x(0) << ", " << D.y(0) << ")" << endl;
	cout << " (" << D.x(0.25) << ", " << D.y(0.25) << ")" << endl;
	cout << " (" << D.x(0.5) << ", " << D.y(0.5) << ")" << endl;
	cout << " (" << D.x(0.75) << ", " << D.y(0.75) << ")" << endl;
	cout << " (" << D.x(1) << ", " << D.y(1) << ")" << "\n" << endl;

	// Create and save a grid (formed by the four straight
	// curves.)
	Domain Grid(A,B,C,D);
	Grid.generate_grid(49,19,false);
	Grid.save2file("task3.bin");

	// Create and save grid with the given exponential
	// function. The computational time for generation of points
	// measured.
	Domain Grid2(E,B,C,D);

	clock_t t;
	t = clock(); // start timing

	Grid2.generate_grid(49,19,false);

	t = clock() - t;

  	printf (" It took %d clicks (%f seconds).\n",(int) t,((float)t)/CLOCKS_PER_SEC);

	Grid2.save2file("task4.bin");

	// Create stretched grid according to task 5
	Grid2.generate_grid(49,19,true);
	Grid2.save2file("task5.bin");


#ifdef _WIN32
	cout << "\n Press any key to quit..." << endl;
	getch();
#endif
	return 0;
}
