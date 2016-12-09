#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

#ifdef _WIN32
#include <conio.h> // loads getch();
#endif

#include "MatrixNEW.h"



using namespace std;


	int main() {


	MatrixNEW A(3,5,1.2);

	//cout << A << endl;

	MatrixNEW B(3,5,7.5);

	//cout << B << endl;

	B = A;

	//cout << B << endl;

	MatrixNEW C(3,5,3.3);

	cout << C << endl;

	A = B + C;

	cout << A << endl;


	MatrixNEW D(3,5,2.3);


	A = A * D;

	cout << A << endl;	

	A *= D;

	cout << A << endl;	

	cout << A.getVal(2,4) << endl;

	#ifdef _WIN32
		cout << "\n Press any key to quit..." << endl;
		getch();
	#endif

	return 0;
}