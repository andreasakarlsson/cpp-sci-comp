#include <iostream>
#include <cmath>
// #include <conio.h> // loads getch();
#include <vector>
#include "matrix.h"
// #include "r8lib.h"
using namespace std;


Matrix matrixExp(Matrix M, double tol){


	// From tol esimate required number of terms n
	Matrix result = M.getSize();
	result.fillIdentityMatrix();
	int n = 0;
	while (true){
		n += 1;
		result.fillIdentityMatrix();					// to avoid pow
		for(int j = n; j >0; j--){ result = M*result/j;}  //(should be better way?)

		if ( result.norm() < tol ){ break; }
	}
	cout << " number of terms = " << n << "\n" << endl;


	//Calculate exponential with Horner's Scheme.
	result.fillMatrix(0);
	Matrix I = M.getSize();
	I.fillIdentityMatrix();

	for(int i = n; i >= 0; i--){
		result = (I + M*result) / ((i == 0) ? 1 : i);
		cout << ((i == 0) ? 1 : i) << endl;
		result.printMatrix();
	}
	return result;


}



int main() {

	Matrix M(4);
	Matrix N(M);
	cout << " matrix N (by default containing only ones): \n";
	N.printMatrix();

	Matrix HE = 3; //uses copy const  works!!

	Matrix G = Matrix(4);
	G.fillMatrix();

	cout << " G (matrix filled with func fillMatrix): \n";
	G.printMatrix();


	G += N;
	cout << " new G ( G += N): \n";
	G.printMatrix();

	Matrix W = 4;

	W = N+N;
	cout << " addition of N and N: \n";
	W.printMatrix();
	cout << " N after addition: \n";
	N.printMatrix();

	cout << "G-norm: " << G.norm() << "\n";
	// Remove argument in norm => enoguh with G.norm();
	// But is this less efficient since const is removed???

	G *= N;  // This gives G*N

	cout << " multiplication of new G and N ( G *= N): \n";
	G.printMatrix();

	Matrix U = 4;
	U = G*N;
	cout << " multiplication of newnew G and N (U = G*N): \n";
	U.printMatrix();
	cout << " N after addition: \n";
	N.printMatrix();
	cout << " G after addition: \n";
	G.printMatrix();

	cout << " G divided by 2.5: \n";
	Matrix A = G/2.5;
	A.printMatrix();

	std::vector<double> v1 = {1,2,3,4};
 	Matrix P(v1);
 	P.printMatrix();

 	cout << " fillMatrix with argument 0 : \n";
 	Matrix V = Matrix(4);
	V.fillMatrix(0);
	V.printMatrix();

	cout << " fillIdentityMatrix:  \n";
 	Matrix  I = Matrix(6);
	I.fillIdentityMatrix();
	I.printMatrix();


	Matrix C = matrixExp(N,0.1);
	cout << " matrixExp of N:  \n";
	C.printMatrix();

	cout << "\n Press any key to quit...";
	// getch();


	return 0;
}
