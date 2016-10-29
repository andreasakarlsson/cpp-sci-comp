#include <iostream>
#include <cmath>
#include <cstdlib>
// #include <conio.h> // loads getch();
#include <vector>
#include "matrix.h"
#include "r8lib.h"
#include "r8mat_expm1.h"
using namespace std;


Matrix matrixExp(Matrix M, double tol){


	// From tol esimate required number of terms n
	Matrix result = M.getSize();
	result.fillIdentityMatrix();
	int n = 0;
	while (true){
		n += 1;
		result.fillIdentityMatrix();	// to avoid pow
		for(int j = n; j >0; j--){
		  result = M * result / j;}
		if ( result.norm() < tol ){ break; }
	}
	cout << "number of terms = " << n << endl;


	//Calculate exponential with Horner's Scheme.
	result.fillMatrix(0);
	Matrix I = M.getSize();
	I.fillIdentityMatrix();

	for(int i = n; i >= 0; i--){
		result = (I + M * result) / ((i == 0) ? 1 : i);
	}
	return result;


}


vector<double> r8matToVec(double* aa, unsigned aaSize) {
  vector<double> vec;
  vec.insert(vec.end(), &aa[0], &aa[aaSize]);
  return vec;
}

int main() {

	Matrix M(4);
	Matrix N(M);
	cout << "matrix N (by default containing only ones): " << N;

	Matrix HE = 3; //uses copy const  works!!

	Matrix G = Matrix(4);
	G.fillMatrix();
	cout << "G (matrix filled with func fillMatrix): " << G;

	G += N;
	cout << "new G ( G += N): " << G;

	Matrix W = 4;

	W = N+N;
	cout << "addition of N and N: " << W;
	cout << "N after addition: " << N;

	// Remove argument in norm => enoguh with G.norm();
	// But is this less efficient since const is removed???
	cout << "G-norm: " << G.norm() << endl;

	G *= N;  // This gives G*N

	cout << "multiplication of new G and N ( G *= N): " << G;

	Matrix U = 4;
	U = G*N;
	cout << "multiplication of newnew G and N (U = G*N): " << U;
	cout << "N after addition: " << N;
	cout << "G after addition: " << G;

	Matrix A = G/2.5;
	cout << "G divided by 2.5: " << A;

	vector<double> v1 = {1,2,3,4};
 	Matrix P(v1);
 	cout << "Matrix from 1-4 vector:  " << P;

 	Matrix V = Matrix(4);
	V.fillMatrix(0);
	cout << "fillMatrix with argument 0 : " << V;

 	Matrix  I = Matrix(6);
	I.fillIdentityMatrix();
	cout << "fillIdentityMatrix: " << I;

	vector<double> vc = {1, 3, 10, 45, 12, 3, 5, 0, 12, 1, 3, 7, 19, 4, 9, 6};
	Matrix VC = Matrix(vc);
	cout << "This is is the VC-matrix."
	  " It will be used for the comparison of the exponential functions: "
	     << VC;

	Matrix C = matrixExp(N,0.1);
	cout << "Results from our matrixExp for matrix filled with ones: " << C;

	Matrix CC = matrixExp(VC,0.001);
	cout << "Results from our matrixExp for VC-matrix (tolerance of 0.001): "
	     << CC;

	double ones[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	double vcArray[16] = {1, 3, 10, 45, 12, 3, 5, 0, 12, 1, 3, 7, 19, 4, 9, 6};
	int ac = 4;

	// From the matlab source code:
	double* ONESE = r8mat_expm1(ac, ones);
	double* VCE = r8mat_expm1(ac, vcArray);

	// Convert matlab results into our Matrix object for comparison
	Matrix ONESEM =  Matrix(r8matToVec(ONESE, 16));
	Matrix VCEM =  Matrix(r8matToVec(VCE, 16));

	cout << "Results from r8mat_expm1 for matrix filled with ones: " << ONESEM;
	cout << "Results from r8mat_expm1 for VC-matrix:" << VCEM;

	Matrix MDIFF = CC - VCEM;
	cout << "Difference for the VC-matrix"
	  " between r8mat_expm1 and our implementation:"
	     << MDIFF;

	// getch();
	return 0;
}
