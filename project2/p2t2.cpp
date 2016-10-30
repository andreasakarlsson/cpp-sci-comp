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
  cout << endl <<
    "======== Testing Constructors =========" << endl;

  Matrix M(4);
  cout << "Matrix from integer: " << M;
  Matrix N(M);
  cout << "Matrix copied from another matrix"
    " (by default containing only ones):" << N;


  cout << endl <<
    "======== Testing Operators & Methods =========" << endl;

  Matrix HE = 3; //uses copy const  works!!
  cout << "Matrix copied with \"=\" operator" << HE;

  Matrix G = Matrix(4);
  G.fillMatrix();
  cout << "Matrix filled using method fillMatrix: " << G;

  G += N;
  cout << "Matrix from the += operator" << G;

  Matrix W = 4;

  W = N+N;
  cout << "Matrix to be doubled using the + operator: " << N;
  cout << "This matrix should be double: " << W;

  // Remove argument in norm => enoguh with G.norm();
  // But is this less efficient since const is removed???
  cout << "Evaluating the norm method: " << G.norm() << endl;

  G *= N;  // This gives G*N

  cout << "Matrix from the *= operator: " << G;

  Matrix U = 4;
  U = G*N;
  cout << "Matrix from the * operator: " << U;

  Matrix A = G/2.5;
  cout << "Matrix from the / operator: " << A;

  vector<double> v1 = {1,2,3,4};
  Matrix P(v1);
  cout << "Matrix from 1-4 vector:  " << P;

  Matrix V = Matrix(4);
  V.fillMatrix(0);
  cout << "fillMatrix with argument 0 : " << V;

  Matrix  I = Matrix(6);
  I.fillIdentityMatrix();
  cout << "Matrix using the fillIdentityMatrix method: " << I;

  cout << endl <<
    "======== Comparing Exponential Functions =========" << endl <<
  "--- Start by testing a matrix filled with zeros ---" << endl;
  double ones[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  int ac = 4;

  Matrix C = matrixExp(N,0.1);
  cout <<
    "Results from our matrixExp for matrix filled with ones"
    "(tolerance of 0.1):" << C;

  // From the matlab source code:
  double* ONESE = r8mat_expm1(ac, ones);

  // Convert matlab results into our Matrix object for comparison
  Matrix ONESEM =  Matrix(r8matToVec(ONESE, 16));

  cout <<
    "Results from r8mat_expm1 for matrix filled with ones: "
       << ONESEM;

Matrix MDIFFO = C - ONESEM;
  cout << "Difference for the matrix filled with ones"
    " when using r8mat_expm1 and our implementation:"
       << MDIFFO;

  cout << "--- We now define a second matrix for the test ---" << endl;
  vector<double> t = {1, 3, 10, 45, 12, 3, 5, 0, 12, 1, 3, 7, 19, 4, 9, 6};
  double tArray[16] = {1, 3, 10, 45, 12, 3, 5, 0, 12, 1, 3, 7, 19, 4, 9, 6};

  Matrix T = Matrix(t);
  cout << "This is the test-matrix."
    " It will be used for the next comparison"
    " of the exponential functions:" << T;

  Matrix CC = matrixExp(T,0.001);
  cout << "Results from our matrixExp for the test-matrix"
    " (tolerance of 0.001):" << CC;

  // From the matlab source code:
  double* TE = r8mat_expm1(ac, tArray);

  // Convert matlab results into our Matrix object for comparison
  Matrix TEM =  Matrix(r8matToVec(TE, 16));

  cout << "Results from r8mat_expm1 for the test-matrix:" << TEM;

  Matrix MDIFF = CC - TEM;
  cout << "Difference for the test-matrix"
    " when using r8mat_expm1 and our implementation:"
       << MDIFF;

  // getch();
  return 0;
}
