#include <iostream>
#include <cmath>
#include <conio.h> // loads getch();
#include <vector>
//#include "r8lib.h"
using namespace std;




class Matrix{

private:
	int m_size;			// size of matirx
	vector<double> Mat;	// vector with values. Length m_size^2

public:

	void printMatrix() const;
	double norm(const Matrix&);

	// CONSTRUCTORS

	Matrix(int m){
		m_size = m;
		Mat.assign(m*m, 1);
	}

	// Copy-constructor?
	Matrix(const Matrix& M) {
		m_size = M.m_size;
		Mat = M.Mat;
	}

	// Constructor that only takes vector.
	Matrix(const vector<double>& V) {
		Mat = V;
		m_size = sqrt((int)V.size());
	}

	// OVERLOADING

	Matrix& operator=(const Matrix& M){
	 	m_size = M.m_size;
	 	Mat = M.Mat;
	 	return *this;
	}

	Matrix operator+(const Matrix& M){
		vector<double> res(m_size*m_size);
		Matrix tempM = Matrix(m_size);
 		for(int i=0; i<m_size*m_size; i++){
 			tempM.Mat[i] = Mat[i] + M.Mat[i];
 		}
	 	return tempM;
 	}

 	Matrix operator/(const double D){
		vector<double> res(m_size*m_size);
		Matrix tempM = Matrix(m_size);
 		for(int i=0; i<m_size*m_size; i++){
 			tempM.Mat[i] = Mat[i]/D;
 		}
	 	return tempM;
 	}

	Matrix& operator+=(const Matrix& M){
		vector<double> res(m_size*m_size);
 		for(int i=0; i<m_size*m_size; i++){
 			Mat[i] = Mat[i] + M.Mat[i];
 		}
	 	return *this;
 	}


	Matrix& operator*=(const Matrix& M){
		vector<double> prodVector(M.m_size*M.m_size);
		for(int i=0; i<M.m_size; i++){
			for(int j=0; j<M.m_size; j++){
				double sum = 0;
				for(int k=0; k<M.m_size; k++){
					sum += Mat[i*m_size+k]*M.Mat[m_size*k+j];
				}
				prodVector[i*m_size+j] = sum;
			}
		}
		Mat = prodVector;
		return *this;
	}

	Matrix operator*(const Matrix& M){
		//vector<int> prodVector(M.m_size*M.m_size);
		Matrix tempM = Matrix(m_size);
		for(int i=0; i<M.m_size; i++){
			for(int j=0; j<M.m_size; j++){
				double sum = 0;
				for(int k=0; k<M.m_size; k++){
					sum += Mat[i*m_size+k]*M.Mat[m_size*k+j];
				}
				//prodVector[i*m_size+j] = sum;
				tempM.Mat[i*m_size+j] = sum;
			}
		}
		//Mat = prodVector;
		return tempM;
	}

	// ADDITIONAL FUNCTIONS

	// Calc. the 1-norm of the matrix. Sum the values in each column and
	// returns the biggest value.
	double norm(){
		double max = 0;
		for(int i=0; i<m_size; i++){
			double sum = 0;
			for(int j=0; j<m_size; j++){
				sum += Mat[m_size*j+i];
			}
			if(sum>max) max = sum;
		}
		return max;
	}

	// prints the matrix 
	void printMatrix(){
		cout << "\n " << m_size << "x" << m_size << " matrix:\n";
		cout << "--------------------\n"; 
		for(int i = 0; i<m_size*m_size; i++){
			cout << " " << Mat[i];
			if((i+1)%m_size==0) cout << "\n";
		}
		cout << "\n";
	}

	// Fill existing matrix with ints between 0 and 10
	void fillMatrix(int max=10){
		for(int i = 0; i<m_size*m_size; i++){
			Mat[i] = 0 + (rand() % (int)(max - 0 + 1));
			//output = min + (rand() % (int)(max - min + 1))
		}
	}

	void fillIdentityMatrix(){
		for(int i = 0; i<m_size*m_size; i++){
			if(i%(m_size+1) == 0)
				Mat[i] = 1;
			else 
				Mat[i] = 0;
		}
	}

	int getSize(){
		return m_size;
	}

};



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
	getch();


	return 0;
}