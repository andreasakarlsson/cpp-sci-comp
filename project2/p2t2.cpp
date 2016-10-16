#include <iostream>
#include <cmath>
//#include <conio.h> // loads getch();
#include <vector>
//#include "r8lib.h"
using namespace std;



class Matrix{

private:
	int m_size;			// size of matirx
	vector<int> Mat;	// vector with values. Length m_size^2

public:

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
	Matrix(const vector<int>& V) {
		Mat = V;
		m_size = sqrt((int)V.size());
	}

	// OVERLOADING

	Matrix& operator=(const Matrix& M){
	 	m_size = M.m_size;
	 	Mat = M.Mat;
	 	return *this;
	}

	Matrix& operator+=(const Matrix& M){
		vector<int> res(m_size*m_size);
 		for(int i=0; i<m_size*m_size; i++){
 			Mat[i] = Mat[i] + M.Mat[i];
 		}
	 	return *this;
 	}


	Matrix& operator*=(const Matrix& M){
		vector<int> prodVector(M.m_size*M.m_size);
		for(int i=0; i<M.m_size; i++){
			for(int j=0; j<M.m_size; j++){
				int sum = 0;
				for(int k=0; k<M.m_size; k++){
					sum += Mat[i*m_size+k]*M.Mat[m_size*k+j];
				}
				prodVector[i*m_size+j] = sum;
			}
		}
		Mat = prodVector;
		return *this;
	}

	// ADDITIONAL FUNCTIONS

	// Calc. the 1-norm of the matrix. Sum the values in each column and
	// returns the biggest value.
	double norm(const Matrix& M){
		double max = 0;
		for(int i=0; i<M.m_size; i++){
			double sum = 0;
			for(int j=0; j<M.m_size; j++){
				sum += M.Mat[M.m_size*j+i];
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
	void fillMatrix(){
		for(int i = 0; i<m_size*m_size; i++){
			Mat[i] = 0 + (rand() % (int)(10 - 0 + 1));
			//output = min + (rand() % (int)(max - min + 1))
		}
	}





};

int main() {

	Matrix M(4);
	Matrix N(M);
	cout << " N: \n";
	N.printMatrix();

	//Matrix O(6);
	//O.fillMatrix();
	//O.printMatrix();

	Matrix HE = 3; //uses copy const  works!!

	Matrix G = Matrix(4);
	G.fillMatrix();

	cout << " G (matrix filled with func fillMatrix): \n";
	G.printMatrix();


	G += N;
	cout << " new G ( G += N): \n";
	G.printMatrix();

	//cout << "G-norm: " << G.norm(G) << "\n";

	G *= N;  // This gives G*N

	cout << " multiplication of new G and N: \n";
	G.printMatrix();



	std::vector<int> v1 = {1,2,3,4};
 	Matrix P(v1);
 	P.printMatrix();


	//cout << "\n Press any key to quit...";
	//getch();


	return 0;
}