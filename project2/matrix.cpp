#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "matrix.h"
using namespace std;


// CONSTRUCTORS

Matrix::Matrix(int m) {
	m_size = m;
	Mat.assign(m*m, 1);
}

// Copy-constructor?
Matrix::Matrix(const Matrix& M) {
	m_size = M.m_size;
	Mat = M.Mat;
}

// Constructor that only takes vector.
Matrix::Matrix(const vector<double>& V) {
	Mat = V;
	m_size = sqrt((int)V.size());
}

// OVERLOADING

Matrix& Matrix::operator=(const Matrix& M) {
 	m_size = M.m_size;
 	Mat = M.Mat;
 	return *this;
}

Matrix Matrix::operator+(const Matrix& M){
	vector<double> res(m_size * m_size);
	Matrix tempM = Matrix(m_size);
	for(int i=0; i < m_size * m_size; i++){
		tempM.Mat[i] = Mat[i] + M.Mat[i];
	}
 	return tempM;
}

Matrix Matrix::operator-(const Matrix& M){
	vector<double> res(m_size * m_size);
	Matrix tempM = Matrix(m_size);
	for(int i=0; i < m_size * m_size; i++){
		tempM.Mat[i] = Mat[i] - M.Mat[i];
	}
 	return tempM;
}

Matrix Matrix::operator/(const double D) {
	vector<double> res(m_size*m_size);
	Matrix tempM = Matrix(m_size);
	for(int i=0; i< m_size * m_size; i++){
		tempM.Mat[i] = Mat[i] / D;
	}
 	return tempM;
}

Matrix& Matrix::operator+=(const Matrix& M) {
	vector<double> res(m_size * m_size);
	for(int i=0; i< m_size * m_size; i++){
		Mat[i] = Mat[i] + M.Mat[i];
	}
 	return *this;
}


Matrix& Matrix::operator*=(const Matrix& M){
	vector<double> prodVector(M.m_size * M.m_size);
	for(int i=0; i < M.m_size; i++){
		for(int j=0; j < M.m_size; j++){
			double sum = 0;
			for(int k=0; k < M.m_size; k++){
				sum += Mat[i * m_size + k] * M.Mat[m_size * k + j];
			}
			prodVector[i * m_size + j] = sum;
		}
	}
	Mat = prodVector;
	return *this;
}

Matrix Matrix::operator*(const Matrix& M){
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

ostream& operator<<( ostream &output, const Matrix &M ) {
  output << endl << M.m_size << "x" << M.m_size << " matrix:" << endl;
  output << "--------------------" << endl;
  for(int i = 0; i < M.m_size * M.m_size; i++){
    output << " " << M.Mat[i];
    if((i + 1) % M.m_size == 0) output << endl;
  }
  output << endl;
  return output;
}

// ADDITIONAL FUNCTIONS

// Calc. the 1-norm of the matrix. Sum the values in each column and
// returns the biggest value.
double Matrix::norm(){
	double max = 0;
	for(int i=0; i<m_size; i++){
		double sum = 0;
		for(int j=0; j < m_size; j++){
			sum += Mat[m_size * j + i];
		}
		if(sum>max) max = sum;
	}
	return max;
}

// Fill existing matrix with ints between 0 and 10
void Matrix::fillMatrix(int max){
	for(int i = 0; i < m_size * m_size; i++){
		Mat[i] = 0 + (rand() % (int)(max - 0 + 1));
		//output = min + (rand() % (int)(max - min + 1))
	}
}

void Matrix::fillIdentityMatrix(){
	for(int i = 0; i < m_size * m_size; i++){
		if(i%(m_size + 1) == 0)
			Mat[i] = 1;
		else
			Mat[i] = 0;
	}
}

int Matrix::getSize(){
	return m_size;
}

double Matrix::getVal(int i){
	return Mat[i];
}
