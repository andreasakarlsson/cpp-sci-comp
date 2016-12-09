#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "MatrixNEW.h"
using namespace std;


// CONSTRUCTORS

MatrixNEW::MatrixNEW(int m = 0, int n = 0, double test = 0.0) : m_(m), n_(n), Mat(nullptr) {
	if (m_*n_ > 0) {
		Mat = new double[m_*n_];
		std::fill(Mat,Mat+m*n,test);
	}
}

// Copy-constructor
MatrixNEW::MatrixNEW(const MatrixNEW& M) : m_(M.m_), n_(M.n_), Mat(nullptr){
	if (m_*n_ > 0){
		Mat = new double[m_*n_];
		std::copy(M.Mat, M.Mat+m_*n_, Mat);
	}
}

// Copy-assignment
MatrixNEW& MatrixNEW::operator=(const MatrixNEW& M) {
	if(this != &M){
		if(m_*n_ != M.m_*M.n_){
			if(Mat != nullptr) delete [] Mat;
			if(M.Mat != nullptr) Mat = new double[M.m_*M.n_];
		}
		m_ = M.m_; n_ = M.n_;
		std::copy(M.Mat, M.Mat+m_*n_, Mat);
	}
 	return *this;
}

// Move
MatrixNEW::MatrixNEW(MatrixNEW&& M) noexcept : m_(M.m_), n_(M.n_), Mat(M.Mat) {
	M.m_ = 0; M.n_ = 0; M.Mat = nullptr;
}

// Move-assignmet
MatrixNEW& MatrixNEW::operator=(MatrixNEW&& M){
	m_ = M.m_;
	n_ = M.n_;
	if(Mat != nullptr) delete [] Mat;
	Mat = M.Mat;
	M.m_ = 0; M.n_ = 0; M.Mat = nullptr;
	return *this;
}

MatrixNEW MatrixNEW::operator+(const MatrixNEW& M) const {
	MatrixNEW tempM(m_, n_);
	if(m_ == M.m_ && n_ == M.n_){
		for(int i = 0; i < m_*n_; i++){
			tempM.Mat[i] = Mat[i] + M.Mat[i];
		}
	}
 	return tempM;
}

MatrixNEW MatrixNEW::operator-(const MatrixNEW& M) const {
	MatrixNEW tempM(m_, n_);
	if(m_ == M.m_ && n_ == M.n_){
		for(int i = 0; i < m_*n_; i++){
			tempM.Mat[i] = Mat[i] - M.Mat[i];
		}
	}
 	return tempM;
}

MatrixNEW MatrixNEW::operator*(const MatrixNEW& M) const {
	MatrixNEW tempM(m_, n_);
	if(m_ == M.m_ && n_ == M.n_){
		for(int i = 0; i < m_*n_; i++){
			tempM.Mat[i] = Mat[i] * M.Mat[i];
		}
	}
 	return tempM;
}

//  Multiplication with scalar (from right side)
MatrixNEW MatrixNEW::operator*(const double& d) const {
	MatrixNEW tempM(m_, n_);
	for(int i = 0; i < m_*n_; i++){
		tempM.Mat[i] = d*Mat[i];
	}
 	return tempM;
}

MatrixNEW MatrixNEW::operator/(const double& d) const {
	MatrixNEW tempM(m_, n_);
	for(int i = 0; i < m_*n_; i++){
		tempM.Mat[i] = Mat[i]/d;
	}
 	return tempM;
}

MatrixNEW& MatrixNEW::operator+=(const MatrixNEW& M){
	for(int i = 0; i < m_*n_; i++){
		Mat[i] = Mat[i] + M.Mat[i];
	}
 	return *this;
}

MatrixNEW& MatrixNEW::operator-=(const MatrixNEW& M){
	for(int i = 0; i < m_*n_; i++){
		Mat[i] = Mat[i] - M.Mat[i];
	}
 	return *this;
}

MatrixNEW& MatrixNEW::operator*=(const MatrixNEW& M){
	for(int i = 0; i < m_*n_; i++){
		Mat[i] = Mat[i] * M.Mat[i];
	}
 	return *this;
}

// Additional methods

ostream& operator<<( ostream &output, const MatrixNEW &M ) {
  output << endl << M.m_ << "x" << M.n_ << " matrix:" << endl;
  output << "--------------------" << endl;
  for(int i = 0; i < M.m_ * M.n_; i++){
    output << " " << M.Mat[i];
    if((i + 1) % M.n_ == 0) output << endl;
  }
  output << endl;
  return output;
}

double MatrixNEW::getVal(const int& i, const int& j){
	return Mat[i+j*m_];
}
