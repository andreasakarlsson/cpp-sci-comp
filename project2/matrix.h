#ifndef MATRIX_H
#define MATRIX_H

using namespace std;

class Matrix{

private:

	int m_size;		// size of matirx
	vector<double> Mat;	// vector with values. Length m_size^2

public:

	Matrix(int m);
	Matrix(const Matrix& M);
	Matrix(const vector<double>& V);
	Matrix& operator=(const Matrix& M);
	Matrix operator+(const Matrix& M);
	Matrix operator-(const Matrix& M);
 	Matrix operator/(const double D);
	Matrix& operator+=(const Matrix& M);
	Matrix& operator*=(const Matrix& M);
	Matrix operator*(const Matrix& M);
	friend ostream& operator<<(ostream& stream, const Matrix& matrix);
	double norm();
	void printMatrix() const;
	void fillMatrix(int max=10);
	void fillIdentityMatrix();
	int getSize();
	double getVal(int i);

};

#endif
