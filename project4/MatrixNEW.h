#ifndef MATRIXNEW_H
#define MATRIXNEW_H

using namespace std;

class MatrixNEW{

private:

	int m_;		// size of matirx
	int n_;
	double *Mat;	// array with values.

public:

	MatrixNEW(int m, int n, double t);			// constr
	MatrixNEW(const MatrixNEW& M);
	MatrixNEW& operator=(const MatrixNEW& M);

	MatrixNEW(MatrixNEW&& M) noexcept;			// move
	MatrixNEW& operator=(MatrixNEW&& M);

	MatrixNEW operator+(const MatrixNEW& M) const;
	MatrixNEW operator-(const MatrixNEW& M) const;
	MatrixNEW operator*(const MatrixNEW& M) const;
	MatrixNEW operator*(const double& d) const;
	MatrixNEW operator/(const double& d) const;

	MatrixNEW& operator+=(const MatrixNEW& M);
	MatrixNEW& operator-=(const MatrixNEW& M);
	MatrixNEW& operator*=(const MatrixNEW& M);

	friend ostream& operator<<(ostream& stream, const MatrixNEW& matrix);
	double& operator()(const int i, const int j) const;
	double* getMatrix();

};

#endif