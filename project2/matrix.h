#include <vector>

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
