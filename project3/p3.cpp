#include <iostream>
#include <cmath>
#include <cstdio>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "curvStraight.h"
#include "curvExp.h"
#include "Curvebase.h"

#ifdef _WIN32
#include <conio.h> // loads getch();
#endif


using namespace std;


class Domain {

public:
	Domain(Curvebase& s1, Curvebase& s2, Curvebase& s3, Curvebase& s4){
		sides[0] = &s1;
		sides[1] = &s2;
		sides[2] = &s3;
		sides[3] = &s4;

		if(!check_consistency(1e-5)){
			sides[0] = sides[1] = sides[2] = sides[3] = nullptr;
			cout << " Not consistent. All curves set to zero.\n" << endl;
		}
	}


	bool check_consistency(double epsilon = 1e-3){


		// Check that the end of a curve is connected to the start of the next.
		if (abs(sides[0]->x(1) - sides[1]->x(0))>epsilon || abs(sides[0]->y(1) - sides[1]->y(0))>epsilon ){
			cout << "\n Curve 0 and 1 not connected \n" << endl;
			return false;
		}
		if ( abs(sides[1]->x(1) - sides[2]->x(1))>epsilon || abs(sides[1]->y(1) - sides[2]->y(1))>epsilon ){
			cout << "\n Curve 1 and 2  not connected \n" << endl;
			return false;
		}
		if ( abs(sides[2]->x(0) - sides[3]->x(1))>epsilon || abs(sides[2]->y(0) - sides[3]->y(1))>epsilon ){
			cout << "\n Curve 2 and 3 not connected \n" << endl;
			return false;
		}
		if ( abs(sides[3]->x(0) - sides[0]->x(0))>epsilon || abs(sides[3]->y(0) - sides[0]->y(0))>epsilon ){
			cout << "\n Curve 3 and 0  not connected \n" << endl;
			return false;
		}

		return true;
	}

  void generate_grid(int n, int m, bool stretch = false){
		if(sides[0]==nullptr){
			cout << " Grid point generation aborted due to nullptr sides.\n" << endl;
			return;
		}
		if((n<1)||(m<1)) {m=20;n=20; cout << "set m and n to 20" << endl;}
		if(n_ != 0) { // if n_0 non-zero => grid allready exist. Must be deleted.
			delete [] x_;
			delete [] y_;
		}
		n_ = n;
		m_ = m;
		x_ = new double[(m_+1)*(n_+1)];
		y_ = new double[(m_+1)*(n_+1)];
		double h1 = 1.0/n_;
		double h2 = 1.0/m_;


		// The eight corner values (x & y). Calculated once
		// before the loop instead of for every iteration to
		// increase speed.
		double s0x0 = sides[0]->x(0);
		double s1x0 = sides[1]->x(0);
		double s3x1 = sides[3]->x(1);
		double s2x1 = sides[2]->x(1);

		double s0y0 = sides[0]->y(0);
		double s1y0 = sides[1]->y(0);
		double s3y1 = sides[3]->y(1);
		double s2y1 = sides[2]->y(1);


		// Values for the sides are calculated before the loop
		// this way we do m+1 calculations instead of
		// (m+1)*(m+1).
		double x_i_s0[(n_+1)];
		double x_i_s2[(n_+1)];
		double y_i_s0[(n_+1)];
		double y_i_s2[(n_+1)];

		for(int i = 0; i<= n_; i++){
			x_i_s0[i] = sides[0]->x(i*h1);
			x_i_s2[i] = sides[2]->x(i*h1);
			y_i_s0[i] = sides[0]->y(i*h1);
			y_i_s2[i] = sides[2]->y(i*h1);
		}

		double x_j_s3[(m_+1)];
		double x_j_s1[(m_+1)];
		double y_j_s3[(m_+1)];
		double y_j_s1[(m_+1)];

		for(int j = 0; j<= m_; j++){
		  double ss; // ss is a stretched version of j*h2
		  if (stretch) ss = (exp(1.5*j*h2)-1)/(exp(1.5)-1); else ss = j*h2;
		  x_j_s3[j] = sides[3]->x(ss);
		  x_j_s1[j] = sides[1]->x(ss);
		  y_j_s3[j] = sides[3]->y(ss);
		  y_j_s1[j] = sides[1]->y(ss);
		}


		for(int i = 0; i <= n_; i++){
			for(int j = 0; j <= m_; j++){
				x_[j+i*(m_+1)] = phi1(i*h1)*x_j_s3[j]
					+ phi2(i*h1)*x_j_s1[j]
					+ phi1(j*h2)*x_i_s0[i]
					+ phi2(j*h2)*x_i_s2[i]
					- phi1(i*h1)*phi1(j*h2)*s0x0
					- phi2(i*h1)*phi1(j*h2)*s1x0
					- phi1(i*h1)*phi2(j*h2)*s3x1
					- phi2(i*h1)*phi2(j*h2)*s2x1;


				y_[j+i*(m_+1)] = phi1(i*h1)*y_j_s3[j]
					+ phi2(i*h1)*y_j_s1[j]
					+ phi1(j*h2)*y_i_s0[i]
					+ phi2(j*h2)*y_i_s2[i]
					- phi1(i*h1)*phi1(j*h2)*s0y0
					- phi2(i*h1)*phi1(j*h2)*s1y0
					- phi1(i*h1)*phi2(j*h2)*s3y1
					- phi2(i*h1)*phi2(j*h2)*s2y1;

			}
		}
	}

	// Saves the x & y-values as a single array in a binary file.
	void save2file(const char* fname = "outfile.bin"){
		// Create new array that contains x_ and y_
		int sizeV = (m_+1)*(n_+1);
		double * result = new double[sizeV + sizeV];
		copy(x_, x_ + sizeV, result);
		copy(y_, y_ + sizeV, result + sizeV);

		FILE *fil;
		fil = fopen(fname,"wb");
		fwrite(result,sizeof(double),2*sizeV,fil);
		fclose(fil);
	}

	~Domain(){
		if (n_ > 0){
			delete [] x_;
			delete [] y_;
		}
	}


private:
	Curvebase *sides[4];
	int m_ = 0; // initialized to zero to enable class to check if
	int n_ = 0; // a grid generation been done.

	double *x_,*y_;

	inline double phi1(double w){
		return 1.0 - w;
	}
	inline double phi2(double w){
		return w;
	}
};


int main() {

	curvExp E(-10,5); // lower boundry

	curvStraight A(-10,5,0,0); // For initial test use lower boundry that is straight
	curvStraight B(0,3,1,5);
	curvStraight C(-10,5,0,3);
	curvStraight D(0,3,1,-10);

	// PRINT OUT TO SHOW THAT THE CURVE OBJECTS GIVE CORRECT VALUES

	cout << "\n Curve E:" << endl;

	cout << " (" << E.x(0) << ", " << E.y(0) << ")" << endl;
	cout << " (" << E.x(0.25) << ", " << E.y(0.25) << ")" << endl;
	cout << " (" << E.x(0.5) << ", " << E.y(0.5) << ")" << endl;
	cout << " (" << E.x(0.75) << ", " << E.y(0.75) << ")" << endl;
	cout << " (" << E.x(1) << ", " << E.y(1) << ")" << endl;

	cout << "\n Curve A:" << endl;

	cout << " (" << A.x(0) << ", " << A.y(0) << ")" << endl;
	cout << " (" << A.x(0.25) << ", " << A.y(0.25) << ")" << endl;
	cout << " (" << A.x(0.5) << ", " << A.y(0.5) << ")" << endl;
	cout << " (" << A.x(0.75) << ", " << A.y(0.75) << ")" << endl;
	cout << " (" << A.x(1) << ", " << A.y(1) << ")" << endl;


	cout << "\n Curve B:" << endl;

	cout << " (" << B.x(0) << ", " << B.y(0) << ")" << endl;
	cout << " (" << B.x(0.25) << ", " << B.y(0.25) << ")" << endl;
	cout << " (" << B.x(0.5) << ", " << B.y(0.5) << ")" << endl;
	cout << " (" << B.x(0.75) << ", " << B.y(0.75) << ")" << endl;
	cout << " (" << B.x(1) << ", " << B.y(1) << ")" << endl;


	cout << "\n Curve C:" << endl;

	cout << " (" << C.x(0) << ", " << C.y(0) << ")" << endl;
	cout << " (" << C.x(0.25) << ", " << C.y(0.25) << ")" << endl;
	cout << " (" << C.x(0.5) << ", " << C.y(0.5) << ")" << endl;
	cout << " (" << C.x(0.75) << ", " << C.y(0.75) << ")" << endl;
	cout << " (" << C.x(1) << ", " << C.y(1) << ")" << endl;


	cout << "\n Durve D:" << endl;

	cout << " (" << D.x(0) << ", " << D.y(0) << ")" << endl;
	cout << " (" << D.x(0.25) << ", " << D.y(0.25) << ")" << endl;
	cout << " (" << D.x(0.5) << ", " << D.y(0.5) << ")" << endl;
	cout << " (" << D.x(0.75) << ", " << D.y(0.75) << ")" << endl;
	cout << " (" << D.x(1) << ", " << D.y(1) << ")" << "\n" << endl;

	// Create and save a grid (formed by the four straight
	// curves.)
	Domain Grid(A,B,C,D);
	Grid.generate_grid(49,19,false);
	Grid.save2file("task3.bin");

	// Create and save grid with the given exponential
	// function. The computational time for generation of points
	// measured.
	Domain Grid2(E,B,C,D);

	clock_t t;
	t = clock(); // start timing

	Grid2.generate_grid(49,19,false);

	t = clock() - t;

  	printf (" It took %d clicks (%f seconds).\n",(int) t,((float)t)/CLOCKS_PER_SEC);

	Grid2.save2file("task4.bin");

	// Create stretched grid according to task 5
	Grid2.generate_grid(49,19,true);
	Grid2.save2file("task5.bin");


#ifdef _WIN32
	cout << "\n Press any key to quit..." << endl;
	getch();
#endif
	return 0;
}
