#include <iostream>
#include <cmath>
#include <cstdlib>
// #include <conio.h> // loads getch();
#include <iomanip> // to get setw
#include <vector>
//#include <stdexcept>
#include <cstdio>
using namespace std;

class Curvebase {

protected:
	double integrate(double (Curvebase::*funcp)(double), double a, double b, double tol1){
		double Int1 = 0.0;
		double Int2 = 0.0;
		double errest = 10;
		int n = 1;

		while (errest > 15*tol1){
			vector<double> lims (n+1);
			double h = (b-a)/(n);
			for(int i = 0; i<n+1; i++) {lims[i] = a+i*h;}
			for(int j = 0; j<n; j++) {Int1 += I(funcp, lims[j],lims[j+1]);}
			errest = abs(Int1-Int2);
			n *= 2;  // update number of intervals used
			Int2 = Int1;
			Int1 = 0;
		}
		return Int2;
	}
	double I(double (Curvebase::*funcp)(double),double a, double b){
		return (b-a)/6*((this->*funcp)(a)+4*(this->*funcp)((a+b)/2)+(this->*funcp)(b));
		//return (b-a)/6*(funcp(a)+4*funcp((a+b)/2)+funcp(b));
	}

	double newton(double (Curvebase::*fp1)(double,double), double (Curvebase::*dfp1)(double), double s, double initG, double tol1){
	double x1;
	double x = initG;
	int it = 0; int maxit = 1000;
	double err = tol1 + 1;
	while( err > tol1 && it < maxit){
		x1 = x-(this->*fp1)(x,s)/(this->*dfp1)(x);
		err = abs(x1-x);
		x = x1;
		it++;
	}
	if(err >= tol1) {cout << s  << initG << " - No convergence Newton. Res = " << x << x1 << endl;}
	return x;
	}


// 	//double pmin; // not necessary according to Hanke.
// 	//double pmax; // not necessary according to Hanke.
	double a_;
	double b_;
	//int rev; // curve orientation. NOT Currently used (determine depending on values of a_ and b_)


	double (Curvebase::*l)(double) = &Curvebase::f; // class member function pointer.
	double (Curvebase::*fpP)(double,double) = &Curvebase::fp;
	double (Curvebase::*dfpP)(double) = &Curvebase::dfp;

	double tol = 1e-6; // tolerance. Set maybe as argument(?)
	//double lb = integrate(l,a_,b_,tol/100);  // length ??

	virtual double xp(double p) = 0;
	virtual double yp(double p) = 0;
	virtual double dxp(double p) = 0;
	virtual double dyp(double p) = 0;

	double f(double p){
		return sqrt(dxp(p)*dxp(p)+dyp(p)*dyp(p)); // see Inheritance slide nr. 10.
	}
	double fp(const double p, const double s){
		return integrate(l,a_,p,tol/100)-s*integrate(l,a_,b_,tol/100);
	}
	double dfp(const double p){
		return f(p); // this is not so nice.
	}


public:
	Curvebase(double a = 0.0, double b = 1.0) : a_(a), b_(b){}
	 // s is the normalized arc length parameter. Interval. [0,1]

	double x(double s){
		if((s<0)||(s>1)){ cout << "Got invalid s" << endl;}
		double p = newton(fpP,dfpP,s,0.5*(a_+b_),tol);
		return xp(p);
	}
	double y(double s){
		if((s<0)||(s>1)){ cout << "Got invalid s" << endl;}
		double p = newton(fpP,dfpP,s,0.5*(a_+b_),tol);
		return yp(p);
	}

	virtual ~Curvebase(){}

};



class Domain {

public:
	Domain(Curvebase& s1, Curvebase& s2, Curvebase& s3, Curvebase& s4){
		sides[0] = &s1;
		sides[1] = &s2;
		sides[2] = &s3;
		sides[3] = &s4;
		
		if(~check_consistency(1e-4)){
			sides[0] = sides[1] = sides[2] = sides[3] = NULL;
			// what more to do?
		}
	}


	bool check_consistency(double epsilon = 1e-3){
		
		bool res = true;
		
		// Check that the end of a curve is connected to the start of the next.
		if (abs(sides[0]->x(1) - sides[1]->x(0))>epsilon || abs(sides[0]->y(1) - sides[1]->y(0))>epsilon ){
			cout << "\n Curve 0 and 1 not connected \n" << endl;
			res = false;
		}
		if ( abs(sides[1]->x(1) - sides[2]->x(1))>epsilon || abs(sides[1]->y(1) - sides[2]->y(1))>epsilon ){
			cout << "\n Curve 1 and 2  not connected \n" << endl;
			res = false;
		}
		if ( abs(sides[2]->x(0) - sides[3]->x(1))>epsilon || abs(sides[2]->y(0) - sides[3]->y(1))>epsilon ){
			cout << "\n Curve 2 and 3 not connected \n" << endl;
			res = false;
		}
		if ( abs(sides[3]->x(0) - sides[0]->x(0))>epsilon || abs(sides[3]->y(0) - sides[0]->y(0))>epsilon ){
			cout << "\n Curve 3 and 0  not connected \n" << endl;
			res = false
		}
		
		// check order and other requirements of curves
		
		return res;
	}

	void generate_grid(int n, int m){
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

		// xi1 in the sildes is i*h1 and xi2 is j*h2

		for(int i = 0; i <= n_; i++){
			for(int j = 0; j <= m_; j++){
				cout << "coordinate: i=" << i << ", j=" << j;
				x_[j+i*(m_+1)] = phi1(i*h1)*sides[3]->x(j*h2)
					+ phi2(i*h1)*sides[1]->x(j*h2)
					+ phi1(j*h2)*sides[0]->x(i*h1)
					+ phi2(j*h2)*sides[2]->x(i*h1)
					- phi1(i*h1)*phi1(j*h2)*sides[0]->x(0)
					- phi2(i*h1)*phi1(j*h2)*sides[1]->x(0)
					- phi1(i*h1)*phi2(j*h2)*sides[3]->x(1)
					- phi2(i*h1)*phi2(j*h2)*sides[2]->x(1);
				cout << "   x-value: " << x_[j+i*(m_+1)];

				y_[j+i*(m_+1)] = phi1(i*h1)*sides[3]->y(j*h2)
					+ phi2(i*h1)*sides[1]->y(j*h2)
					+ phi1(j*h2)*sides[0]->y(i*h1)
					+ phi2(j*h2)*sides[2]->y(i*h1)
					- phi1(i*h1)*phi1(j*h2)*sides[0]->y(0)
					- phi2(i*h1)*phi1(j*h2)*sides[1]->y(0)
					- phi1(i*h1)*phi2(j*h2)*sides[3]->y(1)
					- phi2(i*h1)*phi2(j*h2)*sides[2]->y(1);
				cout <<"   y-value: " << y_[j+i*(m_+1)] << endl;



				// PRINT OUT for testing code.
				// double ttest1 = phi1(i*h1)*sides[3]->x(j*h2);
				// double ttest2 = phi2(i*h1)*sides[1]->x(j*h2);
				// double ttest3 = phi1(j*h2)*sides[0]->x(i*h1);
				// double ttest4 = phi2(j*h2)*sides[2]->x(i*h1);

				// double rtest1 = phi1(i*h1)*phi1(j*h2)*sides[0]->x(0);
				// double rtest2 = phi2(i*h1)*phi1(j*h2)*sides[1]->x(0);
				// double rtest3 = phi1(i*h1)*phi2(j*h2)*sides[3]->x(0);
				// double rtest4 = phi2(i*h1)*phi2(j*h2)*sides[2]->x(0);

				// cout << ttest1 << " " << ttest2 << " " << ttest3 << " " << ttest4 << endl;
				// cout << rtest1 << " " << rtest2 << " " << rtest3 << " " << rtest4 << endl;
			}
		}
	}

	void save2file(){
		// Create new array that contains x_ and y_ that
		// can be written to a binary file.
		int sizeV = (m_+1)*(n_+1);
		double * result = new double[sizeV + sizeV];
		copy(x_, x_ + sizeV, result);
		copy(y_, y_ + sizeV, result + sizeV);

		FILE *fil;
		fil = fopen("outfile.bin","wb");
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

	double phi1(double w){
		return w;
	}
	double phi2(double w){
		return 1.0 - w;
	}
};


class curvStraight: public Curvebase {

public:
	curvStraight(double a1, double b1, int ori, double Sdim) : Curvebase() {
		if ( b1 < a1 ) {
			cout << "\n changed a & b order" << endl;
			a_ = b1;
			b_ = a1;
		} else {
    			a_ = a1;
    			b_ = b1;
		}
    	//rev = 1;
    	o_ = ori;
    	Sdim_ = Sdim;

	}

	~curvStraight(){}

private:

	double xp(double p){
		if(o_ == 0) return p;
		else return Sdim_;
	}
	double yp(double p){
		if(o_ == 0) return Sdim_;
		else return p;
	}
	double dxp(double p){
		if(o_ == 0) return 1.0;
		else return 0.0;
	}
	double dyp(double p){
		if(o_ == 0) return 0.0;
		else return 1.0;
	}


	double lb = integrate(l,a_,b_,tol/100);  // unused. Should be used in abstr. class
	// but since it access a virtual function that is not possible.
	int o_; // orientation parameter. Should be improved.
	double Sdim_;


};

class curvExp: public Curvebase {

public:
	curvExp(double a1, double b1) : Curvebase() {
		if ( b1 < a1 ) {
			cout << "\n changed a & b order" << endl;
			a_ = b1;
			b_ = a1;
		} else {
    			a_ = a1;
    			b_ = b1;
		}
    	//rev = 1;
	}

	~curvExp(){}

private:

	double xp(double p){
		return p;
	}
	double yp(double p){
		if(p<-3)
			return 0.5/(1+exp(-3*(p+6)));
		else
			return 0.5/(1+exp(3*p));
	}
	double dxp(double p){
		return 1.0;
	}
	double dyp(double p){
		if(p<-3)
			return (3/2)*exp(-3*(p+6))/((1+exp(-3*(p+6)))*(1+exp(-3*(p+6))));
		else
			return -(3/2)*exp(3*p)/((1+exp(3*p))*(1+exp(3*p)));
	}

	double lb = integrate(l,a_,b_,tol/100);
	//int o_;
	//double Sdim_;
};

int main() {

	curvExp E(-10,5); // correct lower boundry

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


	// Create and save grid (formed by the four straight curves.)

	Domain Grid(A,B,C,D);

	Grid.generate_grid(4,4);

	// Grid.save2file();

	// All code needed to read in MATLAB: fid = fopen('outfile.bin','r'); c = fread(fid,'double');
	// x = c(1:length(c)/2);  y = c(length(c)/2+1:end); plot(x,y,'*')

	Domain Grid2(E,B,C,D);
	Grid2.generate_grid(10,10);
	Grid2.save2file();

	cout << "\n Press any key to quit..." << endl;
	// getch();
	return 0;
}
