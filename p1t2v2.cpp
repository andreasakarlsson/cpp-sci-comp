#include <iostream>
#include <cmath>
#include <conio.h> // loads getch();
#include <vector>
#include <ctime>
using namespace std;

double f(double x);
double ASI(double (*fp)(double),double a, double b, double tol);
double ASI2(double (*fp)(double),double a, double b, double tol);
double I(double (*fp)(double),double a, double b);


double I(double (*funcp1)(double),double a, double b){
	return (b-a)/6*(funcp1(a)+4*funcp1((a+b)/2)+funcp1(b));
}

double f(double x){
	return 1+sin(exp(3*x));
}

double ASI(double (*funcp)(double),double a, double b, double tol){
	double Int1 = I(funcp,a,b);
	double Int2 = I(funcp,a,(a+b)/2)+I(funcp,(a+b)/2,b);
	double errest = abs(Int1-Int2);
	if (errest < 15*tol)
		return Int2;
	return (ASI(funcp,a,(a+b)/2,tol/2) + ASI(funcp, (a+b)/2,b,tol/2));
}

double ASI2(double (*funcp3)(double),double a, double b, double tol){
	double Int1 = 0.0;
	double Int2 = 0.0;
	double errest = 100;
	int n = 1;

	while (errest > 15*tol){  
		//could set errest calc in the loop condition, 
		//	but it makes it less readable (?).

		vector<double> lims (n+1);
		double h = (b-a)/(n);
		for(int i = 0; i<n+1; i++) {lims[i] = a+i*h;}
		for(int j = 0; j<n; j++) {Int1 += I(funcp3, lims[j],lims[j+1]);}

		// print out to follow loop process.
		//cout << " n = " << n << " " << lims[2*n-2] << "  " << lims[2*n-1] << " Int1 = " << Int1 << endl;
			
		errest = abs(Int1-Int2);
		n *= 2;  // update number of intervals used
		Int2 = Int1;
		Int1 = 0;
	
	}
	return Int2;
}


int main() {

  cout.precision(15);  // gives higher precision

  double (*fp)(double);  // function pointer
  fp = f;

  double MATLAB = 2.500809110336166;
  cout << " correct val (MATLAB): I = 2.500809110336166" << endl;
  cout << " intI    -1 to 1:      I = " << I(fp,-1,1) << endl;

  double tollerance = 0.01;  // 0.001    0.0001

  cout << " input tollerance value: \n" << " tol = ";
  cin >> tollerance; 


  cout << " ASI                   I = " << ASI(fp, -1,1,tollerance) << endl;
  cout << " ASI2                  I = " << ASI2(fp, -1,1,tollerance) << endl;
  

  cout << " ASI-MATLAB            e = " << abs(ASI(fp, -1,1,tollerance)-MATLAB) << endl;
  cout << " ASI2-MATLAB           e = " << abs(ASI2(fp, -1,1,tollerance)-MATLAB) << endl;



  getch();
  return 0;
}