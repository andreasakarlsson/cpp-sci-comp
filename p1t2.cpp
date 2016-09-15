#include <iostream>
#include <cmath>
#include <conio.h> // loads getch();
using namespace std;

double f(double x);
double ASI(double a, double b, double tol);
double I(double a, double b);
double I2(double a, double b);




double ASI(double a, double b, double tol){
	double Int1 = I(a,b);
	double Int2 = I2(a,b);
	double errest = abs(Int1-Int2);
	if (errest < 15*tol)
		return Int2;
	return (ASI(a,(a+b)/2,tol/2) + ASI((a+b)/2,b,tol/2));
}

double I(double a, double b){
	return (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
}

double I2(double a, double b){
	double g = (a+b)/2;
	return I(a,g)+I(g,b);
}

double f(double x){
	return 1+sin(exp(3*x));
}


int main() {
  cout << " correct val (MATLAB): I = 2.500809110336166" << endl;
  cout << " intI    -1 to 1 " << I(-1,1) << endl;
  cout << " intI2   -1 to 1 " << I2(-1,1) << endl;
  cout << " errest          " << abs(I(-1,1)-I2(-1,1)) << endl;

  double tollerance;
  cout << " tolerance" << endl;
  cin >> tolerance;
  cout << " ASI             " << ASI(-1,1,tolerance) << endl;


  getch();
  return 0;
}
