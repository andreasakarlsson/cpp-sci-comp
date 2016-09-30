#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
using namespace std;

double f(double x);
double ASI(double (*fp)(double), double a, double b, double tol);
double ASI2(double (*fp)(double), double a, double b, double tol);
double I(double (*fp)(double), double a, double b);


double I(double (*funcp1)(double), double a, double b){
	return (b - a) / 6 * (funcp1(a) + 4 * funcp1((a + b) / 2) + funcp1(b));
}

double f(double x){
	return 1 + sin(exp(3 * x));
}

double recurisiveASI(double (*funcp)(double), double a, double b, double tol){
	double Int1 = I(funcp, a, b);
	double Int2 = I(funcp, a,(a + b) / 2) + I(funcp, (a + b) / 2, b);
	double errest = abs(Int1 - Int2);
	if (errest < 15 * tol)
		return Int2;
	return (recurisiveASI(funcp, a, (a + b) / 2, tol / 2) +
		recurisiveASI(funcp, (a + b) / 2, b, tol / 2));
}

double ASI(double (*funcp3)(double),double a, double b, double tol){
  double Int1 = 0.0;
  double Int2 = 0.0;
  int n = 1;

  while (true){
    Int1 = 0;
    vector<double> subInterval (n + 1);
    double h = (b - a) / (n);

    // calculating sub-intervals and store in a vector
    for(int i = 0; i < n + 1; i++) {subInterval[i] = a + i * h;}
    // calculating integral from sum of sub-intervals
    for(int j = 0; j < n; j++) {
      Int1 += I(funcp3, subInterval[j], subInterval[j + 1]);
    }

    // checking errest
    if (abs(Int1 - Int2) < 15 * tol) {
      break;
    }
    n *= 2;  // update number of intervals used
    Int2 = Int1;
  }
  return Int1;
}


int main()
{
  cout.precision(15);  // gives higher precision

  double (*fp)(double);  // function pointer
  fp = f;

  double MATLAB = 2.500809110336166;
  vector <double> toleranceVec = {0.01, 0.001, 0.0001};

  for (auto tolerance: toleranceVec) {
    cout << "Setting tolerance = " << tolerance << endl;

    cout << "  ASI recurisive \t I = " << recurisiveASI(fp, -1, 1, tolerance)
	 << "\t error = "
	 << abs(recurisiveASI(fp, -1, 1, tolerance) - MATLAB) << endl;
    cout << "  ASI for-loop   \t I = " << ASI(fp, -1, 1, tolerance)
	 << "\t error = "
	 << abs(ASI(fp, -1, 1, tolerance) - MATLAB) << endl;
  }
  return 0;
}
