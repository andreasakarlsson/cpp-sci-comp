#include <iostream>
#include <cmath>

using namespace std;

double hornerScheme(double * polynomCoef, const int pDegree, double x);



// Horner's Scheme - Evaluated a polynom for a value x. Takes vector of polynomial coefficients (polynomCoef), 
// the degree of the polynomial (pDegree) and an evaluation value (x). 
double hornerScheme(double * polynomCoef, const int pDegree, double x) {
  double hSchemeRes = 0.0;
  for (int i=pDegree; i >=0; i--){
    hSchemeRes = polynomCoef[i]+x*hSchemeRes;
  }
  return hSchemeRes;
}

// recursive factorial function returns double
double factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// taylor expansion for sin
double sinTaylor(int N, double x) {
double sum = 0;
  for (int n = 0; n <= N; n++)
    sum += (pow(-1.0, n) * (pow(x, 2.0 * n + 1) / factorial(2 * n + 1)));
  return sum;
}

// taylor expansion for cos
double cosTaylor(int N, double x) {
double sum = 0;
  for (int n = 0; n <= N; n++)
    sum += (pow(-1.0, n) * (pow(x, 2.0 * n) / factorial(2 * n)));
  return sum;
}

double sinTaylorH(int N, double x) {
double a[2*N]; // number of polyn. coeficients should be the degree + 1
  for (int n = 0; n < N; n++)
    a[2*n+1] = ( (1-2*(n%2)) / factorial(2 * n + 1));
  return hornerScheme(a,2*N-1, x);
}



int main() {
  double factorial(int n);
  int N = 100;
  double x = 3.1416;
  cout << sinTaylor(N, x) << endl;
  cout << cosTaylor(N, x) << endl;
  cout << abs(sin(x) - sinTaylor(N, x)) << endl;
  cout << abs(cos(x) - cosTaylor(N, x)) << endl;
  
  cout << "1 sin (Horner) = "<< sinTaylorH(N, x) << endl;
  cout << "2 sin (Hroner) = "<< sinTaylorH(N, x) << endl;
  
  double poly[3];
  double y = 2;
  b[0] = 5;
  b[1] = 2;
  b[2] = 3;
  cout << hornerScheme(poly, 2, y) << endl;
  
  
  return 0;
}
