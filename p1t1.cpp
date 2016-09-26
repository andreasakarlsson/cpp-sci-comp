#include <iostream>
#include <cmath>
#include <vector>  // vec.clear();
using namespace std;

// recursive factorial function. Returns double.
double factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// factorial function using for-loop. Returns double.
double factorialFor(unsigned int n) {
  unsigned long long fac = 1;
  for(unsigned int j = 1; j<=n; ++j) {fac *=j;}
  return fac;
}

// taylor expansion for sin with Horner's scheme
double sinTaylorH(const int nTerms, double x) {
  double hornerScheme(vector<double> Poly, int degree, double x);
  vector<double> a (2*nTerms); //numb. of polyn. coefficients should be the deg. + 1
  for (int n = 0; n < nTerms; n++)
    a[2*n+1] = ( (1-2*(n%2)) / factorial(2 * n + 1));
  return hornerScheme(a,2*nTerms-1, x);
}

// taylor expansion for sin with Horner's scheme
double cosTaylorH(const int nTerms, double x) {
  double hornerScheme(vector<double> Poly, int degree, double x);
  vector<double> b (2*nTerms); //numb. of polyn. coeficients should be the deg. + 1
  for (int n = 0; n < nTerms; n++)
    b[2*n] = ( (1-2*(n%2)) / factorialFor(2 * n ));
  return hornerScheme(b,2*nTerms-2, x);
}

// Simple implmentation of Horner's Scheme
double hornerScheme(vector<double> PolynomCoef, const int pDegree, double x) {
  double hSchemeRes = 0.0;
  for (int i=pDegree; i >=0; i--){
    hSchemeRes = PolynomCoef[i]+x*hSchemeRes;
  }
  return hSchemeRes;
}

int main() {
  vector<double> x_vec = {-1, 1, 2, 3, 5, 10};
  vector<unsigned int> N_vec = {1, 20, 100, 1000};

  // Learning ranges in C++11
  for (auto N: N_vec) {
    cout << "Setting N = " << N << endl;
    for (auto x: x_vec) {
      cout << "  Setting x = " << x << endl;
      cout << "    sin using Horner, recursive factorial, error:     "
	   << abs(sin(x) - sinTaylorH(N, x)) << endl;
      cout << "    cos using Horner, factorial with for-loop, error: "
	   << abs(cos(x) - cosTaylorH(N, x)) << endl;
    }
  }
  return 0;
}
