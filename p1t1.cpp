#include <iostream>
#include <cmath>
#include <vector>  // vec.clear();
using namespace std;


// taylor expansion for sin with Horner's scheme
double sinTaylorH(const int nTerms, double x) {
  double hornerScheme(vector<double> Poly, int degree, double x);
  vector<double> a (2*nTerms); //numb. of polyn. coefficients should be the deg. + 1
  for (int n = 0; n < nTerms; n++)
    a[ 2 * n + 1] = ( (1-2*(n%2)));
  return hornerScheme(a,2*nTerms-1, x);
}

// taylor expansion for sin with Horner's scheme
double cosTaylorH(const int nTerms, double x) {
  double hornerScheme(vector<double> Poly, int degree, double x);
  vector<double> b (2*nTerms); //numb. of polyn. coeficients should be the deg. + 1
  for (int n = 0; n < nTerms; n++)
    b[2 * n] = ((1 - 2 * (n % 2)));
  return hornerScheme(b, 2 * nTerms - 2, x);
}

// Simple implmentation of Horner's Scheme
double hornerScheme(vector<double> PolynomCoef, const int pDegree, double x) {
  double hSchemeRes = 0.0;
  for (int i=pDegree; i >=0; i--){
    hSchemeRes =  (PolynomCoef[i] + x * hSchemeRes) /  ((i == 0) ? 1 : i);
    // cout << hSchemeRes << endl;
  }
  return hSchemeRes;
}

int main() {
  vector<double> x_vec = {-1, 1, 2, 3, 5, 10};
  vector<unsigned int> N_vec = {1, 10, 20, 100};

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
