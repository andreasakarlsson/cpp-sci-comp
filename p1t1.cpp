#include <iostream>
#include <cmath>
#include <conio.h> // loads getch();
#include <vector>  // vec.clear();
using namespace std;
double hornerScheme(vector<double> Poly, int degree, double x);


// recursive factorial function returns double
double factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// factorial function using for-loop. Returns double.
double factorialFor(unsigned int n)
{
  //return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  unsigned long long fac = 1;
  for(unsigned int j = 1; j<=n; ++j) {fac *=j;}
  return fac;
}

// taylor expansion for sin using pow
double sinTaylor(int N, double x) {
double sum = 0;
  for (int n = 0; n <= N; n++)
    sum += (pow(-1.0, n) * (pow(x, 2.0 * n + 1) / factorial(2 * n + 1)));
  return sum;
}

// taylor expansion for cos using pow
double cosTaylor(int N, double x) {
double sum = 0;
  for (int n = 0; n <= N; n++)
    sum += (pow(-1.0, n) * (pow(x, 2.0 * n) / factorial(2 * n)));
  return sum;
}

// taylor expansion for sin with Horner's scheme
double sinTaylorH(const int nTerms, double x) {
  vector<double> a (2*nTerms); //numb. of polyn. coeficients should be the deg. + 1
  for (int n = 0; n < nTerms; n++)
    a[2*n+1] = ( (1-2*(n%2)) / factorialFor(2 * n + 1));
  return hornerScheme(a,2*nTerms-1, x);
}

// taylor expansion for sin with Horner's scheme
double cosTaylorH(const int nTerms, double x) {
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
  double factorial(int n);
  int N = 10;
  double x = 3.1415; //x = 3.1416;
  cout << "  sin(x) = "<< sinTaylor(N, x) << endl;
  cout << "1 sin(x) = "<< sinTaylorH(N, x) << endl;
  cout << "2 sin(x) = "<< sinTaylorH(N, x) << endl;
  cout << "  cos(x) = "<< cosTaylor(N, x) << endl;
  cout << "1 cos(x) = "<< cosTaylorH(N, x) << endl;
  cout << "2 cos(x) = "<< cosTaylorH(N, x) << endl;
  cout << "error sin " << abs(sin(x) - sinTaylorH(N, x)) << endl;
  cout << "error cos " << abs(cos(x) - cosTaylor(N, x)) << endl;


 // test of Horner's Scheme
  vector<double> testvec (3);
  double y = 2;
  testvec[0] = 5;
  testvec[1] = 2;
  testvec[2] = 3;
  cout << "5+2x+3x^2 for x=2     " << hornerScheme(testvec, 2, y) << endl;

  getch();
  return 0;
}