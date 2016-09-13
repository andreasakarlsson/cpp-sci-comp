#include <iostream>
#include <cmath>

using namespace std;

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

int main() {
  double factorial(int n);
  int N = 100;
  double x = 3.1416;
  cout << sinTaylor(N, x) << endl;
  cout << cosTaylor(N, x) << endl;
  cout << abs(sin(x) - sinTaylor(N, x)) << endl;
  cout << abs(cos(x) - cosTaylor(N, x)) << endl;
  return 0;
}
