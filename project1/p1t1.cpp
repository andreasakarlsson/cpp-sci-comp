#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

// taylor expansion for sin with Horner's scheme
double sinTaylor(int pDegree, double x)
{
  double hSchemeRes = 0.0;
  for (int i = 2 * pDegree + 1; i >= 0; i--){                   // 2*n + 1 for the Taylor expansion
    hSchemeRes =  ((((1 + i) % 2 == 0) ? 2 - (i % 4)  : 0) +	// 0 1 0 -1 0...
		   x * hSchemeRes) /				// Horner's recursion
      ((i == 0) ? 1 : i);					// single factorial factor
  }
  return hSchemeRes;
}

// taylor expansion for cos with Horner's scheme
double cosTaylor(int pDegree, double x)
{
  double hSchemeRes = 0.0;
  for (int i = 2 * pDegree; i >= 0; i--){               // 2*n for the Taylor expansion
    hSchemeRes =  ((( i % 2 == 0) ? 1 - (i % 4)  : 0) + // 1 0 -1 0 1...
		   x * hSchemeRes) /			// Horner's recursion
      ((i == 0) ? 1 : i);				// single factorial factor
  }
  return hSchemeRes;
}

int main()
{
  vector<double> Xvec = {-1, 1, 2, 3, 5, 10};
  vector<unsigned int> Nvec = {1, 10, 100};

  // Learning ranges in C++11
  for (auto N: Nvec) {
    cout << "Setting N = " << N << endl;
    for (auto x: Xvec) {
      cout << "  Setting x = " << x << endl;
      cout << "    sin using Horner: error = "
	   << setw(11)
	   << abs(sin(x) - sinTaylor(N, x))
	   << "    N+1-term = "
	   << abs(sinTaylor(N + 1, x)) << endl;
      cout << "    cos using Horner: error = "
	   << setw(11)
      	   << abs(cos(x) - cosTaylor(N, x))
	   << "    N+1-term = "
	   << abs(cosTaylor(N + 1, x)) << endl;
    }
  }
  return 0;
}
