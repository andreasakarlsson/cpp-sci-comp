#include <iostream>
#include <cmath>
#include <vector>

#ifdef _WIN32
#include <conio.h> // loads getch();
#endif

using namespace std;
double myexp(double x, double tol);




double myexp(double x, double tol){
	int n = 0;

	// The sum of the terms following term n is less than term n.
	// => If term n < tol => Error < tol
	while (true){
		n += 1;
		double res = 1;
		for(int j = n; j > 0; j--){ res = x / j * res;}
			//cout << "n = " << n << " res = "<< res << endl;
		if (abs(res) < tol){ break; }
	}


	//Calculate exponential with Horner's Scheme.
	double result = 0.0;
	for(int i = n; i >= 0; i--){
		result = (1 + x*result) / ((i == 0) ? 1 : i);
	}
	return result;
}


int main() {


  vector<double> x_vec = {-1, 1, 3, 5, 10, 50};
  vector<double> tol_vec = {0.01, 0.001, 1e-8, 1e-10};

  // Learning ranges in C++11
  for (auto tolerance: tol_vec) {
    cout << "tolerance tol = " << tolerance << endl;
    for (auto x: x_vec) {
      cout << "  Setting x = " << x << endl;
      cout << "    myexp and cmath exp, error:     "
	   << abs(myexp(x,tolerance)-exp(x)) << "\n" << endl;

    }
  }

#ifdef _WIN32
	cout << "\n Press any key to quit..." << endl;
	getch();
#endif

	return 0;
}
