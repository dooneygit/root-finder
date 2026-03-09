#include <iostream>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <algorithm>

// Function declarations
int main();
unsigned int roots( double coeffs[], unsigned int degree );
double  horner( double x, double coeffs[], unsigned int degree );
double dhorner( double x, double coeffs[], unsigned int degree );
double newton( double coeffs[], unsigned int degree );
double divide( double r, double poly[], unsigned int degree );
void print( double coeffs[], unsigned int complex_degree, unsigned int degree );

// Function definitions
int main() {
  std::cout.precision( 16 );

  // The output should be:
  //      Remaining polynomial:
  //          3+1x+1x^2
  //      Real roots:
  //          -2
  //          +1
  //          +3
  double p[6]{ 18, -9, -5, -4, -1, 1 };
  unsigned int degree{ roots( p, 5 ) };
  print( p, degree, 5 );

    //////////////////////////////////////////////
   // The sample input and output in Section 9 //
  //////////////////////////////////////////////
  double q0a[1]{ -14.1 };
  unsigned int c0a{ roots( q0a, 0 ) };
  print( q0a, c0a, 0 );
  // Expected output:
  // Remaining polynomial:
  //     -14.1
  // Real roots:

  double q1a[2]{ -5.0, 2.0 };
  unsigned int c1a{ roots( q1a, 1 ) };
  print( q1a, c1a, 1 );
  // Expected output:
  // Remaining polynomial:
  //     2
  // Real roots:
  //     +2.5

  double q2a[3]{ 5.0, 2.0, 2.0 };
  unsigned int c2a{ roots( q2a, 2 ) };
  print( q2a, c2a, 2 );
  // Expected output (no real roots):
  // Remaining polynomial:
  //     5+2x+2x^2
  // Real roots:

  double q2b[3]{ 5.0, 2.0, -2.0 };
  unsigned int c2b{ roots( q2b, 2 ) };
  print( q2b, c2b, 2 );
  // Expected output (two real roots):
  // Remaining polynomial:
  //     -2
  // Real roots:
  //     -1.1583123951777
  //     +2.1583123951777

  double q3a[4]{ -100.8, 8.4, -4.2, 4.2 };
  unsigned int c3a{ roots( q3a, 3 ) };
  print( q3a, c3a, 3 );
  // Expected output (one real and two complex roots):
  // Remaining polynomial:
  //     33.6+8.4x+4.2x^2
  // Real roots:
  //     +3

  return 0;
}

// You must implement the 'roots(...)' function with
// the parameters and return value as specified here.
unsigned int roots( double coeffs[], unsigned int degree ) {
  if(degree == 0) return 0;
  assert(coeffs[degree] != 0.0);

  unsigned int n = degree;
  
  while(n > 0) {
    double root = newton(coeffs, n);

    if(std::isnan(root)) {
      break;
    }
    else {
      divide(root, coeffs, n);
      coeffs[n] = root;
      n--;
    }
  }
  std::sort(&coeffs[n+1], &coeffs[degree + 1]);

  return n;
}

// You must implement Newton's method on polynomials.
// You are welcome to change the signature.
double newton( double coeffs[], unsigned int degree ) {
  if(degree == 0) return NAN;
  assert(coeffs[degree] != 0.0);

  double x0 = 211.59169;
  double eps = 1e-8;

  for(unsigned int k{0}; k < 100; k++) {
    double fx = horner(x0, coeffs, degree);
    double dfx = dhorner(x0, coeffs, degree);

    double x1 = x0 - fx/dfx;

    if(!std::isfinite(x1)) {
      return NAN;
    }

    if((std::abs(x1 - x0) < eps) && (std::abs(horner(x1, coeffs, degree)) < eps)) {
      return x1;
    }

    x0 = x1;
  }

  return NAN;
}

// The polynomial division algorithm is implemented for you. 
double divide( double r, double poly[], unsigned int degree ) {
  double s{ poly[degree] };

  for ( unsigned int k{ degree }; k > 0; --k ) {
    double t{ poly[k - 1] + r*s };
    poly[k - 1] = s;
    s = t;
  }

  return s;
}

// Horner's rule is implemented here for you
double horner( double x, double coeffs[], unsigned int degree ) {
  double result{ coeffs[degree] };

  for ( unsigned int k{ degree }; k > 0; ) {
    --k;
    result = result*x + coeffs[k];
  }

  return result;
}

// You need to implement this
double dhorner( double x, double coeffs[], unsigned int degree ) {
  double result{degree * coeffs[degree]};

  for(unsigned int k{degree}; k > 1;) {
    k--;
    result = result * x + k * coeffs[k];
  }

  return result;
}

void print(
    double coeffs[],
    unsigned int complex_degree,
    unsigned int degree
) {
    // Store the current value of precision
    // and set the precision to 16 digits
    std::streamsize old_precision{ std::cout.precision( 16 ) };
  
    std::cout << "Remaining polynomial:" << std::endl;
    std::cout << "\t" << coeffs[0];

    // Show a "+" sign in front of all positive floating-point numbers
    std::cout << std::showpos;
  
    if ( complex_degree >= 1 ) {
        std::cout << coeffs[1] << "x";
    }

    for ( unsigned int k{ 2 }; k <= complex_degree; ++k ) {
        std::cout << coeffs[k] << "x^" << k;
    }
        
    std::cout << std::endl << "Real roots:" << std::endl;

    for ( unsigned int k{ complex_degree + 1 }; k <= degree; ++k ) {
        std::cout << "\t" << coeffs[k] << std::endl;
    }

    // Stop showing the leading "+"
    std::cout << std::noshowpos << std::endl;

    // Restore the original value of precision
    std::cout.precision( old_precision );
}