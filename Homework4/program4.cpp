/*******************************************************
Program 4.  Finds the least-squares trig poly S_n of 
specified degree n for a given set of data (x_j,y_j), 
j=0...2m-1.  It is assumed 2 <= n < m and that x_j are 
the std nodes in [-pi,pi], i.e. x_j = -pi + j*pi/m.

Inputs:  data (x_j,y_j),  j=0...2m-1
         degree n (entered at runtime)

Outputs: coeffs a_0...a_n, b_1...b_{n-1} of S_n
         error E_n associated with S_n

Note 1: This program is incomplete; you'll need to 
code the a,b-coefficients and fitting error E as 
indicated below.

Note 2: For any given problem, you'll need to set the 
value of m and specify the input file with the xy-data,
and then compile and run.
*******************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Specify name of input data file ***/
const int maxChar = 20;
const char myfile[maxChar]="program4.dat" ;
int main() {

  /*** Input/read (x_j,y_j) data ***/
  int m=50;
  double pi=4.0*atan(1.0); // the number pi
  vector x(2*m), y(2*m);
  ifstream fileRead(myfile);
  for(int j=0; j<2*m; j++) {
    fileRead >> x(j) >> y(j); // read xy-pair 
  }

  /*** Echo input file details ***/
  cout << endl;
  cout << "xy-data read from file: " << myfile << endl;
  cout << endl;
  cout << "Number of data points 2m: " << 2*m << endl;
  cout << endl;
  /*** Compute a,b-coefficients ***/
  int n; 
  cout << "Enter trig poly degree n: " << flush;
  cin >> n;
  if( n<2 || n>= m){
    cerr << "Trig LS: bad degree -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  vector a(n+1), b(n), s(2*m);
  a = 0 ; b = 0 ; s = 0; // initialize for sum

  /*** Computed a vector ***/
  for(int k = 0; k <= n; ++k){
    for(int j = 0; j <= (2*m - 1); ++j){
      a(k) += y(j)*cos(k*x(j));
    }
    a(k) /= m;
  }

  /*** Computed b vector ***/
  for(int k = 1; k <= n-1; ++k){
    for(int j = 0; j <= (2*m - 1); ++j){
      b(k) += y(j)*sin(k*x(j));
    }
    b(k) /= m;
  }

  /*** Computed sn(x) vector for all x ***/
  for(int i = 0; i < 2*m; ++i){
    s(i) = a(0)/2 + a(n)*cos(n*x(i));
    for(int k = 1; k <= n-1; ++k){
      s(i) += a(k)*cos(k*x(i)) + b(k)*sin(k*x(i));
    }
  }

  /*** Print least-squares coeffs to screen ***/
  cout << endl;
  cout << "Trig LS results for n = " << n << ", m = " << m << endl;
  cout << endl;
  cout << "Cosine coeffs: a_0 ... a_n = " << endl;
  cout << a << endl;
  cout << "Sine coeffs: b_1 ... b_{n-1} = " << endl;
  /***  Needed to change this to account for the zero b(0) value ***/
  for(int i = 1; i < dim(b); ++i){
    cout << b(i) << "\n";
  }

  /*** Compute least-squares fitting error ***/
  double E = 0;
  for(int j = 0; j < 2*m; ++j){
    E += pow(y(j) - s(j),2);
  }

  cout << endl;
  cout << "Fitting error: E = " << endl;
  cout << E << endl;
    
  return 0; //terminate main program
}

