/*******************************************************
Program 3.  Finds least-squares polynomial P_n(x) of 
specified degree n for a given set of data (x_k,y_k), 
k = 1...m.  It is assumed that x_k are distinct and 
0 <= n <= m-1.

Inputs: data (x_k,y_k), k=1...m  
        degree n (entered at runtime)

Outputs: coefficients a_0...a_n of poly P_n
         error E_n associated with poly P_n

1) Copy program3.cpp (this file), gauss_elim.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program3 matrix.cpp gauss_elim.cpp program3.cpp"
at the Linux prompt.

3) Type "program3" to run the program.

NOTE: This program is incomplete; you'll need to 
construct the A,b matrices as indicated below.

NOTE: For any given problem, you'll need to set the 
values of m and (x_k,y_k), k=1...m, and then compile
and run.
*******************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare user-defined functions to be used ***/
double poly(vector&, int, double); // defined below
int gauss_elim(matrix&, vector&, vector&);


int main() {
  /*** Input (x,y) data ***/
  int m=5;
  vector x(m), y(m);

  x(0) = 1.1 ;  y(0) = 3.8 ;
  x(1) = 2.0 ;  y(1) = 2.9 ;
  x(2) = 4.3 ;  y(2) = 4.2 ;
  x(3) = 5.2 ;  y(3) = 6.5 ;
  x(4) = 7.9 ;  y(4) = 8.8 ;


  /*** Build least-squares arrays A,b ***/
  int n; 
  cout << "Enter polynomial degree: " << flush;
  cin >> n;
  if( n<0 || n>m-1 ){
    cerr << "Least-squares: bad degree -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  matrix A(n+1,n+1);
  vector b(n+1), alpha(n+1);
  A=0 ; b=0 ; // initialize for sum



           /*(Need to build A,b here)*/



  /*** Print least-squares arrays to screen ***/
  cout << endl;
  cout << "Least-squares degree: n = " << endl;
  cout << n << endl << endl;
  cout << "Least-squares array: A = " << endl;
  cout << A << endl;
  cout << "Least-squares array: b = " << endl;
  cout << b << endl;

  /*** Solve normal equations and print results ***/
  int success_flag;  
  success_flag=gauss_elim(A,b,alpha);
  cout << "Least-squares solution: alpha^(k) = " << endl;
  cout << alpha << endl;

  /*** Compute fitting error and print results ***/
  double E=0;
  for(int k=0; k<m; k++) {
     E = E + pow(y(k)-poly(alpha,n,x(k)),2);
  }
  cout << "Least-squares fitting error: E = " << endl;
  cout << E << endl;

  return 0; //terminate main program
}


/*** Auxiliary function to evalute poly P_n ***/
double poly(vector& alpha, int n, double x) {
  double P=0;
  for(int i=0; i<n+1; i++) {
    P = P + alpha(i)*pow(x,i); 
  }
  return P;
}

