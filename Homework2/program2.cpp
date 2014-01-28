/*******************************************************
Program 2.  Solves Ax=b using Conjugate Gradient method
with pre-conditioning matrix Cinv.   A is assumed to be
symmetric and positive-definite.  

Inputs: A, b, x^(0), Cinv, maxIter, tol
Outputs: x^(k), iteration count k

Here's how to get started:

1) Copy program2.cpp (this file), conjgrad.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program2 matrix.cpp conjgrad.cpp program2.cpp"
at the Linux prompt.

3) Type "program2" to run the program.

4) For any given problem, you'll need to set the 
values of the variables {n,A,b,x,Cinv,maxIter,tol} 
below, and then compile and run as described above.

5) If you add another function for making comparisons, 
e.g. jacobi, gauss_seidel or sor, remember to include 
a header for the function in the main program below, 
and to include the filename when compiling as in 
Step 2 above.  Also, remember to re-initialize x before 
running another method.
*******************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Declare user-defined functions to be used ***/
int conjgrad(matrix&, matrix&, vector&, vector&, int, double);
int jacobi(matrix&, vector&, vector&, int, double);
int gauss_seidel(matrix&, vector&, vector&, int, double);
int sor(matrix&, matrix&, vector&, vector&, int, double);
int sor(matrix&, vector&, vector&, double, int, double);

/*** Main program ***/
int main() {

  /*** Define and input problem data ***/
  int n=100, maxIter=1000, iter;
  double tol=1e-5, w=1.2, r_norm;
  matrix Cinv(n,n), A(n,n);
  vector x(n), b(n);

  // Build A matrix
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      if(j == i - 1) A(i,j) = -1.0;
      else if(j == i) A(i,j) = 2.0 + ((i + 1)/10);
      else if(j == i + 1) A(i,j) = -1.0;
      else A(i,j) = 0.0;
    }
  }

  // Build b vector
  for(int i = 0; i < n; ++i){
    b(i) = 1.0 + ((i + 1)/20);
  }

  x=0; //initialize  x

  /*** Construct pre-conditioning matrix ***/
  Cinv=0;
  for(int i=0; i<n; i++) {
    if(A(i,i) <= 0) {
      cerr << "CG: bad data -- exiting" << endl;
      exit(EXIT_FAILURE);
    } 
    else { 
      Cinv(i,i) = 1/sqrt(A(i,i)); 
      // Cinv(i,i) = 1; //Use this for no pre-cond
    }
  }

  /*** Call conjugate gradient function ***/
  iter=conjgrad(Cinv,A,b,x,maxIter,tol);

  /*** Print results to screen ***/
  cout << endl; 
  cout << "Conjugate Gradient: Iteration index: k = " << iter << endl;
  cout << "Conjugate Gradient: Approximate solution: x^(k) = " << endl;
  cout << x << endl;



  /*** Reset x ***/
  x = 0;

  /*** Call Jacobi function ***/
  iter=jacobi(A,b,x,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl; 
  cout << "Jacobi: Iteration index: k = " << iter << endl;
  cout << "Jacobi: Approximate solution: x^(k) = " << endl;
  cout << x << endl; 
  
  

  /*** Reset x ***/
  x = 0;
  
  /*** Call Gauss Seidel function ***/
  iter=gauss_seidel(A,b,x,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl;
  cout << "Gauss Seidel: Iteration index: k = " << iter << endl;
  cout << "Gauss Seidel: Approximate solution: x^(k) = " << endl;
  cout << x << endl;


  /*** Reset x ***/
  x = 0;
  
  /*** Call Successive Over-relaxation function ***/
  iter=sor(A,b,x,w,maxIter,tol);  

  /*** Print results to screen ***/
  cout << endl;
  cout << "Successive Over-relaxation: Iteration index: k = " << iter << endl;
  cout << "Successive Over-relaxation: Approximate solution: x^(k) = " << endl;
  cout << x << endl;

  return 0;
}

