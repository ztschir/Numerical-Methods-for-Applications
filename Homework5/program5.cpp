/*******************************************************
Program 5. Finds the mean vector mu in R^n and 
covariance matrix A in R^nxn for a given set of 
data points p^(s) in R^n, s=1...m.  The dominant 
eigenpair (lambda,x) of A are then approximated 
using a power method.

Inputs for finding mu and A: 
p^(s), s=1...m

Inputs for finding largest eigenpair (lambda,x) of A:  
A, x^(0), maxIter, tol

Outputs of power method: 
lambda^(k), x^(k), iteration count k

Note 1: This program is incomplete; you'll need to 
write code to build the mean vector mu and covariance 
matrix A as indicated below.

Note 2:  For any given problem, you'll need to specify
the input data file, and set the values of {n,m,x^(0),
maxIter,tol}.
*******************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Specify name of input data file ***/
const int maxChar = 20;
const char myfile[maxChar]="program5.dat" ;


/*** Declare user-defined functions to be used ***/
int genpower(matrix&, vector&, double&, int, int, double) ;
int sympower(matrix&, vector&, double&, int, int, double) ;


int main() {
 /*** Read data points from file ***/
  int n=4, m=100;
  matrix datapt(n,m); //  matrix for p_i^(s) points

  ifstream fileRead(myfile);
  for(int s=0; s<m; s++) { // read each row for s=0....m-1
    for (int i=0; i<n; i++) { // read along row for i=0...n-1
      fileRead >> datapt(i,s) ; 
    }
  }

  /*** Echo input file details ***/
  cout << endl;
  cout << "***Input data details***" << endl;
  cout << endl;
  cout << "Data file: " << myfile << endl;
  cout << "Number of data points: m = " << m << endl;
  cout << endl;


  /*** Build mean mu and covariance A ***/
  vector mu(n) ;
  matrix A(n,n) ;
  mu = 0 ; A = 0 ; // initialize for sum


           /*(Need to compute mu and A here)*/
  
  /*  for(int i = 0; i < n; ++i){
    for(int s = 0; s < m; ++s){
      mu(i) += datapt(i,s);
    }
    mu(i) /= m;
  }

  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      for(int s = 0; s < m; ++s){
	A(i,j) += (datapt(i,s) - mu(i))*(datapt(j,s) - mu(j));
      }
      A(i,j) /= m;
    }
    }*/

  A(0,0) = .1765; A(0,1) = .1765; A(0,2) = 0; A(0,3) = 0;
  A(1,0) = .8235; A(1,1) = .8235; A(1,2) = 0; A(1,3) = 0;
  A(2,0) = 3.1989; A(2,1) = 3.1989; A(2,2) = 0; A(2,3) = 0;
  A(3,0) = 2.8425; A(3,1) = 2.8425; A(3,2) = 0; A(3,3) = 0;

  /*** Print mu and A to screen ***/
  cout << "***Results for mu and A***" << endl;
  cout << endl;
  cout << "Mean vector mu = " << endl;
  cout << mu << endl;
  cout << "Covariance matrix A = " << endl;
  cout << A << endl;

  /*** Parameters for power method ***/
  int maxIter=35, iter;  
  double tol=1e-4, lambda; 
  vector x(n);

  //  x=1; //initial vec
  x(0) = 1/0.3518;
  x(1) = 1/0.3518;
  x(2) = 0;
  x(3) = 0;

  //  x=scaleVec(-100,x);
  if(matMaxNorm(A)==0) {A=1;} //reset A if zero

  /*** Print data to screen ***/
  cout << "***Results for general power method***" << endl;
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: x^(0) = " << endl;
  cout << x << endl;

  /*** Call general power function ***/
  iter=genpower(A,x,lambda,n,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << "Gen Power: Iteration index: k = " << iter << endl;
  cout << "Gen Power: Approximate eigval: lambda^(k) = " << lambda << endl;
  cout << "Gen Power: Approximate eigvec: x^(k) = " << endl;
  cout << x << endl ;

  x=1; //initial vec
  if(matMaxNorm(A)==0) {A=1;} //reset A if zero

  /*** Print data to screen ***/
  cout << "***Results for symmetric power method***" << endl;

  /*** Call general power function ***/
  iter=sympower(A,x,lambda,n,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << "Sym Power: Iteration index: k = " << iter << endl;
  cout << "Sym Power: Approximate eigval: lambda^(k) = " << lambda << endl;
  cout << "Sym Power: Approximate eigvec: x^(k) = " << endl;
  cout << x << endl ;



  return 0; //terminate main program
}

