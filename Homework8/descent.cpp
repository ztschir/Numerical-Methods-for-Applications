/*******************************************************
Function to implement the steepest descent method to 
minimize a function g(x).

Inputs:  
  geval    Function to evaluate g(x), scalar
  dgeval   Function to evaluate dg/dx(x), n vector
  xk       Initial guess for x, n vector  
  maxIter  Maximum number of descent iterations
  tol      Tolerance parameter for ||dg/dx||

  ainit    Initial guess for alpha (step size)
  maxSrch  Maximum number of alpha search steps


Outputs: 
  xk       Approx local minimum point of g(x)
  k        Number of descent steps taken


Note 1: This function is incomplete; you'll need to code 
the quadratic interpolation step as indicated below.

Note 2: The functions geval and dgeval are assumed 
to be defined externally (e.g. by calling program).
*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;

/*** Declare user-defined functions to be used ***/
void geval(vector&, double&) ;  //defined in program8.cpp
void dgeval(vector&, vector&) ; //defined in program8.cpp


int descent(vector& xk, int n, int maxIter, 
                    double tol, int maxSrch, double ainit){
  int k=0, kSrch=0 ;  
  double gk, dgkNorm ;
  double a3, ga3, a0, ga0;
  double y1, y2, y3, k1, k2, k3;
  vector dgk(n), xa3(n), xa0(n) ;
  vector kk1(n), kk2(n), kk3(n);

  geval(xk,gk) ; //evaluate g(xk)
  dgeval(xk,dgk) ;  //evaluate dg(xk)
  dgkNorm = vecL2Norm(dgk) ; 
  
  while(dgkNorm>=tol && k<maxIter) {
    a3 = ainit ;
    xa3 = xk - scaleVec(a3/dgkNorm,dgk) ;
    geval(xa3,ga3) ;

    kSrch = 0 ;
    while( ga3>=gk && kSrch<maxSrch ){ 
      a3 = 0.5*a3 ;
      xa3 = xk - scaleVec(a3/dgkNorm,dgk) ;
      geval(xa3,ga3) ;
      kSrch++ ;
    }
    if(kSrch>=maxSrch){
      cout << "Descent: max iter exceeded in alpha search" << endl;
    }

    a0 = a3 ; ga0 = ga3 ; xa0 = xa3 ;
    // Reset to zero
    y1 = 0; y2 = 0; y3 = 0; k1 = 0; k2 = 0; k3 = 0;
    
    // set k's
    k3 = a3; k2 = a3/2; k1 = 0;
    kk1 = scaleVec(k3, xa3); 
    kk2 = scaleVec(k2, xa3); 
    kk3 = scaleVec(k1, xa3); 

    geval(kk1, y3);
    geval(kk2, y2);
    geval(kk3, y1);
    
    a0 =   (pow(k1,2) * (y2-y3) + 
	    pow(k2,2) * (y3-y1) + 
	    pow(k3,2) * (y1-y2)) /
           (k1 * (y2-y3) + 
	    k2 * (y3-y1) + 
	    k3 * (y1-y2)) / 2;

    xa0 = scaleVec(a0, xa3);
    geval(xa0, ga0);

    if(a0>=0 && a0<=a3 && ga0<ga3){
      xk = xa0 ;
      gk = ga0 ;
    }
    else {
      xk = xa3 ;
      gk = ga3 ;
    }

    dgeval(xk,dgk) ;
    dgkNorm = vecL2Norm(dgk) ; 
    k++ ;
    cout << "Descent: |Gradg_k|_2 = " << dgkNorm << endl ;
  }


  if(dgkNorm < tol) {
    cout << "Descent: solution converged" << endl ;
  } 
  else {
    cout << "Descent: max iterations exceeded" << endl;
  }

  return k;
}

