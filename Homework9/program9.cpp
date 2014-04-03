/*******************************************************
Program 9. Uses the Newton-RK4 shooting method to find 
an approximate solution of a two-point BVP of the form

              y'' = f(x,y,y'),  a<=x<=b  
              y(a) = alpha,  y(b) = beta

Inputs:  
  feval         Function to evaluate f(x,y,y')
  a,b           Interval params
  alpha,beta    Boundary value params
  t             Initial guess for y'(a)
  N             Number of grid points for RK4
  maxIter,tol   Control params for Newton

Outputs: 
  iter    Number of Newton steps taken
  t       Approx value of y'(a) 
  x       Grid point vector: x(j)=a+jh, j=0...N
  y       Approx soln vector: y(j)=soln at x(j), j=0...N

Note 1: For any given problem, the function feval must 
be changed.  This function computes f for any given x, y 
and y'.  (Below we use the symbol yp in place of y').

Note 2: For any given problem the parameters a, b, alpha,
beta, t, N, maxIter and tol must also be specified.

Note 3: To compile this program use the command (all on 
one line)

  c++ -o program9  matrix.cpp shootRK4.cpp program9.cpp

*******************************************************/
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare external function ***/
int shootRK4(int&, double&, double&, double&, double&, 
                double&, int&, double&, vector&, vector&) ;


/*** Define f(x,y,yp) function ***/
void feval(const double& x, const double& y, 
                              const double& yp, double& f){
  double pi=4.0*atan(1.0) ; //the number pi
  double p, q, u, v, w;
  double xp1 = pow(x+1,2);
  double xm1 = pow(x-1,2);
  double yp1 = pow(y+1,2);
  double ym1 = pow(y-1,2);

  p = -3*(x+1)*exp(-xp1-yp1) - 2*(x-1)*exp(-xm1-ym1);
  q = -3*(y+1)*exp(-xp1-yp1) - 2*(y-1)*exp(-xm1-ym1);
  u = 2*(exp(-xm1-ym1) - 2*xm1*exp(-xm1-ym1)) -
      3*(exp(-xp1-yp1) - 2*xp1*exp(-xp1-yp1));
  v =  6*(x+1)*(y+1)*exp(-xp1-yp1) - 4*(x-1)*(y-1)*exp(-xm1-ym1);
  w = 2*(exp(-xm1-ym1) - 2*ym1*exp(-xm1-ym1)) -
      3*(exp(-xp1-yp1) - 2*yp1*exp(-xp1-yp1));

  f = ((p*yp-q)*(u+2*v*yp+w*pow(yp,2))) / (1+pow(p,2)+pow(q,2));
}

double path(vector& x, vector& y) {
  double sum = 0;
  double z, zm1;
  for(int i=1; i<dim(x); i++){
    z = 1.5*exp(-pow(x(i)+1,2)-pow(y(i)+1,2)) + 
      exp(-pow(x(i)-1,2)-pow(y(i)-1,2));
    zm1 = 1.5*exp(-pow(x(i-1)+1,2)-pow(y(i-1)+1,2)) + 
      exp(-pow(x(i-1)-1,2)-pow(y(i-1)-1,2));
    sum += sqrt((x(i)-x(i-1))*(x(i)-x(i-1)) + 
		(y(i)-y(i-1))*(y(i)-y(i-1)) + 
		pow(z-zm1,2));
  }
  return sum;
}

int main() {
  /*** Define problem parameters ***/
  double tol=1e-6;
  int N=40, maxIter=10, iter;  
  double a=-3, b=3, alpha=-2, beta=2, t; 
  vector x(N+1), y(N+1); 
  t=0.45; // initial guess of slope 

  /*** Call Newton-RK4 method ***/
  cout << setprecision(8);
  iter=shootRK4(N,a,b,alpha,beta,t,maxIter,tol,x,y);
  double pathlen = path(x,y);
  /*** Print results to screen ***/
  cout << setprecision(6);
  cout << "Number of RK4 grid points: " << N << endl;
  cout << "Number of Newton iterations: " << iter << endl;
  cout << "Approx solution: t = " << t << endl;
  cout << "Path length: " << pathlen << endl;
  cout << "Approx solution: x_j, y_j =  " << endl;
  for(int j=0; j<N+1; j++){
    cout << "x =" << setw(6) << x(j);
    cout << "   ";
    cout << "y =" << setw(10) << y(j);
    cout << "   "; 
    cout << endl;
  }

  return 0; //terminate main program
}

