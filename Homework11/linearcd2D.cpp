/************************************************************
Function to implement the central-difference method to 
find an approximate solution of an elliptic BVP in a 
rectangular domain of the form 

 P uxx + Q uyy + p ux + q uy + r u = f,  a<=x<=b, c<=y<=d  
 u(a,y) = ga(y),     u(b,y) = gb(y),          c<=y<=d  
 u(x,c) = gc(x),     u(x,d) = gd(x),          a<=x<=b  


Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,f
  BCeval        Function to evaluate ga,gb,gc,gd
  a,b,c,d       Domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln: u(i,j)=soln at x(i),y(j), 
                i=0...N+1, j=0...M+1


Note 1: This function is incomplete; you'll need to code 
the A-matrix and G-vector as indicated below.  The
central-difference equation system is A uint = G.

Note 2: The functions PDEeval and BCeval are assumed 
to be defined externally (e.g. by calling program).

Note 3: The unknowns in the central-difference equation 
system are uint(l) = u(i,j) where l=i+N(M-j), i=1...N,
j=1...M which gives l=1...NM.  The boundary values for 
u(i,j) are obtained from the BCeval function.

Note 4: Gauss elimination with partial pivoting is used to 
solve the system.
************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Ext fun: gauss_elim(A,G,uint) ***/
int gauss_elim(matrix&, vector&, vector&) ;


/*** Ext fun: PDEeval(x,y,P,Q,p,q,r,f) ***/
void PDEeval(const double&, const double&, 
                     double&, double&, double&, 
                          double&, double&, double&) ;

/*** Ext fun: BCeval(x,y,gLeft,gRight,gBottom,gTop) ***/
void BCeval(const double&, const double&, 
                     double&, double&, double&, double&) ;


/*** Aux fun: build A and G for cent-diff system ***/
void AGeval(int N, int M, 
            double a, double b, double c, double d, 
            vector& x, vector& y, matrix& A, vector& G){

  int l, ll ;
  double P, Q, p, q, r, f ;
  double gLeft, gRight, gBottom, gTop ;
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ; //grid parameters
  double al, bl, cl, dl, el ; //A-matrix parameters
  double BCtopl, BCrightl, BCleftl, BCbottoml ; //BC parameters

  A=0 ; G=0 ;
  for(int i=1; i<=N; i++){ //actual i-index on grid
    for(int j=1; j<=M; j++){ //actual j-index on grid
      ll = i + N*(M-j) ; //actual l-label on grid (1...NM)
      l = ll-1 ; //C++ array index (0...NM-1)

      /*** Build the entries of the A-matrix and G-vector here.  
           A few entries of A are shown as an example.  We use 
           the more descriptive notation gLeft,gRight,gBottom,
           gTop in place of ga,gb,gc,gd. ***/

      PDEeval(x(i),y(j),P,Q,p,q,r,f) ; //PDE coeffs at i,j

      al = (dy/dx)*P - (dy/2.0)*p ;
      bl = dx*dy*r - (2.0*dy/dx)*P - (2.0*dx/dy)*Q ;
      cl = (dy/dx)*P + (dy/2.0)*p ;
      dl = (dy/dx)*Q + (dy/2)*q;
      el = (dy/dx)*Q - (dy/2)*q;

      if( i>1 ){ A(l,l-1) = al ; }
      A(l,l) = bl ;
      if( i<N ){ A(l,l+1) = cl ; }
      if(j<M){ A(l,l-N) = dl;}
      if(j>1){ A(l,l+N) = el;}

      G(l) = dy*dx*f;
      BCeval(x(i),y(M+1),gLeft,gRight,gBottom,gTop) ; //BC vals at i,M+1 
      if( j==M ){ G(l) -= dl*gTop ; }
      BCeval(x(N+1),y(j),gLeft,gRight,gBottom,gTop) ; //BC vals at N+1,j 
      if( i==N ){ G(l) -= cl*gRight ; }
      BCeval(x(i),y(0),gLeft,gRight,gBottom,gTop) ; 
      if( j==1 ){ G(l) -= el*gBottom ; }
      BCeval(x(0),y(j),gLeft,gRight,gBottom,gTop) ;  
      if( i==1 ){ G(l) -= al*gLeft ; }
     }
  }
}


/*** Main function: central-difference method ***/
int linearcd2D(int N, int M, 
               double a, double b, double c, double d, 
                            vector& x, vector& y, matrix& u){

  int l, ll, success_flag=0 ;
  matrix A(N*M,N*M) ;
  vector uint(N*M), G(N*M) ;
  double gLeft, gRight, gBottom, gTop ;
 
  /*** Build A-matrix and G-vector ***/
  AGeval(N,M,a,b,c,d,x,y,A,G) ;
  cout << "Linear-CD-2D: arrays assembled" << endl ;

  /*** Solve linear system ***/
  cout << "Linear-CD-2D: solving........." << endl ;
  gauss_elim(A,G,uint) ;
  cout << "Linear-CD-2D: equations solved" << endl ;

  /*** Assemble total solution: interior and boundary ***/
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    for(int j=0; j<=M+1; j++){ //actual j-index on grid
      BCeval(x(i),y(j),gLeft,gRight,gBottom,gTop) ;
      if( j==M+1 ){ u(i,j) = gTop ; }
      if( i==N+1 ){ u(i,j) = gRight ; }
      if( i==0 ){ u(i,j) = gLeft ; }
      if( j==0 ){ u(i,j) = gBottom ; }
      if((i>0)&&(i<N+1)&&(j>0)&&(j<M+1)){
         ll = i + N*(M-j) ; //actual l-label on grid (1...NM)
         l = ll-1 ;         //C++ array index (0...NM-1)
         u(i,j) = uint(l) ; //interior soln
      }
    }
  }
  cout << "Linear-CD-2D: solution assembled" << endl ;

  return success_flag=1 ;
}

