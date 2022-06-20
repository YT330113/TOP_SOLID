/*************************************************************************************************/
/***************************************拓扑变量求解及过滤**************************************/
/*************************************************************************************************/
#include "MMASolver.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

/********************************MMA求解器拓扑变量求解***********************************/


// Math helpers
double Min(double d1, double d2);
double Max(double d1, double d2);
int Min(int d1, int d2);
int Max(int d1, int d2);
double Abs(double d1);


	std::cout<<"\nRun MMA algorithm\n"<<std::endl;


	int n = n_dof;
	int m = m_con;

	double Xmin = 0.0;
	double Xmax = 1.0;
	double movlim = 0.2;

	double *xmma = new double[n];
	double *xold = new double[n];
	double *f = new double[1];
	double *df = new double[n];
	double *g = new double[m];
	double *dg = new double[n*m];
	double *xmin = new double[n];
	double *xmax = new double[n];
  

  	// Initialize MMA
	MMASolver *mma = new MMASolver(n,m);
	
	// 这里是灵明度ObjSens(af,bf,x,n,m,f,df,g,dg);

	// Set outer move limits
	for (int i=0;i<n;i++) {
		xmax[i] = Min(Xmax, xmma[i] + movlim);
		xmin[i] = Max(Xmin, xmma[i] - movlim);
	}

	// Call the update method
	mma->Update(xmma,df,g,dg,xmin,xmax);

	// Compute infnorm on design change
	double ch = 0.0;
	for (int i=0;i<n;i++) {
		ch = Max(ch,Abs(xmma[i]-xold[i]));
		xold[i] = xmma[i];
		x[i]=xmma[i];
	}

	// Print to screen
	printf("设计变量的最大变化量:%f \n", ch);	


//-------------------------------------------------------------------------------------------------------------------------//


/********************************PDE过滤及Heviside过滤***********************************/

/***********************************PDE filter**********************************************/

solve(fvm::laplacian(xp) - fvm::Sp(b, xp) + x * b);
//-------------------------------------------------------------------------------------------------------------------------//

//***********************************Heaviside_filter***************************************//
del = Foam::min(Hdel * opt, 100);
eta0 = 0.0001;
eta1 = 1.0;
y0 = diff(xp, mesh.V(), del, eta0, n);
reduce(y0, sumOp<scalar>());

do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}
//-------------------------------------------------------------------------------------------------------------------------//



double Min(double d1, double d2) {
	return d1<d2 ? d1 : d2;
}

double Max(double d1, double d2) {
	return d1>d2 ? d1 : d2;
}

int Min(int d1, int d2) {
	return d1<d2 ? d1 : d2;
}

int Max(int d1, int d2) {
	return d1>d2 ? d1 : d2;
}

double Abs(double d1) {
	return d1>0 ? d1 : -1.0*d1;
}
