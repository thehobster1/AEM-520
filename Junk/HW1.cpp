#include <iostream>
#include <cmath>
#include <fstream>
//#include "MatlabEngine.hpp"
//#include "MatlabDataArray.hpp"
using namespace std;

long long factorial(int*);
double DDx(double*, double*, double*);
double SINE(double*,int*);
//void writeval(double*,const int*);
int main()
{
 	const double x_low = 0,x_high = 6.283;
	int i,order;
	int n = 1000,gc = 3;
	double domain_step, domain, nodes, dx, x[n + 2*gc] = {0},derivative[n + 2*gc] = {0},A[n + 2*gc] = {0},A1,A2,X,B[n + 2*gc] = {0},B1,B2,C[n + 2*gc] = {0},C1,C2;
	double derivativeB[n + 2*gc] = {0}, derivativeC[n + 2*gc] = {0}, Asin[n + 2*gc] = {0},derivative2[n + 2*gc] = {0};
	//define dx
	dx = (x_high - x_low)/n;
	
	// create domain
	for ( i = 0; i < n + 2*gc; i++) {
	x[i] = x_low - gc*dx + i*dx;
	}
	
	// Analytical Derivative
	for ( i = 1; i < n + 2*gc; i++){
		X = x[i];
		Asin[i] = sin(x[i]);
		//cout << Asin[i] << endl;
	}
	
	for ( i = 1; i < n + 2*gc; i++) {
		derivative[i] = DDx(&Asin[i+1],&Asin[i-1],&dx);		
		cout << derivative[i] << endl;
	}

	for ( i = 1; i < n + 2*gc; i++) {
		derivative2[i] = DDx(&derivative[i+1],&derivative[i-1],&dx);		
		cout << derivative2[i] << endl;
	}
	
	// Taylor Approx

	
	//writeval(&*x,&n);
	//writeval(&*derivative,&n);
	return 0;
	
}



double DDx(double* x2, double* x1, double* dx)
{
	double derivative = (*x2 - *x1)/2.0/ *dx;
	return derivative;
}

//double DDx4(double* x2, double* x1, double* dx)
//{
//	double derivative = (*x2 - *x1)/2.0/ *dx;
//	return derivative;
//}

double SINE (double* x,int* n)
{
	double sum=0;
	int i;
	for (i=*n; i>=0; --i)
	{
		int j=2*i+1;
		long long denom = factorial(&j);
		sum=sum+((pow((-1),i)*pow(*x,((2*i)+1))/denom));
	
	}
	return sum;
}

long long factorial(int* n)
{
	long long fac=1,i;
	for (i=2;i<=*n; i++)
	fac*=i;
	//cout << fac << endl;
	return fac; 
	
}

//void writeval(double* data, const int* n) {
//	std::ofstream myfile;
//	myfile.open ("example.csv");
//	
//	for (int i = 0;i<=*n; i++){
//		myfile << data[i] << endl;
//		myfile << ",";
//	}
//
//}