#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

long long factorial(int*);
double DDx(double*, double*,int*);
double DDx4(double*, double*,int*);
double DDx6(double*, double*,int*);
void writeval(double*,double*,double*,double*,const int*);
int main()
{
 	const double x_low = 0,x_high = 6.283;
	int i,order;
	int n = 10,gc = 3;
	double domain_step, domain, nodes, dx, x[n + 2*gc] = {0},derivative[n + 2*gc] = {0};
	double derivative4[n + 2*gc] = {0}, derivative6[n + 2*gc] = {0}, Asin[n + 2*gc] = {0};
	//define dx
	dx = (x_high - x_low)/n;
	
	// create domain
	for ( i = 0; i < n + 2*gc; i++) {
	x[i] = x_low - gc*dx + i*dx;
	}
	
	
    // Analytical
	for ( i = 1; i < n + 2*gc; i++){
		Asin[i] = sin(x[i]);
	}
	

    // Second Order
    for ( i = 1; i < n + 2*gc; i++) {
		derivative[i] = DDx(&*Asin,&dx,&i);		
	}

    // Fourth Order
     for ( i = 1; i < n + 2*gc; i++) {
		derivative4[i] = DDx4(&*Asin,&dx,&i);		
	}


       for ( i = 1; i < n + 2*gc; i++) {
		derivative6[i] = DDx6(&*Asin,&dx,&i);		
	}
    writeval(&*x,&*derivative,&*derivative4,&*derivative6,&n);
    
    return 0;
	
}



double DDx(double* A, double* dx,int* i)
{
	
    double derivative = (A[*i+1] - A[*i-1])/2.0/ *dx;
    return derivative;
}

double DDx4(double* A, double* dx,int* i)
{
	double derivative = (-A[*i+2] + 8*A[*i+1] - 8*A[*i-1] + A[*i-2])/12.0/ *dx;
    return derivative;
}

double DDx6(double* A, double* dx,int* i)
{
	double derivative = (A[*i+3]/60 - 3*A[*i+2]/20 + 3*A[*i+1]/4 - 3*A[*i-1]/4 + 3*A[*i-2]/20 - A[*i-3]/60)/ *dx;
    return derivative;
}


long long factorial(int* n)
{
	long long fac=1,i;
	for (i=2;i<=*n; i++)
	fac*=i;
	//cout << fac << endl;
	return fac; 
	
}

void writeval(double* data, double* data2, double* data3,double* data4,const int* n) {
	std::ofstream myfile;
	myfile.open ("example.csv");

	for (int i = 0;i<=*n; i++){
		myfile << data[i] << ",";
        myfile << data2[i] << ",";
        myfile << data3[i] << ",";
        myfile << data4[i] << endl;
	}
      
    myfile.close();
}