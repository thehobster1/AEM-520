#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


vector<double> DDx(vector<double>, double *, int *,int *,int);
vector<double> DDxDDx(vector<double>, double *, int *,int *,int);
vector<double> fill_gc(vector<double>, int,int);
vector<double> error(vector<double>,vector<double>);
void create(string,vector<double>,vector<double>,int);

int main()
{
 	const double x_low = 0,x_high = 6.28318531;
	int i,order = 2;
	int n = 10, gc = order/2;
	int sz = n + 2*gc;
	double dx;	
	vector<double> Acos(n + 2*gc),x(n + 2*gc),Asin(n + 2*gc),derivative(n + 2*gc),double_derivative(n + 2*gc),double_derivative_conservative(n + 2*gc),N_Asin(n + 2*gc);
	vector<double> err_d(n + 2*gc), err_dd(n + 2*gc), err_ddc(n + 2*gc);
	string A = "der_1_2_10.csv";
	string B = "der_2_2_10.csv";
	string C = "der_1_4_10.csv";
	string D = "der_2_4_10.csv";
	string E = "der_1_6_10.csv";
	string F = "der_2_6_10.csv";
	string G = "der_1_2_err_10.csv";
	string H = "der_2_2_err_10.csv";
	string I = "der_1_4_err_10.csv";
	string J = "der_2_4_err_10.csv";
	string K = "der_1_6_err_10.csv";
	string L = "der_2_6_err_10.csv";
	string M = "der_2_2_cons_err_10.csv";
	string N = "der_2_4_cons_err_10.csv";
	string O = "der_2_6_cons_err_10.csv";
	string P = "der_2_2_cons_10.csv";
	string Q = "der_2_4_cons_10.csv";
	string R = "der_2_6_cons_10.csv";


// For all orders at 10 cells

	// Second Order
	
	dx = (x_high - x_low)/(n-1);
	
	for ( i = 0; i < n + 2*gc; i++) {
		x[i] = x_low - gc*dx + i*dx; 	// create domain
	}
	

	for ( i = 0; i < n + 2*gc; i++){    // Analytical Solutions for comparison and derivative

		Asin[i]   =  sin(x[i]);
		N_Asin[i] = -sin(x[i]);
		Acos[i]   =  cos(x[i]);
	}
	
	derivative = DDx(Asin,&dx,&gc,&n,order);
	derivative = fill_gc(derivative,gc,n);
	double_derivative_conservative = DDx(derivative,&dx,&gc,&n,order);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(A,x,derivative,sz);
	create(B,x,double_derivative,sz);
	create(P,x,double_derivative_conservative,sz);
	create(G,x,err_d,sz);
	create(H,x,err_dd,sz);
	create(M,x,err_ddc,sz);
	//vector<double>
	


	
	// create(C,x,derivative4_s,sz);
	// create(D,x,derivative4_2_s,sz);
	// create(E,x,derivative6_s,sz);
	// create(F,x,derivative6_2_s,sz);
	
	// create(I,x,err_sm_d1_4,sz);
	// create(J,x,err_sm_d2_4,sz);
	// create(K,x,err_sm_d1_6,sz);
	// create(L,x,err_sm_d2_6,sz);

	// create(N,x,err_sm_d2_4_cons,sz);
	// create(O,x,err_sm_d2_6_cons,sz);

	// create(Q,x,derivative4_2_s_cons,sz);
	// create(R,x,derivative6_2_s_cons,sz);



	return 0;
	
}



vector<double> DDx(vector<double> A, double* dx,int* gc,int* n,int order)
{
    vector<double> derivative(*n + 2 * *gc);

		if (order == 2) {
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+1] - A[i-1])/2.0/ *dx;
	
		} 
		} else if (order == 4) {
		
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (-A[i+2] + 8*A[i+1] - 8*A[i-1] + A[i-2])/12.0/ *dx;
		} 
			
		} else if (order == 6) {

		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+3]/60 - 3*A[i+2]/20 + 3*A[i+1]/4 - 3*A[i-1]/4 + 3*A[i-2]/20 - A[i-3]/60)/ *dx;
		} 
		}

	return derivative;
}

vector<double> DDxDDx(vector<double> A, double* dx,int* gc,int* n,int order)
{
    vector<double> derivative(*n + 2 * *gc);

		if (order == 2) {
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+1] -2*A[i] + A[i-1])/ *dx/ *dx;
	
		} 
		} else if (order == 4) {
		
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (-A[i+2]/12.0 + 4.0*A[i+1]/3.0 - 5.0*A[i]/2.0 + 4*A[i-1]/3 - A[i-2]/12.0)/ *dx/ *dx;
	
		} 
			
		} else if (order == 6) {

		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+3]/90 - 3*A[i+2]/20 + 3*A[i+1]/2 - 49*A[i]/18 + 3*A[i-1]/2 - 3*A[i-2]/20 + A[i-3]/90)/ *dx/ *dx;
	
		} 
		}

	return derivative;
}

vector<double> fill_gc(vector<double> A, int gc,int n){
	
	int sz = n + 2*gc;
	for(int i = 0; i < gc; i++){
		A[i] = A[n - 1 + i];
		A[sz - i - 1] = A[2*gc - i];
	}	 
 
	return A;
}

vector<double> error(vector<double> Numer,vector<double> Exact){
	
	vector<double> error(Numer.size());
	
	for( int i = 0; i < Numer.size(); i++){
		error[i] = abs(Numer[i] - Exact[i]);
	}
	return error;
}


void create(string file, vector<double> x, vector<double> y,int n) 
{ 

	
    // file pointer 
    fstream fout; 

    // opens an existing csv file or creates a new file. 
    fout.open(file, ios::out | ios::app); 
  
    
  
    int i;
  
    // Read the input 
    for (i = 0; i < n; i++) { 
  
        fout << x[i] << ", "
            << y[i] << ", "
            <<"\n"; 
  
    } 
}
