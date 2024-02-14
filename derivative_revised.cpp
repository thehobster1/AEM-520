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
	string A = "der_1_2_10.csv",
	 B = "der_2_2_10.csv",
	 C = "der_1_4_10.csv",
	 D = "der_2_4_10.csv",
	 E = "der_1_6_10.csv",
	 F = "der_2_6_10.csv",
	 G = "der_1_2_err_10.csv",
	 H = "der_2_2_err_10.csv",
	 I = "der_1_4_err_10.csv",
	 J = "der_2_4_err_10.csv",
	 K = "der_1_6_err_10.csv",
	 L = "der_2_6_err_10.csv",
	 M = "der_2_2_cons_err_10.csv",
	 N = "der_2_4_cons_err_10.csv",
	 O = "der_2_6_cons_err_10.csv",
	 P = "der_2_2_cons_10.csv",
	 Q = "der_2_4_cons_10.csv",
	 R = "der_2_6_cons_10.csv",
	 AA = "cos_10.csv",
	 BB = "cos_100.csv",
	 CC = "cos_1000.csv",
	 DD = "sin_10.csv",
	 EE = "sin_100.csv",
	 FF = "sin_1000.csv";

// For all orders at 10 cells

	// Second Order -------------------------------------------------------------------------------------------------------------
	
	dx = (x_high - x_low)/(n-1);
	
	for ( i = 0; i < n + 2*gc; i++) {
		x[i] = x_low - gc*dx + i*dx; 	// create domain
	}
	

	for ( i = 0; i < n + 2*gc; i++){    // Analytical Solutions for comparison and derivative

		Asin[i]   =  sin(x[i]);
		N_Asin[i] = -sin(x[i]);
		Acos[i]   =  cos(x[i]);
	}
	create(AA,x,Acos,sz);
	create(DD,x,Asin,sz);


	derivative = DDx(Asin,&dx,&gc,&n,order);
	derivative = fill_gc(derivative,gc,n);
	double_derivative_conservative = DDx(derivative,&dx,&gc,&n,order);
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(A,x,derivative,sz);
	create(B,x,double_derivative,sz);
	create(P,x,double_derivative_conservative,sz);
	create(G,x,err_d,sz);
	create(H,x,err_dd,sz);
	create(M,x,err_ddc,sz);
	
	// Fourth Order -------------------------------------------------------------------------------------------------------------
	order = 4, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(C,x,derivative,sz);
	create(D,x,double_derivative,sz);
	create(Q,x,double_derivative_conservative,sz);
	create(I,x,err_d,sz);
	create(J,x,err_dd,sz);
	create(N,x,err_ddc,sz);

	// Sixth Order -------------------------------------------------------------------------------------------------------------

	order = 6, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);

	create(E,x,derivative,sz);
	create(F,x,double_derivative,sz);
	create(K,x,double_derivative_conservative,sz);
	create(L,x,err_d,sz);
	create(O,x,err_dd,sz);
	create(R,x,err_ddc,sz);


// For all orders at 100 cells
	n = 100;
	order = 2, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
	
	A = "der_1_2_100.csv";
	B = "der_2_2_100.csv";
	C = "der_1_4_100.csv";
	D = "der_2_4_100.csv";
	E = "der_1_6_100.csv";
	F = "der_2_6_100.csv";
	G = "der_1_2_err_100.csv";
	H = "der_2_2_err_100.csv";
	I = "der_1_4_err_100.csv";
	J = "der_2_4_err_100.csv";
	K = "der_1_6_err_100.csv";
	L = "der_2_6_err_100.csv";
	M = "der_2_2_cons_err_100.csv";
	N = "der_2_4_cons_err_100.csv";
	O = "der_2_6_cons_err_100.csv";
	P = "der_2_2_cons_100.csv";
	Q = "der_2_4_cons_100.csv";
	R = "der_2_6_cons_100.csv";
	
	// Second Order -------------------------------------------------------------------------------------------------------------
	
	dx = (x_high - x_low)/(n-1);
	
	for ( i = 0; i < n + 2*gc; i++) {
		x[i] = x_low - gc*dx + i*dx; 	// create domain
	}
	

	for ( i = 0; i < n + 2*gc; i++){    // Analytical Solutions for comparison and derivative

		Asin[i]   =  sin(x[i]);
		N_Asin[i] = -sin(x[i]);
		Acos[i]   =  cos(x[i]);
	}
	
	create(BB,x,Acos,sz);
	create(EE,x,Asin,sz);


	derivative = DDx(Asin,&dx,&gc,&n,order);
	derivative = fill_gc(derivative,gc,n);
	double_derivative_conservative = DDx(derivative,&dx,&gc,&n,order);
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(A,x,derivative,sz);
	create(B,x,double_derivative,sz);
	create(P,x,double_derivative_conservative,sz);
	create(G,x,err_d,sz);
	create(H,x,err_dd,sz);
	create(M,x,err_ddc,sz);
	
	// Fourth Order -------------------------------------------------------------------------------------------------------------
	order = 4, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(C,x,derivative,sz);
	create(D,x,double_derivative,sz);
	create(Q,x,double_derivative_conservative,sz);
	create(I,x,err_d,sz);
	create(J,x,err_dd,sz);
	create(N,x,err_ddc,sz);

	// Sixth Order -------------------------------------------------------------------------------------------------------------

	order = 6, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);

	create(E,x,derivative,sz);
	create(F,x,double_derivative,sz);
	create(K,x,double_derivative_conservative,sz);
	create(L,x,err_d,sz);
	create(O,x,err_dd,sz);
	create(R,x,err_ddc,sz);


// For all orders at 1000 cells
	n = 1000;
	order = 2, gc = order/2,sz = n + 2*gc;

	 Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	 err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);

	A = "der_1_2_1000.csv";
	B = "der_2_2_1000.csv";
	C = "der_1_4_1000.csv";
	D = "der_2_4_1000.csv";
	E = "der_1_6_1000.csv";
	F = "der_2_6_1000.csv";
	G = "der_1_2_err_1000.csv";
	H = "der_2_2_err_1000.csv";
	I = "der_1_4_err_1000.csv";
	J = "der_2_4_err_1000.csv";
	K = "der_1_6_err_1000.csv";
	L = "der_2_6_err_1000.csv";
	M = "der_2_2_cons_err_1000.csv";
	N = "der_2_4_cons_err_1000.csv";
	O = "der_2_6_cons_err_1000.csv";
	P = "der_2_2_cons_1000.csv";
	Q = "der_2_4_cons_1000.csv";
	R = "der_2_6_cons_1000.csv";
	
	// Second Order -------------------------------------------------------------------------------------------------------------
	
	dx = (x_high - x_low)/(n-1);
	
	for ( i = 0; i < n + 2*gc; i++) {
		x[i] = x_low - gc*dx + i*dx; 	// create domain
	}
	

	for ( i = 0; i < n + 2*gc; i++){    // Analytical Solutions for comparison and derivative

		Asin[i]   =  sin(x[i]);
		N_Asin[i] = -sin(x[i]);
		Acos[i]   =  cos(x[i]);
	}
	
	create(CC,x,Acos,sz);
	create(FF,x,Asin,sz);

	derivative = DDx(Asin,&dx,&gc,&n,order);
	derivative = fill_gc(derivative,gc,n);
	double_derivative_conservative = DDx(derivative,&dx,&gc,&n,order);
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(A,x,derivative,sz);
	create(B,x,double_derivative,sz);
	create(P,x,double_derivative_conservative,sz);
	create(G,x,err_d,sz);
	create(H,x,err_dd,sz);
	create(M,x,err_ddc,sz);
	
	// Fourth Order -------------------------------------------------------------------------------------------------------------
	order = 4, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);


	create(C,x,derivative,sz);
	create(D,x,double_derivative,sz);
	create(Q,x,double_derivative_conservative,sz);
	create(I,x,err_d,sz);
	create(J,x,err_dd,sz);
	create(N,x,err_ddc,sz);

	// Sixth Order -------------------------------------------------------------------------------------------------------------

	order = 6, gc = order/2,sz = n + 2*gc;
	Acos.resize(n + 2*gc),x.resize(n + 2*gc),Asin.resize(n + 2*gc),derivative.resize(n + 2*gc),double_derivative.resize(n + 2*gc),double_derivative_conservative.resize(n + 2*gc),N_Asin.resize(n + 2*gc);
	err_d.resize(n + 2*gc), err_dd.resize(n + 2*gc), err_ddc.resize(n + 2*gc);
	
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
	double_derivative_conservative = fill_gc(double_derivative_conservative,gc,n);
	double_derivative = DDxDDx(Asin,&dx,&gc,&n,order);
	double_derivative = fill_gc(double_derivative,gc,n);

	err_d   = error(derivative,Acos);
	err_dd  = error(double_derivative,N_Asin);
	err_ddc = error(double_derivative_conservative,N_Asin);

	create(E,x,derivative,sz);
	create(F,x,double_derivative,sz);
	create(K,x,double_derivative_conservative,sz);
	create(L,x,err_d,sz);
	create(O,x,err_dd,sz);
	create(R,x,err_ddc,sz);

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
