#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;


double DDx(double *, double *, int *);
double DDx4(double*, double*,int*);
double DDx6(double*, double*,int*);
double DDx_2(double*, double*,int*);
double DDx4_2(double*, double*,int*);
double DDx6_2(double*, double*,int*);
void writeval(double*,double*,double*,double*,const int*);
void writevalmed(double*,double*,double*,double*,const int*);
void writevallarge(double*,double*,double*,double*,const int*);
void create(string,double*,double*,int);

int main()
{
 	const double x_low = 0,x_high = 6.283;
	int i,order;
	int n_sm = 10,n_med = 100,n_lg = 1000,gc = 0,sz_small = n_sm + 2*gc,sz_med = n_med + 2*gc,sz_large = n_lg + 2*gc;
	double domain_step, domain, nodes, dx_sm,dx_med,dx_lg, x_sm[n_sm + 2*gc] = {0},x_med[n_med + 2*gc] = {0},x_lg[n_lg + 2*gc] = {0},derivative[n_sm + 2*gc] = {0},derivative_2[n_sm + 2*gc] = {0};
	double derivative4[n_med + 2*gc] = {0}, derivative6[n_lg + 2*gc] = {0}, derivative6_2[n_lg + 2*gc] = {0}, derivative4_2[n_med + 2*gc] = {0},Asin_sm[n_sm + 2*gc] = {0},Asin_med[n_med + 2*gc] = {0},Asin_lg[n_lg + 2*gc] = {0};
	double Acos_sm[n_sm + 2*gc] = {0}, Acos_med[n_med + 2*gc] = {0},Acos_lg[n_lg + 2*gc] = {0};
	double errsmd1[n_sm + 2*gc] = {0},errmedd1[n_med + 2*gc] = {0},errlgd1[n_lg + 2*gc] = {0}, errsmd2[n_sm + 2*gc] = {0},errmedd2[n_med + 2*gc] = {0},errlgd2[n_lg + 2*gc] = {0};
	double derivative4_s[n_sm + 2*gc] = {0},derivative4_2_s[n_sm + 2*gc] = {0},derivative6_s[n_sm + 2*gc] = {0},derivative6_2_s[n_sm + 2*gc] = {0};
	
	dx_sm = (x_high - x_low)/n_sm;
	dx_med = (x_high - x_low)/n_med;
	dx_lg = (x_high - x_low)/n_lg;
	// create domains
	for ( i = 0; i < n_sm + 2*gc; i++) {
	x_sm[i] = x_low - gc*dx_sm + i*dx_sm;
	}
	
	for ( i = 0; i < n_med + 2*gc; i++) {
	x_med[i] = x_low - gc*dx_med + i*dx_med;
	}
	
	for ( i = 0; i < n_lg + 2*gc; i++) {
	x_lg[i] = x_low - gc*dx_lg + i*dx_lg;
	}

    // Analytical
	for ( i = 1; i < n_sm + 2*gc; i++){
		Asin_sm[i] = sin(x_sm[i]);
		Acos_sm[i] = cos(x_sm[i]);
	}

	for ( i = 1; i < n_med + 2*gc; i++){
		Asin_med[i] = sin(x_med[i]);
		Acos_med[i] = cos(x_med[i]);
	}
	
	for ( i = 1; i < n_lg + 2*gc; i++){
		Asin_lg[i] = sin(x_lg[i]);
		Acos_lg[i] = cos(x_lg[i]);
	}
	
    // For 10 nodes
    for ( i = 1; i < n_sm + 2*gc; i++) {
		derivative[i] = DDx(&*Asin_sm,&dx_sm,&i);
		derivative_2[i] = DDx_2(&*Asin_sm,&dx_sm,&i);
		derivative4_s[i] = DDx4(&*Asin_sm,&dx_sm,&i);
		derivative4_2_s[i] = DDx4_2(&*Asin_sm,&dx_sm,&i);
		derivative6_s[i] = DDx6(&*Asin_lg,&dx_lg,&i);
		derivative6_2_s[i] = DDx6_2(&*Asin_lg,&dx_lg,&i);
		errsmd1[i] = derivative[i] - Acos_sm[i];
		errsmd2[i] = derivative_2[i] + Asin_sm[i];	
	}

	string A = "der_1_2_10.csv";
	string B = "der_2_2_10.csv";
	string C = "der_1_4_10.csv";
	string D = "der_2_4_10.csv";
	string E = "der_1_6_10.csv";
	string F = "der_2_6_10.csv";
	create(A,x_sm,derivative,sz_small);
	create(B,x_sm,derivative_2,sz_small);
	create(C,x_sm,derivative4_s,sz_small);
	create(D,x_sm,derivative4_2_s,sz_small);
	create(E,x_sm,derivative6_s,sz_small);
	create(F,x_sm,derivative6_2_s,sz_small);


    // Fourth Order
    for ( i = 1; i < n_med + 2*gc; i++) {
		derivative4[i] = DDx4(&*Asin_med,&dx_med,&i);
		derivative4_2[i] = DDx4_2(&*Asin_med,&dx_med,&i);
		errmedd1[i] = derivative4[i] - Acos_med[i];		
		errmedd2[i] = derivative4_2[i] + Asin_med[i];
	}

	A = "der_1_4_100.csv";
	B = "der_2_4_100.csv";
	C = "der_1_4_100_error.csv";
	D = "der_2_4_100_error.csv";
	create(A,x_med,derivative4,sz_med);
	create(B,x_med,derivative4_2,sz_med);
	create(C,x_med,derivative4,sz_med);
	create(D,x_med,derivative4_2,sz_med);



    // Sixth Order
	for ( i = 1; i < n_lg + 2*gc; i++) {
		derivative6[i] = DDx6(&*Asin_lg,&dx_lg,&i);
		derivative6_2[i] = DDx6_2(&*Asin_lg,&dx_lg,&i);	
		errlgd1[i] = derivative6[i] - Acos_lg[i];	
		errlgd2[i] = derivative6_2[i] + Asin_lg[i];	
	}
	
	A = "der_1_6_1000.csv";
	B = "der_2_6_1000.csv";
	C = "der_1_6_1000_error.csv";
	D = "der_2_6_1000_error.csv";
	create(A,x_lg,derivative6,sz_large);
	create(B,x_lg,derivative6_2,sz_large);
	create(C,x_lg,derivative6,sz_large);
	create(D,x_lg,derivative6_2,sz_large);
    
	
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

double DDx_2(double* A, double* dx,int* i)
{
	
    double derivative = (A[*i+1] -2*A[*i] + A[*i-1])/ *dx/ *dx;
    return derivative;
}

double DDx4_2(double* A, double* dx,int* i)
{
	double derivative = (-A[*i+2]/12.0 + 4.0*A[*i+1]/3.0 - 5.0*A[*i]/2.0 + 4*A[*i-1]/3 - A[*i-2]/12.0)/ *dx/ *dx;
    return derivative;
}

double DDx6_2(double* A, double* dx,int* i)
{
	double derivative = (A[*i+3]/90 - 3*A[*i+2]/20 + 3*A[*i+1]/2 - 49*A[*i]/18 + 3*A[*i-1]/2 - 3*A[*i-2]/20 + A[*i-3]/90)/ *dx/ *dx;
    return derivative;
}

void writeval(double* data, double* data2, double* data3,double* data4,const int* n) {
	std::ofstream myfile;
	myfile.open ("Small.csv");

	for (int i = 0;i<=*n; i++){
		myfile << data[i] << ",";
        myfile << data2[i] << ",";
        myfile << data3[i] << ",";
        myfile << data4[i] << endl;
	}

    myfile.close();
}

void writevalmed(double* data, double* data2, double* data3,double* data4,const int* n) {
	std::ofstream myfile;
	myfile.open ("Medium.csv");

	for (int i = 0;i<=*n; i++){
		myfile << data[i] << ",";
        myfile << data2[i] << ",";
        myfile << data3[i] << ",";
        myfile << data4[i] << endl;
	}

    myfile.close();
}

void writevallarge(double* data, double* data2, double* data3,double* data4,const int* n) {
	std::ofstream myfile;
	myfile.open ("Large.csv");

	for (int i = 0;i<=*n; i++){
		myfile << data[i] << ",";
        myfile << data2[i] << ",";
        myfile << data3[i] << ",";
        myfile << data4[i] << endl;
	}

    myfile.close();
}



void create(string file, double *x, double *y,int n) 
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
