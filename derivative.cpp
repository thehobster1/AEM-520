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
double* fill_gc(double*, int,int);
void create(string,double*,double*,int);

int main()
{
 	const double x_low = 0,x_high = 6.283;
	int i,order;
	int n_sm = 10,n_med = 100,n_lg = 1000,gc = 0,sz_small = n_sm + 2*gc,sz_med = n_med + 2*gc,sz_large = n_lg + 2*gc;
	double domain_step, domain, nodes, dx_sm,dx_med,dx_lg, x_sm[sz_small] = {0},x_med[sz_med] = {0},x_lg[sz_large] = {0},derivative[sz_small] = {0},derivative_2[sz_small] = {0};
	double derivative4[sz_med] = {0}, derivative6[sz_large] = {0}, derivative6_2[sz_large] = {0}, derivative4_2[sz_med] = {0},Asin_sm[sz_small] = {0},Asin_med[sz_med] = {0},Asin_lg[sz_large] = {0};
	double Acos_sm[sz_small] = {0}, Acos_med[sz_med] = {0},Acos_lg[sz_large] = {0};
	double errsmd1[sz_small] = {0},errmedd1[sz_med] = {0},errlgd1[sz_large] = {0}, errsmd2[sz_small] = {0},errmedd2[sz_med] = {0},errlgd2[sz_large] = {0};
	double derivative4_s[sz_small] = {0},derivative4_2_s[sz_small] = {0},derivative6_s[sz_small] = {0},derivative6_2_s[sz_small] = {0},derivative2_med[sz_med] = {0},derivative2_2_med[sz_med] = {0},derivative6_med[sz_med] = {0},derivative6_2_med[sz_med] = {0};
	double derivative_lg[sz_large] = {0}, derivative_2_lg[sz_large] = {0}, derivative4_lg[sz_large] = {0}, derivative4_2_lg[sz_large] = {0};
	

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
	double err_sm_d1_2[sz_small] = {0},err_sm_d2_2[sz_small] = {0},err_sm_d1_4[sz_small] = {0},err_sm_d2_4[sz_small] = {0},err_sm_d1_6[sz_small] = {0},err_sm_d2_6[sz_small] = {0};
    for ( i = 1; i < n_sm + 2*gc; i++) {
		derivative[i] = DDx(&*Asin_sm,&dx_sm,&i);
		derivative_2[i] = DDx_2(&*Asin_sm,&dx_sm,&i);
		derivative4_s[i] = DDx4(&*Asin_sm,&dx_sm,&i);
		derivative4_2_s[i] = DDx4_2(&*Asin_sm,&dx_sm,&i);
		derivative6_s[i] = DDx6(&*Asin_sm,&dx_sm,&i);
		derivative6_2_s[i] = DDx6_2(&*Asin_sm,&dx_sm,&i);
		err_sm_d1_2[i] = derivative[i] - Acos_sm[i];
		err_sm_d2_2[i] = derivative_2[i] + Asin_sm[i];
		err_sm_d1_4[i] = derivative4[i] - Acos_sm[i];
		err_sm_d2_4[i] = derivative4_2[i] + Asin_sm[i];
		err_sm_d1_6[i] = derivative6[i] - Acos_sm[i];
		err_sm_d2_6[i] = derivative6_2[i] + Asin_sm[i];	
		
	}



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

	create(A,x_sm,derivative,sz_small);
	create(B,x_sm,derivative_2,sz_small);
	create(C,x_sm,derivative4_s,sz_small);
	create(D,x_sm,derivative4_2_s,sz_small);
	create(E,x_sm,derivative6_s,sz_small);
	create(F,x_sm,derivative6_2_s,sz_small);
	create(G,x_sm,err_sm_d1_2,sz_small);
	create(H,x_sm,err_sm_d2_2,sz_small);
	create(I,x_sm,err_sm_d1_4,sz_small);
	create(J,x_sm,err_sm_d2_4,sz_small);
	create(K,x_sm,err_sm_d1_6,sz_small);
	create(L,x_sm,err_sm_d2_6,sz_small);

    // 100 Nodes
	double err_m_d1_2[sz_med] = {0},err_m_d2_2[sz_med] = {0},err_m_d1_4[sz_med] = {0},err_m_d2_4[sz_med] = {0},err_m_d1_6[sz_med] = {0},err_m_d2_6[sz_med] = {0};

    for ( i = 1; i < n_med + 2*gc; i++) {
		derivative2_med[i] = DDx(&*Asin_med,&dx_med,&i);
		derivative2_2_med[i] = DDx_2(&*Asin_med,&dx_med,&i);
		derivative4[i] = DDx4(&*Asin_med,&dx_med,&i);
		derivative4_2[i] = DDx4_2(&*Asin_med,&dx_med,&i);
		derivative6_med[i] = DDx(&*Asin_med,&dx_med,&i);
		derivative6_2_med[i] = DDx_2(&*Asin_med,&dx_med,&i);
		err_m_d1_2[i] = derivative2_med[i] - Acos_med[i];
		err_m_d2_2[i] = derivative2_2_med[i] + Asin_med[i];
		err_m_d1_4[i] = derivative4[i] - Acos_med[i];
		err_m_d2_4[i] = derivative4_2[i] + Asin_med[i];
		err_m_d1_6[i] = derivative6_med[i] - Acos_med[i];
		err_m_d2_6[i] = derivative6_2_med[i] + Asin_med[i];	
	}

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
	
	create(A,x_med,derivative2_med,sz_med);
	create(B,x_med,derivative2_2_med,sz_med);
	create(C,x_med,derivative4,sz_med);
	create(D,x_med,derivative4_2,sz_med);
	create(E,x_med,derivative6_med,sz_med);
	create(F,x_med,derivative6_2_med,sz_med);
	create(G,x_med,err_m_d1_2,sz_med);
	create(H,x_med,err_m_d2_2,sz_med);
	create(I,x_med,err_m_d1_4,sz_med);
	create(J,x_med,err_m_d2_4,sz_med);
	create(K,x_med,err_m_d1_6,sz_med);
	create(L,x_med,err_m_d2_6,sz_med);


    // 1000 nodes
	double err_l_d1_2[sz_large] = {0},err_l_d2_2[sz_large] = {0},err_l_d1_4[sz_large] = {0},err_l_d2_4[sz_med] = {0},err_l_d1_6[sz_large] = {0},err_l_d2_6[sz_large] = {0};

	for ( i = 1; i < n_lg + 2*gc; i++) {
		
		derivative_lg[i] = DDx(&*Asin_lg,&dx_lg,&i);
		derivative_2_lg[i] = DDx_2(&*Asin_lg,&dx_lg,&i);
		derivative4_lg[i] = DDx4(&*Asin_lg,&dx_med,&i);
		derivative4_2_lg[i] = DDx4_2(&*Asin_lg,&dx_med,&i);
		derivative6[i] = DDx6(&*Asin_lg,&dx_lg,&i);
		derivative6_2[i] = DDx6_2(&*Asin_lg,&dx_lg,&i);	
		err_l_d1_2[i] = derivative_lg[i] - Acos_lg[i];
		err_l_d2_2[i] = derivative_2_lg[i] + Asin_lg[i];
		err_l_d1_4[i] = derivative4_lg[i] - Acos_lg[i];
		err_l_d2_4[i] = derivative4_2_lg[i] + Asin_lg[i];
		err_l_d1_6[i] = derivative6[i] - Acos_lg[i];
		err_l_d2_6[i] = derivative6_2[i] + Asin_lg[i];	
	}
	


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
	create(A,x_lg,derivative_lg,sz_large);
	create(B,x_lg,derivative_2_lg,sz_large);
	create(C,x_lg,derivative4_lg,sz_large);
	create(D,x_lg,derivative4_2_lg,sz_large);
	create(E,x_lg,derivative6,sz_large);
	create(F,x_lg,derivative6_2,sz_large);
	create(G,x_lg,err_m_d1_2,sz_large);
	create(H,x_lg,err_m_d2_2,sz_large);
	create(I,x_lg,err_m_d1_4,sz_large);
	create(J,x_lg,err_m_d2_4,sz_large);
	create(K,x_lg,err_m_d1_6,sz_large);
	create(L,x_lg,err_m_d2_6,sz_large);
	
	A = "cos_10.csv";
	B = "cos_100.csv";
	C = "cos_1000.csv";
	D = "sin_10.csv";
	E = "sin_100.csv";
	F = "sin_1000.csv";


	create(A,x_sm,Acos_sm,sz_small);
	create(B,x_med,Acos_med,sz_med);
	create(C,x_lg,Acos_lg,sz_large);
	create(D,x_sm,Asin_sm,sz_small);
	create(E,x_med,Asin_med,sz_med);
	create(F,x_lg,Asin_lg,sz_large);
	
	for ( i = 1; i < n_lg + 2*gc; i++) {
		cout << err_l_d1_6[i] << endl;
	}
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

double* fill_gc(double* A, int gc,int n){
	
	for(int i = 0; i < n; i++){
		A[i] = A[n-gc-i];
		A[n-i] = A[i+gc];
	}	 

	return A;
}

void create(string file, double *x, double *y,int n) 
{ 

	
    // file pointer 
    fstream fout; 
  //	std::string file;
//	fout.remove(file);

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
