#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


vector<double> DDx(vector<double>, double *, int *,int *,int);
vector<double> DDxDDx(vector<double>, double *, int *,int *,int);
vector<double> fill_gc(vector<double>, int,int);
//void create(string,double*,double*,int);

int main()
{
 	const double x_low = 0,x_high = 6.28318531;
	int i,order = 6;
	int n = 10, gc = order/2;
	double dx;	
	vector<double> Acos(n + 2*gc),x(n + 2*gc),Asin(n + 2*gc),derivative(n + 2*gc);
	
	dx = (x_high - x_low)/n;
	// create domains
	for ( i = 0; i < n + 2*gc; i++) {
		x[i] = x_low - gc*dx + i*dx;
	}
	

    // Analytical
	for ( i = 0; i < n + 2*gc; i++){ 
		Asin[i] = sin(x[i]);
		// cout << Asin[i] << "  " << x[i] << endl;
		Acos[i] = cos(x[i]);
	}
	
	derivative = DDx(Asin,&dx,&gc,&n,order);

	derivative = fill_gc(derivative,gc,n);
	// for ( i = 0; i < n + 2*gc; i++) {
	// 	cout << derivative[i] << "  " << x[i] << endl;
	// }
	
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

		if (order = 2) {
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+1] -2*A[i] + A[i-1])/ *dx/ *dx;
	
		} 
		} else if (order = 4) {
		
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (-A[i+2]/12.0 + 4.0*A[i+1]/3.0 - 5.0*A[i]/2.0 + 4*A[i-1]/3 - A[i-2]/12.0)/ *dx/ *dx;
	
		} 
			
		} else if (order = 6) {

		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+3]/90 - 3*A[i+2]/20 + 3*A[i+1]/2 - 49*A[i]/18 + 3*A[i-1]/2 - 3*A[i-2]/20 + A[i-3]/90)/ *dx/ *dx;
	
		} 
		}

	return derivative;
}

vector<double> fill_gc(vector<double> A, int gc,int n){
	
	int sz = n + 2*gc;
	for(int i = 0; i < gc; i++){
		A[i] = A[sz - 2*gc - i];
		//cout << A[i] << "   " << i << endl;
	//	A[sz - i] = A[2*gc - i];
	}	 
	
	for(int j = sz; j > sz - gc; j--){
		A[j] = A[j - 2*gc];
		//cout << A[j] << "    " << j << endl;
	}	 

	return A;
}

// void create(string file, double *x, double *y,int n) 
// { 

	
//     // file pointer 
//     fstream fout; 

//     // opens an existing csv file or creates a new file. 
//     fout.open(file, ios::out | ios::app); 
  
    
  
//     int i;
  
//     // Read the input 
//     for (i = 0; i < n; i++) { 
  
//         fout << x[i] << ", "
//             << y[i] << ", "
//             <<"\n"; 
  
//     } 
// }