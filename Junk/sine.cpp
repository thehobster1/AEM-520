#include <iostream>
#include <cmath>

int main(){
	// Taylor Approx
	order = 6;
	for ( i = 1; i < n + 2*gc; i++){
		X = x[i];
		A[i] = SINE(&X, &order);
		//cout << A[i] << endl;
	}

	for ( i = 1; i < n + 2*gc; i++) {
		A2 = A[i+1];
		A1 = A[i-1];
		derivative[i] = DDx(&A2,&A1,&dx);		
	}
	
	order = 4;
	for ( i = 1; i < n + 2*gc; i++){
		X = x[i];
		B[i] = SINE(&X, &order);
	}

	for ( i = 1; i < n + 2*gc; i++) {
		B2 = B[i+1];
		B1 = B[i-1];
		derivativeB[i] = DDx(&B2,&B1,&dx);		
	}

	order = 2;
	for ( i = 1; i < n + 2*gc; i++){
		X = x[i];
		C[i] = SINE(&X, &order);
	}

	for ( i = 1; i < n + 2*gc; i++) {
		C2 = C[i+1];
		C1 = C[i-1];
		derivativeC[i] = DDx(&C2,&C1,&dx);		
	}
	
}
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