#include<fstream>
#include<cmath>
#include<cstdlib>

using namespace std;

void multiply(double **A, double **B,double **C,int m, int l, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < l; k++)
			{
				C[i][j] +=A[i][k]*B[k][j];
			}
		}
	}
	return;
}

void transpose(double **A, double **B, int m,int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			B[j][i] = A[i][j];
		}
	}
	return;
}

void init(double **A, int n)
{
	for (int i = 0; i < n ; i++)
	{
		A[i][i] = 2;
		if (i > 0)
		{
			A[i][i-1] = A[i-1][i] = -1;
		}
	}
	return;
}

double sgn(double tau)
{
	if (tau > 0)
	{
		return 1;
	}
	else if (tau == 0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

void jcb(double **A, int p, int q, int n)
{
	double **J, **B, **JT;
	double c, s;
	J = new double *[n];
	JT = new double *[n];
	B = new double *[n];
	for (int i = 0; i < n; i++)
	{
		J[i] = new double [n];
		B[i] = new double [n];
		JT[i] = new double [n];
	}
	if (A[p][p] == A[q][q])
	{
		c = pow(1./2,1./2);
		s = c;
	}
	else
	{
		double tau, t;
		tau = (A[p][p]-A[q][q])/(2*A[p][q]);
		t = sgn(tau)/(fabs(tau)+pow(1+tau*tau,1./2));
		c = 1/pow(1+t*t,1./2);
		s = t*c;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if(i == j)
			{
				J[i][j] = 1;
			}
			else
			{
				J[i][j] = 0;
			}
		}
	}
	J[p][p] = c;
	J[q][q] = c;
	J[p][q] = s;
	J[q][p] = -s;
	multiply(J,A,B,n,n,n);
	transpose(J,JT,n,n);
	multiply(B,JT,A,n,n,n);
	
	
	for (int i = 0; i < n; i++)
	{
		delete [] J[i];
		delete [] JT[i];
		delete [] B[i];
	}
	delete [] J;
	delete [] JT;
	delete [] B;
	return;
}

double check(double **A, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				sum += pow(A[i][j],2);
			}
		}
	}
	return sum;
}

void jcbmethod(double **A,int n)
{
	int i,p,q;
	double tol;
	tol = pow(check(A,n)/(n*(n-1)),1./2);
	for (i = 0; tol>10e-10 && i < 20; i++)
	{
		for (p = 0; p < n-1; p++)
		{
			for (q = p+1; q < n; q++)
			{
				if (fabs(A[p][q]) > tol)
				{
					jcb(A,p,q,n);
				}
			}
		}
		tol = check(A,n);
		tol = tol/(n*(n-1));
		tol = pow(tol,1./2);
	}
	return;
}

void beautify(double **A, int n)
{
	double tol = 1e-8;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (fabs(A[i][j]) < tol)
			{
				A[i][j] = 0;
			}
		}
	}
	return;
}


int main()
{
	int n,i,j;
	double **A;
	ofstream outfile;
	outfile.open("losc8-1data");
	for (n = 4; n < 7; n++)
	{
	A = new double *[n];
	for (i = 0; i < n; i++)
	{
		A[i] = new double [n];
	}
	init(A,n);
	jcbmethod(A,n);
	outfile << "n = " << n << endl;
	outfile << "[" << endl;
	beautify(A,n);
	for (i = 0; i < n; i++)
	{
		outfile << "[";
		for (j = 0; j < n; j++)
		{
			outfile << A[i][j] << ",";
		}
		outfile << "]" << endl;
	}
	outfile  << "]" << endl;
	}
	return 0;
}
