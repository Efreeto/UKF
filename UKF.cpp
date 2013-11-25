
#include "UKF.h"

Matrix sigmas(Matrix x, Matrix P, Type c);
Matrix Cholesky(Matrix A);

void UKF::ukf( Matrix& x, const Matrix z)
{
	const int L=2*n+1;
	const Type alpha = 1e-3;                                 //default, tunable
	const Type ki = 0.0;                                     //default, tunable
	const Type beta = 2.0;                                   //default, tunable
	const Type lambda = (alpha*alpha)*(n+ki)-n;              //scaling factor
	Type c = n+lambda;                                 //scaling factor

	/* weight equations are found in the upper part of http://www.cslu.ogi.edu/nsel/ukf/node6.html */
	Matrix Wm(1,L);	//weights for means	
	Matrix Wc = Wm;	//weights for covariance
	Wm(0,0) = lambda/c;
	Wc(0,0) = lambda/c+(1-(alpha*alpha)+beta);
	for (unsigned int k=1; k<L; k++)
	{
		Wm(0,k) = 0.5/c;
		Wc(0,k) = 0.5/c;
	}
	c = sqrt(c);
	Matrix X = sigmas(x, P, c);	//sigma points around x
	
	/* unscented transformation (ut) of process */
	Matrix x1(n,1);
	Matrix X1(n,L);
	for(unsigned int k=0; k<L; k++)
	{
		Matrix Xcol(n,1);
		Matrix X1col(n,1);	/* temp vectors, not used in matlab */

		for (unsigned int i=0; i<n; i++)
		{
			Xcol(i,0) = X(i,k);	// take out a column so that state_function can take it
		}		
		X1col = state_function(Xcol);
		for (unsigned int i=0; i<n; i++)
		{
			x1(i,0) +=  Wm(0,k) * X1col(i,0);
		}
		for (unsigned int i=0; i<n; i++)
		{
			X1(i,k) = X1col(i,0);	// put back the output column
		}	
	}
	Matrix X2(n,L);
	for (unsigned int k=0; k<L; k++)
		for (unsigned int i=0; i<n; i++)
		{
			X2(i,k) = X1(i,k) - x1(i,0);	//X2.Column(k) = X1.Column(k) - x1;
		}
	Matrix diagWm(L,L);
	for(unsigned int k=0; k<L; k++)
		diagWm(k,k) = Wm(0,k);
	Matrix P1 = X2 * diagWm * ~X2 + Q;	/* ~ means transpose */

	/* unscented transformation (ut) of measurements */
	Matrix z1(m,1);
	Matrix Z1(m,L);
	for(unsigned int k=0; k<L; k++)
	{
		Matrix X1col(n,1);
		Matrix Z1col(m,1);	/* temp vectors, not used in matlab */

		for (unsigned int i=0; i<n; i++)
		{
			X1col(i,0) = X1(i,k);	// take out a column so that measurement_function can take it
		}		
		Z1col = measurement_function(X1col);
		for (unsigned int i=0; i<m; i++)
		{
			z1(i,0) += Wm(0,k) * Z1col(i,0);
		}
		for (unsigned int i=0; i<m; i++)
		{
			Z1(i,k) = Z1col(i,0);	// put back the output column
		}	
	}
	Matrix Z2(m,L);
	for (unsigned int k=0; k<L; k++)
		for (unsigned int i=0; i<m; i++)
		{
			Z2(i,k) = Z1(i,k) - z1(i,0);	//Z2.Column(k) = Z1.Column(k) - z1;
		}
	Matrix diagWc(L,L);
	for(unsigned int k=0; k<L; k++)
		diagWc(k,k) = Wc(0,k);
	Matrix P2 = Z2 * diagWc * ~Z2 + R;	

	Matrix P12 = X2 * diagWc * ~Z2;	//transformed cross-covariance
	Matrix K = P12 * !P2;
	x = x1+K*(z-z1);                              //state update
	//cout << x << endl << K << endl;
	P = P1-K*~P12;                                //covariance update
	//cout << P << endl << endl;
}

Matrix sigmas(Matrix x, Matrix P, Type c)
{
	const int n = P.RowNo();
	const int L=2*n+1;
	Matrix Chol = Cholesky(P);
	Matrix A = c * Chol;	/* doesn't need transpose here like matlab b/c Chol gives a lower (not upper) triangle matrix */
			
	Matrix X(n,L);
	unsigned int k=0;
	{
		for (unsigned int i=0; i<n; i++)
		{
			X(i,k) = x(i,0);
		}
	}
	for(k=1; k<n+1; k++)
	{
		for (unsigned int i=0; i<n; i++)
		{
			X(i,k) = x(i,0) + A(i,k-1);
		}
	}
	for(k=n+1; k<L; k++)
	{
		for (unsigned int i=0; i<n; i++)
		{
			X(i,k) = x(i,0) - A(i,k-1-n);
		}
	}

	return X;
}

/* courtesy of http://www.saiensu.co.jp/book_support/4-7819-0855-1/cholesky.txt and http://chips.ncsu.edu/~luw/version2.1.5/HTML_SOURCE/cholesky.c.html */
// returns a lower triangle matrix
Matrix Cholesky(Matrix A)
{
	const int n = A.RowNo();
	Matrix Chol(n,n);
	Type* s = new Type[n];
	Type ss;
	unsigned int i,j,k;
	for (j=0; j<n; j++)
	{
		if (j==0) 
		{
			for (i=j; i<n; i++)
			{
				s[i]=A(i,j);
			} 
		}
		if (j!=0)
		{
			for (i=j; i<n; i++)
			{
				ss=0.;
				for (k=0; k<=j-1; k++)
				{
					ss+=Chol(i,k)*Chol(j,k);
				}
				s[i]=A(i,j)-ss;
			}
		}
		if (fabs(s[j])<0.000001)
		{
			cout << "ERROR: ukf.cpp Cholesky - Matrix not positive definite\n";
			delete [] s;
			return A;	// ERROR
		}
		for (i=j; i<n; i++)
		{
			Chol(i,j)=s[i]/sqrt(s[j]);
		}
	}
	delete [] s;
	return Chol;
}

Matrix UKF::state_function (Matrix x)	// override this function with inherited UKF class
{
	cout << "UKF: Please form your own state function!" << endl;
	return x;
}

Matrix UKF::measurement_function (Matrix x)	// override this function with inherited UKF class
{
	cout << "UKF: Please form your own measurement function!" << endl;
	return x;
}