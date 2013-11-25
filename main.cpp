
#include "UKF.h"

class UKF2DPoint : public UKF {
public:
	Matrix state_function (Matrix x)
	{
		Matrix state(4,1);
		state(0,0) = x(0,0)+x(2,0);
		state(1,0) = x(1,0)+x(3,0);
		state(2,0) = x(2,0);
		state(3,0) = x(3,0);
		return state;
	}

	Matrix measurement_function (Matrix x)
	{
		Matrix measurement(2,1);
		measurement(0,0) = x(0,0);
		measurement(1,0) = x(1,0);
		return measurement;
	}
};

int main ()
{
	cout.precision(5);

	unsigned int n=4;
	unsigned int m=2;
	unsigned int N = 20;	// total dynamic steps

	UKF2DPoint Example;
	Example.n = n;
	Example.m = m;
	
	Matrix I4(n,n);	// 4x4 Identity Matrix
	I4(0,0) = 1;	I4(1,1) = 1;	I4(2,2) = 1;	I4(3,3) = 1;
	Matrix I2(m,m);	// 2x2 Identity Matrix
	I2(0,0) = 1;	I2(1,1) = 1;
	
	Matrix s(n,1);	// initial state
	s(0,0) = 1;	s(1,0) = 1;	s(2,0) = 0;	s(3,0) = 0;

	Matrix x = s; /*x=s+q*randn(2,0);*/	// initial state with noise
	const double q=0.1;	//std of process. "smoothness". lower the value, smoother the curve
	const double r=0.1;	//std of measurement. "tracking". lower the value, faster the track
	Example.P = I4;	// state covriance
	Example.Q = (q*q) * I4;	// covariance of process	(size must be nxn)
	Example.R = (r*r) * I2;	// covariance of measurement (size must be mxm)
	
	Matrix xV(n,N);	//estmate        // allocate memory to show outputs
	Matrix sV(n,N);	//actual
	Matrix zV(m,N);

	for (unsigned int k=0; k<N; k++)
	{
		Matrix z = Example.measurement_function(s);  // measurments
		z += 0.01;	// artificial noise (remove this line in your application!)
		for (unsigned int i=0; i<n; i++)
		{
			sV(i,k) = s(i,0);	// save actual state
		}
		for (unsigned int i=0; i<m; i++)
		{
			zV(i,k) = z(i,0);	// save measurement
		}
		Example.ukf(x, z);
		for (unsigned int i=0; i<n; i++)
		{
			xV(i,k) = x(i,0);	// save estimate
		}
		s = Example.state_function(s);	// update process
		s += 0.01;	// artificial noise (remove this line in your application!)
	}

	cout << "sV:\n" << sV << endl;
	cout << "\nxV:\n" << xV << endl;
	system("pause");

	/* Your application should look like this! */
	//for (unsigned int k=0; k<N; k++)
	//{
	//	Matrix z = [measurements];  // measurment
	//	YourUKF.ukf(x, z);
	//	Matrix w = YourUKF.measurement_function(x);	// estimated measurement
	//}

	return 1;
}