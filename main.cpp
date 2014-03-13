
#include <iostream>
#include <random>
#include "UKF.h"

class UKF2DPoint : public UKF {
public:
	Matrix state_function (Matrix s)
	{
		Matrix state(4,1);
		state(0,0) = s(0,0)+s(2,0);	// x position in 2D point
		state(1,0) = s(1,0)+s(3,0);	// y position in 2D point
		state(2,0) = s(2,0);	// velocity in x
		state(3,0) = s(3,0);	// velocity in y
		return state;
	}

	Matrix measurement_function (Matrix m)
	{
		Matrix measurement(2,1);
		measurement(0,0) = m(0,0);	// measured x position in 2D point
		measurement(1,0) = m(1,0);	// measured y position in 2D point
		return measurement;
	}
};

int main ()
{
	cout.precision(5);
	cout << fixed;
	default_random_engine generator;	// for artificial noise
	normal_distribution<double> distribution(0.0,1.0);	// for artificial noise

	unsigned int n = 4;
	unsigned int m = 2;
	unsigned int N = 20;	// total dynamic steps

	UKF2DPoint tracked_point;
	tracked_point.n = n;
	tracked_point.m = m;
	
	Matrix I4(n,n);	// 4x4 Identity Matrix
	I4(0,0) = 1;	I4(1,1) = 1;	I4(2,2) = 1;	I4(3,3) = 1;
	Matrix I2(m,m);	// 2x2 Identity Matrix
	I2(0,0) = 1;	I2(1,1) = 1;
	
	Matrix s(n,1);	// initial state
	s(0,0) = 1;	s(1,0) = 1;	s(2,0) = 0;	s(3,0) = 0;

	Matrix x = s; // initial state
	const double q=0.1;	//std of process. "smoothness". lower the value, smoother the curve
	const double r=0.1;	//std of measurement. "tracking". lower the value, faster the track
	tracked_point.P = I4;	// state covriance
	tracked_point.Q = (q*q) * I4;	// covariance of process	(size must be nxn)
	tracked_point.R = (r*r) * I2;	// covariance of measurement (size must be mxm)
	
	Matrix xV(n,N);	//estmate        // allocate memory to show outputs
	Matrix sV(n,N);	//actual
	Matrix zV(m,N);

	for (unsigned int k=0; k<N; k++)
	{
		double noise = distribution(generator);
		Matrix z = tracked_point.measurement_function(s);  // make measurments
		z += noise;	// add artificial noise in the measurement

		for (unsigned int i=0; i<n; i++)
		{
			sV(i,k) = s(i,0);	// save actual state
		}
		for (unsigned int i=0; i<m; i++)
		{
			zV(i,k) = z(i,0);	// save measurement
		}
		tracked_point.ukf(x, z);
		for (unsigned int i=0; i<n; i++)
		{
			xV(i,k) = x(i,0);	// save estimate
		}

		s = tracked_point.state_function(s);	// update process with artificial increment
		s += 1.0;
	}

	//cout << "sV:\n" << sV << endl;
	//cout << "\nxV:\n" << xV << endl;
	cout << "\nDifference between actual states and tracked measurements:\n";
	for (unsigned int k=0; k<N; k++)
	{		
		cout << k << ": " << abs( sV(0,k) - xV(0,k) ) << "  (" << sV(0,k) << '-' << xV(0,k) << ")\n";
	}
	system("pause");

	/* In actual practice, your application should look like this! */
	//for (unsigned int k=0; k<N; k++)
	//{
	//	Matrix z = [measurements];  // noised measurment
	//	YourUKF.ukf(x, z);
	//	Matrix w = YourUKF.measurement_function(x);	// estimated/corrected measurement
	//}

	return 1;
}