
#include "matrix.h"

using namespace std;
using namespace math;	// needed for matrix types

#ifndef Type
typedef double Type;	// You can change <double> to <float>
#endif
typedef matrix<Type> Matrix;

class UKF {
public:		
	void ukf( Matrix& x, const Matrix z);
	virtual Matrix state_function (Matrix x);
	virtual Matrix measurement_function (Matrix x);

	int n;      //number of state
	int m;      //number of measurement
	Matrix Q;	//noise covariance of process	(size must be nxn)
	Matrix R;	//noise covariance of measurement (size must be mxm)
	Matrix P;	//state covariance
};
