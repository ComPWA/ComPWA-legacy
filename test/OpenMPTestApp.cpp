/*
 * OpenMPTestApp.cpp
 *
 *  Created on: Aug 29, 2016
 *      Author: weidenka
 */

#include "Physics/AmplitudeSum/AmpFlatteRes.cpp"
//#include <omp.h>

using namespace ComPWA;
using namespace ComPWA::Physics::AmplitudeSum;

int main()
{
	std::vector<double> m13 = {1.0,1.1,1.2,1.3};
	std::vector<double> m23 = {1.0,1.1,1.2,1.3};
	double m0 = 1.2;
	double ma = 0.4;
	double mb = 0.1;
	double g1 = 3;
	double massB1 = 0.3;
	double massB2 = 0.3;
	double g2 = 2;
	int i,n;
	int spin = 0;
	double d = 1.5;
#pragma omp parallel for
	for (auto it = m23.begin(); it != m23.end(); it++) /* i is private by default */{
		double mSq = (*it);
		std::complex<double> res = AmpFlatteRes::dynamicalFunction(
				mSq,m0,ma,mb,g1,massB1,massB2,g2,0,0,0,
				spin,d, formFactorType(0));
		std::cout<<res<<std::endl;
	}
}
