//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//     Peter Weidenkaff -  assignment of final state particle masses
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpWigner.hpp"

#include "qft++.h"

AmpWigner::AmpWigner()
{
	toEvaluate=false;
}
AmpWigner::AmpWigner(const char *name,
		unsigned int spin, unsigned int m,  unsigned int n,unsigned int subSys) :
						  _inSpin(spin),
						  _outSpin1(m),
						  _outSpin2(n),
						  _subSys(subSys)
{
	toEvaluate=true;
	initialise();
}

AmpWigner::AmpWigner(const AmpWigner& other, const char* newname) :
						  _inSpin(other._inSpin),
						  _outSpin1(other._outSpin1),
						  _outSpin2(other._outSpin2),
						  _subSys(other._subSys),
						  toEvaluate(other.toEvaluate)
{
	initialise();
}

AmpWigner::~AmpWigner() 
{
}

void AmpWigner::initialise() 
{
	static dataPoint* point = dataPoint::instance();
	_M=point->DPKin.M;
	//if(_subSys==5){
		_m1=point->DPKin.m1;
		_m2=point->DPKin.m2;
		_m3=point->DPKin.m3;//}
	//if(_subSys==4){
		//_m1=point->DPKin.m2;
		//_m2=point->DPKin.m3;
		//_m3=point->DPKin.m1;//}
	//if(_subSys==3){
		//_m1=point->DPKin.m3;
		//_m2=point->DPKin.m2;
		//_m3=point->DPKin.m1;//}
//	cout<<"AmpWigner DEBUG set masses to m1="<<_m1<< " m2="<<_m2<<" m3=" <<_m3<<endl;

}    
void AmpWigner::setDecayMasses(double m1, double m2, double m3, double M){
	_M=M; _m1=m1; _m2=m2; _m3=m3;
	return;
}

double AmpWigner::evaluate() const {
	if(!toEvaluate) return 1.;

	double locmin_sq, locmax_sq, beta;

//	cout<<"==== "<<_subSys<< " " <<_m1<< " "<<_m2<<" " <<_m3<<endl;
	double invM1 = dataPoint::instance()->getM(_subSys);
	/*
	 * For WignerD functions we need one more invariant mass:
	 */
	int mod=0;
	if(_subSys==5) mod=4; //5->3 work also without beta=nan, what is correct?
	if(_subSys==4) mod=5;
	if(_subSys==3) mod=5;
	double invM2 = dataPoint::instance()->getM(mod);

//	dataPoint* point = dataPoint::instance();
	//	cout<<point->getM(3)<<" " <<point->getM(4)<< " " << point->getM(5)<<endl;
	//	double locmin_sq2 = s2min(_m23*_m23,_M,_m1,_m2,_m3);
	//	double locmax_sq2 = s2max(_m23*_m23,_M,_m1,_m2,_m3);
	//	double beta2=acos((2.*_m13*_m13-locmax_sq2-locmin_sq2)/(locmax_sq2-locmin_sq2));

	  switch(_subSys){
		case 5:{ //reso in m23
		  locmin_sq = s2min(invM1*invM1,_M,_m1,_m2,_m3);
		  locmax_sq = s2max(invM1*invM1,_M,_m1,_m2,_m3);
		  break;
		}
		case 4:{ //reso in m13
		  locmin_sq = s1min(invM1*invM1,_M,_m1,_m2,_m3);
		  locmax_sq = s1max(invM1*invM1,_M,_m1,_m2,_m3);
		  break;
		}
		case 3:{ //reso in m12
		  //return 1;
		  locmin_sq = s1min(invM1*invM1,_M,_m1,_m3,_m2);
		  locmax_sq = s1max(invM1*invM1,_M,_m1,_m3,_m2);
		  //if(beta!=beta) return 1.;
		  break;
		}
	  }
	beta=acos((2.*invM2*invM2-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
	//if(_subSys!=5) beta=acos(1);

	//	cout<<"==== "<<_m23<< " "<<_m13<< " "<<invM1<<" " <<invM2<<endl;
	//	cout<< "wwww " <<_subSys<< " "<<mod<<" " <<beta<< " "<<beta2<<endl;
	//	if(beta!=beta){
	//		std::cout<<beta<<std::endl;
	//		std::cout<<_M<< " "<<_m1<<" " <<_m2<<" " <<_m3<<std::endl;
	//		std::cout<<(2.*_m13*_m13-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq)<<std::endl;
	//		std::cout<<_m13<< " " <<_m23<<" " <<locmin_sq << " "<<locmax_sq<<std::endl;
	//	}
	//double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
	//if( _x*_x>locmax_sq || _x*_x<locmin_sq )
	//  return 0.;

	Spin j((double) _inSpin), m((double) _outSpin1), n((double)_outSpin2);
	//  cout<<Wigner_d(j,m,n,beta)<<endl;
	double result = Wigner_d(j,m,n,beta);
	if( ( result!=result ) || (beta!=beta)) {
		std::cout<< "NAN! J="<< _inSpin<<" M="<<_outSpin1<<" N="<<_outSpin2<<" beta="<<beta<<std::endl;
		std::cout<< "subSys: "<<_subSys<<" ("<<mod<<") "<<invM1*invM1 << " " <<invM2*invM2<< " cos(beta)="<<(2.*invM2*invM2-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq)<<std::endl;
		return 0;
	}
	//	cout<<"result: "<<result<<endl;
	return result;
}


double AmpWigner::lambda(double x, double y, double z)const{
	return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
}

double AmpWigner::s2min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);

	return result;
}

double AmpWigner::s2max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double AmpWigner::s3min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);

	return result;
}

double AmpWigner::s3max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double AmpWigner::s1min(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);

	return result;
}

double AmpWigner::s1max(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);

	return result;
}
