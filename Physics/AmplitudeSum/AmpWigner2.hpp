//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel
//     Peter Weidenkaff
//-------------------------------------------------------------------------------

//! Angular distribution based on WignerD functions
/*!
 * @file AmpWigner2.hpp
 *\class AmpWigner2
 *The helicity angle for sub system \_subSys is calculated and the value of the WignerD function is returned
 */

#ifndef AMPWIGNER2
#define AMPWIGNER2

#include <vector>

#include "qft++.h"

#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Functions.hpp"

using namespace std;

class WignerDStrategy : public Strategy {
public:
  WignerDStrategy(const std::string resonanceName):name(resonanceName){
    //name = +resonanceName;
  }

  virtual const std::string to_str() const {
    return ("WignerD of "+name);
  }

  virtual std::shared_ptr<AbsParameter> execute(ParameterList& paras) {

    double Gamma0, GammaV;
    double _m23 = double(paras.GetParameterValue("m23"));
    double _m13 = double(paras.GetParameterValue("m13"));
    double _m12 = double(paras.GetParameterValue("m12"));
    double _M  = double(paras.GetParameterValue("ParOfNode_M"));
    double _m1 = double(paras.GetParameterValue("ParOfNode_m1"));
    double _m2 = double(paras.GetParameterValue("ParOfNode_m2"));
    double _m3 = double(paras.GetParameterValue("ParOfNode_m3"));
    //double locmax_sq = double(paras.GetParameterValue("mb_"+name));
    //double locmin_sq = double(paras.GetParameterValue("mb_"+name));
    unsigned int _subSysFlag = (unsigned int)(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
    double _inSpin = double(paras.GetParameterValue("ParOfNode_spin_"+name));
    double _outSpin1 = double(paras.GetParameterValue("ParOfNode_m_"+name));
    double _outSpin2 = double(paras.GetParameterValue("ParOfNode_n_"+name));

    double locmin_sq, locmax_sq, beta;

    switch(_subSysFlag){
      case 5:{ //reso in m23
        locmin_sq = s2min(_m23*_m23,_M,_m1,_m2,_m3);
        locmax_sq = s2max(_m23*_m23,_M,_m1,_m2,_m3);
        beta=acos((2.*_m13*_m13-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        break;
      }
      case 4:{ //reso in m13
        locmin_sq = s1min(_m13*_m13,_M,_m1,_m2,_m3);
        locmax_sq = s1max(_m13*_m13,_M,_m1,_m2,_m3);
        beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        break;
      }
      case 3:{ //reso in m12
        //return 1;
        locmin_sq = s1min(_m12*_m12,_M,_m1,_m3,_m2);
        locmax_sq = s1max(_m12*_m12,_M,_m1,_m3,_m2);
        beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        if(beta!=beta) return NULL;
        break;
      }
    }

    //double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
    //if( _x*_x>locmax_sq || _x*_x<locmin_sq )
    //  return 0.;

    Spin j(_inSpin), m(_outSpin1), n(_outSpin2);
    std::shared_ptr<DoubleParameter> result = std::shared_ptr<DoubleParameter>(new DoubleParameter("WignerD of "+name+" result",Wigner_d(j,m,n,beta)));
    return result;
  }

protected:
  std::string name;

  double lambda(double x, double y, double z)const{
    return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
  }

  double s2min(double s1, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	  double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);

  	return result;
  }

  double s2max(double s1, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	  double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);

  	return result;
  }

  double s3min(double s1, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	  double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);

  	return result;
  }

  double s3max(double s1, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	  double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);

  	return result;
  }

  double s1min(double s2, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );

	  double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);

  	return result;
  }

  double s1max(double s2, double m0, double m1, double m2, double m3)const
  {
	  double s      = m0*m0;
	  double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );

	  double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);

  	return result;
  }

};

class AmpWigner2{
public:
  AmpWigner2( unsigned int subSys, unsigned int resSpin);

  AmpWigner2(const AmpWigner2&, const char*);

  virtual ~AmpWigner2() {};

  virtual double evaluate() const;
  //double evaluateTree(const ParameterList& paras, const std::string name) const;
  void setDecayMasses(double, double, double, double);

protected:

  virtual void initialise();

  unsigned int _resSpin;
  unsigned int _subSys;

  double _M;
  double _m1;
  double _m2;
  double _m3;

 unsigned int _spinM;
 unsigned int _spin1;
 unsigned int _spin2;
 unsigned int _spin3;
private:

};

#endif
