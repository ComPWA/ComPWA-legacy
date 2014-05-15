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
#include <memory>

#include "qft++.h"

#include "Core/ParameterList.hpp"
#include "Core/Functions.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/DataPoint.hpp"

//using namespace std;
class AmpWigner2{
public:
	AmpWigner2( unsigned int subSys, unsigned int resSpin);

	AmpWigner2(const AmpWigner2&, const char*);

	virtual ~AmpWigner2() {};

	virtual double evaluate(dataPoint& point) ;
	void setDecayMasses(double, double, double, double);
	//  static double theta(Spin J, double m23sq, double m13sq, double m12sq, unsigned int subSys, const double& M, const double& m1, const double& m2, const double& m3){
	//    double cosTheta=-999;
	//    DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//    dataPoint p; p.setVal(0,m23sq); p.setVal(1,m13sq);
	//    cosTheta = kin->calcHelicityAngle(subSys,p);
	//    switch(subSys){
	//    case 3:
	////      cosTheta = kin->calcHelicityAngle(point->getMsq(3),point->getMsq(4),_M,_m3,_m1,_m2);
	//        cosTheta = kin->calcHelicityAngle(m12sq,m23sq,M,m3,m1,m2);
	//        break;
	//    case 4:
	//        cosTheta = kin->calcHelicityAngle(m13sq,m23sq,M,m2,m1,m3); //angle between 2 and 1
	////        cosTheta = kin->calcHelicityAngle(m13sq,m23sq,M,m2,m3,m1); //angle between 2 and 3
	//        break;
	//    case 5:
	//        cosTheta = kin->calcHelicityAngle(m23sq,m13sq,M,m1,m2,m3);
	//        break;
	//    default:
	//        BOOST_LOG_TRIVIAL(fatal)<<"AmpWigner2: wrong subSystem! Exit!"; exit(1);
	//    }
	//    if(cosTheta>1.) cosTheta=1.;
	//    if(cosTheta<-1.) cosTheta=-1.;
	//    return std::acos(cosTheta);
	//  }

protected:
	bool massIdsSet;
	unsigned int id23;
	unsigned int id13;

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

class WignerDStrategy : public Strategy {
public:
	WignerDStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){
		//name = +resonanceName;
	}

	virtual const std::string to_str() const {
		return ("WignerD of "+name);
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" Wigner strat"));
			return false;
		}

		//		double _M  = double(paras.GetParameterValue("ParOfNode_M"));
		//		double _m1 = double(paras.GetParameterValue("ParOfNode_m1"));
		//		double _m2 = double(paras.GetParameterValue("ParOfNode_m2"));
		//		double _m3 = double(paras.GetParameterValue("ParOfNode_m3"));
		//double locmax_sq = double(paras.GetParameterValue("mb_"+name));
		//double locmin_sq = double(paras.GetParameterValue("mb_"+name));
		unsigned int _subSysFlag = (unsigned int)(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		double _inSpin = double(paras.GetParameterValue("ParOfNode_spin_"+name));
		double _outSpin1 = double(paras.GetParameterValue("ParOfNode_m_"+name));
		double _outSpin2 = double(paras.GetParameterValue("ParOfNode_n_"+name));

		double beta;
		double _m23,_m13,_m12;

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MDOUBLE){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
				std::shared_ptr<MultiDouble> _pm23 = paras.GetMultiDouble("m23sq");
				std::shared_ptr<MultiDouble> _pm13 = paras.GetMultiDouble("m13sq");
				std::shared_ptr<MultiDouble> _pm12 = paras.GetMultiDouble("m12sq");

				std::vector<double> results(nElements, 0.);

				for(unsigned int ele=0; ele<nElements; ele++){
					_m23 = double(_pm23->GetValue(ele));
					_m13 = double(_pm13->GetValue(ele));
					_m12 = double(_pm12->GetValue(ele));
					Spin j(_inSpin), m(_outSpin1), n(_outSpin2);

					/*switch(_subSysFlag){
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
              //if(beta!=beta) return false;
              break;
            }
          }*/
					//          beta = AmpWigner2::theta(j, _m23, _m13, _m12, _subSysFlag, _M, _m1, _m2, _m3);
					dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
					DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
					beta = std::acos(kin->helicityAngle(_subSysFlag, point));
					if(beta!=beta) beta=1.;

					results[ele]=(2*j+1)*Wigner_d(j,m,n,beta);
				}//end element loop

				out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
				return true;
			}
		}else if(checkType == ParType::DOUBLE){ //one dim output
			_m23 = double(paras.GetParameterValue("m23sq"));
			_m13 = double(paras.GetParameterValue("m13sq"));
			_m12 = double(paras.GetParameterValue("m12sq"));
			Spin j(_inSpin), m(_outSpin1), n(_outSpin2);

			/*switch(_subSysFlag){
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
          if(beta!=beta) return false;
          break;
        }
      }*/

			//double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
			//if( _x*_x>locmax_sq || _x*_x<locmin_sq )
			//  return 0.;
			//      beta = AmpWigner2::theta(j, _m23, _m13, _m12, _subSysFlag, _M, _m1, _m2, _m3);
			dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
			DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
			beta = std::acos(kin->helicityAngle(_subSysFlag, point));
			if(beta!=beta) beta=1.;

			out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),(2*j+1)*Wigner_d(j,m,n,beta)));
			return true;
		}else{
			return false;
		}
	}

protected:
	std::string name;



	//  double lambda(double x, double y, double z)const{
	//    return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
	//  }
	//
	//  double s2min(double s1, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
	//
	//	  double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);
	//
	//  	return result;
	//  }
	//
	//  double s2max(double s1, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
	//
	//	  double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);
	//
	//  	return result;
	//  }
	//
	//  double s3min(double s1, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
	//
	//	  double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);
	//
	//  	return result;
	//  }
	//
	//  double s3max(double s1, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
	//
	//	  double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);
	//
	//  	return result;
	//  }
	//
	//  double s1min(double s2, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );
	//
	//	  double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);
	//
	//  	return result;
	//  }
	//
	//  double s1max(double s2, double m0, double m1, double m2, double m3)const
	//  {
	//	  double s      = m0*m0;
	//	  double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );
	//
	//	  double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);
	//
	//  	return result;
	//  }

};
class WignerDphspStrategy : public Strategy {
public:
	WignerDphspStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){
		//name = +resonanceName;
	}

	virtual const std::string to_str() const {
		return ("WignerD of "+name);
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" Wigner strat"));
			return false;
		}

//		double _M  = double(paras.GetParameterValue("ParOfNode_M"));
//		double _m1 = double(paras.GetParameterValue("ParOfNode_m1"));
//		double _m2 = double(paras.GetParameterValue("ParOfNode_m2"));
//		double _m3 = double(paras.GetParameterValue("ParOfNode_m3"));
		unsigned int _subSysFlag = (unsigned int)(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		double _inSpin = double(paras.GetParameterValue("ParOfNode_spin_"+name));
		double _outSpin1 = double(paras.GetParameterValue("ParOfNode_m_"+name));
		double _outSpin2 = double(paras.GetParameterValue("ParOfNode_n_"+name));
		double _m23,_m13,_m12;

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MDOUBLE){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
				std::shared_ptr<MultiDouble> _pm23 = paras.GetMultiDouble("m23sq_phsp");
				std::shared_ptr<MultiDouble> _pm13 = paras.GetMultiDouble("m13sq_phsp");
				std::shared_ptr<MultiDouble> _pm12 = paras.GetMultiDouble("m12sq_phsp");

				std::vector<double> results(nElements, 0.);

				for(unsigned int ele=0; ele<nElements; ele++){
					_m23 = double(_pm23->GetValue(ele));
					_m13 = double(_pm13->GetValue(ele));
					_m12 = double(_pm12->GetValue(ele));
					Spin j(_inSpin), m(_outSpin1), n(_outSpin2);

					dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
					DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
					double beta = std::acos(kin->helicityAngle(_subSysFlag, point));
					//          double beta = AmpWigner2::theta(j, _m23, _m13, _m12, _subSysFlag, _M, _m1, _m2, _m3);
					if(beta!=beta) beta=1.;

					results[ele]=Wigner_d(j,m,n,beta);
				}//end element loop
				out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
				return true;
			}
		}else if(checkType == ParType::DOUBLE){ //one dim output
			_m23 = double(paras.GetParameterValue("m23sq_phsp"));
			_m13 = double(paras.GetParameterValue("m13sq_phsp"));
			_m12 = double(paras.GetParameterValue("m12sq_phsp"));
			Spin j(_inSpin), m(_outSpin1), n(_outSpin2);
			dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
			DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
			double beta = std::acos(kin->helicityAngle(_subSysFlag, point));
			//      double beta = AmpWigner2::theta(j, _m23, _m13, _m12, _subSysFlag, _M, _m1, _m2, _m3);
			if(beta!=beta) beta=1.;
			out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),Wigner_d(j,m,n,beta)));
			return true;
		}else{
			return false;
		}
	}

protected:
	std::string name;
};


#endif
