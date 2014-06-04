//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

//#include "TObject.h"
//#include "TString.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

	AmpRelBreitWignerRes(const char *name,
			DoubleParameter& _resMass, DoubleParameter& _resWidth,
			double& _radius,
			int _subsys,
			int resSpin, int m, int n
	) ;
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

	virtual ~AmpRelBreitWignerRes();

	//  double operator() (double *x, size_t dim, void*);

	static std::complex<double> dynamicalFunction(double mSq, double mR, double ma, double mb, double gamma0, unsigned int J, double mesonRadius);
	virtual void initialise();
	virtual std::complex<double> evaluate(dataPoint& point) { return _norm*evaluateAmp(point)*evaluateWignerD(point); }
	virtual std::complex<double> evaluateAmp(dataPoint& point);
	virtual double evaluateWignerD(dataPoint& point) {
		if(_spin==0) return 1.0;//save some computing time
		return _wignerD.evaluate(point);
	};

	void setDecayMasses(double, double, double, double);
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };

protected:
	DoubleParameter& _resWidth;
	AmpWigner2 _wignerD;
	bool foundMasses;
	unsigned int id23, id13;
	//	AmpWigner _wignerD;

private:

};

class BreitWignerStrategy : public Strategy {
public:
	BreitWignerStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}

	virtual const std::string to_str() const {
		return ("relativistic BreitWigner of "+name);
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double Gamma0, m0, d;
		unsigned int spin, subSys;
		try{
			m0 = double(paras.GetParameterValue("mass_"+name));
		}catch(BadParameter& e){
			m0 = double(paras.GetParameterValue("ParOfNode_m0_"+name));
		}
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		d = double(paras.GetParameterValue("ParOfNode_d_"+name));
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));

		//m  = double(paras.GetParameterValue("mym"));
		double ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		double mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));

		try{
			Gamma0 = double(paras.GetParameterValue("width_"+name));
		}catch(BadParameter& e){
			Gamma0 = double(paras.GetParameterValue("ParOfNode_resWidth_"+name));
		}
		/*GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);

    std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
    std::complex<double> res(m0 * Gamma0);
    res = res / denom;*/

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ //reso in sys of particles 1&2
					mp  = (paras.GetMultiDouble("m12sq"));
					//map  = (paras.GetMultiDouble("m23"));
					// mbp  = (paras.GetMultiDouble("m13"));
					break;
				}
				case 4:{ //reso in sys of particles 1&3
					mp  = (paras.GetMultiDouble("m13sq"));
					// map  = (paras.GetMultiDouble("m12"));
					// mbp  = (paras.GetMultiDouble("m23"));
					break;
				}
				case 5:{ //reso in sys of particles 2&3
					mp  = (paras.GetMultiDouble("m23sq"));
					// map  = (paras.GetMultiDouble("m13"));
					// mbp  = (paras.GetMultiDouble("m12"));
					break;
				}
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
//					double norm =1.0;
//					double BLWeiss2 = BLprime2(ma,mb,m0,sqrt(mSq),d,spin);
//					double qTerm = std::pow( ( q(ma,mb,sqrt(mSq)) / q(ma,mb,m0) ) , (2.*spin + 1.) );
//					//Gamma0 = 1;
//					//if(ele==0) std::cout << " DEBUG  " << q(map->GetValue(ele),mbp->GetValue(ele),mp->GetValue(ele)) << " " << q0(map->GetValue(ele),mbp->GetValue(ele),m0) << std::endl;
//					GammaV = Gamma0 * qTerm * (m0 / sqrt(mSq)) * BLWeiss2;
//					std::complex<double> denom(m0*m0 - mSq, -m0 * GammaV);
//					results[ele] = (std::complex<double>(norm*(2*spin+1))) / denom; //Laura++ (old) definition*/

					results[ele] = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
//					if(ele<10) std::cout<<"Strategy BWrel "<<results[ele]<<std::endl;
				}

				//std::vector<std::complex<double> > resultsTMP(nElements, std::complex<double>(1.));
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output


		//Only StandardDim Paras in input
		//  double spinTerm = evaluateWignerD(); //spinTerm =1;
		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			mSq  = (double(paras.GetParameterValue("m12sq")));
			//ma  = (double(paras.GetParameterValue("m23")));
			// mb  = (double(paras.GetParameterValue("m13")));
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			mSq  = (double(paras.GetParameterValue("m13sq")));
			// ma  = (double(paras.GetParameterValue("m12")));
			// mb  = (double(paras.GetParameterValue("m23")));
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			mSq  = (double(paras.GetParameterValue("m23sq")));
			// ma  = (double(paras.GetParameterValue("m13")));
			// mb  = (double(paras.GetParameterValue("m12")));
			break;
		}
		}
		//		BLWeiss2 = BLprime2(ma,mb,m0,m,d,spin);
		//		qTerm = std::pow(q(ma,mb,m) / q(ma,mb,m0), (2.*spin + 1.));
		//		//double Gamma0 = _resWidth.GetValue();
		//		GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
		//		std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
		//
		//		std::complex<double> result = std::complex<double>(norm) / denom; //Laura++ (old) definition*/

		//std::complex<double> result (res.re(),res.im());
		std::complex<double> result = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;

//	double q(const double& ma, const double& mb, const double& x) const {
//		double mapb = ma + mb;
//		double mamb = ma - mb;
//
//		if( (x*x - mapb*mapb) < 0 ) {
//			//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
//			return 1; //below threshold
//		}
//
//		return std::sqrt ( (x*x - mapb*mapb) * (x*x - mamb*mamb) ) / (2. * x );
//	}

	// compute part of the Blatt-Weisskopf barrier factor
	//   BLprime = sqrt (F(q0)/F(q))
	/* double F(const double& p, const double& d, unsigned int& spin) const {
    double retVal = 1;

    if (spin == 0)
      retVal = 1;
    else if (spin == 1)
      retVal = 1 + p*p * d*d;
    else if (spin == 2) {
      double z = p*p * d*d;
      retVal = (z-3.)*(z-3.) + 9*z;
    }
    return retVal;
  }*/

	// compute square of Blatt-Weisskopf barrier factor
//	double BLprime2(const double& ma, const double& mb, const double& m0, const double& x, const double& d, unsigned int& spin) const {
//		double t0= q(ma, mb, m0)*q(ma, mb, m0) * d*d;
//		double t= q(ma, mb, x)*q(ma, mb, x) * d*d;
//		return FormFactor(t0,t,spin);
//	}
//
//	double FormFactor(double& z0, double& z, unsigned int& spin) const{
//		double nom=0, denom=0;
//		switch(spin){
//		case 0:{
//			return 1.;
//		}
//		case 1:{
//			//if(_type==barrierType::BWPrime){
//			nom = 1 + z0;
//			denom = 1 + z;
//			//} else if(_type==barrierType::BW){
//			//   nom = 2*z;
//			//   denom = 1 + z;
//			//} else {
//			//   std::cout<<"Wrong BLW factor definition: "<<_type<<std::endl;
//			//   return 1;
//			//}
//			break;
//		}
//		case 2:{
//			//if(_type==barrierType::BWPrime){
//			nom = (z0-3)*(z0-3)+9*z0;
//			denom = (z-3)*(z-3)+9*z;
//			// } else if(_type==barrierType::BW){
//			//  nom = 13*z*z;
//			//   denom = (z-3)*(z-3)+9*z;
//			// } else {
//			//    std::cout<<"Wrong BLW factor definition: "<<_type<<std::endl;
//			//    return 1;
//			// }
//			break;
//		}
//		default:{
//			std::cout<<"Wrong spin value! BLW factors only implemented for spin 0,1 and 2! "<<std::endl;
//			return 0;
//		}
//		}
//		return nom/denom;
//	}

};
class BreitWignerPhspStrategy : public BreitWignerStrategy {
public:
	BreitWignerPhspStrategy(const std::string resonanceName, ParType in):BreitWignerStrategy(resonanceName,in){}
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double Gamma0, m0, d;
		unsigned int spin, subSys;
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			m0 = double(paras.GetParameterValue("ParOfNode_m0_"+name));
		}
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		d = double(paras.GetParameterValue("ParOfNode_d_"+name));
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));

		//m  = double(paras.GetParameterValue("mym"));
		double ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		double mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));

		try{
			Gamma0 = double(paras.GetParameterValue("resWidth_"+name));
		}catch(BadParameter& e){
			Gamma0 = double(paras.GetParameterValue("ParOfNode_resWidth_"+name));
		}
		/*GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);

    std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
    std::complex<double> res(m0 * Gamma0);
    res = res / denom;*/

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ //reso in sys of particles 1&2
					mp  = (paras.GetMultiDouble("m12sq_phsp"));
					//map  = (paras.GetMultiDouble("m23"));
					// mbp  = (paras.GetMultiDouble("m13"));
					break;
				}
				case 4:{ //reso in sys of particles 1&3
					mp  = (paras.GetMultiDouble("m13sq_phsp"));
					// map  = (paras.GetMultiDouble("m12"));
					// mbp  = (paras.GetMultiDouble("m23"));
					break;
				}
				case 5:{ //reso in sys of particles 2&3
					mp  = (paras.GetMultiDouble("m23sq_phsp"));
					// map  = (paras.GetMultiDouble("m13"));
					// mbp  = (paras.GetMultiDouble("m12"));
					break;
				}
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
					//          BLWeiss2 = BLprime2(ma,mb,m0,m,d,spin);
					//          qTerm = std::pow( ( q(ma,mb,m) / q(ma,mb,m0) ) , (2.*spin + 1.) );
					//          //Gamma0 = 1;
					//          //if(ele==0) std::cout << " DEBUG  " << q(map->GetValue(ele),mbp->GetValue(ele),mp->GetValue(ele)) << " " << q0(map->GetValue(ele),mbp->GetValue(ele),m0) << std::endl;
					//          GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
					//
					//          std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
					//
					//          results[ele] = (std::complex<double>(norm*(2*spin+1))) / denom; //Laura++ (old) definition*/
					results[ele] = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
//					if(ele<10) std::cout<<"Strategy BWrelPHSP ("<<ele<<"/"<<nElements<<") "<<results[ele]<<std::endl;
				}
				//std::vector<std::complex<double> > resultsTMP(nElements, std::complex<double>(1.));
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output


		//Only StandardDim Paras in input
		//  double spinTerm = evaluateWignerD(); //spinTerm =1;
		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			mSq  = (double(paras.GetParameterValue("m12sq_phsp")));
			//ma  = (double(paras.GetParameterValue("m23")));
			// mb  = (double(paras.GetParameterValue("m13")));
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			mSq  = (double(paras.GetParameterValue("m13sq_phsp")));
			// ma  = (double(paras.GetParameterValue("m12")));
			// mb  = (double(paras.GetParameterValue("m23")));
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			mSq  = (double(paras.GetParameterValue("m23sq_phsp")));
			// ma  = (double(paras.GetParameterValue("m13")));
			// mb  = (double(paras.GetParameterValue("m12")));
			break;
		}
		}
		//    BLWeiss2 = BLprime2(ma,mb,m0,m,d,spin);
		//    qTerm = std::pow(q(ma,mb,m) / q(ma,mb,m0), (2.*spin + 1.));
		//    //double Gamma0 = _resWidth.GetValue();
		//    GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
		//    std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
		//
		//    std::complex<double> result = std::complex<double>(norm) / denom; //Laura++ (old) definition*/
		//std::complex<double> result (res.re(),res.im());
		std::complex<double> result = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
		std::cout<<"Strategy BWrelPHSP single value "<<result<<std::endl;

		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));

		return true;
	}
}; // end BreitwignerPhspStrategy


#endif
