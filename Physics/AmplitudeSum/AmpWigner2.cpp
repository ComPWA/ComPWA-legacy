//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff  -
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

#include <cmath>
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include "qft++.h"
#include <boost/math/special_functions/legendre.hpp>

AmpWigner2::AmpWigner2(unsigned int subSys, unsigned int spin) : _subSys(subSys), _spin(spin) { }

double AmpWigner2::evaluate(dataPoint& point) {
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	return dynamicalFunction(_spin,0,0,kin->helicityAngle(_subSys, point));
}

double AmpWigner2::dynamicalFunction(int J, int mu, int muPrime, double cosTheta){
	/* We assume that we have spin 0 particles only and the Wigner_d functions simplifies to
	 * ordinary Legendre polynomials. We normalize the square of these to one by the pre factor
	 * sqrt(2J+1). The factor was obtained by trial and error. No idea for why thats the
	 * normalization.  */
	//	double norm = sqrt(2*J+1);
	double norm = 1;
	if(J==0) return norm; //assure that angular term is properly normalized
	if(cosTheta>1 ||cosTheta<-1)
		throw std::runtime_error( "AmpWigner2::dynamicalFunction() | "
				"scattering angle out of range! Datapoint beyond phsp?");
	// Calling WignerD function (boost and qft++ version give same results)
	double result = Wigner_d( Spin(J), Spin(mu), Spin(muPrime), acos(cosTheta) );
	//	double result = boost::math::legendre_p<double>(J,cosTheta);
	if( ( result!=result ) )
		throw std::runtime_error("AmpWigner2::evaluate() | Result is NaN!");

	return (norm*(2*J+1)*result);
}

bool WignerDStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw( WrongParType( std::string("Output Type ")
		+ParNames[out->type()]+std::string(" conflicts expected type ")
		+ParNames[checkType]+std::string(" of ")+name+" Wigner strat") );
		return false;
	}

	unsigned int _subSysFlag = (unsigned int)(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
	double _inSpin = double(paras.GetParameterValue("ParOfNode_spin_"+name));
	double _outSpin1 = double(paras.GetParameterValue("ParOfNode_m_"+name));
	double _outSpin2 = double(paras.GetParameterValue("ParOfNode_n_"+name));

	double _m23,_m13,_m12;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MDOUBLE){
		if(paras.GetNMultiDouble()){
			unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
			std::shared_ptr<MultiDouble> _pm23 = paras.GetMultiDouble("m23sq");
			std::shared_ptr<MultiDouble> _pm13 = paras.GetMultiDouble("m13sq");

			std::vector<double> results(nElements, 0.);
			for(unsigned int ele=0; ele<nElements; ele++){
				_m23 = double(_pm23->GetValue(ele));
				_m13 = double(_pm13->GetValue(ele));
				dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
				double theta = acos(kin->helicityAngle(_subSysFlag, point));
				//results[ele]=(2*_inSpin+1)*
				//		Wigner_d(Spin(_inSpin),Spin(_outSpin1),Spin(_outSpin2),theta);
				results[ele]=AmpWigner2::dynamicalFunction(
						_inSpin,_outSpin1,_outSpin2,kin->helicityAngle(_subSysFlag, point));
			}//end element loop
			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
			return true;
		}
	}else if(checkType == ParType::DOUBLE){ //one dim output
		_m23 = double(paras.GetParameterValue("m23sq"));
		_m13 = double(paras.GetParameterValue("m13sq"));
		dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
		double result = AmpWigner2::dynamicalFunction(
				_inSpin,_outSpin1,_outSpin2,kin->helicityAngle(_subSysFlag, point));
		out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),result));
		return true;
	}
	return false;
}
bool WignerDPhspStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw( WrongParType( std::string("Output Type ")
		+ParNames[out->type()]+std::string(" conflicts expected type ")
		+ParNames[checkType]+std::string(" of ")+name+" Wigner strat") );
		return false;
	}

	unsigned int _subSysFlag = (unsigned int)(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
	double _inSpin = double(paras.GetParameterValue("ParOfNode_spin_"+name));
	double _outSpin1 = double(paras.GetParameterValue("ParOfNode_m_"+name));
	double _outSpin2 = double(paras.GetParameterValue("ParOfNode_n_"+name));

	double _m23,_m13,_m12;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MDOUBLE){
		if(paras.GetNMultiDouble()){
			unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
			std::shared_ptr<MultiDouble> _pm23 = paras.GetMultiDouble("m23sq_phsp");
			std::shared_ptr<MultiDouble> _pm13 = paras.GetMultiDouble("m13sq_phsp");

			std::vector<double> results(nElements, 0.);
			for(unsigned int ele=0; ele<nElements; ele++){
				_m23 = double(_pm23->GetValue(ele));
				_m13 = double(_pm13->GetValue(ele));
				dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
				double theta = acos(kin->helicityAngle(_subSysFlag, point));
				//results[ele]=(2*_inSpin+1)*
				//		Wigner_d(Spin(_inSpin),Spin(_outSpin1),Spin(_outSpin2),theta);
				results[ele]=AmpWigner2::dynamicalFunction(
						_inSpin,_outSpin1,_outSpin2,kin->helicityAngle(_subSysFlag, point));
			}//end element loop
			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
			return true;
		}
	}else if(checkType == ParType::DOUBLE){ //one dim output
		_m23 = double(paras.GetParameterValue("m23sq_phsp"));
		_m13 = double(paras.GetParameterValue("m13sq_phsp"));
		dataPoint point; point.setVal(0,_m23); point.setVal(1,_m13);
		double result = AmpWigner2::dynamicalFunction(
				_inSpin,_outSpin1,_outSpin2,kin->helicityAngle(_subSysFlag, point));
		out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),result));
		return true;
	}
	return false;
}
