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

AmpWigner2::AmpWigner2(unsigned int subSys, unsigned int spin,
		unsigned int mu, unsigned int muPrime) :
		_subSys(subSys), _spin(spin), _mu(mu), _muPrime(muPrime)
{

}

double AmpWigner2::evaluate(dataPoint& point) {
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(
			Kinematics::instance());
	return dynamicalFunction(
			_spin, _mu, _muPrime, point.getVal(_subSys)
			);
}

double AmpWigner2::dynamicalFunction(int J, int mu, int muPrime, double cosTheta){
	/* We assume that we have spin 0 particles only and the Wigner_d functions simplifies to
	 * ordinary Legendre polynomials. We normalize the square of these to one by the pre factor
	 * sqrt(2J+1). The factor was obtained by trial and error. No idea for why thats the
	 * normalization.  */
//	double norm = 1/sqrt(2*J+1);
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

std::shared_ptr<FunctionTree> AmpWigner2::SetupTree(
			allMasses& theMasses, std::string suffix){

	std::shared_ptr<MultiDouble> m23sq(
			new MultiDouble("m23sq",theMasses.masses_sq.at( std::make_pair(2,3) )) );
	std::shared_ptr<MultiDouble> m13sq(
			new MultiDouble("m13sq",theMasses.masses_sq.at( std::make_pair(1,3) )) );
	std::shared_ptr<MultiDouble> m12sq(
			new MultiDouble("m12sq",theMasses.masses_sq.at( std::make_pair(1,2) )) );


	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

	//----Strategies needed
	std::shared_ptr<WignerDStrategy> angdStrat(
			new WignerDStrategy("AngD"+suffix,ParType::MDOUBLE) );
	newTree->createHead("AngD_"+suffix, angdStrat, theMasses.nEvents);

	newTree->createLeaf("subSysFlag", _subSys, "AngD_"+suffix); //subSysFlag
	newTree->createLeaf("spin",_spin, "AngD_"+suffix); //spin
	newTree->createLeaf("m", _mu, "AngD_"+suffix); //OutSpin 1
	newTree->createLeaf("n", _muPrime, "AngD_"+suffix); //OutSpin 2
	newTree->createLeaf("m12sq", m12sq, "AngD_"+suffix); //mc
	newTree->createLeaf("m13sq", m13sq, "AngD_"+suffix); //mb
	newTree->createLeaf("m23sq", m23sq, "AngD_"+suffix); //ma

	return newTree;
}


bool WignerDStrategy::execute(ParameterList& paras,
		std::shared_ptr<AbsParameter>& out)
{
	if( checkType != out->type() ) {
		throw( WrongParType( std::string("Output Type ")
		+ParNames[out->type()]+std::string(" conflicts expected type ")
		+ParNames[checkType]+std::string(" of ")+name+" Wigner strat") );
		return false;
	}

	unsigned int _subSysFlag = (unsigned int)paras.GetDoubleParameter(0)->GetValue();
	double _inSpin = paras.GetDoubleParameter(1)->GetValue();
	double _outSpin1 = paras.GetDoubleParameter(2)->GetValue();
	double _outSpin2 = paras.GetDoubleParameter(3)->GetValue();

	double _m23,_m13,_m12;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(
			Kinematics::instance());
	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MDOUBLE){
		if(paras.GetNMultiDouble()){
			std::shared_ptr<MultiDouble> _pm12 = paras.GetMultiDouble(0);
			std::shared_ptr<MultiDouble> _pm13 = paras.GetMultiDouble(1);
			std::shared_ptr<MultiDouble> _pm23 = paras.GetMultiDouble(2);

			std::vector<double> results(_pm23->GetNValues(), 0.);
			for(unsigned int ele=0; ele<_pm23->GetNValues(); ele++){
				_m23 = double(_pm23->GetValue(ele));
				_m13 = double(_pm13->GetValue(ele));
				dataPoint point;
				try{
					Kinematics::instance()->FillDataPoint(0,1,_m23,_m13,point);
				} catch (BeyondPhsp& ex){
					BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | Data Beyond Phsp!";
					throw;
				}

				double cosTheta;
				try{
					cosTheta = kin->helicityAngle(_subSysFlag, point);
				} catch (std::exception &ex) {
					BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | "
							<<ex.what();
					throw std::runtime_error("WignerDStrategy::execute() | "
							"calculation of helicity angle failed!");
				}
				try{
					results[ele]=AmpWigner2::dynamicalFunction(
							_inSpin,_outSpin1,_outSpin2,cosTheta
					);
				} catch (std::exception &ex) {
					BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | "
							<<ex.what();
					throw std::runtime_error("WignerDStrategy::execute() | "
							"Evaluation of dynamical function failed!");
				}
			}//end element loop
			out = std::shared_ptr<AbsParameter>(
					new MultiDouble(out->GetName(),results));
			return true;
		}
	}else if(checkType == ParType::DOUBLE){ //one dim output
		_m23 = paras.GetDoubleParameter(4)->GetValue();
		_m13 = paras.GetDoubleParameter(5)->GetValue();
		dataPoint point;
		try{
			Kinematics::instance()->FillDataPoint(0,1,_m23,_m13,point);
		} catch (BeyondPhsp& ex){
			BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | Data Beyond Phsp!";
			throw;
		}

		double cosTheta;
		try{
			cosTheta = acos(kin->helicityAngle(_subSysFlag, point));
		} catch (std::exception &ex) {
			BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | "
					<<ex.what();
			throw std::runtime_error("WignerDStrategy::execute() | "
					"calculation of helicity angle failed!");
		}

		double result;
		try{
			result = AmpWigner2::dynamicalFunction(
				_inSpin,_outSpin1,_outSpin2,cosTheta
				);
		} catch (std::exception &ex) {
			BOOST_LOG_TRIVIAL(error) << "WignerDStrategy::execute() | "
					<<ex.what();
			throw std::runtime_error("WignerDStrategy::execute() | "
					"Evaluation of dynamical function failed!");
		}
		out = std::shared_ptr<AbsParameter>(
				new DoubleParameter(out->GetName(),result));
		return true;
	}
	return false;
}
