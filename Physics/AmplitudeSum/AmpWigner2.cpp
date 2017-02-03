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

#include <boost/math/special_functions/legendre.hpp>

//#include "qft++.h"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {
AmpWigner2::AmpWigner2(unsigned int varId, ComPWA::Spin spin,
		unsigned int mu, unsigned int muPrime) :
		_varId(varId), _spin(spin), _mu(mu), _muPrime(muPrime)
{

}

double AmpWigner2::evaluate(dataPoint& point) {
	ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
			dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(
					Kinematics::instance()
	);
	return dynamicalFunction(
			_spin.Val(), _mu, _muPrime, point.getVal(_varId)
			);
}

double AmpWigner2::dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu, ComPWA::Spin muPrime, double cosTheta)
{
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

    
	double result = 1.0;
	//TODO: this is the only point in which we use qft++, we should be rid of it
//	::Spin _J(J.GetNumerator(), J.GetDenominator());
//	::Spin _mu(mu.GetNumerator(), mu.GetDenominator());
//	::Spin _muPrime(muPrime.GetNumerator(), muPrime.GetDenominator());
//    	result = Wigner_d( _J, _mu, _muPrime, acos(cosTheta) );

#ifdef DEBUG
	if( ( result!=result ) )
		throw
		std::runtime_error("AmpWigner2::dynamicalFunction() | Result is NaN!");
#endif

	return (norm*(2*J.GetSpin()+1)*result);
}

std::complex<double> AmpWigner2::dynamicalFunction(double cosAlpha,
		double cosBeta, double cosGamma, ComPWA::Spin J, ComPWA::Spin mu,
		ComPWA::Spin muPrime )
{
	if(J==0) return 1.0;

	std::complex<double> i(0,1);

	double tmp = AmpWigner2::dynamicalFunction(J, mu, muPrime, cosBeta);
	std::complex<double> result =
			tmp * std::exp( -i*(mu.GetSpin()*acos(cosAlpha)
					+ muPrime.GetSpin()*acos(cosGamma)) );

#ifdef DEBUG
	if( ( result.real != result.real || result.imag != result.imag ) )
		throw
		std::runtime_error("AmpWigner2::dynamicalFunction() | Result is NaN!");
#endif

	return result;
}

std::shared_ptr<FunctionTree> AmpWigner2::SetupTree(
			ParameterList& sample, std::string suffix)
{
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

	int sampleSize = sample.GetMultiDouble(0)->GetNValues();
	//----Strategies needed
	std::shared_ptr<WignerDStrategy> angdStrat(
			new WignerDStrategy("AngD"+suffix) );
	newTree->createHead("AngD_"+suffix, angdStrat, sampleSize);

	newTree->createLeaf("spin",_spin.Val(), "AngD_"+suffix); //spin
	newTree->createLeaf("m", _mu, "AngD_"+suffix); //OutSpin 1
	newTree->createLeaf("n", _muPrime, "AngD_"+suffix); //OutSpin 2
	newTree->createLeaf("AngD_sample", sample.GetMultiDouble(_varId), "AngD_"+suffix);

	return newTree;
}

bool WignerDStrategy::execute(ParameterList& paras,
		std::shared_ptr<AbsParameter>& out)
{
#ifdef DEBUG
	if( checkType != out->type() ) {
		throw( WrongParType( std::string("Output Type ")
		+ParNames[out->type()]+std::string(" conflicts expected type ")
		+ParNames[checkType]+std::string(" of ")+name+" Wigner strat") );
		return false;
	}
#endif

	double _inSpin = paras.GetDoubleParameter(0)->GetValue();
	double _outSpin1 = paras.GetDoubleParameter(1)->GetValue();
	double _outSpin2 = paras.GetDoubleParameter(2)->GetValue();

	ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
			dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(
			Kinematics::instance()
	);

	std::shared_ptr<MultiDouble> _angle = paras.GetMultiDouble(0);

	std::vector<double> results(_angle->GetNValues(), 0.);
	for(unsigned int ele=0; ele<_angle->GetNValues(); ele++){
		try{
			results.at(ele)=AmpWigner2::dynamicalFunction(
					_inSpin,_outSpin1,_outSpin2,_angle->GetValue(ele)
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

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
