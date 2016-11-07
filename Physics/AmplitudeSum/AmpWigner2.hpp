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
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/Utility.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

class AmpWigner2 {
public:
	AmpWigner2( unsigned int varId=0, ComPWA::Spin spin = ComPWA::Spin(0),
			unsigned int mu=0, unsigned int muPrime=0);

	virtual ~AmpWigner2() {};

	virtual double evaluate(dataPoint& point) ;

	static double dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu,
			ComPWA::Spin muPrime, double cosTheta);

	static std::complex<double> dynamicalFunction(double cosAlpha,
			double cosBeta, double cosGamma, ComPWA::Spin J, ComPWA::Spin mu,
			ComPWA::Spin muPrime);

	virtual std::shared_ptr<FunctionTree> SetupTree(
			ParameterList& sample, std::string suffix="");

	virtual unsigned int GetVarId() const { return _varId; };

	virtual void SetVarId(unsigned int id) { _varId = id; };


protected:
	unsigned int _varId;
	ComPWA::Spin _spin;
	unsigned int _mu;
	unsigned int _muPrime;
};

class WignerDStrategy: public Strategy {
public:
	WignerDStrategy(const std::string resonanceName) :
		Strategy(ParType::MDOUBLE), name(resonanceName) { }

	virtual const std::string to_str() const { return ("WignerD of "+name);	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
