//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Wrapper of the Minuit2 Optimizer library.
/*! \class MinuitIF
 * @file MinuitIF.hpp
 * This class provides a wrapper around the Minuit2 library. It fulfills the
 * Optimizer interface to be easily adapted to other modules. The data needs to
 * be provided with the ControlParameter interface.
 */

#ifndef _OIFMINUIT_HPP
#define _OIFMINUIT_HPP

#include <vector>
//#include <boost/shared_ptr.hpp>
#include <memory>

#include "Optimizer/ControlParameter.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Minuit2/MnStrategy.h"
#include <fstream>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>


using namespace ROOT::Minuit2;

class MinuitIF : public Optimizer {

public:
	/// Default Constructor (0x0)
	MinuitIF(std::shared_ptr<ControlParameter> esti, ParameterList& par);
	virtual std::shared_ptr<FitResult> exec(ParameterList& par);

	/** Destructor */
	virtual ~MinuitIF();

protected:

private:
	MinuitFcn _myFcn;
	std::shared_ptr<ControlParameter> estimator;
	// vector<string> paramNames;
};
class MinuitStrategy : public MnStrategy
{
public:
	MinuitStrategy(unsigned int i=1) : MnStrategy(i) {
		fGradNCyc = GradientNCycles();
		fGradTlrStp = GradientStepTolerance();
		fGradTlr = GradientTolerance();
		fHessNCyc = HessianNCycles();
		fHessTlrStp = HessianStepTolerance();
		fHessTlrG2 = HessianG2Tolerance();
		fHessGradNCyc = HessianGradientNCycles();
	};
	void init(){
		SetGradientNCycles(fGradNCyc);
		SetGradientStepTolerance(fGradTlrStp);
		SetGradientTolerance(fGradTlr);
		SetHessianNCycles(fHessNCyc);
		SetHessianStepTolerance(fHessTlrStp);
		SetHessianG2Tolerance(fHessTlrG2);
		SetHessianGradientNCycles(fHessGradNCyc);
	}

private:
	unsigned int fGradNCyc;
	double fGradTlrStp;
	double fGradTlr;
	unsigned int fHessNCyc;
	double fHessTlrStp;
	double fHessTlrG2;
	unsigned int fHessGradNCyc;

	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		//		ar & BOOST_SERIALIZATION_NVP(fStrategy);
		ar & BOOST_SERIALIZATION_NVP(fGradNCyc);
		ar & BOOST_SERIALIZATION_NVP(fGradTlrStp);
		ar & BOOST_SERIALIZATION_NVP(fGradTlr);
		ar & BOOST_SERIALIZATION_NVP(fHessNCyc);
		ar & BOOST_SERIALIZATION_NVP(fHessTlrStp);
		ar & BOOST_SERIALIZATION_NVP(fHessTlrG2);
		ar & BOOST_SERIALIZATION_NVP(fHessGradNCyc);
	}
};

#endif /* _OIFMinuit_HPP */
